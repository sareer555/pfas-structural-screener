"""
pfas_screener.py

Core PFAS (per- and polyfluoroalkyl substance) structural screening logic.

Implements the OECD (2021) structural definition:
    "PFASs are defined as fluorinated substances that contain at least one
    fully fluorinated methyl (CF3-) or methylene (-CF2-) carbon atom
    (without any H/Cl/Br/I atom attached to it)."
    -- OECD, "Reconciling Terminology of the Universe of Per- and
       Polyfluoroalkyl Substances", Series on Risk Management No. 61,
       ENV/CBC/MONO(2021)25.

This is a *structural* screen, not a lookup against a curated regulatory
list. It will flag any molecule containing a qualifying -CF3 or -CF2-
carbon, whether or not that specific substance has been individually
named by a regulator. It deliberately does NOT flag:
  - aromatic C-F bonds (e.g. fluorobenzene) -- these are not methyl/
    methylene carbons in the OECD sense,
  - partially fluorinated carbons that still carry an H (e.g. -CHF2,
    -CH2F) -- not "fully fluorinated",
  - carbons that mix F with Cl/Br/I (e.g. CFC/HCFC refrigerants like
    CCl2F2) -- explicitly excluded by the "without any H/Cl/Br/I" clause.

Use resolve_to_smiles() to turn a name/CAS number into a SMILES string via
PubChem, then classify_pfas() to run the structural screen on that SMILES.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Literal

import requests
from rdkit import Chem

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

IdentifierType = Literal["name", "cas", "smiles"]


class ResolutionError(Exception):
    """Raised when a name/CAS/SMILES identifier could not be resolved to a valid structure."""


@dataclass
class MatchedGroup:
    group_type: Literal["CF3", "CF2"]
    carbon_atom_idx: int
    fluorine_atom_idxs: list[int]


@dataclass
class PfasResult:
    input_identifier: str
    identifier_type: IdentifierType
    resolved_smiles: str | None
    is_pfas: bool
    matched_groups: list[MatchedGroup] = field(default_factory=list)
    n_cf3_groups: int = 0
    n_cf2_groups: int = 0
    highlight_atoms: list[int] = field(default_factory=list)
    highlight_bonds: list[int] = field(default_factory=list)
    note: str = ""
    error: str | None = None

    @property
    def flag_label(self) -> str:
        if self.error:
            return "ERROR"
        return "PFAS (structural match)" if self.is_pfas else "Not PFAS"


def resolve_to_smiles(identifier: str, identifier_type: IdentifierType) -> str:
    """Resolve a chemical name, CAS number, or SMILES string to a canonical SMILES.

    Name and CAS lookups go through PubChem's PUG REST API (CAS numbers are
    queried via PubChem's "name" endpoint, which accepts registry numbers).
    Raises ResolutionError if the identifier cannot be resolved or does not
    parse as a valid structure.
    """
    identifier = identifier.strip()
    if not identifier:
        raise ResolutionError("Empty identifier.")

    if identifier_type == "smiles":
        mol = Chem.MolFromSmiles(identifier)
        if mol is None:
            raise ResolutionError(f"'{identifier}' is not a valid SMILES string.")
        return Chem.MolToSmiles(mol)

    # name or cas -> PubChem PUG REST lookup
    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(identifier)}/property/CanonicalSMILES/TXT"
    try:
        resp = requests.get(url, timeout=15)
    except requests.RequestException as exc:
        raise ResolutionError(f"Network error contacting PubChem for '{identifier}': {exc}") from exc

    if resp.status_code == 404:
        raise ResolutionError(f"PubChem has no record for '{identifier}'.")
    if not resp.ok:
        raise ResolutionError(f"PubChem lookup failed for '{identifier}' (HTTP {resp.status_code}).")

    smiles = resp.text.strip().splitlines()[0].strip()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ResolutionError(f"PubChem returned an unparseable structure for '{identifier}'.")
    return Chem.MolToSmiles(mol)


def _is_fully_fluorinated_no_other_halogens(atom) -> tuple[int, int]:
    """For a carbon atom, count attached F atoms and check no H/Cl/Br/I are attached.

    Returns (fluorine_count, other_halogen_or_H_count). The atom only
    qualifies for CF3/CF2 purposes if other_halogen_or_H_count == 0.
    """
    f_count = 0
    bad_count = 0  # H, Cl, Br, I attached directly to this carbon
    bad_count += atom.GetTotalNumHs()  # implicit + explicit H
    for nbr in atom.GetNeighbors():
        sym = nbr.GetSymbol()
        if sym == "F":
            f_count += 1
        elif sym in ("Cl", "Br", "I"):
            bad_count += 1
    return f_count, bad_count


def classify_pfas(smiles: str) -> tuple[bool, list[MatchedGroup]]:
    """Apply the OECD structural rule to a SMILES string.

    Returns (is_pfas, matched_groups). Only saturated (sp3) carbon atoms
    are considered, matching the OECD definition's "methyl"/"methylene"
    language (aromatic and sp2 C-F bonds are not counted).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ResolutionError(f"'{smiles}' is not a valid SMILES string.")

    matches: list[MatchedGroup] = []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        if atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        if atom.GetIsAromatic():
            continue

        f_count, bad_count = _is_fully_fluorinated_no_other_halogens(atom)
        if bad_count > 0:
            continue  # has an H, Cl, Br, or I attached -> not "fully fluorinated"

        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != "F"]
        n_heavy = len(heavy_neighbors)

        f_atom_idxs = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "F"]

        # Methyl carbon fully fluorinated: -CF3 attached to exactly one other heavy atom
        if f_count == 3 and n_heavy == 1:
            matches.append(MatchedGroup("CF3", atom.GetIdx(), f_atom_idxs))
        # Methylene carbon fully fluorinated: -CF2- attached to exactly two other heavy atoms
        elif f_count == 2 and n_heavy == 2:
            matches.append(MatchedGroup("CF2", atom.GetIdx(), f_atom_idxs))
        # CF4 (carbon tetrafluoride): 4 F, 0 heavy neighbours -- degenerate case,
        # still a "fully fluorinated methyl-like" carbon; count it as CF3-type.
        elif f_count == 4 and n_heavy == 0:
            matches.append(MatchedGroup("CF3", atom.GetIdx(), f_atom_idxs))

    return (len(matches) > 0, matches)


def screen(identifier: str, identifier_type: IdentifierType) -> PfasResult:
    """End-to-end: resolve an identifier to a structure and run the PFAS screen."""
    try:
        smiles = resolve_to_smiles(identifier, identifier_type)
    except ResolutionError as exc:
        return PfasResult(
            input_identifier=identifier,
            identifier_type=identifier_type,
            resolved_smiles=None,
            is_pfas=False,
            error=str(exc),
        )

    try:
        is_pfas, matches = classify_pfas(smiles)
    except ResolutionError as exc:
        return PfasResult(
            input_identifier=identifier,
            identifier_type=identifier_type,
            resolved_smiles=smiles,
            is_pfas=False,
            error=str(exc),
        )

    n_cf3 = sum(1 for m in matches if m.group_type == "CF3")
    n_cf2 = sum(1 for m in matches if m.group_type == "CF2")

    highlight_atoms = [m.carbon_atom_idx for m in matches] + [
        idx for m in matches for idx in m.fluorine_atom_idxs
    ]

    if is_pfas:
        parts = []
        if n_cf3:
            parts.append(f"{n_cf3} fully-fluorinated -CF3 group{'s' if n_cf3 != 1 else ''}")
        if n_cf2:
            parts.append(f"{n_cf2} fully-fluorinated -CF2- group{'s' if n_cf2 != 1 else ''}")
        note = "Structural PFAS match: " + " and ".join(parts) + " (OECD 2021 definition)."
    else:
        note = "No fully-fluorinated -CF3 or -CF2- carbon found; not a structural PFAS match under the OECD definition."

    return PfasResult(
        input_identifier=identifier,
        identifier_type=identifier_type,
        resolved_smiles=smiles,
        is_pfas=is_pfas,
        matched_groups=matches,
        n_cf3_groups=n_cf3,
        n_cf2_groups=n_cf2,
        highlight_atoms=highlight_atoms,
        note=note,
    )


def screen_batch(
    identifiers: list[tuple[str, IdentifierType]],
    pubchem_delay_s: float = 0.2,
) -> list[PfasResult]:
    """Screen a batch of (identifier, identifier_type) pairs.

    A small delay is inserted between PubChem lookups to stay well under
    PubChem's public rate limits (max ~5 requests/second).
    """
    results = []
    for identifier, id_type in identifiers:
        results.append(screen(identifier, id_type))
        if id_type != "smiles":
            time.sleep(pubchem_delay_s)
    return results
