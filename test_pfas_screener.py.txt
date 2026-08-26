"""
Unit tests for the OECD-definition PFAS structural screen.

Run with:  python3 -m pytest test_pfas_screener.py -v
(or just:  python3 test_pfas_screener.py   -- runs as a plain script too)
"""

from pfas_screener import classify_pfas

CASES = [
    # (name, SMILES, expected_is_pfas, reason)
    ("Trifluoroacetic acid (TFA)", "OC(=O)C(F)(F)F", True,
     "classic edge case -- has one -CF3 group, IS considered PFAS under OECD def"),
    ("Perfluorooctanoic acid (PFOA)", "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", True,
     "long perfluorinated chain: CF3 terminus + many CF2 backbone carbons"),
    ("Perfluorooctanesulfonic acid (PFOS)",
     "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", True,
     "long perfluorinated sulfonic acid chain"),
    ("6:2 Fluorotelomer sulfonic acid",
     "OS(=O)(=O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", True,
     "polyfluorinated (not fully-fluorinated throughout) but still has CF2/CF3 segments"),
    ("Perfluorodecalin", "FC1(F)C(F)(F)C(F)(F)C(F)(F)C2(F)C1(F)C(F)(F)C(F)(F)C(F)(F)C2(F)F", True,
     "fully fluorinated bicyclic ring system, all-CF2 ring carbons"),
    ("Benzene", "c1ccccc1", False, "no fluorine at all"),
    ("Fluorobenzene", "Fc1ccccc1", False, "aromatic C-F, not an sp3 methyl/methylene carbon"),
    ("1,1-Difluoroethane", "CC(F)F", False, "CHF2 group still carries an H -- not fully fluorinated"),
    ("Dichlorodifluoromethane (CFC-12)", "FC(F)(Cl)Cl", False,
     "carbon has F AND Cl attached -- excluded by the 'without any Cl/Br/I' clause"),
    ("Vinyl fluoride", "C=CF", False, "sp2 (alkene) carbon, not sp3 methyl/methylene"),
    ("Ethanol", "CCO", False, "no fluorine"),
    ("Perfluoromethane (CF4)", "FC(F)(F)F", True, "degenerate fully-fluorinated methane, counted as CF3-type"),
    ("Difluoromethane", "FCF", False, "CH2F2 -- carbon still carries 2 H atoms, not fully fluorinated"),
]


def run():
    failures = []
    for name, smiles, expected, reason in CASES:
        is_pfas, matches = classify_pfas(smiles)
        status = "PASS" if is_pfas == expected else "FAIL"
        if status == "FAIL":
            failures.append(name)
        print(f"[{status}] {name:45s} expected={expected!s:5s} got={is_pfas!s:5s}  ({reason})")
        if matches:
            for m in matches:
                print(f"        matched: {m.group_type} at atom {m.carbon_atom_idx}")

    print()
    if failures:
        print(f"{len(failures)} FAILURE(S): {failures}")
        raise SystemExit(1)
    else:
        print(f"All {len(CASES)} test cases passed.")


if __name__ == "__main__":
    run()
