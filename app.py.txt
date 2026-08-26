"""
PFAS Structural Screener -- Streamlit app.

Screens chemicals (by name, CAS number, or SMILES) for PFAS status using
the OECD (2021) structural definition, either one at a time or as a batch
(CSV upload) -- built for compliance teams working through TSCA 8(a)(7)
reporting, Minnesota's PFAS reporting law, or general REACH/SVHC-style
watchlist screening.
"""

import io

import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from pfas_screener import IdentifierType, PfasResult, screen, screen_batch

st.set_page_config(page_title="PFAS Structural Screener", page_icon="🧪", layout="wide")


def mol_image_png(smiles: str, highlight_atoms: list[int], size=(360, 300)) -> bytes:
    mol = Chem.MolFromSmiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors={idx: (1.0, 0.55, 0.55) for idx in highlight_atoms},
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def render_result(res: PfasResult):
    if res.error:
        st.error(f"Could not screen '{res.input_identifier}': {res.error}")
        return

    col1, col2 = st.columns([1, 1.4])
    with col1:
        st.markdown("**Structure**")
        png_bytes = mol_image_png(res.resolved_smiles, res.highlight_atoms)
        st.image(io.BytesIO(png_bytes), width="stretch")
        st.caption(f"SMILES: `{res.resolved_smiles}`")

    with col2:
        st.markdown("**PFAS screening result**")
        if res.is_pfas:
            st.error(f"🚩 {res.flag_label}")
        else:
            st.success(f"✅ {res.flag_label}")
        st.write(res.note)
        if res.matched_groups:
            st.markdown("**Matched groups**")
            for m in res.matched_groups:
                st.write(f"- `{m.group_type}` group at carbon atom index {m.carbon_atom_idx}")
        st.caption(
            "Structural screen only, per the OECD (2021) definition — this is not a lookup "
            "against EPA's TSCA PFAS list, ECHA's SVHC list, or any other curated regulatory "
            "list. Treat a flag as 'this molecule's structure meets the scientific PFAS "
            "definition,' not as 'this exact substance is named on a specific regulation.'"
        )


st.title("🧪 PFAS Structural Screener")
st.markdown(
    "Screens chemicals for PFAS status using the **OECD (2021) structural definition**: "
    "any molecule containing at least one fully-fluorinated methyl (`-CF3`) or methylene "
    "(`-CF2-`) carbon, with no H/Cl/Br/I attached to that carbon. Built to help with "
    "**TSCA 8(a)(7)** historical-use reporting, **Minnesota's PFAS reporting law**, and "
    "general formulation/BOM screening."
)

tab_single, tab_batch, tab_about = st.tabs(["Single chemical", "Batch upload (CSV)", "About / methodology"])

with tab_single:
    st.subheader("1. Enter a chemical")
    id_type_label = st.radio(
        "How would you like to identify the chemical?",
        options=["Chemical name", "CAS number", "SMILES string"],
        horizontal=True,
    )
    id_type_map: dict[str, IdentifierType] = {
        "Chemical name": "name",
        "CAS number": "cas",
        "SMILES string": "smiles",
    }
    identifier = st.text_input(
        "Identifier",
        placeholder="e.g. perfluorooctanoic acid, 335-67-1, or OC(=O)C(F)(F)C(F)(F)F",
    )
    st.caption(
        "Try one of these: perfluorooctanoic acid · trifluoroacetic acid · "
        "dichlorodifluoromethane · fluorobenzene · benzene"
    )

    if st.button("Screen chemical", type="primary"):
        if not identifier.strip():
            st.warning("Enter a chemical name, CAS number, or SMILES string first.")
        else:
            with st.spinner("Resolving structure and screening..."):
                result = screen(identifier, id_type_map[id_type_label])
            st.divider()
            render_result(result)

with tab_batch:
    st.subheader("Screen a list of chemicals at once")
    st.write(
        "Upload a CSV with one chemical per row. The file needs a column called "
        "`identifier` and a column called `type` (values: `name`, `cas`, or `smiles`). "
        "Mixed types in the same file are fine."
    )
    st.download_button(
        "Download a CSV template",
        data="identifier,type\nperfluorooctanoic acid,name\n335-67-1,cas\nOC(=O)C(F)(F)F,smiles\nbenzene,name\n",
        file_name="pfas_screener_template.csv",
        mime="text/csv",
    )

    uploaded = st.file_uploader("Upload CSV", type=["csv"])
    if uploaded is not None:
        try:
            df_in = pd.read_csv(uploaded)
        except Exception as exc:  # noqa: BLE001
            st.error(f"Could not read that CSV: {exc}")
            df_in = None

        if df_in is not None:
            missing_cols = {"identifier", "type"} - set(df_in.columns.str.lower())
            df_in.columns = [c.lower() for c in df_in.columns]
            if missing_cols:
                st.error(f"CSV is missing required column(s): {', '.join(missing_cols)}")
            else:
                bad_types = set(df_in["type"].str.lower().unique()) - {"name", "cas", "smiles"}
                if bad_types:
                    st.error(f"Unrecognized 'type' value(s): {', '.join(bad_types)} (use name / cas / smiles)")
                else:
                    st.write(f"Loaded {len(df_in)} chemical(s). Preview:")
                    st.dataframe(df_in.head(10), width="stretch")

                    if st.button("Run batch screen", type="primary"):
                        pairs = list(zip(df_in["identifier"].astype(str), df_in["type"].str.lower()))
                        progress = st.progress(0.0, text="Screening...")
                        results: list[PfasResult] = []
                        for i, (identifier, id_type) in enumerate(pairs):
                            results.append(screen(identifier, id_type))
                            progress.progress((i + 1) / len(pairs), text=f"Screening {i + 1}/{len(pairs)}...")
                        progress.empty()

                        out_rows = []
                        for r in results:
                            out_rows.append(
                                {
                                    "identifier": r.input_identifier,
                                    "type": r.identifier_type,
                                    "resolved_smiles": r.resolved_smiles,
                                    "pfas_flag": r.flag_label,
                                    "n_cf3_groups": r.n_cf3_groups,
                                    "n_cf2_groups": r.n_cf2_groups,
                                    "note": r.error if r.error else r.note,
                                }
                            )
                        df_out = pd.DataFrame(out_rows)

                        n_flagged = (df_out["pfas_flag"] == "PFAS (structural match)").sum()
                        n_errors = (df_out["pfas_flag"] == "ERROR").sum()
                        c1, c2, c3 = st.columns(3)
                        c1.metric("Chemicals screened", len(df_out))
                        c2.metric("Flagged as PFAS", int(n_flagged))
                        c3.metric("Could not resolve", int(n_errors))

                        st.dataframe(df_out, width="stretch")
                        st.download_button(
                            "Download full results (CSV)",
                            data=df_out.to_csv(index=False),
                            file_name="pfas_screening_results.csv",
                            mime="text/csv",
                            type="primary",
                        )

with tab_about:
    st.markdown(
        """
### Methodology

This tool applies the OECD's 2021 structural definition of PFAS
(*"Reconciling Terminology of the Universe of Per- and Polyfluoroalkyl
Substances"*, Series on Risk Management No. 61, ENV/CBC/MONO(2021)25):

> PFASs are defined as fluorinated substances that contain at least one
> fully fluorinated methyl (CF3-) or methylene (-CF2-) carbon atom
> (without any H/Cl/Br/I atom attached to it).

Concretely, for every carbon atom in the molecule, the screen checks:

1. Is the carbon sp3 (saturated) and non-aromatic?
2. Does it carry **zero** hydrogens, chlorines, bromines, or iodines?
3. Does it carry exactly **3 fluorines** (a `-CF3` methyl group) or exactly
   **2 fluorines** with its two remaining bonds to other heavy atoms
   (a `-CF2-` methylene group)?

If any carbon in the molecule meets these criteria, the molecule is
flagged as a structural PFAS match.

**What this deliberately excludes**, matching the OECD text:
- Aromatic C-F bonds (e.g. fluorobenzene) — not a methyl/methylene carbon.
- Partially fluorinated carbons that still carry an H (e.g. `-CHF2`, `-CH2F`).
- Carbons that mix fluorine with chlorine/bromine/iodine (e.g. CFC and HCFC
  refrigerants) — explicitly excluded by the "without any Cl/Br/I" clause.

**What this is not:** a lookup against any specific regulatory list (EPA's
TSCA PFAS reporting list, ECHA's SVHC Candidate List, Minnesota's
Chapter 116 list, etc.). Those lists name specific substances; this tool
tells you whether a structure *meets the scientific definition* those
regulations are increasingly built around. For a compliance filing, treat
a flag here as a strong signal to look the specific substance up against
the applicable list, not as the filing determination itself.

### Data source

Chemical names and CAS numbers are resolved to structures via the
[PubChem PUG REST API](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest).
        """
    )
