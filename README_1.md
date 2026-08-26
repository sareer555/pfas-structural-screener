# PFAS Structural Screener

Screens chemicals (by name, CAS number, or SMILES) for PFAS status using the
**OECD (2021) structural definition** — any molecule containing at least one
fully-fluorinated methyl (`-CF3`) or methylene (`-CF2-`) carbon, with no
H/Cl/Br/I attached to that carbon. Single-chemical lookup or batch CSV
upload, both with a downloadable results report.

Built for the same use case as the GHS Hazard Screening tool: a live demo
that also doubles as the delivery mechanism for a paid Project Catalog
service (PFAS/TSCA 8(a)(7) screening, Minnesota PFAS reporting, or general
REACH/SVHC-style watchlist checks).

## Files

- `pfas_screener.py` — core logic: PubChem name/CAS → SMILES resolution,
  and the OECD structural rule itself (pure RDKit, no network needed once
  you have a SMILES string).
- `app.py` — the Streamlit UI (single lookup, batch CSV, methodology tab).
- `test_pfas_screener.py` — 13 unit tests against known PFAS and
  known-tricky non-PFAS cases (CFCs, aromatic C-F, partially-fluorinated
  carbons). Run with `python3 test_pfas_screener.py` — all should pass.
- `requirements.txt` — pinned dependencies.

## Running locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then open the local URL it prints (usually `http://localhost:8501`).

## Deploying (same flow as your GHS tool)

1. Push this folder to a new GitHub repo (public or private both work).
2. Go to [share.streamlit.io](https://share.streamlit.io), sign in with the
   same GitHub account you used for `ghs-hazard-classification`.
3. "New app" → point it at this repo, branch `main`, main file `app.py`.
4. Deploy. First load will be slow while it installs RDKit; subsequent
   loads are fast (same free-tier sleep-after-inactivity behavior as your
   GHS app, so it may need a "wake up" click after idle periods).

## What this tool is (and isn't)

It's a **structural** screen, not a lookup against any specific regulatory
list. It tells you a molecule's structure meets the scientific PFAS
definition that regulators are increasingly building policy around
(EPA TSCA 8(a)(7), Minnesota Chapter 116, the EU's PFAS restriction
proposal) — not that the exact substance is individually named on a given
list. For a real compliance filing, a flag here is a strong signal to check
the specific list that filing requires, not the filing determination
itself. The app says this explicitly on its "About / methodology" tab and
under every result, so a client can't mistake it for legal/regulatory
sign-off.

## Extending this for a Project Catalog listing

Ideas for a v2, if this gets traction:
- Add curated-list cross-referencing (EPA's TSCA PFAS list, ECHA's SVHC
  Candidate List) as a second, separate flag column — so a report can show
  "structural match" vs. "named on EPA's list" vs. "named on ECHA's list"
  side by side.
- PDF report export (client-ready, like the GHS tool's PDF deliverable) —
  the underlying `PfasResult` objects already have everything a report
  template would need.
- A "flag changes over time" mode for the SVHC/REACH OSOA use case,
  since those lists update periodically and clients would want to know
  when a previously-clear substance gets added.
