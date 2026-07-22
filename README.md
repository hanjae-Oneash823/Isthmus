# Isthmus

Flask web app for browsing `master_surveyor` hits: filter/select transcripts
and drugs, inspect predicted structures, and export selections to BIOVIA
Discovery Studio. Extracted out of the Viscacha_pipeline repo into its own
history on 2026-07-22.

## Layout

- `dossier_server/` — the Flask app itself (this is what's branded "Isthmus"
  in the UI). Entry point: `python -m dossier_server.app`.
- `master_surveyor/` — the hit-selection / structure-prediction / docking-export
  engine that `dossier_server` calls into directly (not just reads output
  files from) via `m0_select`, `m1_ligands`, `m2_structures`, `m2b_structure_qc`,
  `m3_export`.
- `dossier/` — standalone per-candidate report generator (`generate_dossier.py`,
  `generate_index.py`, `generate_all.py`). `dossier_server` reuses its
  `data.py` / `manifest.py` modules for the live `/api/hits` endpoint, but the
  CLI report generators are otherwise independent of the Flask app.

All three packages were originally three sibling directories inside a larger
pipeline repo's `02_SURVEYOR/` folder; that nesting is preserved here as three
sibling top-level packages, so their existing `sys.path` setup and internal
imports (`from master_surveyor import ...`, `from dossier_server import ...`)
work unchanged.

`dossier/sequence_diff.py` used to import a private function
(`_align_proteins`) from a fourth package, `junior_surveyor`, that stayed
behind in the pipeline repo. That function is vendored verbatim into
`dossier/_protein_align.py` here so this repo has no cross-repo Python import.

## Data dependency

None of these three packages generate their own primary data — they all read
from a pipeline's `outputs/` tree (hit tables, pseudobulk CSVs, structure
caches, etc.) that lives in the Viscacha_pipeline repo (or wherever that
pipeline's outputs get written).

Every `config.py` in this repo resolves data paths as:

```python
REPO_ROOT = Path(os.environ.get("ISTHMUS_DATA_ROOT", <this repo's root>))
```

Set `ISTHMUS_DATA_ROOT` to the pipeline checkout that has the real
`outputs/` directory, e.g.:

```bash
export ISTHMUS_DATA_ROOT=/path/to/Viscacha_pipeline
```

Leaving it unset defaults to this repo's own root (i.e. expects an
`outputs/` directory here), which is only useful for a from-scratch local
setup.

Some machine-specific absolute paths (ColabFold binary, ESMFold conda env,
headless Chrome binary) are still hardcoded in `master_surveyor/config.py`
and `dossier/config.py` — update those for whatever host actually runs
structure prediction / dossier PDF rendering.

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# optional: pip install -r requirements-optional.txt for offline
# dossier-generation scripts / local ESMFold

export ISTHMUS_DATA_ROOT=/path/to/Viscacha_pipeline
python -m dossier_server.app   # serves on 0.0.0.0:5057
```
