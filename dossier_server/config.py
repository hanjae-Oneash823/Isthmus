"""DOSSIER_SERVER — central configuration."""

import os
from pathlib import Path

# Root of the pipeline's outputs/ tree (e.g. a Viscacha_pipeline checkout).
# Defaults to this repo's own root; override to point at wherever the
# pipeline actually wrote its outputs/ directory.
REPO_ROOT = Path(os.environ.get("ISTHMUS_DATA_ROOT", Path(__file__).resolve().parent.parent)).resolve()

HOST = "0.0.0.0"   # server-hosted, reachable off this box -- single-user
                    # internal tool, no auth layer (see plan doc's Context)
PORT = 5057

BASE_DIR   = REPO_ROOT / "outputs/master_surveyor"
CART_JSON  = BASE_DIR / "cart.json"
JOBS_DIR   = BASE_DIR / "jobs"
DOCKING_DIR = BASE_DIR / "docking"
