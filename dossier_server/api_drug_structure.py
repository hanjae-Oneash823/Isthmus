"""GET /api/drug_structure/<name> -- 2D chemical structure image for a drug
name.

Resolved via PubChem's PUG REST API, which looks compounds up by name
directly -- no ChEMBL/DGIdb molecule ID plumbing needed, so this works
uniformly across all three drug-evidence sources (ChEMBL, Open Targets,
DGIdb) even though only ChEMBL's own pipeline stage ever sees a ChEMBL ID.

Every name is resolved at most once: both hits and misses are cached to
disk keyed by a hash of the name, so repeat page loads/re-renders never
re-hit PubChem for a name already known to work or not resolve.
"""

from __future__ import annotations

import hashlib
import sys
import urllib.parse
from pathlib import Path

import requests
from flask import Blueprint, Response, abort

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dossier_server.config import DRUG_STRUCTURE_CACHE_DIR

bp = Blueprint("api_drug_structure", __name__)

_PUBCHEM_IMG_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/PNG?image_size=150x150"
_TIMEOUT = 10
# PubChem's usage policy asks requesting tools to identify themselves.
_HEADERS = {"User-Agent": "Isthmus-dossier-viewer/1.0 (single-user internal research tool)"}
_MISS_SUFFIX = ".miss"


def _cache_key(name: str) -> str:
    return hashlib.sha1(name.encode("utf-8")).hexdigest()


@bp.route("/api/drug_structure/<path:name>")
def drug_structure(name: str):
    DRUG_STRUCTURE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    key = _cache_key(name)
    img_path = DRUG_STRUCTURE_CACHE_DIR / f"{key}.png"
    miss_path = DRUG_STRUCTURE_CACHE_DIR / f"{key}{_MISS_SUFFIX}"

    if img_path.exists():
        return Response(img_path.read_bytes(), mimetype="image/png")
    if miss_path.exists():
        abort(404)

    url = _PUBCHEM_IMG_URL.format(name=urllib.parse.quote(name, safe=""))
    try:
        resp = requests.get(url, timeout=_TIMEOUT, headers=_HEADERS)
    except requests.RequestException:
        # Network hiccup, not a confirmed "PubChem doesn't know this name"
        # -- don't cache a miss, so a later request can retry cleanly.
        abort(502)

    if resp.status_code == 404:
        # Only a confirmed "no compound by this name" is worth caching
        # forever. Other non-200s (503 ServerBusy under PubChem's rate
        # limit being the common one) are transient -- caching those as a
        # permanent miss would silently hide the image after every retry,
        # even once PubChem is no longer busy.
        miss_path.write_bytes(b"")
        abort(404)

    if resp.status_code != 200 or not resp.content:
        abort(502)

    img_path.write_bytes(resp.content)
    return Response(resp.content, mimetype="image/png")
