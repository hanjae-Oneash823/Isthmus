"""Vendored copy of junior_surveyor.j2_protein_diff._align_proteins.

Isthmus is a separate repo from the Viscacha_pipeline that produces
junior_surveyor's outputs, so this function (and its two small private
helpers) is copied here verbatim rather than imported across repos.
master_surveyor/align_utils.py made the same call for the same reason.
"""

from __future__ import annotations

import re

import edlib

_CIGAR_RE = re.compile(r"(\d+)([=XID])")

_N_TERM_TRUNC = 5
_N_EXT_ABS    = 10


def _parse_cigar(cigar: str) -> list[tuple[int, str]]:
    return [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]


def _align_proteins(canonical: str, alt: str) -> dict:
    """Align alt to canonical with edlib NW and classify the protein change.

    CIGAR conventions (edlib, query=canonical):
        = : match            — both advance
        X : mismatch         — both advance; canonical position recorded
        I : insertion in query (canonical has it, alt doesn't)
              → deletion from alt's perspective; canonical position recorded
        D : deletion in query (alt has it, canonical doesn't)
              → insertion in alt; anchor = current canonical pos (clamped ≥ 1)
    """
    empty = {
        "protein_change_type": "no_sequence",
        "protein_length_diff": 0,
        "pct_identity":        0.0,
        "changed_aa_start":    0,
        "changed_aa_end":      0,
        "mismatch_aa_count":   0,
        "indel_aa_count":      0,
        "changed_aa_fraction": 0.0,
        "premature_stop":      False,
    }
    if not canonical or not alt:
        return empty

    premature_stop = "*" in alt[:-1]

    if canonical == alt:
        return {**empty, "protein_change_type": "identical",
                "pct_identity": 1.0, "premature_stop": False}

    length_diff = len(alt) - len(canonical)
    canon_len   = len(canonical)

    result = edlib.align(canonical, alt, mode="NW", task="path")
    ops = _parse_cigar(result["cigar"])

    can_pos        = 0
    match_count    = 0
    mismatch_count = 0
    indel_count    = 0
    changed_positions: list[int] = []

    for count, op in ops:
        if op == "=":
            can_pos     += count
            match_count += count
        elif op == "X":
            for _ in range(count):
                can_pos += 1
                mismatch_count += 1
                changed_positions.append(can_pos)
        elif op == "I":
            # Canonical residues absent from alt (deletion in alt)
            for _ in range(count):
                can_pos += 1
                indel_count += 1
                changed_positions.append(can_pos)
        elif op == "D":
            # Alt residues absent from canonical (insertion in alt)
            for _ in range(count):
                indel_count += 1
                # Anchor to current canonical position (clamped to ≥ 1 so
                # an N-terminal insertion before the first residue is not
                # recorded as position 0, which would spuriously trigger
                # the N_truncation threshold check of ≤ 5).
                changed_positions.append(max(can_pos, 1))

    unique_positions    = sorted(set(changed_positions))
    changed_aa_start    = unique_positions[0]  if unique_positions else 0
    changed_aa_end      = unique_positions[-1] if unique_positions else 0
    changed_aa_fraction = round(len(unique_positions) / canon_len, 4) if canon_len else 0.0

    aligned_pairs = match_count + mismatch_count   # residue–residue pairs only
    pct_identity  = round(match_count / aligned_pairs, 4) if aligned_pairs else 0.0

    # ── Relative threshold for N/C extension boundary (5% of protein, min 10) ──
    n_ext_thresh = max(_N_EXT_ABS, int(canon_len * 0.05))
    c_ext_thresh = n_ext_thresh

    # ── Classifier ────────────────────────────────────────────────────────────
    if premature_stop:
        change_type = "frameshift_stop"

    elif length_diff < 0:
        if not unique_positions:
            change_type = "C_truncation"
        elif changed_aa_start <= _N_TERM_TRUNC:
            change_type = "N_truncation"
        elif changed_aa_end >= canon_len - _N_TERM_TRUNC:
            change_type = "C_truncation"
        else:
            change_type = "internal_indel"

    elif length_diff > 0:
        if not unique_positions:
            # No divergence in shared region → pure C-terminal addition
            change_type = "C_extension"
        elif changed_aa_start <= n_ext_thresh:
            change_type = "N_extension"
        elif changed_aa_end >= canon_len - c_ext_thresh:
            change_type = "C_extension"
        else:
            change_type = "internal_insertion"

    else:
        # length_diff == 0: with edlib unit-cost, same-length swaps always
        # produce X ops (not I/D pairs), so indel_count is 0 here in practice.
        change_type = "substitution" if mismatch_count > 0 else "identical"

    return {
        "protein_change_type": change_type,
        "protein_length_diff": length_diff,
        "pct_identity":        pct_identity,
        "changed_aa_start":    changed_aa_start,
        "changed_aa_end":      changed_aa_end,
        "mismatch_aa_count":   mismatch_count,
        "indel_aa_count":      indel_count,
        "changed_aa_fraction": changed_aa_fraction,
        "premature_stop":      premature_stop,
    }
