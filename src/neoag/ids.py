"""Stable, short variant identifiers.

The original code used `hashids` to turn (transcript, chromosome, position) into
a compact token. We keep the same *purpose* — a deterministic ID that is unique
per variant and stable across the separate per-k-mer netMHC runs so results can
be joined back together — but implement it with the standard library (no extra
dependency).

netMHC truncates the FASTA identifier, so IDs are capped at 10 characters.
"""
from __future__ import annotations

import hashlib

ID_LENGTH = 10


def variant_id(transcript: str, chrom: str, position: str | int) -> str:
    """Deterministic 10-char alphanumeric id for a variant."""
    key = f"{transcript}_{chrom}_{position}".encode()
    return hashlib.sha1(key).hexdigest()[:ID_LENGTH]
