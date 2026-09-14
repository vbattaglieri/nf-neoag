"""Thin wrapper around the netMHC-4.0 executable.

netMHC is license-restricted and provided by the user (env var ``NETMHC`` or an
explicit path). It is invoked once per peptide length and its output is parsed
for strong binders (rows flagged ``SB``), mirroring the original ``grep -w SB``
plus column-extraction step.
"""
from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Binder:
    """One strong-binder hit parsed from netMHC output."""
    ident: str      # FASTA identifier (variant id), truncated to 10 chars
    allele: str
    peptide: str
    affinity: str   # predicted IC50 (nM)
    rank: str       # %rank


def resolve_netmhc(explicit: str | None = None) -> str:
    """Locate the netMHC executable from an explicit path or the NETMHC env var."""
    candidate = explicit or os.environ.get("NETMHC")
    if not candidate:
        raise RuntimeError(
            "netMHC not found: pass --netmhc or set the NETMHC environment variable."
        )
    p = Path(candidate)
    if p.is_dir():
        p = p / "netMHC"
    if not p.exists():
        raise RuntimeError(f"netMHC executable not found at {p}")
    return str(p)


def run_netmhc(netmhc: str, fasta: str, alleles: str, length: int) -> str:
    """Run netMHC for one peptide length; return raw stdout text.

    Uses ``tcsh`` because the netMHC 4.0 launcher is a csh script.
    """
    cmd = ["tcsh", netmhc, "-a", alleles, "-s", "-l", str(length), "-f", fasta]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return proc.stdout


def parse_strong_binders(netmhc_stdout: str) -> list[Binder]:
    """Extract strong-binder (``SB``) rows.

    Column layout follows netMHC 4.0 ``-s`` output as used by the original
    pipeline: whitespace-split fields [1]=allele, [2]=peptide, [10]=identifier,
    [12]=affinity, [13]=%rank; kept only when the row is flagged ``SB``.
    """
    binders: list[Binder] = []
    for line in netmhc_stdout.splitlines():
        fields = line.split()
        if "SB" not in fields:
            continue
        if len(fields) <= 13:
            continue
        binders.append(
            Binder(
                ident=fields[10][:10],
                allele=fields[1],
                peptide=fields[2],
                affinity=fields[12],
                rank=fields[13],
            )
        )
    return binders
