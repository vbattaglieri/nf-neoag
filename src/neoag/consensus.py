"""Collapse per-peptide netMHC strong-binder hits into a per-variant consensus.

Replaces the original ``Pep_consensus.py`` (SNV) and ``Pep_consensus_indel.py``
(Indel, with Levenshtein "family" grouping). Behaviour is preserved:

* the "best" hit per variant is the one with the lowest %rank;
* "others" columns list the set of alleles / peptides seen for that variant;
* for indels, N.FAMILY counts groups of peptides within an edit distance.

Output rows are written tab-separated. Row widths intentionally match the
originals (the downstream merge step relies on them).
"""
from __future__ import annotations

import sys
from collections import OrderedDict
from typing import Iterable, Sequence

SNV_HEADER = [
    "GENE", "ACCESSION", "COORDINATES", "NCHANGE", "AACHANGE", "ALLELE FREQ",
    "MHC HAPLOTYPE", "PEPTIDE (best)", "IC50 (best)", "RANK (best)",
    "MHC HAPLOTYPE (others)", "PEPTIDE (others)",
]

INDEL_HEADER = [
    "GENE", "ACCESSION", "COORDINATES", "TYPE", "FRAMESHIFT", "ALLELE FREQ",
    "MHC HAPLOTYPE", "PEPTIDE (best)", "IC50 (best)", "RANK (best)",
    "MHC HAPLOTYPE (others)", "PEPTIDE (others)", "N.FAMILY",
]

# A record: (ident, info_fields, (allele, peptide, affinity, rank))
Record = tuple


def _group(records: Iterable[Record]) -> "OrderedDict[str, tuple]":
    grouped: "OrderedDict[str, tuple]" = OrderedDict()
    for ident, info, hit in records:
        if ident not in grouped:
            grouped[ident] = (list(info), [list(hit)])
        else:
            grouped[ident][1].append(list(hit))
    return grouped


def _best(hits: Sequence[Sequence[str]]) -> list[str]:
    """Hit with the minimum %rank (4th field)."""
    return min(hits, key=lambda h: float(h[3]))


def _others(hits: Sequence[Sequence[str]]) -> tuple[str, str]:
    alleles = sorted({h[0] for h in hits})
    peptides = sorted({h[1] for h in hits})
    return ",".join(alleles), ",".join(peptides)


def levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a):
        cur = [i + 1]
        for j, cb in enumerate(b):
            cost = 0 if ca == cb else 1
            cur.append(min(cur[j] + 1, prev[j + 1] + 1, prev[j] + cost))
        prev = cur
    return prev[-1]


def count_families(peptides: Sequence[str], dist: int) -> int:
    """Number of connected components where edit distance < ``dist``."""
    peps = list(dict.fromkeys(peptides))  # unique, order-preserving
    seen: set[int] = set()
    families = 0
    for i in range(len(peps)):
        if i in seen:
            continue
        families += 1
        stack = [i]
        seen.add(i)
        while stack:
            k = stack.pop()
            for j in range(len(peps)):
                if j not in seen and levenshtein(peps[k], peps[j]) < dist:
                    seen.add(j)
                    stack.append(j)
    return families


def write_snv_consensus(records: Iterable[Record], out=sys.stdout) -> None:
    print("\t".join(SNV_HEADER), file=out)
    for info, hits in _group(records).values():
        if len(hits) == 1:
            print("\t".join(info + hits[0]), file=out)
        else:
            best = _best(hits)
            hapl, pep = _others(hits)
            print("\t".join(info + best + [hapl, pep]), file=out)


def write_indel_consensus(records: Iterable[Record], dist: int, out=sys.stdout) -> None:
    print("\t".join(INDEL_HEADER), file=out)
    for info, hits in _group(records).values():
        if len(hits) == 1:
            print("\t".join(info + hits[0] + ["-", "-", "1"]), file=out)
        else:
            best = _best(hits)
            hapl, pep = _others(hits)
            nfam = count_families([h[1] for h in hits], dist)
            print("\t".join(info + best + [hapl, pep, str(nfam)]), file=out)
