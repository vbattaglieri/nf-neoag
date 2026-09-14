"""Neoantigen peptide generation from SNVs (pipeline step 10).

Clean re-implementation of the original ``Neoantigens.v4.sh`` filter + awk
parsing + ``trad_yar.py``. For each qualifying non-synonymous SNV it builds the
mutated coding sequence and translates a window of (2*k - 1) residues centred on
the mutated codon, so that every k-mer overlapping the mutation is represented.

Input files
-----------
* IDEA table  : ``*.snp.str-bias.Idea.OK`` (the step-5 output).
* cDNA table  : ``refFlat_mRNA.<build>.noalt.parsed.tab`` =
                ``accession <tab> coord <tab> sequence`` (lowercase = non-coding).
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass

from .genetic_code import translate_codon
from .ids import variant_id

# 0-based column indices in the IDEA .OK table
_COL_AACHANGE = 3     # e.g. p.S553A
_COL_EFFECT = 4       # e.g. nonsynonymous
_COL_PCT_MUT = 5      # % mutant reads
_COL_NCHANGE = 6      # e.g. c.T1657G
_COL_ACCESSION = 7    # e.g. NM_005807
_COL_COORD = 12       # e.g. chr1:186307376
_COL_SNP = 35         # dbSNP id or '--'

_DIGITS = re.compile(r"\d+")
_NONCODING = re.compile(r"[acgt]")


@dataclass
class SnvVariant:
    accession: str
    gene: str
    coord: str
    aa_pos: int
    aa_ref: str
    aa_alt: str
    nc_pos: int
    nc_ref: str
    nc_alt: str
    freq: str

    @property
    def ident(self) -> str:
        chrom, _, pos = self.coord.partition(":")
        return variant_id(self.accession, chrom, pos)

    @property
    def aachange(self) -> str:
        return f"{self.aa_ref}{self.aa_pos}{self.aa_alt}"

    @property
    def nchange(self) -> str:
        return f"{self.nc_ref}{self.nc_pos}{self.nc_alt}"

    @property
    def info(self) -> list[str]:
        """Per-variant columns for the consensus output."""
        return [self.gene, self.accession, self.coord,
                self.nchange, self.aachange, self.freq]


def _split_change(text: str) -> tuple[str, int, str]:
    """'p.S553A' / 'c.T1657G' -> (first_char, number, last_char)."""
    body = text.split(".", 1)[-1]
    m = _DIGITS.search(body)
    if not m:
        raise ValueError(f"no position in change string: {text!r}")
    return body[0], int(m.group()), body[-1]


def parse_idea(path: str) -> list[SnvVariant]:
    """Read and filter the IDEA table.

    Keeps rows that are non-synonymous, supported by >1% mutant reads and not a
    known dbSNP polymorphism (same criteria as the original awk filter).
    """
    variants: list[SnvVariant] = []
    with open(path) as fh:
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if len(row) <= _COL_SNP:
                continue
            snp = row[_COL_SNP]
            if re.search(r"rs[0-9]", snp):
                continue
            try:
                if float(row[_COL_PCT_MUT]) <= 1:
                    continue
            except ValueError:
                continue
            if row[_COL_EFFECT] == "synonymous":
                continue
            accession = row[_COL_ACCESSION]
            if not accession or accession.startswith("?"):
                continue
            try:
                aa_ref, aa_pos, aa_alt = _split_change(row[_COL_AACHANGE])
                nc_ref, nc_pos, nc_alt = _split_change(row[_COL_NCHANGE])
            except (ValueError, IndexError):
                continue
            variants.append(SnvVariant(
                accession=accession, gene=row[1], coord=row[_COL_COORD],
                aa_pos=aa_pos, aa_ref=aa_ref, aa_alt=aa_alt,
                nc_pos=nc_pos, nc_ref=nc_ref, nc_alt=nc_alt,
                freq=row[_COL_PCT_MUT],
            ))
    return variants


def load_cdna(path: str) -> dict[str, str]:
    """accession -> coding sequence (non-coding lowercase bases removed)."""
    cdna: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            cdna[parts[0]] = _NONCODING.sub("", parts[2])
    return cdna


def _codon_bounds(pos0: int, rest: int) -> tuple[int, int]:
    """Slice bounds of the codon containing 0-based position ``pos0``."""
    if rest == 0:
        return pos0 - 2, pos0 + 1
    if rest == 1:
        return pos0, pos0 + 3
    return pos0 - 1, pos0 + 2  # rest == 2


def _upstream(cdna: str, start: int, kmer: int) -> list[str]:
    k = start // 3 if start <= kmer * 3 else kmer
    return [translate_codon(cdna[start - 3 * i: start - 3 * i + 3])
            for i in range(1, k)]


def _downstream(cdna: str, end: int, kmer: int) -> list[str]:
    k = (len(cdna) - end) // 3 if end + kmer * 3 >= len(cdna) else kmer
    return [translate_codon(cdna[end - 3 + 3 * i: end - 3 + 3 * i + 3])
            for i in range(1, k)]


def mutant_peptide(cdna: str, v: SnvVariant, kmer: int) -> str | None:
    """Translate the (2k-1)-residue window around the mutation, or None."""
    cdna_new = cdna[:v.nc_pos - 1] + v.nc_alt + cdna[v.nc_pos:]
    pos0 = v.nc_pos - 1
    rest = v.nc_pos % 3
    start, end = _codon_bounds(pos0, rest)
    if end > len(cdna_new) or start < 0:
        print(f"WARNING: codon out of transcript bounds: {v.nchange} "
              f"{v.accession} {v.gene}", file=sys.stderr)
        return None
    before = _upstream(cdna_new, start, kmer)
    after = _downstream(cdna_new, end, kmer)
    return "".join(before[::-1] + [v.aa_alt] + after)
