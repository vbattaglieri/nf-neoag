"""Neoantigen peptide generation from Indels (pipeline step 11).

Faithful Python 3 port of the original ``indel_frameshift.py`` (flag=0 path)
plus the ``oneline_sequence.sh`` YAR reader. The frameshift / in-frame
translation logic — including reverse-strand handling and complex indels — is
preserved on purpose; only the language and I/O plumbing were modernized.

Two deliberate deviations from the legacy code, both flagged in the project docs:

1. The YAR reference is parsed with its refFlat columns anchored from the RIGHT,
   so an optional leading ``>id`` column (present in the current reference files)
   is tolerated. The legacy script assumed fixed left indices.
2. For each peptide length k we emit the peptide built specifically for that
   length (WT tail of k-1 residues + novel residues). The legacy bash used
   ``$1=kmer`` (assignment, always-true) instead of ``$1==kmer``, which pooled
   peptides across lengths; we use the length-matched peptide, which is the
   scientifically intended behaviour.
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass

from .genetic_code import translate as create_prot
from .ids import variant_id
from .genetic_code import reverse_complement as reverse

_NONCODING = re.compile(r"[acgtn]")

# 0-based columns in the indel (cosmic.Indel.altX) table
_COL_COORD = 1
_COL_GENE = 2
_COL_TYPE = 3
_COL_EFFECT = 5
_COL_FREQ = 6
_COL_SEQUENCE = 12


@dataclass
class Exon:
    ex_start: int      # genomic exon start (0-based, refFlat)
    ex_end: int        # genomic exon end
    accession: str
    tx_len: int        # transcript length (for longest-isoform tie-break)
    seq_index: int     # index into the sequence list
    mrna_start: int    # 0-based start of this exon within the spliced mRNA
    mrna_end: int
    strand: str


# ---------------------------------------------------------------------------
# YAR reference (spliced-mRNA FASTA with refFlat header) -> transcript model
# ---------------------------------------------------------------------------
def load_transcripts(yar_path: str) -> tuple[dict[str, list[Exon]], list[str]]:
    """Parse the YAR FASTA. Returns (exons_by_chrom, sequences)."""
    transcripts: dict[str, list[Exon]] = {}
    sequences: list[str] = []

    header: list[str] | None = None
    seq_parts: list[str] = []

    def flush(hdr: list[str], seq: str) -> None:
        # refFlat fields anchored from the right (tolerates a leading >id column)
        accession = hdr[-10]
        chrom = hdr[-9]
        strand = hdr[-8]
        tx_start = int(hdr[-7])
        tx_end = int(hdr[-6])
        exon_starts = [int(x) for x in hdr[-2].strip(",").split(",") if x]
        exon_ends = [int(x) for x in hdr[-1].strip(",").split(",") if x]
        seq_index = len(sequences)
        sequences.append(seq)
        tx_len = tx_end - tx_start
        nx = -1
        for es, ee in zip(exon_starts, exon_ends):
            m_start = nx + 1
            m_end = m_start + (ee - es) - 1
            nx = m_end
            transcripts.setdefault(chrom, []).append(
                Exon(es, ee, accession, tx_len, seq_index, m_start, m_end, strand)
            )

    with open(yar_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    flush(header, "".join(seq_parts))
                header = line.split("\t")
                seq_parts = []
            elif line:
                seq_parts.append(line)
        if header is not None:
            flush(header, "".join(seq_parts))

    # order exons by genomic start, then longest transcript first
    for chrom in transcripts:
        transcripts[chrom].sort(key=lambda e: (e.ex_start, -e.tx_len))
    return transcripts, sequences


def _exons_at(chrom: str, pos: int, transcripts: dict[str, list[Exon]]) -> list[Exon]:
    """Exons (across isoforms) overlapping a genomic position."""
    exons = transcripts.get(chrom, [])
    matches: list[Exon] = []
    for n, e in enumerate(exons):
        if e.ex_start <= pos <= e.ex_end:
            matches.append(e)
            for j in exons[n + 1:]:
                if e.ex_end < j.ex_start:
                    break
                if j.ex_start <= pos <= j.ex_end:
                    matches.append(j)
            break
    return matches


# ---------------------------------------------------------------------------
# Frameshift / in-frame translation (verbatim logic from the original)
# ---------------------------------------------------------------------------
def _clean_len(s: str) -> int:
    return len(_NONCODING.sub("", s))


def calculate_fs(coord: list, exon: Exon, seq: str) -> tuple[str, str, str]:
    """Return (indel_bases, novel_protein, wt_5prime_protein).

    ``coord`` = [chrom, start(1-based), end, type, inserted_or_complex_seq].
    """
    indel_start = coord[1] - exon.ex_start + exon.mrna_start
    indel_end = coord[2] - exon.ex_start + exon.mrna_start - 1

    if exon.strand == "+":
        if coord[3] == "Deletion":
            indel = seq[indel_start:indel_end]
            only_fs = seq[indel_end:]
            if coord[4]:                      # complex deletion
                only_fs = coord[4] + only_fs
        elif coord[3] == "Insertion":
            indel = coord[4]
            only_fs = coord[4] + seq[indel_start:]
        else:
            raise ValueError(f"bad indel type at {coord[0]}:{coord[1]}-{coord[2]}")

        rest = _clean_len(seq[:indel_start]) % 3
        if rest:
            only_fs = seq[:indel_start][-rest:] + only_fs
            five_prime = seq[:indel_start][:-rest]
        else:
            five_prime = seq[:indel_start]

    elif exon.strand == "-":
        seq = reverse(seq)
        if coord[3] == "Deletion":
            indel = seq[indel_start:indel_end]
            only_fs = seq[:indel_start]
            if coord[4]:
                only_fs += coord[4]
        elif coord[3] == "Insertion":
            indel = coord[4]
            only_fs = seq[:indel_start] + coord[4]
        else:
            raise ValueError(f"bad indel type at {coord[0]}:{coord[1]}-{coord[2]}")

        rest = _clean_len(seq[indel_end:]) % 3
        if rest:
            only_fs = only_fs + seq[indel_end:][:rest]
            five_prime = seq[indel_end:][rest:]
        else:
            five_prime = seq[indel_end:]

        only_fs = reverse(only_fs)
        five_prime = reverse(five_prime)
    else:
        raise ValueError(f"bad strand {exon.strand!r}")

    frameshift = create_prot(_NONCODING.sub("", only_fs))
    wt = create_prot(_NONCODING.sub("", five_prime))
    if "*" in frameshift:
        frameshift = frameshift[:frameshift.index("*") + 1]
    return indel, frameshift, wt


# ---------------------------------------------------------------------------
# Per-length peptide records
# ---------------------------------------------------------------------------
@dataclass
class IndelPeptide:
    kmer: int
    peptide: str
    ident: str
    gene: str
    accession: str
    coord: str
    label: str          # e.g. "Ins_I-18_fs12*"  (consensus "TYPE" column)
    frameshift: str     # novel protein (consensus "FRAMESHIFT" column)
    freq: str

    @property
    def info(self) -> list[str]:
        return [self.gene, self.accession, self.coord,
                self.label, self.frameshift, self.freq]


def generate(indel_path: str, transcripts: dict[str, list[Exon]],
             sequences: list[str], min_kmer: int, max_kmer: int) -> list[IndelPeptide]:
    """Yield per-(variant, k) peptides for netMHC input."""
    out: list[IndelPeptide] = []
    with open(indel_path) as fh:
        next(fh, None)  # header
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if len(row) <= _COL_SEQUENCE:
                continue
            chrom, _, span = row[_COL_COORD].partition(":")
            c0, _, c1 = span.partition("-")
            type_indel = row[_COL_TYPE]
            sequence = row[_COL_SEQUENCE].strip('"')
            gene = row[_COL_GENE]
            effect = row[_COL_EFFECT]
            freq = row[_COL_FREQ]
            coord = [chrom, int(c0), int(c1), type_indel, sequence]

            if chrom not in transcripts:
                print(f"WARNING: {chrom} not in annotated transcripts", file=sys.stderr)
                continue
            if "N" in sequence:
                print(f"WARNING: indel sequence contains N at {chrom}:{c0}-{c1} in {gene}",
                      file=sys.stderr)
                continue

            candidates = _exons_at(chrom, coord[1], transcripts)
            if not candidates:
                print(f"WARNING: transcript not found at {chrom}:{c0}-{c1} in {gene}",
                      file=sys.stderr)
                continue

            # longest isoform whose exon fully contains the deletion end
            candidates = sorted(candidates,
                                key=lambda e: -_clean_len(sequences[e.seq_index]))
            exon = None
            seq = ""
            for e in candidates:
                if coord[2] > e.ex_end:
                    print(f"WARNING: deletion spans two exons, skipping "
                          f"{chrom}:{c0}-{c1} in {gene} {e.accession}", file=sys.stderr)
                    continue
                exon = e
                seq = sequences[e.seq_index]
                break
            if not seq or seq.islower():
                print(f"WARNING: skipped deletion {chrom}:{c0}-{c1}", file=sys.stderr)
                continue

            indel, frameshift, wt = calculate_fs(coord, exon, seq)

            complex_indel = False
            if type_indel == "Deletion" and sequence:
                complex_indel = True
                if not abs(len(sequence) - len(indel)) % 3:
                    effect = "in-frame"
            coord[4] = indel

            if frameshift[:-1].count("*") or wt.count("*"):
                print(f"WARNING: wrong YAR annotation, skipped {chrom}:{c0}-{c1} "
                      f"in {gene} {exon.accession}", file=sys.stderr)
                continue
            if indel.islower():
                print(f"WARNING: indel not in exonic region {chrom}:{c0}-{c1} "
                      f"in {gene} {exon.accession}", file=sys.stderr)
                continue
            if not frameshift:
                print(f"WARNING: no frameshift created at {chrom}:{c0}-{c1} in {gene}",
                      file=sys.stderr)
                continue

            ident = variant_id(exon.accession, chrom, c0 + c1)
            typ = type_indel[:3]

            for k in range(min_kmer, max_kmer + 1):
                kwt = wt[-(k - 1):]
                if effect == "frameshift":
                    fs = frameshift.replace("*", "")
                    aachange = "fs" + str(len(fs)) + ("*" if "*" in frameshift else "")
                    out.append(IndelPeptide(
                        kmer=k, peptide=kwt + fs, ident=ident, gene=gene,
                        accession=exon.accession, coord=row[_COL_COORD],
                        label=f"{typ}_{indel}_{aachange}",
                        frameshift=frameshift, freq=freq,
                    ))
                elif effect == "in-frame":
                    lifs = k + 1
                    if type_indel == "Insertion":
                        lifs += len(indel) // 3
                    elif type_indel == "Deletion" and complex_indel:
                        lifs += len(sequence) // 3 + (1 if len(sequence) % 3 else 0)
                    ifs = frameshift.replace("*", "")[:lifs]
                    out.append(IndelPeptide(
                        kmer=k, peptide=kwt + ifs, ident=ident, gene=gene,
                        accession=exon.accession, coord=row[_COL_COORD],
                        label=f"{typ}_{indel}_{effect}",
                        frameshift=frameshift[:(lifs - k)], freq=freq,
                    ))
                else:
                    print(f"WARNING: unexpected effect {effect!r}, skipped "
                          f"{chrom}:{c0}-{c1} in {gene}", file=sys.stderr)
    return out
