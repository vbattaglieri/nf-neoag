"""Merge the SNV and Indel per-variant consensus tables into one.

Replaces the original ``Neoantigens.merge_SNV_FS.sh``. The legacy script assumed
the two consensus files had already been expression-annotated (extra
``expected_counts`` / ``FPKM`` columns). That annotation is a separate step
outside the scope of pipeline points 10 & 11, so this reimplementation operates
directly on the raw consensus outputs produced by :mod:`neoag.consensus`.

Two files are produced:

* ``<out>.neoantigens.consensus.all`` — SNV rows and Indel rows harmonised into
  a single schema, tagged by ``SOURCE``;
* ``<out>.neoantigens.consensus.all.expanded`` — the same, but rows with
  ``N.FAMILY > 1`` are duplicated once per family (coordinate suffixed ``.i``),
  mirroring the legacy expansion.
"""
from __future__ import annotations

MERGED_HEADER = [
    "GENE", "ACCESSION", "COORDINATES", "VARIANT", "ALLELE FREQ",
    "MHC HAPLOTYPE (best)", "PEPTIDE (best)", "IC50 (best)", "RANK (best)",
    "MHC HAPLOTYPE (others)", "PEPTIDE (others)", "N.FAMILY", "SOURCE",
]


def _read_rows(path: str) -> list[list[str]]:
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    return rows[1:] if rows else []   # drop header


def _norm_snv(row: list[str]) -> list[str]:
    """SNV consensus row -> merged schema. 10 cols (single) or 12 (multi)."""
    gene, acc, coord, nchange, aachange, freq = row[:6]
    hapl, pep, ic50, rank = row[6:10]
    others_h = row[10] if len(row) > 10 else "-"
    others_p = row[11] if len(row) > 11 else "-"
    variant = f"{nchange}/{aachange}"
    return [gene, acc, coord, variant, freq, hapl, pep, ic50, rank,
            others_h, others_p, "1", "SNV"]


def _norm_indel(row: list[str]) -> list[str]:
    """Indel consensus row -> merged schema. Always 13 cols."""
    gene, acc, coord, typ, frameshift, freq = row[:6]
    hapl, pep, ic50, rank = row[6:10]
    others_h = row[10] if len(row) > 10 else "-"
    others_p = row[11] if len(row) > 11 else "-"
    nfam = row[12] if len(row) > 12 else "1"
    variant = f"{typ}/{frameshift}"
    return [gene, acc, coord, variant, freq, hapl, pep, ic50, rank,
            others_h, others_p, nfam, "Indel"]


def merge(snv_path: str, indel_path: str, out_basename: str) -> tuple[str, str]:
    merged = [_norm_snv(r) for r in _read_rows(snv_path) if r and r[0]]
    merged += [_norm_indel(r) for r in _read_rows(indel_path) if r and r[0]]

    all_path = f"{out_basename}.neoantigens.consensus.all"
    with open(all_path, "w") as fh:
        fh.write("\t".join(MERGED_HEADER) + "\n")
        for row in merged:
            fh.write("\t".join(row) + "\n")

    expanded_path = f"{out_basename}.neoantigens.consensus.all.expanded"
    with open(expanded_path, "w") as fh:
        fh.write("\t".join(MERGED_HEADER) + "\n")
        for row in merged:
            fh.write("\t".join(row) + "\n")
            try:
                nfam = int(row[11])
            except (ValueError, IndexError):
                nfam = 1
            for i in range(2, nfam + 1):
                copy = list(row)
                copy[2] = f"{row[2]}.{i}"   # suffix the coordinate
                fh.write("\t".join(copy) + "\n")

    return all_path, expanded_path
