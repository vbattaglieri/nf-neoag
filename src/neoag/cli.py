"""Command-line entrypoints: ``neo-snv`` (step 10) and ``neo-indel`` (step 11).

Each entrypoint reproduces the orchestration of the original driver scripts
(``Neoantigens.v4.sh`` / ``Neoantigens.indel.v3.sh``):

1. parse the variant table into per-(variant, k) peptides;
2. write one FASTA per peptide length k (8..11 by default);
3. run netMHC once per length and collect strong-binder (SB) hits;
4. pool hits per variant and write the per-variant consensus.

The identifier written into the FASTA (:mod:`neoag.ids`) is stable across the
separate per-length netMHC runs, so hits are joined back to the originating
variant regardless of which length produced them.
"""
from __future__ import annotations

import argparse
import os
import sys

from . import consensus as consensus_mod
from . import indel as indel_mod
from . import merge as merge_mod
from . import snv as snv_mod
from .netmhc import parse_strong_binders, resolve_netmhc, run_netmhc

DEFAULT_MIN_KMER = 8
DEFAULT_MAX_KMER = 11
DEFAULT_FAMILY_DIST = 3  # indel peptide "family" edit-distance threshold


def _write_fasta(entries: list[tuple[str, str]], path: str) -> None:
    with open(path, "w") as fh:
        for ident, peptide in entries:
            fh.write(f">{ident}\n{peptide}\n")


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("output", help="output basename (results go to "
                                  "<basename>.neoantigens.results.consensus)")
    p.add_argument("hla", help="comma-separated netMHC allele list, e.g. "
                               "HLA-A0201,HLA-B0702")
    p.add_argument("--netmhc", default=None,
                   help="path to the netMHC executable or its directory "
                        "(defaults to the NETMHC environment variable)")
    p.add_argument("--min-kmer", type=int, default=DEFAULT_MIN_KMER,
                   help="shortest peptide length (default: 8)")
    p.add_argument("--max-kmer", type=int, default=DEFAULT_MAX_KMER,
                   help="longest peptide length (default: 11)")
    p.add_argument("--keep-intermediate", action="store_true",
                   help="keep the per-length FASTA and netMHC files")


# ---------------------------------------------------------------------------
# SNV (step 10)
# ---------------------------------------------------------------------------
def snv_main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="neo-snv",
        description="Predict neoantigens from non-synonymous SNVs (pipeline step 10).",
    )
    p.add_argument("idea", help="IDEA table (*.snp.str-bias.Idea.OK)")
    p.add_argument("cdna", help="cDNA table (refFlat_mRNA.<build>.noalt.parsed.tab)")
    _add_common(p)
    args = p.parse_args(argv)

    netmhc = resolve_netmhc(args.netmhc)
    variants = snv_mod.parse_idea(args.idea)
    cdna = snv_mod.load_cdna(args.cdna)

    info_by_id: dict[str, list[str]] = {}
    records: list[tuple] = []

    for k in range(args.min_kmer, args.max_kmer + 1):
        entries: list[tuple[str, str]] = []
        for v in variants:
            seq = cdna.get(v.accession)
            if seq is None:
                continue
            peptide = snv_mod.mutant_peptide(seq, v, k)
            if peptide is None or "*" in peptide or len(peptide) < k:
                continue
            info_by_id.setdefault(v.ident, v.info)
            entries.append((v.ident, peptide))
        if not entries:
            continue
        fasta = f"{args.output}.{k}.fa"
        _write_fasta(entries, fasta)
        stdout = run_netmhc(netmhc, fasta, args.hla, k)
        for b in parse_strong_binders(stdout):
            info = info_by_id.get(b.ident)
            if info is not None:
                records.append((b.ident, info, [b.allele, b.peptide, b.affinity, b.rank]))
        if not args.keep_intermediate:
            os.remove(fasta)
        print(f"completed partial results for {k}..", file=sys.stderr)

    out_path = f"{args.output}.neoantigens.results.consensus"
    with open(out_path, "w") as fh:
        consensus_mod.write_snv_consensus(records, out=fh)
    print(f"wrote {out_path}", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# Indel (step 11)
# ---------------------------------------------------------------------------
def indel_main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="neo-indel",
        description="Predict neoantigens from indels/frameshifts (pipeline step 11).",
    )
    p.add_argument("indel", help="indel table (cosmic.Indel.altX)")
    p.add_argument("yar", help="YAR reference (refFlat_mRNA.<build>.yar)")
    _add_common(p)
    p.add_argument("--dist", type=int, default=DEFAULT_FAMILY_DIST,
                   help="edit-distance threshold for peptide family grouping "
                        "(default: 3)")
    args = p.parse_args(argv)

    netmhc = resolve_netmhc(args.netmhc)
    transcripts, sequences = indel_mod.load_transcripts(args.yar)
    peptides = indel_mod.generate(
        args.indel, transcripts, sequences, args.min_kmer, args.max_kmer)

    info_by_id: dict[str, list[str]] = {}
    by_kmer: dict[int, list[indel_mod.IndelPeptide]] = {}
    for pep in peptides:
        info_by_id.setdefault(pep.ident, pep.info)
        by_kmer.setdefault(pep.kmer, []).append(pep)

    records: list[tuple] = []
    for k in range(args.min_kmer, args.max_kmer + 1):
        entries = [(pep.ident, pep.peptide) for pep in by_kmer.get(k, [])
                   if len(pep.peptide) >= k]
        if not entries:
            continue
        fasta = f"{args.output}.{k}.fa"
        _write_fasta(entries, fasta)
        stdout = run_netmhc(netmhc, fasta, args.hla, k)
        for b in parse_strong_binders(stdout):
            info = info_by_id.get(b.ident)
            if info is not None:
                records.append((b.ident, info, [b.allele, b.peptide, b.affinity, b.rank]))
        if not args.keep_intermediate:
            os.remove(fasta)
        print(f"completed partial results for {k}..", file=sys.stderr)

    out_path = f"{args.output}.neoantigens.results.consensus"
    with open(out_path, "w") as fh:
        consensus_mod.write_indel_consensus(records, args.dist, out=fh)
    print(f"wrote {out_path}", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# Merge SNV + Indel consensus (optional)
# ---------------------------------------------------------------------------
def merge_main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="neo-merge",
        description="Merge SNV and Indel consensus tables into one.",
    )
    p.add_argument("indel_consensus", help="Indel *.neoantigens.results.consensus")
    p.add_argument("snv_consensus", help="SNV *.neoantigens.results.consensus")
    p.add_argument("output", help="output basename")
    args = p.parse_args(argv)

    all_path, expanded_path = merge_mod.merge(
        args.snv_consensus, args.indel_consensus, args.output)
    print(f"wrote {all_path}\nwrote {expanded_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(snv_main())
