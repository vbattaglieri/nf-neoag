"""Genetic code utilities: translation and reverse-complement.

Identical codon table and semantics to the original scripts; '*' marks a stop.
"""
from __future__ import annotations

CODON_TABLE = {
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "AAT": "N", "AAC": "N",
    "GAT": "D", "GAC": "D",
    "TGT": "C", "TGC": "C",
    "CAA": "Q", "CAG": "Q",
    "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "AAA": "K", "AAG": "K",
    "TTT": "F", "TTC": "F",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TAA": "*", "TGA": "*", "TAG": "*",
}

_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "N": "N",
    "a": "t", "c": "g", "g": "c", "t": "a", "n": "n",
}


def translate_codon(codon: str, context: str = "") -> str:
    """Translate a single 3-nt codon to one amino-acid letter."""
    if len(codon) != 3:
        raise ValueError(f"wrong codon length: {codon!r} {context}".strip())
    return CODON_TABLE[codon.upper()]


def translate(dna: str) -> str:
    """Translate a DNA string in frame; trailing 1-2 nt are ignored."""
    return "".join(
        translate_codon(dna[i:i + 3]) for i in range(0, len(dna) - len(dna) % 3, 3)
    )


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA string (preserves case; allows N)."""
    try:
        return "".join(_COMPLEMENT[b] for b in reversed(seq))
    except KeyError as exc:  # pragma: no cover - mirrors original hard failure
        raise ValueError(f"illegal base in sequence: {seq!r}") from exc
