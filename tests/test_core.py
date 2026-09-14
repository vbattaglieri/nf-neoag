"""Unit tests for the pure-Python analysis logic (no netMHC / Nextflow needed).

Run with:  python -m unittest discover -s tests
"""
import io
import unittest

from neoag import consensus, genetic_code as gc, ids, merge, snv


class GeneticCode(unittest.TestCase):
    def test_translate_in_frame(self):
        self.assertEqual(gc.translate("ATGGCC"), "MA")

    def test_translate_ignores_trailing(self):
        self.assertEqual(gc.translate("ATGGCCA"), "MA")

    def test_stop_codon(self):
        self.assertEqual(gc.translate("ATGTAA"), "M*")

    def test_reverse_complement(self):
        self.assertEqual(gc.reverse_complement("ATGC"), "GCAT")
        self.assertEqual(gc.reverse_complement("atgcN"), "Ngcat")


class VariantIds(unittest.TestCase):
    def test_stable_and_capped(self):
        a = ids.variant_id("NM_1", "chr1", 100)
        b = ids.variant_id("NM_1", "chr1", 100)
        c = ids.variant_id("NM_1", "chr1", 101)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertEqual(len(a), ids.ID_LENGTH)


class SnvHelpers(unittest.TestCase):
    def test_split_change(self):
        self.assertEqual(snv._split_change("p.S553A"), ("S", 553, "A"))
        self.assertEqual(snv._split_change("c.T1657G"), ("T", 1657, "G"))

    def test_codon_bounds(self):
        self.assertEqual(snv._codon_bounds(0, 1), (0, 3))
        self.assertEqual(snv._codon_bounds(2, 0), (0, 3))
        self.assertEqual(snv._codon_bounds(1, 2), (0, 3))

    def test_mutant_peptide_centres_mutation(self):
        # Long cDNA (40 codons) so the (2k-1) window is not truncated at the ends.
        # Mutate a codon well inside the transcript.
        cdna = "ATG" * 40
        v = snv.SnvVariant(
            accession="NM_x", gene="G", coord="chr1:1", aa_pos=20,
            aa_ref="M", aa_alt="A", nc_pos=58, nc_ref="A", nc_alt="G",
            freq="10",
        )
        pep = snv.mutant_peptide(cdna, v, kmer=8)
        self.assertIsNotNone(pep)
        # 2k-1 window => 15 residues, mutated residue at the centre (index k-1=7)
        self.assertEqual(len(pep), 15)
        self.assertEqual(pep[7], "A")


class Consensus(unittest.TestCase):
    def test_levenshtein(self):
        self.assertEqual(consensus.levenshtein("AAAA", "AAAB"), 1)
        self.assertEqual(consensus.levenshtein("", "ABC"), 3)

    def test_families(self):
        self.assertEqual(
            consensus.count_families(["AAAA", "AAAB", "ZZZZ"], dist=3), 2)

    def test_snv_consensus_header_only_when_empty(self):
        buf = io.StringIO()
        consensus.write_snv_consensus([], out=buf)
        self.assertEqual(buf.getvalue().strip(), "\t".join(consensus.SNV_HEADER))

    def test_snv_consensus_best_by_rank(self):
        info = ["G", "NM_1", "chr1:100", "c.A1T", "p.M1L", "5.0"]
        recs = [
            ("id1", info, ["HLA-A0201", "PEP1", "120", "0.5"]),
            ("id1", info, ["HLA-B0702", "PEP2", "300", "0.9"]),
        ]
        buf = io.StringIO()
        consensus.write_snv_consensus(recs, out=buf)
        rows = buf.getvalue().strip().splitlines()
        self.assertEqual(len(rows), 2)  # header + 1 variant
        fields = rows[1].split("\t")
        self.assertEqual(fields[7], "PEP1")  # best peptide = lowest rank
        self.assertEqual(fields[9], "0.5")


class Merge(unittest.TestCase):
    def test_merge_roundtrip(self, ):
        import tempfile
        import os
        d = tempfile.mkdtemp()
        snv_path = os.path.join(d, "s.consensus")
        indel_path = os.path.join(d, "i.consensus")
        with open(snv_path, "w") as fh:
            fh.write("\t".join(consensus.SNV_HEADER) + "\n")
            fh.write("\t".join(["G", "NM_1", "chr1:1", "c.A1T", "p.M1L",
                                "5.0", "HLA-A0201", "PEP", "100", "0.4"]) + "\n")
        with open(indel_path, "w") as fh:
            fh.write("\t".join(consensus.INDEL_HEADER) + "\n")
            fh.write("\t".join(["G2", "NM_2", "chr2:2", "Ins_T_fs5*", "AB*",
                                "9.0", "HLA-B0702", "PEP2", "50", "0.1",
                                "-", "-", "2"]) + "\n")
        all_path, exp_path = merge.merge(snv_path, indel_path, os.path.join(d, "out"))
        with open(all_path) as fh:
            all_rows = fh.read().strip().splitlines()
        with open(exp_path) as fh:
            exp_rows = fh.read().strip().splitlines()
        self.assertEqual(len(all_rows), 3)         # header + SNV + Indel
        self.assertEqual(len(exp_rows), 4)         # indel N.FAMILY=2 -> +1 row
        self.assertTrue(all_rows[1].endswith("SNV"))
        self.assertTrue(all_rows[2].endswith("Indel"))


if __name__ == "__main__":
    unittest.main()
