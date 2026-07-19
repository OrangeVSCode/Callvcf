import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vcf_service import classify_variant, genotype_alleles, normalize_genotype, parse_loci
from vcf_service import PurePythonVCFService
from advanced_analysis import AdvancedAnalyzer, genotype_dosage, pairwise_r2


class CoreTests(unittest.TestCase):
    def test_genotype_categories(self):
        self.assertEqual(normalize_genotype("0/0"), "HOM_REF")
        self.assertEqual(normalize_genotype("0|1"), "HET")
        self.assertEqual(normalize_genotype("1/1"), "HOM_ALT")
        self.assertEqual(normalize_genotype("./."), "MISSING")
        self.assertEqual(normalize_genotype("1/2"), "OTHER")

    def test_allele_rendering(self):
        self.assertEqual(genotype_alleles("0/1", "A", "G"), "A/G")
        self.assertEqual(genotype_alleles("1|1", "TC", "T"), "T|T")
        self.assertEqual(genotype_alleles("./.", "A", "G"), ".")

    def test_variant_types(self):
        self.assertEqual(classify_variant("A", "G", "SNP"), "SNP")
        self.assertEqual(classify_variant("AT", "A", "INDEL"), "INDEL")
        self.assertEqual(classify_variant("N", "<DEL>", "OTHER"), "SV")
        self.assertEqual(classify_variant("N", "N]3:100]", "BND"), "SV")
        self.assertEqual(classify_variant("A", "T", "SNP", "SVTYPE=DEL;SVLEN=120"), "SV")

    def test_locus_parsing(self):
        self.assertEqual(parse_loci("5:10\n3 20, D07:30"), [("5", 10), ("3", 20), ("D07", 30)])
        self.assertEqual(parse_loci(["chr1:5", "chr1:5"]), [("chr1", 5)])

    def test_ld_math(self):
        self.assertEqual(genotype_dosage("0|1:99"), 1)
        self.assertIsNone(genotype_dosage("1/2"))
        r2, n = pairwise_r2([0, 1, 2, None], [0, 1, 2, 0], min_samples=3)
        self.assertAlmostEqual(r2, 1.0)
        self.assertEqual(n, 3)

    def test_advanced_analysis(self):
        fixtures = Path(__file__).resolve().parent / "fixtures"
        analyzer = AdvancedAnalyzer(PurePythonVCFService())
        result = analyzer.analyze({
            "path": str(fixtures / "tiny.vcf"), "lead_locus": "1:100",
            "window_kb": 1, "r2_threshold": 0.8, "min_samples": 3,
            "options": {"ld": True, "gene": True, "function": True,
                        "domain": True, "phenotype": True, "ldblockshow": False},
            "gff_path": str(fixtures / "tiny.gff3"),
            "annotation_path": str(fixtures / "annotation.tsv"),
            "domain_path": str(fixtures / "domains.tsv"),
            "phenotype_path": str(fixtures / "tiny.ps"),
        })
        self.assertEqual(result["ld"]["linked_variant_count"], 2)
        self.assertEqual(result["ld"]["linkage_region"]["end"], 300)
        self.assertEqual(result["gene"]["matches"][0]["relation"], "CDS")
        self.assertEqual(result["function"]["vcf_info"][0]["effect"], "missense_variant")
        self.assertEqual(result["domain"]["matches"][0]["domain"], "Kinase")
        self.assertEqual(result["phenotype"]["traits"][0]["min_p"], 0.001)


if __name__ == "__main__":
    unittest.main()
