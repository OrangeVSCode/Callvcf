import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vcf_service import classify_variant, genotype_alleles, normalize_genotype, parse_loci
from vcf_service import PurePythonVCFService
from advanced_analysis import AdvancedAnalyzer, genotype_dosage, pairwise_r2, pairwise_dprime, parse_region, parse_trait_directions
from tool_manager import tools_status


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

    def test_optional_tool_status(self):
        status = tools_status()
        self.assertIn("plink", status)
        self.assertIn("ldblockshow", status)
        self.assertIn("tool_root", status)

    def test_locus_parsing(self):
        self.assertEqual(parse_loci("5:10\n3 20, D07:30"), [("5", 10), ("3", 20), ("D07", 30)])
        self.assertEqual(parse_loci(["chr1:5", "chr1:5"]), [("chr1", 5)])

    def test_ld_math(self):
        self.assertEqual(genotype_dosage("0|1:99"), 1)
        self.assertIsNone(genotype_dosage("1/2"))
        r2, n = pairwise_r2([0, 1, 2, None], [0, 1, 2, 0], min_samples=3)
        self.assertAlmostEqual(r2, 1.0)
        self.assertEqual(n, 3)
        dprime, dprime_n = pairwise_dprime([0, 1, 2, 0], [0, 1, 2, 0], min_samples=3)
        self.assertAlmostEqual(dprime, 1.0)
        self.assertEqual(dprime_n, 4)
        self.assertEqual(parse_region("1:10-500"), ("1", 10, 500))
        self.assertEqual(parse_trait_directions("Yield=HIGH\nDisease=LOW"), {"Yield": 1, "Disease": -1})

    def test_advanced_analysis(self):
        fixtures = Path(__file__).resolve().parent / "fixtures"
        analyzer = AdvancedAnalyzer(PurePythonVCFService())
        result = analyzer.analyze({
            "path": str(fixtures / "tiny.vcf"), "lead_locus": "1:100",
            "window_kb": 1, "r2_threshold": 0.8, "min_samples": 3,
            "options": {"ld": True, "gene_track": True, "gene": True, "function": True,
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
        self.assertEqual(result["phenotype"]["traits"][0]["min_p"], 0.000001)
        self.assertEqual(result["ld"]["heatmap"]["plotted_count"], 2)
        self.assertEqual(result["ld"]["heatmap"]["matrix_dprime"][0][0], 1.0)
        self.assertEqual(result["gene_track"]["models"][0]["name"], "GeneA")

    def test_multi_lead_and_sample_profile(self):
        fixtures = Path(__file__).resolve().parent / "fixtures"
        analyzer = AdvancedAnalyzer(PurePythonVCFService())
        result = analyzer.analyze({
            "path": str(fixtures / "tiny.vcf"), "lead_locus": "1:100",
            "lead_loci": "1:100\n1:300", "ld_mode": "multi_lead_region",
            "region": "1:1-500", "window_kb": 1, "r2_threshold": 0.8,
            "min_samples": 3, "heatmap_max_variants": 120,
            "options": {"ld": True, "phenotype": True},
            "phenotype_path": str(fixtures / "tiny.ps"),
        })
        self.assertEqual(result["ld"]["mode"], "multi_lead_region")
        self.assertEqual(len(result["ld"]["lead_results"]), 2)
        self.assertEqual(result["ld"]["linked_variant_count"], 2)
        self.assertEqual(result["ld"]["heatmap"]["plotted_count"], 2)
        self.assertEqual(result["ld"]["heatmap"]["matrix"][0][1], result["ld"]["heatmap"]["matrix"][1][0])
        self.assertEqual(result["ld"]["heatmap"]["matrix_dprime"][0][1], result["ld"]["heatmap"]["matrix_dprime"][1][0])

        profile = analyzer.sample_lead_profile({
            "path": str(fixtures / "tiny.vcf"), "samples": ["S2"],
            "lead_loci": "1:100\n1:300", "phenotype_path": str(fixtures / "tiny.ps"),
            "trait_directions": "tiny=HIGH", "significance_threshold": 0.01,
        })
        self.assertEqual(profile["summaries"][0]["distinct_leads"], 2)
        self.assertEqual(profile["summaries"][0]["significant_records"], 2)
        self.assertIsNotNone(profile["summaries"][0]["trend_index"])
        no_direction = analyzer.sample_lead_profile({
            "path": str(fixtures / "tiny.vcf"), "samples": ["S2"],
            "lead_loci": "1:100", "phenotype_path": str(fixtures / "tiny.ps"),
            "trait_directions": "", "significance_threshold": 0.01,
        })
        self.assertIsNone(no_direction["summaries"][0]["trend_index"])


if __name__ == "__main__":
    unittest.main()
