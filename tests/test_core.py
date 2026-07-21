import gzip
import struct
import sys
import tempfile
import unittest
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vcf_service import classify_variant, detect_compression, genotype_alleles, normalize_genotype, parse_loci
from vcf_service import PurePythonVCFService
from advanced_analysis import AdvancedAnalyzer, genotype_dosage, pairwise_r2, pairwise_dprime, parse_region, parse_trait_directions
from tool_manager import tools_status
from quality_engine import QualityEvaluator, QualityJobManager, render_report, _svg_ld_decay
from quality_profiles import resolve_profile
from population_analysis import PopulationAnalyzer
from repair_engine import RepairExecutor


def bgzf_block(data):
    compressor = zlib.compressobj(level=6, wbits=-15)
    payload = compressor.compress(data) + compressor.flush()
    block_size = 18 + len(payload) + 8
    header = (
        b"\x1f\x8b\x08\x04" + b"\x00\x00\x00\x00" + b"\x00\xff" +
        struct.pack("<H", 6) + b"BC" + struct.pack("<H", 2) + struct.pack("<H", block_size - 1)
    )
    footer = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data) & 0xFFFFFFFF)
    return header + payload + footer


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

    def test_gzip_vcf_all_core_paths(self):
        fixtures = Path(__file__).resolve().parent / "fixtures"
        source = fixtures / "tiny.vcf"
        with tempfile.TemporaryDirectory(prefix="callvcf-gzip-test-") as temp_name:
            compressed = Path(temp_name) / "tiny.vcf.gz"
            with gzip.open(compressed, "wb", compresslevel=6) as handle:
                handle.write(source.read_bytes())

            service = PurePythonVCFService()
            metadata = service.inspect(str(compressed))
            self.assertEqual(detect_compression(compressed), "gzip")
            self.assertTrue(metadata["compressed"])
            self.assertEqual(metadata["compression"], "gzip")
            self.assertFalse(metadata["index_usable"])
            self.assertIn("no decompressed", metadata["space_mode"])
            self.assertEqual(metadata["sample_count"], 4)

            existence = service.check_loci(str(compressed), "1:100\n1:999")
            self.assertTrue(existence["results"][0]["exists"])
            self.assertFalse(existence["results"][1]["exists"])

            distribution = service.distributions(str(compressed), "1:100")
            self.assertEqual(distribution["records"][0]["counts"]["HET"], 1)
            matrix = service.sample_locus_matrix(str(compressed), "1:100\n1:300", ["S2"])
            self.assertEqual([x["gt"] for x in matrix["rows"][0]["values"]], ["0/1", "0/1"])
            stats = service.sample_stats(str(compressed), ["S2"])
            self.assertEqual(stats["processed_records"], 3)
            self.assertEqual(stats["results"][0]["counts"]["HET"], 2)
            counted = service.count_records(str(compressed))
            self.assertEqual(counted["record_count"], 3)
            self.assertFalse(counted["cached"])
            self.assertEqual(service.inspect(str(compressed))["record_count"], 3)
            self.assertTrue(service.count_records(str(compressed))["cached"])

            analyzer = AdvancedAnalyzer(service)
            result = analyzer.analyze({
                "path": str(compressed), "lead_locus": "1:100", "window_kb": 1,
                "r2_threshold": 0.8, "min_samples": 3,
                "options": {"ld": True, "gene_track": True, "gene": True, "function": True},
                "gff_path": str(fixtures / "tiny.gff3"),
            })
            self.assertEqual(result["ld"]["linked_variant_count"], 2)
            self.assertEqual(result["gene"]["matches"][0]["relation"], "CDS")
            self.assertEqual(result["function"]["vcf_info"][0]["effect"], "missense_variant")

    def test_bgzf_content_detection_and_direct_read(self):
        fixtures = Path(__file__).resolve().parent / "fixtures"
        with tempfile.TemporaryDirectory(prefix="callvcf-bgzf-test-") as temp_name:
            compressed = Path(temp_name) / "renamed.data"
            compressed.write_bytes(bgzf_block((fixtures / "tiny.vcf").read_bytes()))
            service = PurePythonVCFService()
            metadata = service.inspect(str(compressed))
            self.assertEqual(metadata["compression"], "bgzf")
            self.assertEqual(metadata["storage"], "BGZF VCF（压缩直读）")
            self.assertEqual(metadata["sample_count"], 4)
            self.assertTrue(service.check_loci(str(compressed), "1:300")["results"][0]["exists"])

    def test_quality_profiles_and_report(self):
        profile = resolve_profile({
            "profile_id": "cotton_inbred", "species_name": "陆地棉",
            "thresholds": {"sample_missing_warn": 0.03},
        })
        self.assertEqual(profile["ploidy"], 4)
        self.assertEqual(profile["subgenomes"], 2)
        self.assertEqual(profile["thresholds"]["sample_missing_warn"], 0.03)
        self.assertEqual(profile["threshold_sources"]["sample_missing_warn"], "用户自定义")
        self.assertFalse(profile["interpretation"]["hwe_enabled"])

        fixtures = Path(__file__).resolve().parent / "fixtures"
        evaluator = QualityEvaluator(PurePythonVCFService())
        result = evaluator.evaluate(str(fixtures / "tiny.vcf"), {
            "profile": {"profile_id": "generic_diploid_plant"},
            "scan_mode": "full",
        })
        self.assertEqual(result["summary"]["record_count"], 3)
        self.assertEqual(result["summary"]["sample_count"], 4)
        self.assertEqual(result["scan"]["evaluated_records"], 3)
        self.assertTrue(result["repair_policy"]["dangerous_requires_second_confirmation"])
        self.assertEqual(result["samples"][0]["missing_rate"], 0.0)
        report = render_report(result)
        self.assertIn("VCF质量评估报告", report)
        self.assertIn("打印/另存为PDF", report)
        self.assertNotIn("&quot;schema_version&quot;", report)

    def test_quality_job_artifacts(self):
        fixtures = Path(__file__).resolve().parent / "fixtures"
        with tempfile.TemporaryDirectory(prefix="callvcf-qc-") as temp_name:
            manager = QualityJobManager(PurePythonVCFService())
            job = manager.start({
                "path": str(fixtures / "tiny.vcf"), "output_dir": temp_name,
                "profile": {"profile_id": "plant_inbred"}, "scan_mode": "full",
            })
            for _ in range(100):
                current = manager.status(job["id"])
                if current["status"] in {"complete", "failed", "cancelled"}:
                    break
                import time
                time.sleep(0.02)
            self.assertEqual(current["status"], "complete", current.get("error"))
            names = {x["name"] for x in current["artifacts"]}
            self.assertTrue({"report.html", "report_summary.json", "sample_metrics.tsv", "warnings.tsv", "run_manifest.json", "CallVCF_QC_report.zip"}.issubset(names))
            self.assertTrue(manager.artifact(job["id"], "report.html").is_file())

    def test_sv_heterozygosity_uses_cohort_outliers(self):
        samples = ["S{:02d}".format(i) for i in range(30)]
        header = [
            "##fileformat=VCFv4.2",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
        ]
        rows = []
        for pos in range(1, 101):
            genotypes = []
            for index in range(30):
                # Every sample has the same 20% cohort-background SV het rate.
                gt = "0/1" if (pos + index) % 5 == 0 else "0/0"
                # One genuinely extreme sample is heterozygous at every site.
                if index == 29:
                    gt = "0/1"
                genotypes.append(gt)
            rows.append("1\t{}\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END={}\tGT\t{}".format(pos, pos + 10, "\t".join(genotypes)))
        with tempfile.TemporaryDirectory(prefix="callvcf-sv-het-") as temp_name:
            path = Path(temp_name) / "uniform_sv.vcf"
            path.write_text("\n".join(header + rows) + "\n", encoding="utf-8")
            result = QualityEvaluator(PurePythonVCFService()).evaluate(str(path), {
                "profile": {"profile_id": "cotton_inbred"}, "scan_mode": "full",
            })
        per_sample = [x for x in result["warnings"] if x["scope"] == "sample" and x["code"] in {"HIGH_HET", "HET_OUTLIER"}]
        self.assertLessEqual(len(per_sample), 1)
        self.assertEqual(per_sample[0]["target"], "S29")
        self.assertEqual(per_sample[0]["level"], "warning")
        self.assertTrue(any(x["code"] == "SV_HET_RELATIVE_ONLY" for x in result["warnings"]))
        self.assertEqual(result["cohort_metrics"]["het_rule_mode"], "cohort_relative_sv")

    def test_population_analysis_unavailable_is_nonfatal(self):
        analyzer = PopulationAnalyzer(plink=str(Path("definitely-missing-plink")))
        analyzer.plink = None
        with tempfile.TemporaryDirectory(prefix="callvcf-pop-unavailable-") as temp_name:
            result = analyzer.run(
                str(Path(__file__).resolve().parent / "fixtures" / "tiny.vcf"),
                resolve_profile({"profile_id": "cotton_inbred"}), temp_name,
                {"hwe": True, "pca": False, "kinship": False, "ld": False, "roh": False},
            )
        self.assertEqual(result["status"], "unavailable")
        self.assertEqual(result["modules"]["hwe"]["status"], "unavailable")
        self.assertEqual(result["modules"]["hwe"]["interpretation"], "descriptive_only_not_scored_for_profile")

    def test_repair_plan_blocks_overwrite_and_requires_unique_phrase(self):
        class Service:
            bcftools = "bcftools-test-double"
        executor = RepairExecutor(Service())
        fixture = Path(__file__).resolve().parent / "fixtures" / "tiny.vcf"
        with self.assertRaises(Exception):
            executor.plan({"action": "sort_copy", "path": str(fixture), "output_path": str(fixture)})
        with tempfile.TemporaryDirectory(prefix="callvcf-repair-plan-") as temp_name:
            plan = executor.plan({
                "action": "filter_copy", "path": str(fixture),
                "output_path": str(Path(temp_name) / "filtered.vcf.gz"),
                "expression": "QUAL>=30",
            })
        self.assertEqual(plan["risk"], "dangerous")
        self.assertTrue(plan["confirmation_phrase"].startswith("确认执行-"))
        self.assertIn("<temporary-output>", plan["command_preview"])

    def test_ld_decay_chart_auto_scales_small_r2_values(self):
        svg = _svg_ld_decay([
            {"distance_bin_kb": "0-10", "pair_count": 426, "mean_r2": 0.0231},
            {"distance_bin_kb": "500-1000", "pair_count": 8073, "mean_r2": 0.0185},
        ])
        self.assertIn("纵轴自动缩放", svg)
        self.assertIn("0.0231", svg)
        self.assertIn("8,073 pairs", svg)
        self.assertNotIn("0–1</text>", svg)


if __name__ == "__main__":
    unittest.main()
