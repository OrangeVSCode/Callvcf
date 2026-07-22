import gzip
import os
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
from quality_engine import QualityEvaluator, QualityJobManager, render_report, _svg_ld_decay, _svg_sample_qc, build_analysis_readiness
from quality_profiles import crop_catalog, profile_catalog, resolve_profile
from population_analysis import PopulationAnalyzer
from repair_engine import RepairExecutor, _quality_comparison
from phenotype_engine import PhenotypeAnalyzer, phenotype_catalog


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
    def test_phenotype_qc_blue_blup_outliers_and_artifacts(self):
        rows = ["sample,trait,value,year,location,latitude,longitude,replicate,group"]
        for sample_index, sample in enumerate(("C1", "C2", "C3", "C4")):
            for year_index, year in enumerate((2022, 2023, 2024)):
                for rep in (1, 2):
                    ph = 90 + sample_index * 3 + year_index * 2 + rep * .2
                    if sample == "C4" and year == 2024 and rep == 2:
                        ph = 400
                    location = "A" if year_index % 2 == 0 else "B"
                    rows.append(f"{sample},PH,{ph},{year},{location},{35+year_index},{110+year_index},{rep},G{sample_index%2+1}")
                    rows.append(f"{sample},FL,{33+sample_index*.4+year_index*.2},{year},{location},{35+year_index},{110+year_index},{rep},G{sample_index%2+1}")
        with tempfile.TemporaryDirectory(prefix="gpa-phenotype-") as temp_name:
            source = Path(temp_name) / "phenotype.csv"
            source.write_text("\n".join(rows) + "\n", encoding="utf-8")
            analyzer = PhenotypeAnalyzer()
            output = analyzer.analyze({
                "path": str(source), "output_dir": str(Path(temp_name) / "reports"),
                "format": "long", "cotton_species": "cotton_barbadense", "density_mode": "overlay",
                "outlier": {"iqr_factor": 1.5, "sigma": 3, "mad_z": 3.5, "tail_fraction": .005, "consensus": 2},
            })
            result = output["result"]
            self.assertEqual(result["summary"]["traits"], 2)
            self.assertEqual(result["summary"]["samples"], 4)
            self.assertEqual(result["summary"]["years"], 3)
            self.assertTrue(any(x["trait"] == "PH" for x in result["outliers"]))
            ph = next(x for x in result["traits"] if x["trait"] == "PH")
            self.assertEqual(ph["threshold"]["mean"], [100, 180])
            self.assertEqual(len(ph["year_summaries"]), 3)
            self.assertIsNotNone(ph["variance_components"]["broad_sense_h2_entry_mean"])
            names = {x["name"] for x in output["artifacts"]}
            self.assertTrue({"phenotype_report.html", "trait_statistics.tsv", "replicate_averages.tsv", "blue_blup_estimates.tsv", "outlier_candidates.tsv", "cotton_thresholds.tsv", "GPA_Accelerator_phenotype_QC.zip"}.issubset(names))
            report = analyzer.artifact(output["run_id"], "phenotype_report.html").read_text(encoding="utf-8")
            self.assertIn("时间动态", report)
            self.assertIn("地理分布动态", report)
            self.assertIn("BLUE/BLUP", report)

            from openpyxl import Workbook
            workbook = Workbook(); sheet = workbook.active; sheet.title = "Raw phenotype"
            sheet.append(rows[0].split(","))
            for line in rows[1:]: sheet.append(line.split(","))
            xlsx = Path(temp_name) / "phenotype.xlsx"; workbook.save(xlsx)
            xlsx_output = analyzer.analyze({"path": str(xlsx), "output_dir": str(Path(temp_name) / "xlsx-reports"), "format": "auto", "sheet": "Raw phenotype"})
            self.assertEqual(xlsx_output["result"]["input"]["source"]["format"], "xlsx")
            self.assertEqual(xlsx_output["result"]["summary"]["traits"], 2)

        catalog = phenotype_catalog()
        self.assertEqual({x["id"] for x in catalog["cotton_species"]}, {"cotton_upland", "cotton_barbadense", "cotton_herbaceum"})
        self.assertTrue(any(x["code"] == "PD" and x["single_range"] is None for x in catalog["cotton_traits"]))

    def test_crop_presets_are_editable_and_traceable(self):
        crops = crop_catalog()
        crop_ids = {item["id"] for item in crops}
        self.assertTrue({
            "rice", "wheat", "wheat_durum", "wheat_wild_emmer", "wheat_einkorn",
            "cotton", "cotton_barbadense", "cotton_herbaceum",
            "soybean", "maize", "arabidopsis", "tomato", "potato", "tobacco",
        }.issubset(crop_ids))
        by_id = {item["id"]: item for item in crops}
        self.assertEqual((by_id["wheat"]["ploidy"], by_id["wheat"]["subgenomes"]), (6, 3))
        self.assertEqual((by_id["wheat_durum"]["ploidy"], by_id["wheat_durum"]["subgenomes"]), (4, 2))
        self.assertEqual((by_id["wheat_einkorn"]["ploidy"], by_id["wheat_einkorn"]["subgenomes"]), (2, 1))
        self.assertEqual(by_id["cotton_herbaceum"]["scientific_name"], "Gossypium herbaceum")
        self.assertEqual((by_id["cotton"]["ploidy"], by_id["cotton_barbadense"]["ploidy"]), (4, 4))
        self.assertEqual(by_id["cotton_herbaceum"]["ploidy"], 2)
        rice = resolve_profile({"crop_id": "rice", "profile_id": "plant_inbred"})
        self.assertEqual(rice["analysis_scope"], "plant_crop_preset")
        self.assertEqual(rice["thresholds"]["median_dp_min_warn"], 10)
        self.assertEqual(rice["thresholds"]["expected_titv_min"], 2.3)
        self.assertEqual(rice["threshold_sources"]["median_dp_min_warn"], "作物预设：水稻")
        edited = resolve_profile({
            "crop_id": "rice", "profile_id": "plant_inbred",
            "species_name": "本地水稻群体", "thresholds": {"median_dp_min_warn": 8, "site_missing_warn": .10, "expected_titv_min": None},
        })
        self.assertEqual(edited["thresholds"]["median_dp_min_warn"], 8)
        self.assertEqual(edited["threshold_sources"]["median_dp_min_warn"], "用户自定义（基于水稻）")
        self.assertEqual(edited["threshold_sources"]["site_missing_warn"], "作物预设：水稻")
        self.assertIsNone(edited["thresholds"]["expected_titv_min"])
        self.assertEqual(edited["threshold_sources"]["expected_titv_min"], "用户自定义（基于水稻）")
        self.assertEqual(edited["field_sources"]["species_name"], "用户自定义（基于水稻）")
        maize = resolve_profile({"crop_id": "maize", "profile_id": "generic_diploid_plant"})
        self.assertIsNone(maize["thresholds"]["het_rate_warn"])
        grass_cotton = resolve_profile({"crop_id": "cotton_herbaceum", "profile_id": "plant_inbred"})
        self.assertEqual(grass_cotton["species_name"], "草棉（二倍体）（Gossypium herbaceum）")
        self.assertEqual(grass_cotton["analysis_scope"], "plant_crop_preset")
        self.assertEqual(grass_cotton["subgenomes"], 1)
        with self.assertRaisesRegex(Exception, "作物预设与Profile不一致"):
            resolve_profile({"crop_id": "wheat", "profile_id": "plant_inbred"})

    def test_quality_profiles_are_plant_focused_and_animal_requires_explicit_custom_parameters(self):
        catalog = profile_catalog()
        self.assertFalse(any(x["id"] == "generic_diploid_animal" for x in catalog))
        self.assertTrue(all(x["kingdom"] == "plant" for x in catalog if x["id"] != "custom"))
        with self.assertRaisesRegex(Exception, "内置Profile仅用于植物"):
            resolve_profile({"profile_id": "generic_diploid_plant", "kingdom": "animal"})

        animal = {
            "profile_id": "custom", "kingdom": "animal", "species_name": "测试动物",
            "ploidy": 2, "genotype_ploidy": 2, "subgenomes": 1, "mating_system": "outcrossing",
            "thresholds": {
                "sample_missing_warn": .05, "sample_missing_critical": .10,
                "site_missing_warn": .10, "site_missing_critical": .20,
                "median_dp_min_warn": 5, "median_dp_max_warn": 100,
                "median_gq_warn": 20, "median_gq_critical": 10,
                "ab_lower": .25, "ab_upper": .75, "maf_structure": .05, "hwe_p_warn": 1e-6,
            },
        }
        with self.assertRaisesRegex(Exception, "必须确认"):
            resolve_profile(animal)
        animal["animal_parameters_confirmed"] = True
        resolved = resolve_profile(animal)
        self.assertEqual(resolved["analysis_scope"], "non_plant_custom")
        self.assertEqual(resolved["parameter_responsibility"], "user")
        self.assertTrue(all(resolved["threshold_sources"][key] == "用户自定义" for key in animal["thresholds"]))

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
        self.assertIn("quality_field_summaries", result["site_metrics"])
        self.assertIn("phase_rate", result["samples"][0])
        self.assertIn("module_availability", result)
        self.assertEqual(result["qc_recommendation"]["parameters"]["max_site_missing"], .10)
        self.assertFalse(result["qc_recommendation"]["parameters"]["remove_samples"])
        self.assertIn(result["analysis_readiness"]["code"], {"ready", "caution", "hold"})
        self.assertIn("score_dimensions", result["analysis_readiness"])
        report = render_report(result)
        self.assertIn("VCF质量评估报告", report)
        self.assertIn("打印/另存为PDF", report)
        self.assertIn("report-nav", report)
        self.assertIn("filterWarnings", report)
        self.assertIn("分析", result["analysis_readiness"]["title"])
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
            self.assertTrue({"report.html", "report_summary.json", "sample_metrics.tsv", "variant_metrics.tsv", "site_metric_summary.tsv", "density_windows.tsv", "module_availability.tsv", "recommend_filters.tsv", "sv_metrics.tsv", "group_batch_metrics.tsv", "subgenome_metrics.tsv", "fake_heterozygosity_windows.tsv", "annotation_consequences.tsv", "analysis_priorities.tsv", "score_dimensions.tsv", "warnings.tsv", "run_manifest.json", "GPA_Accelerator_VCF_QC_report.zip"}.issubset(names))
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

    def test_quality_conditional_reference_metadata_and_bed(self):
        fixture = Path(__file__).resolve().parent / "fixtures" / "tiny.vcf"
        with tempfile.TemporaryDirectory(prefix="callvcf-qc-inputs-") as temp_name:
            root = Path(temp_name)
            fasta = root / "ref.fa"
            sequence = list("A" * 1000000)
            sequence[200] = "T"; sequence[299] = "N"
            fasta.write_bytes(b">1\n" + "".join(sequence).encode("ascii") + b"\n")
            Path(str(fasta) + ".fai").write_text("1\t1000000\t3\t1000000\t1000001\n", encoding="ascii")
            metadata = root / "samples.tsv"
            metadata.write_text("sample_id\tgroup\tbatch\nS1\tG1\tB1\nS2\tG1\tB1\nS3\tG2\tB2\nS4\tG2\tB2\n", encoding="utf-8")
            bed = root / "regions.bed"
            bed.write_text("1\t99\t101\trepeat\n", encoding="utf-8")
            result = QualityEvaluator(PurePythonVCFService()).evaluate(str(fixture), {
                "profile": {"profile_id": "generic_diploid_plant"}, "scan_mode": "full",
                "reference_path": str(fasta), "sample_meta_path": str(metadata), "region_bed_path": str(bed),
            })
        self.assertEqual(result["reference_audit"]["status"], "complete")
        self.assertEqual(result["reference_audit"]["mismatch_records"], 0)
        self.assertEqual(result["metadata_audit"]["matched_vcf_samples"], 4)
        self.assertEqual(result["region_audit"]["overlap_records"], 1)
        self.assertEqual(result["samples"][0]["group"], "G1")
        self.assertIn("sample_score", result["samples"][0])

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
        with tempfile.TemporaryDirectory(prefix="callvcf-repair-tags-") as temp_name:
            tags_plan = executor.plan({
                "action": "fill_tags_copy", "path": str(fixture),
                "output_path": str(Path(temp_name) / "tagged.vcf.gz"),
            })
        self.assertIn("+fill-tags", tags_plan["command_preview"])
        self.assertEqual(tags_plan["risk"], "dangerous")
        with tempfile.TemporaryDirectory(prefix="callvcf-repair-more-") as temp_name:
            dedup = executor.plan({"action": "deduplicate_copy", "path": str(fixture), "output_path": str(Path(temp_name) / "dedup.vcf.gz")})
            subset = executor.plan({"action": "subset_samples_copy", "path": str(fixture), "output_path": str(Path(temp_name) / "subset.vcf.gz"), "samples": "S1\nS2"})
            mask = executor.plan({"action": "mask_genotypes_copy", "path": str(fixture), "output_path": str(Path(temp_name) / "mask.vcf.gz"), "expression": "FMT/DP<5 || FMT/GQ<20"})
        self.assertIn("exact", dedup["command_preview"])
        self.assertIn("S1,S2", subset["command_preview"])
        self.assertIn("+setGT", mask["command_preview"])
        with tempfile.TemporaryDirectory(prefix="callvcf-qc-plan-") as temp_name:
            qc = executor.plan({
                "action": "quality_control_copy", "path": str(fixture),
                "output_path": str(Path(temp_name) / "qc.vcf.gz"),
                "quality_profile": {"profile_id": "generic_diploid_plant"},
                "qc_parameters": {"mode": "profile_adaptive", "max_site_missing": .1, "min_dp": 5, "min_gq": 20},
            })
        self.assertEqual(qc["risk"], "dangerous")
        self.assertEqual(qc["qc_parameters"]["max_site_missing"], .1)
        self.assertIn("前后复评", qc["command_preview"])

    def test_ld_decay_chart_auto_scales_small_r2_values(self):
        svg = _svg_ld_decay([
            {"distance_bin_kb": "0-10", "pair_count": 426, "mean_r2": 0.0231},
            {"distance_bin_kb": "500-1000", "pair_count": 8073, "mean_r2": 0.0185},
        ])
        self.assertIn("纵轴自动缩放", svg)
        self.assertIn("0.0231", svg)
        self.assertIn("8,073 pairs", svg)
        self.assertNotIn("0–1</text>", svg)

    def test_sample_qc_chart_uses_robust_axis_and_clipped_outlier(self):
        rows = [{"sample_id": "S{}".format(i), "missing_rate": .002 + i / 100000, "het_rate": .03, "sample_status": "pass"} for i in range(30)]
        rows.append({"sample_id": "OUT", "missing_rate": .4, "het_rate": .6, "sample_status": "warning"})
        svg = _svg_sample_qc(rows, {"sample_missing_warn": .05, "het_rate_warn": .05})
        self.assertIn("超出P98主区间", svg)
        self.assertIn("OUT", svg)
        self.assertIn("缺失提醒线", svg)

    def test_quality_comparison_reports_delta_and_next_steps(self):
        before = {"summary": {"score": 70, "status": "warning", "record_count": 100, "sample_count": 2, "critical_count": 0, "warning_count": 1}, "samples": [{"sample_status": "warning"}, {"sample_status": "pass"}], "warnings": [{"code": "A", "scope": "site", "target": "x", "advice": "复核A"}]}
        after = {"summary": {"score": 85, "status": "pass", "record_count": 80, "sample_count": 2, "critical_count": 0, "warning_count": 0}, "samples": [{"sample_status": "pass"}, {"sample_status": "pass"}], "warnings": []}
        result = _quality_comparison(before, after, {"max_site_missing": .1})
        self.assertEqual(result["change"]["score_delta"], 15)
        self.assertEqual(result["change"]["records_removed"], 20)
        self.assertEqual(result["change"]["resolved_warning_count"], 1)

    def test_readiness_and_page_workflow_connect_analysis_steps(self):
        readiness = build_analysis_readiness({
            "summary": {"score": 55},
            "warnings": [{"level": "critical", "scope": "reference", "target": "vcf", "code": "REF_MISMATCH", "message": "参考不一致", "evidence": "1/10", "advice": "核对FASTA"}],
            "samples": [], "header_audit": {"qc_evidence_score": 50},
        })
        self.assertEqual(readiness["code"], "hold")
        self.assertEqual(readiness["top_priorities"][0]["anchor"], "section-input")
        root = Path(__file__).resolve().parents[1]
        html = (root / "static" / "index.html").read_text(encoding="utf-8")
        js = (root / "static" / "app.js").read_text(encoding="utf-8")
        self.assertIn('id="workflowPanel"', html)
        self.assertIn('id="animalCustomGuard"', html)
        self.assertIn("动物（仅自定义）", html)
        self.assertIn("function activateTab", js)
        self.assertIn("function validateSpeciesFocus", js)
        self.assertIn("quality-samples", js)

    def test_quality_control_copy_runs_full_audit_with_simulated_backend(self):
        class Service(PurePythonVCFService):
            bcftools = "simulated-bcftools"

        class SimulatedExecutor(RepairExecutor):
            def _bcftools_command(self, arguments, path_indexes=()):
                return [str(x) for x in arguments]

            def _run_process(self, command, log, cancel):
                if command[0] == "index":
                    Path(command[command.index("-o") + 1]).write_bytes(b"simulated-csi")
                    return
                output = Path(command[command.index("-o") + 1])
                source = Path(command[-1] if command[0] in {"view", "norm", "sort"} else command[1])
                raw = source.read_bytes()
                if raw.startswith(b"\x1f\x8b"):
                    raw = gzip.decompress(raw)
                with gzip.open(output, "wb") as handle:
                    handle.write(raw)

            def _run_to_file(self, command, destination, log, cancel):
                Path(destination).write_text("SN\t0\tnumber of records:\t3\n", encoding="utf-8")

        fixture = Path(__file__).resolve().parent / "fixtures" / "tiny.vcf"
        with tempfile.TemporaryDirectory(prefix="callvcf-qc-execute-") as temp_name:
            output = Path(temp_name) / "qc.vcf.gz"
            previous_local = os.environ.get("LOCALAPPDATA")
            os.environ["LOCALAPPDATA"] = temp_name
            service = Service()
            service.bcftools = "simulated-bcftools"
            executor = SimulatedExecutor(service)
            plan = executor.plan({
                "action": "quality_control_copy", "path": str(fixture), "output_path": str(output),
                "quality_profile": {"profile_id": "generic_diploid_plant"},
                "qc_parameters": {"max_site_missing": .1, "min_dp": 5, "min_gq": 20},
            })
            executor.execute({"plan_id": plan["id"], "confirmation": plan["confirmation_phrase"]})
            for _ in range(100):
                status = executor.status(plan["id"])
                if status["status"] in {"complete", "failed", "cancelled"}:
                    break
                import time
                time.sleep(.02)
            self.assertEqual(status["status"], "complete", status.get("error"))
            self.assertTrue(output.is_file())
            self.assertIn("qc_comparison", status["result"])
            self.assertTrue(Path(status["result"]["qc_comparison_report"]).is_file())
            if previous_local is None:
                os.environ.pop("LOCALAPPDATA", None)
            else:
                os.environ["LOCALAPPDATA"] = previous_local

    def test_population_pair_evidence_is_merged_conservatively(self):
        result = {
            "summary": {"score": 92, "status": "pass", "critical_count": 0, "warning_count": 0},
            "warnings": [],
            "population_analysis": {"modules": {
                "kinship": {"status": "complete", "summary": {"critical_pairs": 1, "warning_pairs": 2}},
                "pca": {"status": "complete", "summary": {"robust_pc_outliers": 1}},
                "hwe": {"status": "complete", "summary": {"scoring_enabled": False, "tested_sites": 100, "p_below_profile_threshold": 100}},
            }},
        }
        QualityJobManager._merge_population_evidence(result)
        codes = {item["code"] for item in result["warnings"]}
        self.assertIn("POSSIBLE_DUPLICATE_PAIRS", codes)
        self.assertIn("RELATED_SAMPLE_PAIRS", codes)
        self.assertIn("PCA_ROBUST_OUTLIERS", codes)
        self.assertNotIn("HWE_EXCESS_DEVIATION", codes)
        self.assertEqual(result["summary"]["status"], "critical")
        self.assertLess(result["summary"]["score"], 92)


if __name__ == "__main__":
    unittest.main()
