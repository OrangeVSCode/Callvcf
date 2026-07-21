#!/usr/bin/env python3
"""Read-only smoke test for large real SNP/INDEL/SV .vcf.gz files.

The source files are never modified or fully decompressed.  A small first-contig
excerpt is written directly as gzip inside a temporary directory and deleted at
the end of the run.
"""
import argparse
import gzip
import json
import math
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from advanced_analysis import AdvancedAnalyzer, genotype_dosage
from vcf_service import PurePythonVCFService


def make_excerpt(service, source, destination, max_records=400, max_span_bp=2_000_000):
    first_chrom = None
    first_pos = None
    positions = []
    with service._open_vcf(source) as reader, gzip.open(destination, "wt", encoding="utf-8", newline="") as writer:
        for line in reader:
            if line.startswith("#"):
                writer.write(line)
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos = parts[0], int(parts[1])
            if first_chrom is None:
                first_chrom, first_pos = chrom, pos
            if chrom != first_chrom or pos > first_pos + max_span_bp or len(positions) >= max_records:
                break
            writer.write(line)
            positions.append(pos)
    if not positions:
        raise RuntimeError("VCF 没有可测试的变异记录：{}".format(source))
    return first_chrom, positions


def choose_lead(records, min_samples=20):
    for record in records:
        if "," in record["alt"]:
            continue
        valid = sum(genotype_dosage(gt) is not None for gt in record["genotypes"].values())
        if valid >= min_samples:
            return record
    raise RuntimeError("测试片段中没有足够完整的二等位 Lead 位点")


def smoke_one(source, expected_type):
    service = PurePythonVCFService()
    source_meta = service.inspect(source)
    if source_meta["compression"] != "bgzf":
        raise AssertionError("真实输入应为 BGZF：{}".format(source_meta["storage"]))
    if source_meta["sample_count"] < 20:
        raise AssertionError("样本数不足，无法测试 LD")

    with tempfile.TemporaryDirectory(prefix="callvcf-real-vcfgz-") as temp_name:
        excerpt = Path(temp_name) / (Path(source).name + ".excerpt.vcf.gz")
        chrom, positions = make_excerpt(service, source, excerpt)
        excerpt_meta = service.inspect(str(excerpt), force=True)
        samples = excerpt_meta["samples"][:3]
        loci = "\n".join("{}:{}".format(chrom, pos) for pos in positions[:3])

        # Every VCF-dependent page path is exercised against a compressed excerpt.
        existence = service.check_loci(str(excerpt), loci)
        distribution = service.distributions(str(excerpt), loci)
        matrix = service.sample_locus_matrix(str(excerpt), loci, samples)
        stats = service.sample_stats(str(excerpt), samples)
        _, cohort_samples, records = service.query_region(
            str(excerpt), chrom, min(positions), max(positions), max_records=max(1000, len(positions) + 1)
        )
        lead = choose_lead(records)
        window_kb = max(1, math.ceil(max(abs(pos - lead["pos"]) for pos in positions) / 1000))
        analysis = AdvancedAnalyzer(service).analyze({
            "path": str(excerpt),
            "lead_locus": "{}:{}".format(chrom, lead["pos"]),
            "window_kb": min(window_kb, 50_000),
            "r2_threshold": 0.6,
            "min_samples": 20,
            "heatmap_max_variants": 120,
            "options": {"ld": True, "function": True},
        })

        observed = source_meta.get("observed_types", {})
        if observed.get(expected_type, 0) == 0:
            raise AssertionError("没有识别到预期类型 {}：{}".format(expected_type, observed))
        if not all(item["exists"] for item in existence["results"]):
            raise AssertionError("真实片段位点存在性检查失败")
        if not distribution["records"] or len(matrix["rows"]) != len(samples):
            raise AssertionError("基因型分布或样本矩阵检查失败")

        return {
            "file": source_meta["name"],
            "source_size": source_meta["file_size"],
            "source_storage": source_meta["storage"],
            "samples": source_meta["sample_count"],
            "contigs": source_meta["contig_count"],
            "expected_type": expected_type,
            "observed_types_first_500": observed,
            "excerpt_chrom": chrom,
            "excerpt_records": len(positions),
            "tested_samples": samples,
            "distribution_records": len(distribution["records"]),
            "sample_stats_records": stats["processed_records"],
            "ld_lead": analysis["lead_locus"],
            "ld_tested_variants": analysis["ld"]["tested_variant_count"],
            "ld_linked_variants": analysis["ld"]["linked_variant_count"],
            "heatmap_variants": analysis["ld"]["heatmap"]["plotted_count"],
            "temporary_excerpt_deleted_after_test": True,
            "cohort_samples_in_excerpt": len(cohort_samples),
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snp", required=True)
    parser.add_argument("--indel", required=True)
    parser.add_argument("--sv", required=True)
    args = parser.parse_args()
    results = [
        smoke_one(args.snp, "SNP"),
        smoke_one(args.indel, "INDEL"),
        smoke_one(args.sv, "SV"),
    ]
    print(json.dumps({"ok": True, "results": results}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
