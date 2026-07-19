#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from advanced_analysis import AdvancedAnalyzer
from vcf_service import create_service


def main():
    parser = argparse.ArgumentParser(description="Run one real-VCF Lead LD smoke test")
    parser.add_argument("vcf")
    parser.add_argument("locus")
    parser.add_argument("--window-kb", type=float, default=100)
    parser.add_argument("--r2", type=float, default=0.6)
    parser.add_argument("--min-samples", type=int, default=20)
    args = parser.parse_args()
    result = AdvancedAnalyzer(create_service()).analyze({
        "path": args.vcf,
        "lead_locus": args.locus,
        "window_kb": args.window_kb,
        "r2_threshold": args.r2,
        "min_samples": args.min_samples,
        "options": {"ld": True},
    })
    ld = result["ld"]
    print(json.dumps({
        "lead": result["lead_record"],
        "tested_variant_count": ld["tested_variant_count"],
        "linked_variant_count": ld["linked_variant_count"],
        "linked_snp_count": ld["linked_snp_count"],
        "linkage_region": ld["linkage_region"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
