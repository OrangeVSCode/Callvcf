#!/usr/bin/env python3
"""Read-only quality smoke test for real compressed SNP/INDEL/SV files."""

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quality_engine import QualityEvaluator
from vcf_service import PurePythonVCFService


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snp", required=True)
    parser.add_argument("--indel", required=True)
    parser.add_argument("--sv", required=True)
    parser.add_argument("--target-records", type=int, default=20000)
    args = parser.parse_args()
    output = []
    for kind, path in (("SNP", args.snp), ("INDEL", args.indel), ("SV", args.sv)):
        started = time.time()
        evaluator = QualityEvaluator(PurePythonVCFService())
        result = evaluator.evaluate(path, {
            "profile": {"crop_id": "cotton", "profile_id": "cotton_inbred"},
            "scan_mode": "smart",
            "target_records": args.target_records,
        })
        output.append({
            "kind": kind,
            "summary": result["summary"],
            "scan": result["scan"],
            "variant_types": result["site_metrics"]["variant_types"],
            "median_het_rate": result["cohort_metrics"]["median_het_rate"],
            "elapsed_wall_seconds": round(time.time() - started, 3),
        })
    print(json.dumps({"ok": True, "results": output}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
