#!/usr/bin/env python3
import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from quality_engine import QualityJobManager
from vcf_service import create_service


def main():
    parser = argparse.ArgumentParser(description="Generate a persistent GPA-Accelerator quality report for manual validation")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--profile", default="cotton_inbred")
    parser.add_argument("--target-records", type=int, default=20000)
    args = parser.parse_args()

    manager = QualityJobManager(create_service())
    job = manager.start({
        "path": args.vcf,
        "profile": {"profile_id": args.profile},
        "scan_mode": "smart",
        "target_records": args.target_records,
        "output_dir": args.output,
    })
    while True:
        status = manager.status(job["id"])
        if status["status"] in {"complete", "failed", "cancelled"}:
            print(json.dumps(status, ensure_ascii=False, indent=2))
            if status["status"] != "complete":
                raise SystemExit(1)
            return
        time.sleep(0.5)


if __name__ == "__main__":
    main()
