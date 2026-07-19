#!/usr/bin/env python3
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vcf_service import VCFService


def main(paths):
    service = VCFService()
    summaries = []
    for path in paths:
        meta = service.inspect(path, force=True)
        first = service._observe_types(Path(path), limit=1)
        if not first:
            raise RuntimeError("No records in {}".format(path))
        locus = "{}:{}".format(first[0]["chrom"], first[0]["pos"])
        checked = service.check_loci(path, locus)
        dist = service.distributions(path, locus)
        if not checked["results"][0]["exists"] or not dist["records"]:
            raise RuntimeError("Query failed for {} at {}".format(path, locus))
        summaries.append({
            "file": meta["name"],
            "storage": meta["storage"],
            "indexed": meta["indexed"],
            "samples": meta["sample_count"],
            "observed_types": meta["observed_types"],
            "test_locus": locus,
            "record_type": dist["records"][0]["record"]["variant_type"],
            "counts": dist["records"][0]["counts"],
        })
    print(json.dumps(summaries, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("Usage: integration_smoke.py VCF [VCF ...]")
    main(sys.argv[1:])

