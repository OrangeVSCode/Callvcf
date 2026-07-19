#!/usr/bin/env python3
import json
import sys
from urllib.request import Request, urlopen


def post(base, route, payload):
    request = Request(
        base.rstrip("/") + route,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urlopen(request, timeout=120) as response:
        result = json.loads(response.read().decode("utf-8"))
    if not result.get("ok"):
        raise RuntimeError(result.get("error"))
    return result["data"]


def main(base, vcf):
    with urlopen(base.rstrip("/") + "/api/health", timeout=10) as response:
        health = json.loads(response.read().decode("utf-8"))
    assert health["ok"]

    meta = post(base, "/api/inspect", {"path": vcf})
    samples = meta["samples"][:2]
    loci = "5:2539274\n3:1394517"
    checked = post(base, "/api/check", {"path": vcf, "loci": loci})
    distribution = post(base, "/api/distribution", {"path": vcf, "loci": loci})
    matrix = post(base, "/api/sample-loci", {"path": vcf, "loci": loci, "samples": samples})
    stats = post(base, "/api/sample-stats", {"path": vcf, "loci": loci, "samples": samples})

    assert len(checked["results"]) == 2
    assert distribution["records"]
    assert len(matrix["rows"]) == 2
    assert len(stats["results"]) == 2
    print(json.dumps({
        "health": health["service"],
        "file": meta["name"],
        "sample_count": meta["sample_count"],
        "check_results": len(checked["results"]),
        "distribution_records": len(distribution["records"]),
        "matrix_rows": len(matrix["rows"]),
        "stats_rows": len(stats["results"]),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: api_smoke.py BASE_URL VCF")
    main(sys.argv[1], sys.argv[2])

