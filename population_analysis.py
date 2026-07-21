"""Optional PLINK-backed population analyses for CallVCF quality reports."""

import csv
import bisect
import math
import os
import shutil
import statistics
import subprocess
import time
from pathlib import Path

from tool_manager import resolve_plink
from vcf_service import VCFError


DEFAULT_OPTIONS = {
    "hwe": True,
    "pca": True,
    "kinship": True,
    "ld": True,
    "roh": True,
    "site_missing": 0.10,
    "maf": None,
    "prune_window": 50,
    "prune_step": 5,
    "prune_r2": 0.20,
    "pca_components": 10,
    "ld_max_markers": 5000,
    "ld_window_kb": 1000,
}


def _float(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _read_space_table(path):
    path = Path(path)
    if not path.is_file() or not path.stat().st_size:
        return []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().strip().lstrip("#").split()
        return [dict(zip(header, line.split())) for line in handle if line.strip()]


def _write_tsv(path, rows, fields=None):
    path = Path(path)
    rows = list(rows)
    fields = list(fields or (rows[0].keys() if rows else []))
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


class PopulationAnalyzer:
    """Run bounded, auditable PLINK 1.9 analyses without changing the input VCF."""

    def __init__(self, plink=None):
        self.plink = resolve_plink(plink)
        self.log_path = None

    def _run(self, args, cancel, label):
        command = [self.plink] + [str(x) for x in args]
        creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
        started = time.time()
        with self.log_path.open("a", encoding="utf-8", errors="replace") as log:
            log.write("\n=== {} ===\n$ {}\n".format(label, " ".join(command)))
            log.flush()
            proc = subprocess.Popen(
                command, stdout=log, stderr=subprocess.STDOUT,
                creationflags=creationflags,
            )
            while proc.poll() is None:
                if cancel is not None and cancel.is_set():
                    proc.kill()
                    proc.wait()
                    raise VCFError("群体分析已取消")
                time.sleep(0.15)
            log.write("exit={} elapsed={:.2f}s\n".format(proc.returncode, time.time() - started))
        if proc.returncode:
            raise VCFError("{}失败（PLINK退出码{}；详见population_analysis.log）".format(label, proc.returncode))

    @staticmethod
    def _module(enabled, interpretation=None):
        return {
            "enabled": bool(enabled),
            "status": "pending" if enabled else "disabled",
            "interpretation": interpretation or "exploratory_non_scoring",
            "summary": {}, "preview": [], "artifacts": [], "error": None,
        }

    def run(self, vcf_path, profile, run_dir, requested=None, progress=None, cancel=None):
        options = dict(DEFAULT_OPTIONS)
        options.update(requested or {})
        for key in ("hwe", "pca", "kinship", "ld", "roh"):
            options[key] = bool(options.get(key))
        options["site_missing"] = min(1.0, max(0.0, float(options.get("site_missing", 0.10))))
        options["maf"] = float(options.get("maf") if options.get("maf") not in (None, "") else profile["thresholds"].get("maf_structure", 0.05))
        options["maf"] = min(0.5, max(0.0, options["maf"]))
        options["hwe_p_warn"] = profile["thresholds"].get("hwe_p_warn")

        hwe_interpretation = "quality_context" if profile.get("interpretation", {}).get("hwe_enabled") else "descriptive_only_not_scored_for_profile"
        modules = {
            "hwe": self._module(options["hwe"], hwe_interpretation),
            "pca": self._module(options["pca"]),
            "kinship": self._module(options["kinship"], "IBD_PI_HAT_not_KING"),
            "ld": self._module(options["ld"], "pruned_panel_ld_decay"),
            "roh": self._module(options["roh"], "exploratory_non_scoring"),
        }
        result = {
            "backend": "PLINK 1.9", "backend_path": self.plink,
            "status": "disabled", "options": options, "panel": {}, "modules": modules,
            "artifacts": [],
        }
        if not any(options[x] for x in modules):
            return result
        if not self.plink:
            result["status"] = "unavailable"
            for module in modules.values():
                if module["enabled"]:
                    module["status"] = "unavailable"
                    module["error"] = "未找到PLINK 1.9，请先在“本地软件资源”安装"
            return result

        run_dir = Path(run_dir)
        work = run_dir / ".population_work"
        work.mkdir(parents=True, exist_ok=True)
        self.log_path = run_dir / "population_analysis.log"
        self.log_path.write_text("CallVCF population analysis audit log\n", encoding="utf-8")
        result["artifacts"].append(self.log_path.name)
        panel = work / "panel"
        input_flag = "--bcf" if str(vcf_path).lower().endswith(".bcf") else "--vcf"
        try:
            if progress:
                progress(58, "正在建立群体分析标记面板")
            self._run([
                input_flag, vcf_path, "--allow-extra-chr", "--double-id",
                "--biallelic-only", "strict", "--set-missing-var-ids", "@:#:$1:$2",
                "--geno", options["site_missing"], "--maf", options["maf"],
                "--make-bed", "--out", panel,
            ], cancel, "建立二等位标记面板")
            with panel.with_suffix(".fam").open("r", encoding="utf-8", errors="replace") as handle:
                sample_count = sum(1 for _ in handle)
            with panel.with_suffix(".bim").open("r", encoding="utf-8", errors="replace") as handle:
                variant_count = sum(1 for _ in handle)
            result["panel"] = {"sample_count": sample_count, "variant_count": variant_count}
            if sample_count < 2 or variant_count < 1:
                raise VCFError("过滤后的群体分析面板样本或位点不足")
        except Exception as exc:
            result["status"] = "failed"
            result["panel"]["error"] = str(exc)
            for module in modules.values():
                if module["enabled"]:
                    module["status"] = "failed"
                    module["error"] = str(exc)
            shutil.rmtree(work, ignore_errors=True)
            return result

        prune_prefix = work / "prune"
        prune_ids = []
        need_prune = any(options[x] for x in ("pca", "kinship", "ld"))
        if need_prune:
            try:
                if progress:
                    progress(64, "正在进行LD剪枝")
                self._run([
                    "--bfile", panel, "--allow-extra-chr", "--indep-pairwise",
                    int(options["prune_window"]), int(options["prune_step"]), float(options["prune_r2"]),
                    "--out", prune_prefix,
                ], cancel, "LD剪枝")
                prune_file = Path(str(prune_prefix) + ".prune.in")
                prune_ids = [x.strip() for x in prune_file.read_text(encoding="utf-8", errors="replace").splitlines() if x.strip()]
                result["panel"]["pruned_variant_count"] = len(prune_ids)
            except Exception as exc:
                result["panel"]["pruning_error"] = str(exc)

        steps = [("hwe", 69, self._hwe), ("pca", 75, self._pca), ("kinship", 81, self._kinship), ("ld", 87, self._ld), ("roh", 93, self._roh)]
        for name, percent, method in steps:
            module = modules[name]
            if not module["enabled"]:
                continue
            if cancel is not None and cancel.is_set():
                raise VCFError("群体分析已取消")
            try:
                if name in {"pca", "kinship", "ld"} and not prune_ids:
                    raise VCFError("LD剪枝后没有足够标记")
                if progress:
                    progress(percent, "正在计算{}".format({"hwe":"HWE", "pca":"PCA", "kinship":"亲缘关系/IBD", "ld":"LD衰减", "roh":"ROH"}[name]))
                method(panel, prune_prefix, prune_ids, options, module, run_dir, cancel)
                module["status"] = "complete"
                result["artifacts"].extend(module["artifacts"])
            except Exception as exc:
                module["status"] = "failed"
                module["error"] = str(exc)

        complete = sum(x["status"] == "complete" for x in modules.values())
        failed = sum(x["status"] == "failed" for x in modules.values())
        result["status"] = "complete" if complete and not failed else "partial" if complete else "failed"
        result["module_counts"] = {"complete": complete, "failed": failed}
        shutil.rmtree(work, ignore_errors=True)
        return result

    def _hwe(self, panel, prune_prefix, prune_ids, options, module, run_dir, cancel):
        prefix = Path(panel).parent / "hwe"
        self._run(["--bfile", panel, "--allow-extra-chr", "--hardy", "midp", "--out", prefix], cancel, "HWE")
        rows = _read_space_table(str(prefix) + ".hwe")
        all_rows = [x for x in rows if x.get("TEST") == "ALL"] or rows
        threshold = options.get("hwe_p_warn")
        parsed = []
        for row in all_rows:
            p = _float(row.get("P"))
            parsed.append({"chrom": row.get("CHR"), "variant_id": row.get("SNP"), "test": row.get("TEST"), "obs_het": row.get("O(HET)"), "exp_het": row.get("E(HET)"), "p": p})
        out = _write_tsv(run_dir / "hwe_sites.tsv", parsed, ["chrom", "variant_id", "test", "obs_het", "exp_het", "p"])
        p_values = [x["p"] for x in parsed if x["p"] is not None]
        module["summary"] = {
            "tested_sites": len(p_values),
            "p_below_profile_threshold": sum(p < threshold for p in p_values) if threshold is not None else None,
            "profile_threshold": threshold,
            "scoring_enabled": module["interpretation"] == "quality_context",
        }
        module["preview"] = sorted(parsed, key=lambda x: x["p"] if x["p"] is not None else 2)[:20]
        module["artifacts"] = [out.name]

    def _pca(self, panel, prune_prefix, prune_ids, options, module, run_dir, cancel):
        with Path(str(panel) + ".fam").open("r", encoding="utf-8", errors="replace") as handle:
            sample_count = sum(1 for _ in handle)
        components = max(1, min(int(options["pca_components"]), sample_count - 1, len(prune_ids)))
        prefix = Path(panel).parent / "pca"
        self._run(["--bfile", panel, "--allow-extra-chr", "--extract", str(prune_prefix) + ".prune.in", "--pca", components, "tabs", "--out", prefix], cancel, "PCA")
        parsed = []
        with Path(str(prefix) + ".eigenvec").open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                values = line.strip().split()
                if len(values) < components + 2:
                    continue
                item = {"sample_id": values[1]}
                for index in range(1, components + 1):
                    item["PC{}".format(index)] = _float(values[index + 1])
                parsed.append(item)
        eigen = [_float(x) for x in Path(str(prefix) + ".eigenval").read_text(encoding="utf-8", errors="replace").splitlines()]
        eigen = [x for x in eigen if x is not None]
        total = sum(eigen)
        explained = [x / total if total else None for x in eigen]
        fields = ["sample_id"] + ["PC{}".format(x) for x in range(1, components + 1)]
        out = _write_tsv(run_dir / "pca_scores.tsv", parsed, fields)
        eigen_out = _write_tsv(run_dir / "pca_eigenvalues.tsv", [{"component": "PC{}".format(i + 1), "eigenvalue": value, "explained_fraction": explained[i]} for i, value in enumerate(eigen)], ["component", "eigenvalue", "explained_fraction"])
        module["summary"] = {"samples": len(parsed), "components": components, "pc1_relative_to_reported": explained[0] if explained else None, "pc2_relative_to_reported": explained[1] if len(explained) > 1 else None}
        module["preview"] = parsed[:1000]
        module["artifacts"] = [out.name, eigen_out.name]

    def _kinship(self, panel, prune_prefix, prune_ids, options, module, run_dir, cancel):
        prefix = Path(panel).parent / "kinship"
        self._run(["--bfile", panel, "--allow-extra-chr", "--extract", str(prune_prefix) + ".prune.in", "--genome", "full", "--out", prefix], cancel, "IBD/PI_HAT")
        rows = _read_space_table(str(prefix) + ".genome")
        parsed = []
        for row in rows:
            ibs0, ibs1, ibs2 = _int(row.get("IBS0")), _int(row.get("IBS1")), _int(row.get("IBS2"))
            compared = sum(x or 0 for x in (ibs0, ibs1, ibs2))
            discordance = ((ibs0 or 0) + (ibs1 or 0)) / compared if compared else None
            ibs_similarity = _float(row.get("DST"))
            pi_hat = _float(row.get("PI_HAT"))
            parsed.append({
                "sample_1": row.get("IID1"), "sample_2": row.get("IID2"),
                "z0": _float(row.get("Z0")), "z1": _float(row.get("Z1")), "z2": _float(row.get("Z2")),
                "pi_hat": pi_hat, "ibs_similarity": ibs_similarity, "genotype_discordance": discordance,
                "ibs0": ibs0, "ibs1": ibs1, "ibs2": ibs2,
            })
        raw_differences = [1 - x["ibs_similarity"] for x in parsed if x["ibs_similarity"] is not None]
        ranked = sorted(raw_differences)
        for item in parsed:
            similarity = item["ibs_similarity"]
            difference = 1 - similarity if similarity is not None else None
            if difference is None or not ranked:
                item["difference_score"] = None
            else:
                below = bisect.bisect_left(ranked, difference)
                equal = bisect.bisect_right(ranked, difference) - below
                item["difference_score"] = round(100 * (below + .5 * equal) / len(ranked), 3)
            pi_hat = item["pi_hat"] or 0
            discordance = item["genotype_discordance"]
            duplicate_evidence = similarity is not None and similarity >= .98 and discordance is not None and discordance <= .005
            if pi_hat >= .90 or duplicate_evidence:
                item["warning_level"] = "critical"
                item["relationship_hint"] = "疑似重复/同一样本"
            elif pi_hat >= .375:
                item["warning_level"] = "warning"
                item["relationship_hint"] = "约一级亲缘；需结合群体频率复核"
            elif pi_hat >= .177:
                item["warning_level"] = "info"
                item["relationship_hint"] = "约二级亲缘；探索性提示"
            else:
                item["warning_level"] = ""
                item["relationship_hint"] = ""
        parsed.sort(key=lambda x: x["pi_hat"] if x["pi_hat"] is not None else -1, reverse=True)
        pair_fields = ["sample_1", "sample_2", "z0", "z1", "z2", "pi_hat", "ibs_similarity", "genotype_discordance", "difference_score", "warning_level", "relationship_hint", "ibs0", "ibs1", "ibs2"]
        out = _write_tsv(run_dir / "pairwise_similarity.tsv", parsed, pair_fields)
        suspicious = _write_tsv(run_dir / "suspicious_pairs.tsv", [x for x in parsed if x["warning_level"]], pair_fields)
        nearest = []
        sample_ids = sorted({x["sample_1"] for x in parsed} | {x["sample_2"] for x in parsed})
        for sample_id in sample_ids:
            candidates = []
            for pair in parsed:
                if sample_id not in {pair["sample_1"], pair["sample_2"]}:
                    continue
                other = pair["sample_2"] if pair["sample_1"] == sample_id else pair["sample_1"]
                candidates.append((pair["ibs_similarity"] if pair["ibs_similarity"] is not None else -1, other, pair))
            for rank, (_, other, pair) in enumerate(sorted(candidates, reverse=True)[:10], 1):
                nearest.append({"sample_id": sample_id, "rank": rank, "neighbor": other, "ibs_similarity": pair["ibs_similarity"], "pi_hat": pair["pi_hat"], "genotype_discordance": pair["genotype_discordance"], "difference_score": pair["difference_score"], "warning_level": pair["warning_level"]})
        nearest_out = _write_tsv(run_dir / "sample_nearest_neighbors.tsv", nearest, ["sample_id", "rank", "neighbor", "ibs_similarity", "pi_hat", "genotype_discordance", "difference_score", "warning_level"])
        het_prefix = Path(panel).parent / "inbreeding"
        self._run(["--bfile", panel, "--allow-extra-chr", "--extract", str(prune_prefix) + ".prune.in", "--het", "--out", het_prefix], cancel, "样本近交系数F")
        het_rows = [{"sample_id": x.get("IID"), "observed_homozygotes": _int(x.get("O(HOM)")), "expected_homozygotes": _float(x.get("E(HOM)")), "nonmissing_autosomal": _int(x.get("N(NM)")), "inbreeding_f": _float(x.get("F"))} for x in _read_space_table(str(het_prefix) + ".het")]
        het_out = _write_tsv(run_dir / "sample_inbreeding.tsv", het_rows, ["sample_id", "observed_homozygotes", "expected_homozygotes", "nonmissing_autosomal", "inbreeding_f"])
        module["summary"] = {"pairs": len(parsed), "critical_pairs": sum(x["warning_level"] == "critical" for x in parsed), "warning_pairs": sum(x["warning_level"] == "warning" for x in parsed), "pi_hat_ge_0_9": sum((x["pi_hat"] or 0) >= .9 for x in parsed), "samples_with_f": len(het_rows), "estimator": "PLINK1.9 IBD PI_HAT + IBS DST; not KING"}
        module["preview"] = parsed[:20]
        module["artifacts"] = [out.name, suspicious.name, nearest_out.name, het_out.name]

    def _ld(self, panel, prune_prefix, prune_ids, options, module, run_dir, cancel):
        maximum = max(2, min(20000, int(options["ld_max_markers"])))
        if len(prune_ids) <= maximum:
            selected = prune_ids
        else:
            step = len(prune_ids) / maximum
            selected = [prune_ids[min(len(prune_ids) - 1, int(i * step))] for i in range(maximum)]
        if len(selected) < 2:
            raise VCFError("用于LD衰减的剪枝标记少于2个")
        marker_file = Path(panel).parent / "ld_markers.txt"
        marker_file.write_text("\n".join(selected) + "\n", encoding="utf-8")
        prefix = Path(panel).parent / "ld"
        self._run(["--bfile", panel, "--allow-extra-chr", "--extract", marker_file, "--r2", "gz", "--ld-window", 100, "--ld-window-kb", int(options["ld_window_kb"]), "--ld-window-r2", 0, "--out", prefix], cancel, "LD衰减")
        ld_path = Path(str(prefix) + ".ld.gz")
        if not ld_path.is_file():
            ld_path = Path(str(prefix) + ".ld")
        if ld_path.suffix == ".gz":
            import gzip
            handle = gzip.open(ld_path, "rt", encoding="utf-8", errors="replace")
        else:
            handle = ld_path.open("r", encoding="utf-8", errors="replace")
        bins = [(0, 10), (10, 50), (50, 100), (100, 250), (250, 500), (500, 1000)]
        values = {bounds: [] for bounds in bins}
        with handle:
            header = handle.readline().strip().split()
            for line in handle:
                row = dict(zip(header, line.split()))
                if row.get("CHR_A") != row.get("CHR_B"):
                    continue
                distance = abs((_int(row.get("BP_B")) or 0) - (_int(row.get("BP_A")) or 0)) / 1000.0
                r2 = _float(row.get("R2"))
                if r2 is None:
                    continue
                for bounds in bins:
                    if bounds[0] <= distance < bounds[1]:
                        values[bounds].append(r2)
                        break
        parsed = []
        for bounds in bins:
            current = values[bounds]
            parsed.append({"distance_bin_kb": "{}-{}".format(*bounds), "pair_count": len(current), "mean_r2": statistics.fmean(current) if current else None, "median_r2": statistics.median(current) if current else None})
        out = _write_tsv(run_dir / "ld_decay.tsv", parsed, ["distance_bin_kb", "pair_count", "mean_r2", "median_r2"])
        module["summary"] = {"markers": len(selected), "pairs": sum(x["pair_count"] for x in parsed), "window_kb": int(options["ld_window_kb"]), "panel": "LD-pruned uniformly capped panel"}
        module["preview"] = parsed
        module["artifacts"] = [out.name]

    def _roh(self, panel, prune_prefix, prune_ids, options, module, run_dir, cancel):
        prefix = Path(panel).parent / "roh"
        self._run(["--bfile", panel, "--allow-extra-chr", "--homozyg", "--homozyg-kb", 1000, "--homozyg-snp", 50, "--homozyg-window-snp", 50, "--out", prefix], cancel, "ROH")
        segments_raw = _read_space_table(str(prefix) + ".hom")
        samples_raw = _read_space_table(str(prefix) + ".hom.indiv")
        segments = [{"sample_id": x.get("IID"), "chrom": x.get("CHR"), "start_bp": _int(x.get("POS1")), "end_bp": _int(x.get("POS2")), "length_kb": _float(x.get("KB")), "snp_count": _int(x.get("NSNP"))} for x in segments_raw]
        samples = [{"sample_id": x.get("IID"), "roh_count": _int(x.get("NSEG")), "total_roh_kb": _float(x.get("KB")), "average_roh_kb": _float(x.get("KBAVG"))} for x in samples_raw]
        samples.sort(key=lambda x: x["total_roh_kb"] if x["total_roh_kb"] is not None else -1, reverse=True)
        seg_out = _write_tsv(run_dir / "roh_segments.tsv", segments, ["sample_id", "chrom", "start_bp", "end_bp", "length_kb", "snp_count"])
        sample_out = _write_tsv(run_dir / "roh_samples.tsv", samples, ["sample_id", "roh_count", "total_roh_kb", "average_roh_kb"])
        module["summary"] = {"samples": len(samples), "segments": len(segments), "total_roh_mb": sum((x["length_kb"] or 0) for x in segments) / 1000.0, "minimum_roh_kb": 1000}
        module["preview"] = samples[:20]
        module["artifacts"] = [seg_out.name, sample_out.name]
