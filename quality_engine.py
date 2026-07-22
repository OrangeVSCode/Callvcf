"""Streaming VCF quality evaluation and self-contained report generation."""

import csv
import bisect
import gzip
import hashlib
import html
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import threading
import time
import uuid
import zipfile
from collections import Counter
from pathlib import Path

from quality_profiles import crop_catalog, profile_catalog, resolve_profile
from population_analysis import PopulationAnalyzer
from vcf_service import VCFError, classify_variant, detect_compression, parse_info


def _safe_float(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _safe_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _median_counter(counter):
    total = sum(counter.values())
    if not total:
        return None
    left = (total - 1) // 2
    right = total // 2
    found = []
    seen = 0
    for value, count in sorted(counter.items()):
        next_seen = seen + count
        if seen <= left < next_seen:
            found.append(value)
        if seen <= right < next_seen:
            found.append(value)
        if len(found) == 2:
            break
        seen = next_seen
    return round(sum(found) / len(found), 3) if found else None


def _median(values):
    values = [x for x in values if x is not None]
    return statistics.median(values) if values else None


def _mad(values, center=None):
    values = [x for x in values if x is not None]
    if not values:
        return None
    center = statistics.median(values) if center is None else center
    return statistics.median(abs(x - center) for x in values)


def _quantile(values, probability):
    values = sorted(x for x in values if x is not None)
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    position = max(0.0, min(1.0, probability)) * (len(values) - 1)
    left = int(math.floor(position))
    right = int(math.ceil(position))
    if left == right:
        return values[left]
    fraction = position - left
    return values[left] * (1 - fraction) + values[right] * fraction


def _fmt(value, digits=3):
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "是" if value else "否"
    if isinstance(value, int):
        return "{:,}".format(value)
    if isinstance(value, float):
        return ("{:,.%df}" % digits).format(value)
    return str(value)


def _ratio(value):
    return "—" if value is None else "{:.2%}".format(value)


def _level_rank(level):
    return {"info": 0, "warning": 1, "critical": 2}.get(level, 0)


def _warning(level, scope, target, code, message, evidence=None, advice=None):
    return {
        "level": level,
        "scope": scope,
        "target": target,
        "code": code,
        "message": message,
        "evidence": evidence or "",
        "advice": advice or "",
    }


def build_analysis_readiness(result):
    """Summarize whether the VCF is ready for downstream work and where the user should drill down next."""
    warnings = list(result.get("warnings") or [])
    summary = result.get("summary") or {}
    hard_codes = {
        "HEADER_REQUIRED", "UNSORTED", "MALFORMED", "CONTIG_LENGTH_MISMATCH", "REF_MISMATCH",
    }
    hard_blockers = [item for item in warnings if item.get("level") == "critical" and item.get("code") in hard_codes]
    critical = [item for item in warnings if item.get("level") == "critical"]
    score = int(summary.get("score") or 0)
    if hard_blockers or score < 60:
        code, title = "hold", "暂停下游分析，先处理关键问题"
        description = "检测到结构、参考一致性或其他高风险证据；修复或复核后应重新评估。"
    elif critical or score < 80:
        code, title = "caution", "可继续探索，但主分析前建议复核"
        description = "当前文件可用于定位问题和试算；正式关联分析前建议处理首页优先事项。"
    else:
        code, title = "ready", "质量证据支持继续下游分析"
        description = "未发现阻断性问题；仍建议保留报告、参数和原始VCF以保证可追溯。"

    anchor_by_scope = {
        "file": "section-input", "reference": "section-input", "metadata": "section-input",
        "site": "section-site", "region": "section-site", "annotation": "section-annotation",
        "sample": "section-samples", "sample_pair": "section-population", "cohort": "section-population",
        "batch": "section-population",
    }
    priorities = []
    for item in sorted(warnings, key=lambda row: (-_level_rank(row.get("level")), row.get("scope", ""), row.get("target", ""))):
        if item.get("level") not in {"critical", "warning"}:
            continue
        priorities.append({
            "level": item.get("level"), "code": item.get("code"), "scope": item.get("scope"),
            "target": item.get("target"), "message": item.get("message"), "evidence": item.get("evidence"),
            "advice": item.get("advice"), "anchor": anchor_by_scope.get(item.get("scope"), "section-warnings"),
        })
        if len(priorities) >= 5:
            break

    scope_counts = Counter((item.get("scope"), item.get("level")) for item in warnings)
    samples = result.get("samples") or []
    sample_score = _median([item.get("sample_score") for item in samples])
    file_reference_score = max(0, 100 - 18 * sum(scope_counts[(scope, "critical")] for scope in ("file", "reference")) - 6 * sum(scope_counts[(scope, "warning")] for scope in ("file", "reference")))
    site_score = max(0, 100 - 10 * scope_counts[("site", "critical")] - 3 * scope_counts[("site", "warning")])
    population_score = max(0, 100 - 15 * scope_counts[("sample_pair", "critical")] - 5 * scope_counts[("sample_pair", "warning")] - 3 * scope_counts[("cohort", "warning")])
    modules = (result.get("population_analysis") or {}).get("modules") or {}
    requested_modules = [module for module in modules.values() if module.get("enabled")]
    module_score = round(100 * sum(module.get("status") == "complete" for module in requested_modules) / len(requested_modules)) if requested_modules else None
    dimensions = {
        "file_reference": file_reference_score,
        "sample_median": round(sample_score) if sample_score is not None else None,
        "site": site_score,
        "population_integrity": population_score if requested_modules else None,
        "provenance": (result.get("header_audit") or {}).get("qc_evidence_score"),
        "requested_module_completion": module_score,
    }

    next_steps = []
    if hard_blockers:
        next_steps.append({"label": "核对输入与参考", "tab": "quality", "anchor": "section-input"})
    if any(item.get("scope") == "sample" and item.get("level") in {"critical", "warning"} for item in warnings):
        next_steps.append({"label": "查看异常样本", "tab": "stats", "anchor": "section-samples"})
    if any(item.get("scope") in {"sample_pair", "cohort", "batch"} and item.get("level") in {"critical", "warning"} for item in warnings):
        next_steps.append({"label": "查看群体与亲缘证据", "tab": "quality", "anchor": "section-population"})
    if result.get("qc_recommendation"):
        next_steps.append({"label": "核对智能质控参数", "tab": "quality", "anchor": "section-recommendations"})
    if code == "ready":
        next_steps.append({"label": "进入位点或Lead分析", "tab": "existence", "anchor": "section-overview"})
    return {
        "code": code, "title": title, "description": description, "score": score,
        "hard_blocker_count": len(hard_blockers), "top_priorities": priorities,
        "score_dimensions": dimensions, "next_steps": next_steps[:4],
    }


def build_qc_recommendation(result):
    """Build a conservative, auditable QC preset from the selected profile and observed VCF fields."""
    profile = result.get("profile") or {}
    thresholds = profile.get("thresholds") or {}
    header = result.get("header_audit") or {}
    site = result.get("site_metrics") or {}
    format_ids = set(header.get("format_ids") or [])
    metric_summary = site.get("quality_field_summaries") or {}
    reference_ready = (result.get("reference_audit") or {}).get("status") == "complete"
    parameters = {
        "max_site_missing": float(thresholds.get("site_missing_warn", .10)),
        "max_sample_missing": float(thresholds.get("sample_missing_warn", .05)),
        "min_dp": float(thresholds.get("median_dp_min_warn", 5)) if "DP" in format_ids else None,
        "max_dp": None,
        "min_gq": float(thresholds.get("median_gq_warn", 20)) if "GQ" in format_ids else None,
        "min_qual": None,
        "min_maf": 0.0,
        "min_qd": None,
        "min_mq": None,
        "max_fs": None,
        "max_sor": None,
        "pass_only": False,
        "biallelic_only": False,
        "normalize": bool(reference_ready and site.get("nonminimal_indel_count")),
        "deduplicate": bool(site.get("adjacent_duplicate_count")),
        "remove_samples": False,
    }
    broad_types = site.get("broad_variant_types") or {}
    dominant_type = max(broad_types, key=broad_types.get) if broad_types else None
    dominant_fraction = broad_types.get(dominant_type, 0) / max(sum(broad_types.values()), 1)
    qual_reference = thresholds.get("site_qual_warn")
    if qual_reference is not None and (metric_summary.get("QUAL") or {}).get("coverage", 0) >= .80:
        parameters["min_qual"] = float(qual_reference)
    if dominant_fraction >= .90 and (metric_summary.get("QD") or {}).get("coverage", 0) >= .80:
        parameters["min_qd"] = 2.0
    if dominant_type == "SNP" and dominant_fraction >= .90:
        if (metric_summary.get("MQ") or {}).get("coverage", 0) >= .80:
            parameters["min_mq"] = 40.0
        if (metric_summary.get("FS") or {}).get("coverage", 0) >= .80:
            parameters["max_fs"] = 60.0
        if (metric_summary.get("SOR") or {}).get("coverage", 0) >= .80:
            parameters["max_sor"] = 3.0
    elif dominant_type == "INDEL" and dominant_fraction >= .90 and (metric_summary.get("FS") or {}).get("coverage", 0) >= .80:
        parameters["max_fs"] = 200.0
    reasons = [
        "位点缺失上限采用“{}”Profile的提醒线 {}".format(profile.get("name") or profile.get("id") or "当前", _ratio(parameters["max_site_missing"])),
        "低质量GT仅在VCF实际含DP/GQ时掩蔽；缺失字段不会被当作0",
        "默认保留稀有、多等位和非PASS记录，避免在未知caller/研究目的下过度过滤",
        "样本删除默认关闭；候选样本只列出，必须由用户显式启用",
    ]
    if any(parameters[key] is not None for key in ("min_qd", "min_mq", "max_fs", "max_sor")):
        reasons.append("检测到{}占比高且相应质量字段覆盖充分，启用caller-aware位点质量起始线".format(dominant_type))
    else:
        reasons.append("位点质量字段覆盖不足或VCF类型混合，默认不机械套用QD/MQ/FS/SOR阈值")
    if parameters["min_qual"] is not None:
        reasons.append("当前作物预设提供QUAL建议起始线且QUAL覆盖率≥80%；该值仍需结合caller和原始过滤流程复核")
    if reference_ready:
        reasons.append("已完成参考FASTA核验，可在存在非最简INDEL时启用标准化")
    else:
        reasons.append("未完成参考FASTA核验，自动方案不执行左对齐/标准化")
    candidates = [
        item.get("sample_id") for item in result.get("samples", [])
        if item.get("missing_rate") is not None and item["missing_rate"] >= parameters["max_sample_missing"]
    ]
    sampled_rows = result.get("_variant_metrics_rows") or []
    evaluable = [item for item in sampled_rows if item.get("missing_rate") is not None]
    estimated_removed = sum(item["missing_rate"] > parameters["max_site_missing"] for item in evaluable)
    mask_parts = []
    if parameters["min_dp"] is not None:
        mask_parts.append("FMT/DP<{}".format(_fmt(parameters["min_dp"])))
    if parameters["min_gq"] is not None:
        mask_parts.append("FMT/GQ<{}".format(_fmt(parameters["min_gq"])))
    return {
        "mode": "profile_adaptive",
        "profile_id": profile.get("id"),
        "profile_name": profile.get("name"),
        "parameters": parameters,
        "available_fields": {
            "format": sorted(format_ids),
            "site_quality": {key: value.get("coverage") for key, value in metric_summary.items()},
        },
        "sample_exclusion_candidates": [x for x in candidates if x],
        "estimated_site_removal_fraction": estimated_removed / len(evaluable) if evaluable else None,
        "estimated_from_records": len(evaluable),
        "mask_expression_preview": " || ".join(mask_parts) or None,
        "site_expression_preview": "F_MISSING<={}".format(_fmt(parameters["max_site_missing"])),
        "reasons": reasons,
    }


def parse_gt(gt_text):
    gt = str(gt_text or ".").replace("|", "/")
    parts = gt.split("/")
    if not parts or any(x in {"", "."} for x in parts):
        return None
    try:
        alleles = [int(x) for x in parts]
    except ValueError:
        return None
    return alleles


def genotype_class(alleles):
    if not alleles:
        return "missing"
    unique = set(alleles)
    if unique == {0}:
        return "hom_ref"
    if len(unique) > 1:
        return "het"
    if len(unique) == 1 and next(iter(unique)) > 0:
        return "hom_alt"
    return "other"


def _metric_summary(values, evaluated, available_n=None):
    values = [x for x in values if x is not None]
    available_n = len(values) if available_n is None else available_n
    return {
        "available": bool(available_n), "n": available_n,
        "coverage": available_n / evaluated if evaluated else None,
        "quantile_sample_n": len(values),
        "min": min(values) if values else None,
        "q05": _quantile(values, .05), "median": _median(values),
        "q95": _quantile(values, .95), "max": max(values) if values else None,
    }


def _detail_variant_type(ref, alt, info, broad_type):
    svtype = str(info.get("SVTYPE") or "").upper()
    if svtype:
        return svtype
    if broad_type == "SV":
        if "[" in alt or "]" in alt:
            return "BND"
        return "SV"
    alleles = alt.split(",")
    if alleles and all(len(ref) == len(value) > 1 for value in alleles):
        return "MNP"
    return broad_type


def _sv_length(info, pos):
    raw = str(info.get("SVLEN") or "").split(",", 1)[0]
    value = _safe_int(raw)
    if value is not None:
        return abs(value)
    end = _safe_int(info.get("END"))
    return abs(end - pos) + 1 if end is not None else None


class IndexedFasta:
    """Minimal read-only FASTA/.fai accessor; never creates or edits an index."""

    def __init__(self, path_text):
        self.path = Path(path_text).expanduser().resolve()
        if not self.path.is_file():
            raise VCFError("参考基因组FASTA不存在：{}".format(self.path))
        if self.path.name.lower().endswith((".gz", ".bgz")):
            raise VCFError("压缩FASTA的随机REF核验需要bcftools/samtools；请提供未压缩FASTA及同名.fai")
        fai = Path(str(self.path) + ".fai")
        if not fai.is_file():
            raise VCFError("未找到FASTA索引：{}；请先生成.fai，软件不会静默修改参考文件".format(fai))
        self.entries = {}
        with fai.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                values = line.rstrip("\r\n").split("\t")
                if len(values) < 5:
                    continue
                try:
                    self.entries[values[0]] = tuple(int(x) for x in values[1:5])
                except ValueError:
                    continue
        if not self.entries:
            raise VCFError("FASTA .fai中没有可用contig记录")
        self.handle = self.path.open("rb")

    def fetch(self, chrom, pos, length):
        entry = self.entries.get(chrom)
        if not entry or pos < 1 or length < 1 or pos + length - 1 > entry[0]:
            return None
        seq_len, offset, line_bases, line_width = entry
        zero = pos - 1
        byte_offset = offset + (zero // line_bases) * line_width + (zero % line_bases)
        self.handle.seek(byte_offset)
        chunks = []
        remaining = length
        while remaining > 0:
            chunk = self.handle.read(remaining + 8)
            if not chunk:
                break
            clean = chunk.replace(b"\r", b"").replace(b"\n", b"")
            take = clean[:remaining]
            chunks.append(take)
            remaining -= len(take)
        return b"".join(chunks).decode("ascii", errors="replace").upper() if remaining == 0 else None

    def close(self):
        self.handle.close()


class BedRegions:
    def __init__(self, path_text):
        path = Path(path_text).expanduser().resolve()
        if not path.is_file():
            raise VCFError("区域BED不存在：{}".format(path))
        opener = gzip.open if path.name.lower().endswith(".gz") else open
        raw = {}
        with opener(str(path), "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                values = line.rstrip("\r\n").split("\t")
                if len(values) < 3:
                    values = line.split()
                try:
                    chrom, start, end = values[0], int(values[1]), int(values[2])
                except (IndexError, ValueError):
                    continue
                if end > start >= 0:
                    raw.setdefault(chrom, []).append((start + 1, end))
        self.regions = {}
        for chrom, intervals in raw.items():
            merged = []
            for start, end in sorted(intervals):
                if merged and start <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self.regions[chrom] = (merged, [x[0] for x in merged])
        self.interval_count = sum(len(x[0]) for x in self.regions.values())

    def contains(self, chrom, pos):
        data = self.regions.get(chrom)
        if not data:
            return False
        intervals, starts = data
        index = bisect.bisect_right(starts, pos) - 1
        return index >= 0 and intervals[index][0] <= pos <= intervals[index][1]


def load_sample_metadata(path_text):
    path = Path(path_text).expanduser().resolve()
    if not path.is_file():
        raise VCFError("样本元数据不存在：{}".format(path))
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    first = next((line for line in text.splitlines() if line.strip()), "")
    delimiter = "\t" if "\t" in first else ","
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    if not reader.fieldnames:
        raise VCFError("样本元数据缺少表头")
    aliases = {str(x).strip().lower(): x for x in reader.fieldnames}
    id_key = next((aliases[x] for x in ("sample_id", "sample", "iid", "id", "材料", "样本") if x in aliases), None)
    if not id_key:
        raise VCFError("样本元数据必须包含sample_id/sample/IID/id列")
    group_key = next((aliases[x] for x in ("group", "variety", "population", "period4", "品种", "群体") if x in aliases), None)
    batch_key = next((aliases[x] for x in ("batch", "year", "lane", "批次", "年份") if x in aliases), None)
    rows = {}
    for row in reader:
        sample_id = str(row.get(id_key) or "").strip()
        if sample_id:
            rows[sample_id] = {
                "group": str(row.get(group_key) or "").strip() if group_key else "",
                "batch": str(row.get(batch_key) or "").strip() if batch_key else "",
            }
    return rows, {"path": str(path), "rows": len(rows), "id_column": id_key, "group_column": group_key, "batch_column": batch_key}


class QualityEvaluator:
    def __init__(self, service):
        self.service = service

    def _iter_lines(self, path):
        path = Path(path)
        if path.name.lower().endswith(".bcf"):
            if not self.service.bcftools:
                raise VCFError("BCF质量评估需要bcftools")
            proc = subprocess.Popen(
                [str(self.service.bcftools), "view", "-Ov", str(path)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                encoding="utf-8", errors="replace",
            )
            try:
                for line in proc.stdout:
                    yield line
                code = proc.wait()
                error = proc.stderr.read() if proc.stderr else ""
                if code:
                    raise VCFError((error or "bcftools读取BCF失败")[-3000:])
            finally:
                if proc.poll() is None:
                    proc.terminate()
            return
        compression = detect_compression(path)
        opener = gzip.open if compression in {"gzip", "bgzf"} else open
        with opener(str(path), "rt", encoding="utf-8", errors="replace") as handle:
            yield from handle

    @staticmethod
    def _smart_stride(metadata, path, target_records):
        total = metadata.get("record_count")
        if total:
            return max(1, int(math.ceil(total / max(target_records, 1))))
        size = Path(path).stat().st_size
        estimated_records = max(1, int(size / (90 if metadata.get("compressed") else 350)))
        return max(1, int(math.ceil(estimated_records / max(target_records, 1))))

    @staticmethod
    def _subgenome_label(chrom, profile):
        text = str(chrom)
        match = re.match(r"^(?:chr)?([AD])(?:t)?0*([0-9]+)$", text, re.IGNORECASE)
        if match:
            return match.group(1).upper()
        if str(profile.get("id") or "").startswith("cotton_") and text.isdigit():
            number = int(text)
            if 1 <= number <= 13:
                return "A（按1–13推断）"
            if 14 <= number <= 26:
                return "D（按14–26推断）"
        return "未分配"

    def evaluate(self, path_text, config, progress=None, cancel_event=None):
        metadata = self.service.inspect(path_text)
        path = metadata["path"]
        profile = resolve_profile(config.get("profile"))
        scan_mode = str(config.get("scan_mode") or "smart")
        if scan_mode not in {"smart", "full"}:
            raise VCFError("未知扫描模式")
        target_records = int(config.get("target_records") or 200000)
        target_records = max(1000, min(target_records, 1000000))
        stride = 1 if scan_mode == "full" else self._smart_stride(metadata, path, target_records)
        samples = list(metadata.get("samples") or [])
        sample_acc = {
            sample: {
                "sample_id": sample, "evaluated": 0, "called": 0, "missing": 0,
                "het": 0, "hom_ref": 0, "hom_alt": 0, "other": 0,
                "ploidy_counts": Counter(), "dp": Counter(), "gq": Counter(),
                "ab_n": 0, "ab_out": 0, "phased": 0,
            } for sample in samples
        }
        header = {
            "fileformat": None, "reference": None, "sources": [], "contigs": [],
            "filters": [], "info_ids": [], "format_ids": [], "has_chrom_header": False,
            "contig_lengths": {}, "info_descriptions": {},
        }
        type_counts = Counter()
        broad_type_counts = Counter()
        filter_counts = Counter()
        svtype_counts = Counter()
        sv_length_bins = Counter()
        sv_imprecise = 0
        sv_interval_uncertainty = 0
        sv_support_tagged = 0
        contig_counts = Counter()
        density_windows = Counter()
        region_overlap_count = 0
        annotation_tag_counts = Counter()
        consequence_counts = Counter()
        fake_het_windows = {}
        subgenome_acc = {}
        maf_bins = Counter()
        missing_bins = Counter()
        site_warn_count = 0
        site_critical_count = 0
        multiallelic = 0
        malformed = 0
        transitions = 0
        transversions = 0
        nonminimal_indels = 0
        adjacent_duplicates = 0
        previous_key = None
        metric_values = {key: [] for key in ("QUAL", "QD", "MQ", "FS", "SOR", "MQRankSum", "ReadPosRankSum")}
        metric_available_counts = Counter()
        metric_cap = max(10000, min(target_records, 100000))
        variant_export_cap = int(config.get("variant_export_cap") or 100000)
        variant_export_cap = max(1000, min(variant_export_cap, 200000))
        variant_rows = []
        sampled_records = 0
        gt_ploidy_counts = Counter()
        total_records = 0
        unsorted_records = 0
        last_chrom = None
        last_pos = -1
        seen_chroms = []
        seen_chrom_set = set()
        started = time.time()
        thresholds = profile["thresholds"]
        conditional_warnings = []
        reference_accessor = None
        reference_audit = {"requested": False, "status": "not_provided", "path": None, "checked_records": 0, "mismatch_records": 0, "missing_contig_records": 0, "contig_length_mismatches": []}
        reference_path = str(config.get("reference_path") or "").strip()
        if reference_path:
            reference_audit["requested"] = True
            reference_audit["path"] = str(Path(reference_path).expanduser().resolve())
            try:
                reference_accessor = IndexedFasta(reference_path)
                reference_audit["status"] = "complete"
                reference_audit["fai_contigs"] = len(reference_accessor.entries)
            except VCFError as exc:
                reference_audit["status"] = "unavailable"
                reference_audit["error"] = str(exc)
                conditional_warnings.append(_warning("warning", "reference", metadata["name"], "REFERENCE_AUDIT_UNAVAILABLE", "已提供参考基因组，但REF随机核验不可用", str(exc), "提供未压缩FASTA及同名.fai，或安装bcftools后重试"))
        sample_metadata = {}
        metadata_audit = {"requested": False, "status": "not_provided"}
        sample_meta_path = str(config.get("sample_meta_path") or "").strip()
        if sample_meta_path:
            sample_metadata, metadata_audit = load_sample_metadata(sample_meta_path)
            metadata_audit.update({"requested": True, "status": "complete"})
        bed_regions = None
        bed_audit = {"requested": False, "status": "not_provided", "path": None, "interval_count": 0}
        region_bed_path = str(config.get("region_bed_path") or "").strip()
        if region_bed_path:
            bed_regions = BedRegions(region_bed_path)
            bed_audit = {"requested": True, "status": "complete", "path": str(Path(region_bed_path).expanduser().resolve()), "interval_count": bed_regions.interval_count}

        for line in self._iter_lines(path):
            if cancel_event and cancel_event.is_set():
                raise VCFError("用户已取消质量评估")
            if line.startswith("##"):
                text = line.rstrip("\r\n")
                if text.startswith("##fileformat="):
                    header["fileformat"] = text.split("=", 1)[1]
                elif text.startswith("##reference="):
                    header["reference"] = text.split("=", 1)[1]
                elif text.startswith("##source="):
                    header["sources"].append(text.split("=", 1)[1])
                elif text.startswith("##contig=<ID="):
                    contig_id = text.split("##contig=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0]
                    header["contigs"].append(contig_id)
                    length_match = re.search(r"(?:^|,)length=([0-9]+)", text, re.IGNORECASE)
                    if length_match:
                        header["contig_lengths"][contig_id] = int(length_match.group(1))
                elif text.startswith("##FILTER=<ID="):
                    header["filters"].append(text.split("##FILTER=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
                elif text.startswith("##INFO=<ID="):
                    info_id = text.split("##INFO=<ID=", 1)[1].split(",", 1)[0]
                    header["info_ids"].append(info_id)
                    description = re.search(r'Description="([^"]*)"', text)
                    if description:
                        header["info_descriptions"][info_id] = description.group(1)
                elif text.startswith("##FORMAT=<ID="):
                    header["format_ids"].append(text.split("##FORMAT=<ID=", 1)[1].split(",", 1)[0])
                continue
            if line.startswith("#CHROM"):
                header["has_chrom_header"] = True
                continue
            if line.startswith("#"):
                continue
            total_records += 1
            # Avoid splitting hundreds/thousands of sample columns for records
            # which are not part of the genotype sample panel.
            parts = line.rstrip("\r\n").split("\t", 8)
            if len(parts) < 8:
                malformed += 1
                continue
            chrom = parts[0]
            pos = _safe_int(parts[1])
            if pos is None:
                malformed += 1
                continue
            if chrom == last_chrom and pos < last_pos:
                unsorted_records += 1
            elif chrom != last_chrom:
                if chrom in seen_chrom_set:
                    unsorted_records += 1
                else:
                    seen_chroms.append(chrom)
                    seen_chrom_set.add(chrom)
            last_chrom, last_pos = chrom, pos
            ref, alt, qual_text, filt, info_text = parts[3], parts[4], parts[5], parts[6], parts[7]
            info = parse_info(info_text)
            broad_type = classify_variant(ref, alt, "", info_text)
            variant_type = _detail_variant_type(ref, alt, info, broad_type)
            type_counts[variant_type] += 1
            broad_type_counts[broad_type] += 1
            contig_counts[chrom] += 1
            density_windows[(chrom, (pos - 1) // 1000000)] += 1
            if bed_regions and bed_regions.contains(chrom, pos):
                region_overlap_count += 1
            record_key = (chrom, pos, ref, alt)
            if record_key == previous_key:
                adjacent_duplicates += 1
            previous_key = record_key
            if "," in alt:
                multiallelic += 1
            if broad_type == "INDEL":
                alleles = alt.split(",")
                reducible = False
                for value in alleles:
                    if not value or value.startswith("<"):
                        continue
                    left, right = ref, value
                    while len(left) > 1 and len(right) > 1 and left[-1] == right[-1]:
                        left, right, reducible = left[:-1], right[:-1], True
                    while len(left) > 1 and len(right) > 1 and left[0] == right[0]:
                        left, right, reducible = left[1:], right[1:], True
                if reducible:
                    nonminimal_indels += 1
            for value in filt.split(";"):
                filter_counts[value or "."] += 1
            if broad_type == "SV":
                svtype_counts[info.get("SVTYPE") or variant_type or "未标注"] += 1
                length = _sv_length(info, pos)
                length_bucket = "未知" if length is None else "<100 bp" if length < 100 else "100 bp–1 kb" if length < 1000 else "1–10 kb" if length < 10000 else "10–100 kb" if length < 100000 else ">=100 kb"
                sv_length_bins[length_bucket] += 1
                sv_imprecise += int("IMPRECISE" in info_text.split(";"))
                sv_interval_uncertainty += int("CIPOS" in info or "CIEND" in info)
                sv_support_tagged += int(any(key in info for key in ("SU", "SUPPORT", "RE", "PE", "SR", "DV", "RV")))
            for tag in ("ANN", "CSQ", "BCSQ"):
                if info.get(tag):
                    annotation_tag_counts[tag] += 1
                    entries = str(info[tag]).split(",")
                    for entry in entries[:20]:
                        fields = entry.split("|")
                        if tag == "ANN" and len(fields) > 1:
                            labels = fields[1].split("&")
                        elif tag == "CSQ":
                            description = header["info_descriptions"].get("CSQ", "")
                            format_match = re.search(r"Format:\s*([^\"]+)", description)
                            names = [x.strip() for x in format_match.group(1).split("|")] if format_match else []
                            consequence_index = names.index("Consequence") if "Consequence" in names else 1
                            labels = fields[consequence_index].split("&") if len(fields) > consequence_index else []
                        else:
                            labels = fields[0].split("&") if fields else []
                        for label in labels:
                            if label.strip():
                                consequence_counts[label.strip()] += 1

            sampled = stride == 1 or ((sum(ord(x) for x in chrom) * 131 + pos) % stride == 0)
            if sampled:
                sampled_records += 1
                current_metrics = {"QUAL": _safe_float(qual_text)}
                for metric in ("QD", "MQ", "FS", "SOR", "MQRankSum", "ReadPosRankSum"):
                    current_metrics[metric] = _safe_float(str(info.get(metric) or "").split(",", 1)[0])
                for metric, value in current_metrics.items():
                    if value is not None:
                        metric_available_counts[metric] += 1
                        if len(metric_values[metric]) < metric_cap:
                            metric_values[metric].append(value)
                if broad_type == "SNP" and len(ref) == 1 and len(alt) == 1:
                    pair = (ref.upper(), alt.upper())
                    if pair in {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}:
                        transitions += 1
                    else:
                        transversions += 1
                sample_tail = parts[8].split("\t") if len(parts) > 8 else []
                format_keys = sample_tail[0].split(":") if sample_tail else []
                index = {key: i for i, key in enumerate(format_keys)}
                missing = 0
                missing_rate = None
                maf = None
                allele_count = Counter()
                allele_number = 0
                window_key = (chrom, (pos - 1) // 1000000)
                window_acc = fake_het_windows.setdefault(window_key, {"sites": 0, "called": 0, "het": 0, "ab_n": 0, "quarter_ab": 0, "dp_sum": 0, "dp_n": 0})
                window_acc["sites"] += 1
                subgenome = self._subgenome_label(chrom, profile)
                sub_acc = subgenome_acc.setdefault(subgenome, {"sites": 0, "called": 0, "het": 0})
                sub_acc["sites"] += 1
                for sample_index, sample in enumerate(samples):
                    field = sample_tail[1 + sample_index] if 1 + sample_index < len(sample_tail) else "."
                    values = field.split(":")
                    gt = values[index["GT"]] if "GT" in index and index["GT"] < len(values) else "."
                    alleles = parse_gt(gt)
                    acc = sample_acc[sample]
                    acc["evaluated"] += 1
                    category = genotype_class(alleles)
                    acc[category] += 1
                    if alleles is None:
                        missing += 1
                        continue
                    acc["called"] += 1
                    window_acc["called"] += 1
                    sub_acc["called"] += 1
                    if category == "het":
                        window_acc["het"] += 1
                        sub_acc["het"] += 1
                    if "|" in gt:
                        acc["phased"] += 1
                    acc["ploidy_counts"][len(alleles)] += 1
                    gt_ploidy_counts[len(alleles)] += 1
                    for allele in alleles:
                        allele_count[allele] += 1
                        allele_number += 1
                    if "DP" in index and index["DP"] < len(values):
                        value = _safe_int(values[index["DP"]])
                        if value is not None and value >= 0:
                            acc["dp"][min(value, 10000)] += 1
                            window_acc["dp_sum"] += value
                            window_acc["dp_n"] += 1
                    if "GQ" in index and index["GQ"] < len(values):
                        value = _safe_int(values[index["GQ"]])
                        if value is not None and value >= 0:
                            acc["gq"][min(value, 999)] += 1
                    if category == "het" and "AD" in index and index["AD"] < len(values):
                        ads = [_safe_float(x) for x in values[index["AD"]].split(",")]
                        if len(ads) >= 2 and ads[0] is not None and ads[1] is not None and ads[0] + ads[1] > 0:
                            ab = ads[1] / (ads[0] + ads[1])
                            acc["ab_n"] += 1
                            window_acc["ab_n"] += 1
                            if min(abs(ab - .25), abs(ab - .75)) <= .10:
                                window_acc["quarter_ab"] += 1
                            if ab < thresholds["ab_lower"] or ab > thresholds["ab_upper"]:
                                acc["ab_out"] += 1
                if samples:
                    missing_rate = missing / len(samples)
                    if missing_rate >= thresholds["site_missing_critical"]:
                        site_critical_count += 1
                    elif missing_rate >= thresholds["site_missing_warn"]:
                        site_warn_count += 1
                    bucket = "0–1%" if missing_rate <= .01 else "1–5%" if missing_rate <= .05 else "5–10%" if missing_rate <= .10 else "10–20%" if missing_rate <= .20 else ">20%"
                    missing_bins[bucket] += 1
                if allele_number and len(allele_count) >= 1:
                    ref_freq = allele_count.get(0, 0) / allele_number
                    alt_freq = max(0.0, 1.0 - ref_freq)
                    maf = min(ref_freq, alt_freq)
                    bucket = "0" if maf == 0 else "(0,1%]" if maf <= .01 else "(1%,5%]" if maf <= .05 else "(5%,10%]" if maf <= .10 else "(10%,30%]" if maf <= .30 else "(30%,50%]"
                    maf_bins[bucket] += 1
                if len(variant_rows) < variant_export_cap:
                    variant_rows.append({
                        "chrom": chrom, "pos": pos, "id": parts[2], "ref": ref, "alt": alt,
                        "type": variant_type, "qual": _safe_float(qual_text), "filter": filt,
                        "missing_rate": missing_rate, "maf": maf,
                        "qd": _safe_float(info.get("QD")), "mq": _safe_float(info.get("MQ")),
                        "fs": _safe_float(info.get("FS")), "sor": _safe_float(info.get("SOR")),
                        "mq_rank_sum": _safe_float(info.get("MQRankSum")),
                        "read_pos_rank_sum": _safe_float(info.get("ReadPosRankSum")),
                        "svtype": info.get("SVTYPE"), "svlen": _sv_length(info, pos) if broad_type == "SV" else None,
                        "imprecise": "IMPRECISE" in info_text.split(";"),
                    })
                if reference_accessor and ref and len(ref) <= 1000 and re.fullmatch(r"[ACGTNacgtn]+", ref):
                    observed = reference_accessor.fetch(chrom, pos, len(ref))
                    if observed is None:
                        reference_audit["missing_contig_records"] += 1
                    else:
                        reference_audit["checked_records"] += 1
                        if observed != ref.upper():
                            reference_audit["mismatch_records"] += 1

            if progress and total_records % 5000 == 0:
                expected = metadata.get("record_count")
                pct = min(99.0, total_records * 100 / expected) if expected else None
                progress(total_records, pct, "正在流式扫描VCF")

        self.service._remember_record_count(path, total_records)
        if reference_accessor:
            for chrom, length in header["contig_lengths"].items():
                ref_entry = reference_accessor.entries.get(chrom)
                if ref_entry and ref_entry[0] != length:
                    reference_audit["contig_length_mismatches"].append({"chrom": chrom, "vcf_length": length, "fasta_length": ref_entry[0]})
            reference_accessor.close()
        configured_gt_ploidy = profile.get("genotype_ploidy", "auto")
        detected_gt_ploidy = gt_ploidy_counts.most_common(1)[0][0] if gt_ploidy_counts else None
        effective_gt_ploidy = detected_gt_ploidy if configured_gt_ploidy == "auto" else configured_gt_ploidy
        profile["detected_gt_ploidy"] = detected_gt_ploidy
        profile["effective_gt_ploidy"] = effective_gt_ploidy
        profile["interpretation"]["hwe_enabled"] = bool(
            thresholds.get("hwe_p_warn") is not None and effective_gt_ploidy == 2
        )
        profile["interpretation"]["diploid_ab_model"] = effective_gt_ploidy == 2
        report_samples = []
        for sample, acc in sample_acc.items():
            evaluated = acc["evaluated"]
            called = acc["called"]
            report_samples.append({
                "sample_id": sample,
                "evaluated_records": evaluated,
                "call_rate": called / evaluated if evaluated else None,
                "missing_rate": acc["missing"] / evaluated if evaluated else None,
                "het_rate": acc["het"] / called if called else None,
                "hom_ref_rate": acc["hom_ref"] / called if called else None,
                "hom_alt_rate": acc["hom_alt"] / called if called else None,
                "median_dp": _median_counter(acc["dp"]),
                "median_gq": _median_counter(acc["gq"]),
                "ab_outlier_rate": acc["ab_out"] / acc["ab_n"] if acc["ab_n"] else None,
                "phase_rate": acc["phased"] / called if called else None,
                "ploidy_mismatch_rate": (
                    sum(count for value, count in acc["ploidy_counts"].items() if value != effective_gt_ploidy) / called
                    if called and effective_gt_ploidy else None
                ),
            })

        matched_meta = 0
        for item in report_samples:
            meta = sample_metadata.get(item["sample_id"], {})
            item["group"] = meta.get("group") or ""
            item["batch"] = meta.get("batch") or ""
            matched_meta += int(bool(meta))
        metadata_audit["matched_vcf_samples"] = matched_meta
        metadata_audit["unmatched_vcf_samples"] = len(samples) - matched_meta if sample_metadata else None
        metadata_audit["metadata_only_samples"] = len(set(sample_metadata) - set(samples)) if sample_metadata else None

        group_summaries = []
        for field in ("group", "batch"):
            labels = sorted({x[field] for x in report_samples if x.get(field)})
            for label in labels:
                rows = [x for x in report_samples if x.get(field) == label]
                group_summaries.append({
                    "dimension": field, "label": label, "sample_count": len(rows),
                    "median_missing_rate": _median([x["missing_rate"] for x in rows]),
                    "median_het_rate": _median([x["het_rate"] for x in rows]),
                    "median_dp": _median([x["median_dp"] for x in rows]),
                    "median_gq": _median([x["median_gq"] for x in rows]),
                })

        subgenome_summary = []
        for label, values in subgenome_acc.items():
            subgenome_summary.append({
                "subgenome": label, "evaluated_sites": values["sites"], "called_genotypes": values["called"],
                "heterozygous_genotypes": values["het"],
                "heterozygosity_rate": values["het"] / values["called"] if values["called"] else None,
            })
        subgenome_summary.sort(key=lambda x: x["subgenome"])

        window_rows = []
        for (chrom, window), values in fake_het_windows.items():
            window_rows.append({
                "chrom": chrom, "start": window * 1000000 + 1, "end": (window + 1) * 1000000,
                "evaluated_sites": values["sites"],
                "het_rate": values["het"] / values["called"] if values["called"] else None,
                "quarter_ab_rate": values["quarter_ab"] / values["ab_n"] if values["ab_n"] else None,
                "mean_dp": values["dp_sum"] / values["dp_n"] if values["dp_n"] else None,
            })
        eligible_windows = [x for x in window_rows if x["evaluated_sites"] >= 3 and x["het_rate"] is not None]
        window_het_q99 = _quantile([x["het_rate"] for x in eligible_windows], .99)
        window_dp_q95 = _quantile([x["mean_dp"] for x in eligible_windows], .95)
        for item in window_rows:
            item["fake_het_candidate"] = bool(
                profile.get("subgenomes", 1) > 1 and item["evaluated_sites"] >= 3 and
                window_het_q99 is not None and item["het_rate"] is not None and item["het_rate"] >= window_het_q99 and
                item["quarter_ab_rate"] is not None and item["quarter_ab_rate"] >= .30 and
                (window_dp_q95 is None or item["mean_dp"] is None or item["mean_dp"] >= window_dp_q95)
            )
        window_rows.sort(key=lambda x: (not x["fake_het_candidate"], -(x["het_rate"] or 0)))

        warnings = list(conditional_warnings)
        if not header["fileformat"] or not header["has_chrom_header"]:
            warnings.append(_warning("critical", "file", metadata["name"], "HEADER_REQUIRED", "VCF关键Header不完整", advice="修复Header后再进入下游分析"))
        if not header["reference"]:
            warnings.append(_warning("warning", "file", metadata["name"], "REFERENCE_UNKNOWN", "未声明参考基因组", advice="补充参考版本后再进行功能或数据库注释"))
        if not header["sources"]:
            warnings.append(_warning("info", "file", metadata["name"], "SOURCE_UNKNOWN", "未在Header中识别到来源程序"))
        if not header["filters"] or set(filter_counts) <= {".", "PASS"}:
            warnings.append(_warning("warning", "file", metadata["name"], "QC_EVIDENCE_WEAK", "过滤与质控痕迹较弱", evidence="FILTER定义或记录过滤状态不足", advice="核对上游calling和过滤流程"))
        if unsorted_records:
            warnings.append(_warning("critical", "file", metadata["name"], "UNSORTED", "检测到坐标顺序异常", evidence="{}条记录发生倒序或染色体回跳".format(unsorted_records), advice="生成新文件进行排序，保留原文件"))
        if malformed:
            warnings.append(_warning("critical", "file", metadata["name"], "MALFORMED", "存在无法解析的记录", evidence="{}条".format(malformed)))
        if nonminimal_indels:
            warnings.append(_warning(
                "warning", "site", metadata["name"], "NONMINIMAL_ALLELES",
                "检测到可继续最简化表达的INDEL", "{}条记录".format(nonminimal_indels),
                "用参考FASTA生成bcftools norm新副本；不要覆盖原文件",
            ))
        if adjacent_duplicates:
            warnings.append(_warning(
                "warning", "site", metadata["name"], "ADJACENT_DUPLICATES",
                "检测到相邻重复变异记录", "{}条CHROM/POS/REF/ALT重复".format(adjacent_duplicates),
                "标准化后再次审计，再决定是否去重",
            ))
        if reference_audit["status"] == "complete":
            mismatch_count = reference_audit["mismatch_records"]
            checked_count = reference_audit["checked_records"]
            if reference_audit["contig_length_mismatches"]:
                warnings.append(_warning("critical", "reference", metadata["name"], "CONTIG_LENGTH_MISMATCH", "VCF Header与FASTA的contig长度不一致", "{}个contig".format(len(reference_audit["contig_length_mismatches"])), "确认参考版本后再做标准化或功能注释"))
            if mismatch_count:
                rate = mismatch_count / checked_count if checked_count else None
                level = "critical" if mismatch_count >= 10 or (rate is not None and rate >= .001) else "warning"
                warnings.append(_warning(level, "reference", metadata["name"], "REF_MISMATCH", "VCF REF与参考FASTA不一致", "{}/{}（{}）".format(mismatch_count, checked_count, _ratio(rate)), "停止自动标准化，先核对build、染色体命名与FASTA版本"))
        if sample_metadata and matched_meta < len(samples):
            rate = (len(samples) - matched_meta) / max(len(samples), 1)
            warnings.append(_warning("warning" if rate > .05 else "info", "metadata", metadata["name"], "SAMPLE_METADATA_UNMATCHED", "部分VCF样本未匹配到元数据", "{}/{}（{}）".format(len(samples) - matched_meta, len(samples), _ratio(rate)), "检查样本ID的大小写、前后缀和隐藏空格"))
        batch_rows = [x for x in group_summaries if x["dimension"] == "batch" and x["sample_count"] >= 3]
        if len(batch_rows) >= 2:
            missing_values = [x["median_missing_rate"] for x in batch_rows if x["median_missing_rate"] is not None]
            dp_values = [x["median_dp"] for x in batch_rows if x["median_dp"] is not None and x["median_dp"] > 0]
            if missing_values and max(missing_values) - min(missing_values) >= .05:
                warnings.append(_warning("warning", "batch", metadata["name"], "BATCH_MISSING_SHIFT", "不同批次的中位缺失率差异较大", "范围{}–{}".format(_ratio(min(missing_values)), _ratio(max(missing_values))), "结合测序批次、建库与过滤流程复核"))
            if len(dp_values) >= 2 and max(dp_values) / min(dp_values) >= 2:
                warnings.append(_warning("warning", "batch", metadata["name"], "BATCH_DEPTH_SHIFT", "不同批次的中位深度相差两倍以上", "范围{}–{}".format(_fmt(min(dp_values)), _fmt(max(dp_values))), "检查批次覆盖度与统一过滤阈值"))
        fake_candidates = [x for x in window_rows if x["fake_het_candidate"]]
        if fake_candidates:
            warnings.append(_warning("warning", "region", metadata["name"], "POLYPLOID_FAKE_HET_WINDOWS", "检测到高深度且AB偏向0.25/0.75的假杂合候选窗口", "{}个1 Mb窗口".format(len(fake_candidates)), "优先检查homeologous错配、collapsed repeats和亚基因组注释"))

        qual_reference = thresholds.get("site_qual_warn")
        qual_values = metric_values.get("QUAL") or []
        qual_coverage = metric_available_counts["QUAL"] / sampled_records if sampled_records else 0
        if qual_reference is not None and qual_coverage >= .80 and qual_values:
            low_qual_fraction = sum(value < qual_reference for value in qual_values) / len(qual_values)
            if low_qual_fraction >= .10:
                warnings.append(_warning(
                    "warning", "site", metadata["name"], "LOW_QUAL_FRACTION",
                    "较多位点低于所选作物的QUAL建议起始线",
                    "抽样中{:.2%}低于QUAL {}（覆盖率{:.2%}）".format(low_qual_fraction, _fmt(qual_reference), qual_coverage),
                    "QUAL依赖caller；请结合FILTER、QD/MQ/FS/SOR和上游流程复核后再过滤",
                ))
        titv_reference = thresholds.get("expected_titv_min")
        titv_observed = transitions / transversions if transversions else None
        titv_snp_n = transitions + transversions
        if titv_reference is not None and titv_observed is not None and titv_snp_n >= 1000 and titv_observed < titv_reference:
            warnings.append(_warning(
                "warning", "site", metadata["name"], "TITV_BELOW_CROP_REFERENCE",
                "Ti/Tv低于所选作物的大规模SNP集经验参考线",
                "Ti/Tv={:.3f}，参考线={}，抽样双等位SNP={}条".format(titv_observed, _fmt(titv_reference), titv_snp_n),
                "该指标不适用于SV/INDEL或小位点集；请检查参考版本、变异过滤和测序错误，不要单独据此删位点",
            ))

        het_values = [x["het_rate"] for x in report_samples if x["het_rate"] is not None]
        het_center = _median(het_values)
        het_mad = _mad(het_values, het_center)
        het_q01 = _quantile(het_values, 0.01)
        het_q99 = _quantile(het_values, 0.99)
        robust_scale = 1.4826 * het_mad if het_mad else None
        relative_upper = None
        relative_critical = None
        if len(het_values) >= 20 and het_center is not None:
            if robust_scale:
                relative_upper = max(het_q99, het_center + 5 * robust_scale)
                relative_critical = max(_quantile(het_values, 0.995), het_center + 8 * robust_scale)
            elif het_q99 is not None and het_q99 > het_center:
                relative_upper = het_q99
                relative_critical = _quantile(het_values, 0.995)
        variant_classes = set(broad_type_counts)
        sv_only = bool(variant_classes) and variant_classes == {"SV"}
        het_rule_mode = "cohort_relative_sv" if sv_only else "profile_absolute_plus_relative"

        if sv_only and thresholds.get("het_rate_warn") is not None:
            warnings.append(_warning(
                "info", "cohort", metadata["name"], "SV_HET_RELATIVE_ONLY",
                "SV文件不使用SNP/INDEL的绝对杂合率阈值逐样本报警",
                "队列中位杂合率{}；改用队列相对极端值".format(_ratio(het_center)),
                "结合SV caller、变异类型和群体分组解释",
            ))
        elif thresholds.get("het_rate_warn") is not None and het_center is not None:
            if thresholds.get("het_rate_critical") is not None and het_center >= thresholds["het_rate_critical"]:
                warnings.append(_warning(
                    "warning", "cohort", metadata["name"], "COHORT_HIGH_HET",
                    "队列整体杂合率超过当前Profile的严重参考线", _ratio(het_center),
                    "先核对材料类型、变异过滤和亚基因组错配；不会因此给全部样本重复报警",
                ))
            elif het_center >= thresholds["het_rate_warn"]:
                warnings.append(_warning(
                    "warning", "cohort", metadata["name"], "COHORT_HIGH_HET",
                    "队列整体杂合率超过当前Profile提醒线", _ratio(het_center),
                    "结合材料世代和变异类型解释；逐样本仅报告相对极端值",
                ))
        for item in report_samples:
            sid = item["sample_id"]
            missing = item["missing_rate"]
            if missing is not None and missing >= thresholds["sample_missing_critical"]:
                warnings.append(_warning("critical", "sample", sid, "SAMPLE_MISSING", "样本缺失率过高", _ratio(missing), "优先复核或排除该样本"))
            elif missing is not None and missing >= thresholds["sample_missing_warn"]:
                warnings.append(_warning("warning", "sample", sid, "SAMPLE_MISSING", "样本缺失率偏高", _ratio(missing), "检查测序深度、批次与过滤条件"))
            gq = item["median_gq"]
            if gq is not None and gq < thresholds["median_gq_critical"]:
                warnings.append(_warning("critical", "sample", sid, "LOW_GQ", "样本中位GQ过低", str(gq)))
            elif gq is not None and gq < thresholds["median_gq_warn"]:
                warnings.append(_warning("warning", "sample", sid, "LOW_GQ", "样本中位GQ偏低", str(gq)))
            dp = item["median_dp"]
            if dp is not None and dp < thresholds["median_dp_min_warn"]:
                warnings.append(_warning("warning", "sample", sid, "LOW_DP", "样本中位DP低于当前阈值", str(dp)))
            elif dp is not None and dp > thresholds["median_dp_max_warn"]:
                warnings.append(_warning("warning", "sample", sid, "HIGH_DP", "样本中位DP高于当前阈值", str(dp), "检查重复区、拷贝数和比对偏倚"))
            het = item["het_rate"]
            if het is not None and relative_upper is not None and het > relative_upper:
                if relative_critical is not None and het > relative_critical:
                    message = "杂合率为队列中的极端高值"
                else:
                    message = "杂合率为队列中的相对高值"
                corroborated = (
                    missing is not None and missing >= thresholds["sample_missing_critical"]
                ) or (
                    gq is not None and gq < thresholds["median_gq_critical"]
                )
                level = "critical" if corroborated else "warning"
                warnings.append(_warning(
                    level, "sample", sid, "HET_OUTLIER", message, _ratio(het),
                    "单一杂合率离群不判定样本失败；请结合品种分组、批次、缺失率和GQ复核",
                ))
            mismatch = item["ploidy_mismatch_rate"]
            if mismatch is not None and mismatch > .05:
                warnings.append(_warning("warning", "sample", sid, "PLOIDY_MISMATCH", "GT倍性与所选Profile不一致", _ratio(mismatch), "确认VCF GT编码与倍性设置"))

        sample_levels = {}
        for warning in warnings:
            if warning["scope"] == "sample":
                old = sample_levels.get(warning["target"], "info")
                if _level_rank(warning["level"]) > _level_rank(old):
                    sample_levels[warning["target"]] = warning["level"]
        score_weights = {"callrate": .25, "depth": .15, "gq": .15, "ab": .10, "het": .15, "ploidy": .10}
        for item in report_samples:
            components = {}
            if item["call_rate"] is not None:
                components["callrate"] = max(0, min(1, item["call_rate"]))
            if item["median_dp"] is not None:
                dp = item["median_dp"]
                if dp < thresholds["median_dp_min_warn"]:
                    components["depth"] = max(0, dp / max(thresholds["median_dp_min_warn"], 1e-9))
                elif dp > thresholds["median_dp_max_warn"]:
                    components["depth"] = max(0, thresholds["median_dp_max_warn"] / max(dp, 1e-9))
                else:
                    components["depth"] = 1.0
            if item["median_gq"] is not None:
                components["gq"] = max(0, min(1, item["median_gq"] / max(thresholds["median_gq_warn"], 1e-9)))
            if item["ab_outlier_rate"] is not None:
                components["ab"] = max(0, 1 - item["ab_outlier_rate"])
            if item["ploidy_mismatch_rate"] is not None:
                components["ploidy"] = max(0, 1 - item["ploidy_mismatch_rate"])
            het = item["het_rate"]
            if het is not None:
                warn_line = relative_upper
                critical_line = relative_critical
                if profile.get("mating_system") == "selfing" and not sv_only and thresholds.get("het_rate_warn") is not None:
                    warn_line = thresholds["het_rate_warn"]
                    critical_line = thresholds.get("het_rate_critical") or max(warn_line * 2, warn_line + .01)
                if warn_line is None or het <= warn_line:
                    components["het"] = 1.0
                elif critical_line and critical_line > warn_line:
                    components["het"] = max(0, 1 - (het - warn_line) / (critical_line - warn_line))
                else:
                    components["het"] = .5
            denominator = sum(score_weights[key] for key in components)
            item["score_components"] = {key: round(value * 100, 2) for key, value in components.items()}
            item["sample_score"] = round(100 * sum(score_weights[key] * value for key, value in components.items()) / denominator) if denominator else None
            level = sample_levels.get(item["sample_id"])
            item["sample_status"] = "critical" if level == "critical" else "warning" if level == "warning" else "pass"

        warnings.sort(key=lambda x: (-_level_rank(x["level"]), x["scope"], x["target"], x["code"]))
        critical_n = sum(x["level"] == "critical" for x in warnings)
        warning_n = sum(x["level"] == "warning" for x in warnings)
        file_critical = sum(x["level"] == "critical" and x["scope"] == "file" for x in warnings)
        file_warning = sum(x["level"] == "warning" and x["scope"] == "file" for x in warnings)
        critical_samples = len({x["target"] for x in warnings if x["level"] == "critical" and x["scope"] == "sample"})
        warning_samples = len({x["target"] for x in warnings if x["level"] == "warning" and x["scope"] == "sample"})
        sample_denominator = max(len(samples), 1)
        site_denominator = max(sampled_records, 1)
        score = 100.0
        score -= min(30.0, file_critical * 12.0)
        score -= min(12.0, file_warning * 4.0)
        score -= 30.0 * critical_samples / sample_denominator
        score -= 15.0 * warning_samples / sample_denominator
        score -= 20.0 * site_critical_count / site_denominator
        score -= 8.0 * site_warn_count / site_denominator
        score = max(0, round(score))
        if any(x["code"] in {"HEADER_REQUIRED", "UNSORTED", "MALFORMED"} and x["level"] == "critical" for x in warnings):
            score = min(score, 59)
        status = "critical" if critical_n else "warning" if warning_n else "pass"
        evidence_items = {
            "fileformat": bool(header["fileformat"]), "chrom_header": header["has_chrom_header"],
            "reference": bool(header["reference"]), "source": bool(header["sources"]),
            "filter_definition": bool(header["filters"]), "format_gt": "GT" in header["format_ids"],
            "info_definitions": bool(header["info_ids"]), "sorted": not unsorted_records,
        }
        qc_evidence_score = round(100 * sum(evidence_items.values()) / len(evidence_items))
        quality_summaries = {key: _metric_summary(values, sampled_records, metric_available_counts[key]) for key, values in metric_values.items()}
        density_all = []
        for (chrom, window), count in density_windows.most_common():
            density_all.append({"chrom": chrom, "start": window * 1000000 + 1, "end": (window + 1) * 1000000, "variant_count": count})
        sv_total = sum(svtype_counts.values())
        phase_values = [x["phase_rate"] for x in report_samples if x.get("phase_rate") is not None]
        sample_scores = [x["sample_score"] for x in report_samples if x.get("sample_score") is not None]
        result = {
            "schema_version": "callvcf-qc-1.2",
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "input": {
                "path": path, "name": metadata["name"], "file_size": metadata["file_size"],
                "storage": metadata["storage"], "backend": metadata.get("backend"),
            },
            "profile": profile,
            "scan": {
                "mode": scan_mode, "stride": stride, "total_records": total_records,
                "evaluated_records": sampled_records,
                "sampling_fraction": sampled_records / total_records if total_records else None,
                "elapsed_seconds": round(time.time() - started, 3),
                "method_note": "所有记录均用于总数、类型和FILTER统计；样本/基因型指标按跨全基因组确定性抽样计算。" if stride > 1 else "全部记录均进入位点与样本指标计算。",
            },
            "summary": {
                "score": score, "status": status, "sample_count": len(samples),
                "record_count": total_records, "contig_count": len(contig_counts),
                "critical_count": critical_n, "warning_count": warning_n,
                "median_sample_score": _median(sample_scores),
            },
            "header_audit": {**header, "qc_evidence_items": evidence_items, "qc_evidence_score": qc_evidence_score},
            "site_metrics": {
                "variant_types": dict(type_counts), "broad_variant_types": dict(broad_type_counts), "filters": dict(filter_counts),
                "svtypes": dict(svtype_counts), "contigs": dict(contig_counts),
                "multiallelic_count": multiallelic, "malformed_count": malformed,
                "unsorted_count": unsorted_records, "site_warning_count": site_warn_count,
                "site_critical_count": site_critical_count,
                "missing_rate_bins": dict(missing_bins), "maf_bins": dict(maf_bins),
                "transitions": transitions, "transversions": transversions,
                "titv": transitions / transversions if transversions else None,
                "quality_field_summaries": quality_summaries,
                "nonminimal_indel_count": nonminimal_indels,
                "adjacent_duplicate_count": adjacent_duplicates,
                "density_window_bp": 1000000,
                "density_hotspots": density_all[:100], "density_windows": density_all,
                "sv_length_bins": dict(sv_length_bins),
                "sv_imprecise_rate": sv_imprecise / sv_total if sv_total else None,
                "sv_interval_uncertainty_rate": sv_interval_uncertainty / sv_total if sv_total else None,
                "sv_support_tag_coverage": sv_support_tagged / sv_total if sv_total else None,
            },
            "cohort_metrics": {
                "median_het_rate": het_center, "het_mad": het_mad,
                "het_q01": het_q01, "het_q99": het_q99,
                "het_relative_upper": relative_upper,
                "het_relative_critical": relative_critical,
                "het_rule_mode": het_rule_mode,
                "median_phase_rate": _median(phase_values),
                "group_summaries": group_summaries,
                "subgenome_summary": subgenome_summary,
            },
            "reference_audit": reference_audit,
            "metadata_audit": metadata_audit,
            "region_audit": {**bed_audit, "overlap_records": region_overlap_count, "overlap_rate": region_overlap_count / total_records if total_records else None},
            "annotation_audit": {"tag_record_counts": dict(annotation_tag_counts), "consequence_counts": dict(consequence_counts), "annotated_record_count": sum(annotation_tag_counts.values())},
            "fake_heterozygosity_windows": window_rows,
            "samples": report_samples,
            "module_availability": [
                {"module": "VCF/Header/位点/样本指标", "status": "complete", "reason": "VCF内生数据"},
                {"module": "参考REF与contig长度核验", "status": reference_audit["status"] if config.get("reference_path") else "conditional", "reason": reference_audit.get("error") or "同版本FASTA及.fai"},
                {"module": "批次/品种内比较", "status": metadata_audit["status"] if config.get("sample_meta_path") else "conditional", "reason": "已匹配{}个VCF样本".format(matched_meta) if sample_metadata else "需要含sample_id、group/variety、batch的样本表"},
                {"module": "重复区/低复杂度重叠", "status": bed_audit["status"] if config.get("region_bed_path") else "conditional", "reason": "已载入{}个合并区间".format(bed_audit["interval_count"]) if bed_regions else "需要BED注释"},
                {"module": "功能后果与结构域", "status": "conditional", "reason": "需要ANN/CSQ/BCSQ或GFF3/外部注释表"},
                {"module": "污染定量", "status": "not_applicable", "reason": "VCF只能给代理信号；可靠估计需要BAM/CRAM与专用模型"},
            ],
            "_variant_metrics_rows": variant_rows,
            "variant_export": {"rows": len(variant_rows), "cap": variant_export_cap, "truncated": sampled_records > len(variant_rows)},
            "warnings": warnings,
            "repair_policy": {
                "safe_actions": ["建立缺失索引", "生成新的排序副本"],
                "dangerous_actions": ["自适应/自定义组合质控", "标准化/拆分", "补全INFO统计标签", "去重", "样本子集", "按条件掩蔽GT", "位点过滤", "覆盖原VCF", "改写REF/ALT", "坐标转换", "染色体批量重命名", "删除原文件"],
                "dangerous_requires_second_confirmation": True,
                "original_is_never_silently_overwritten": True,
            },
        }
        result["qc_recommendation"] = build_qc_recommendation(result)
        result["analysis_readiness"] = build_analysis_readiness(result)
        if progress:
            progress(total_records, 100.0, "质量评估完成，正在生成报告")
        return result


def _svg_bars(items, title, color="#2f725f", width=760, height=250):
    items = [(str(k), int(v)) for k, v in items if int(v) >= 0]
    if not items:
        return '<div class="empty">无可绘制数据</div>'
    margin_l, margin_r, margin_t, margin_b = 60, 20, 38, 58
    plot_w, plot_h = width - margin_l - margin_r, height - margin_t - margin_b
    max_value = max(v for _, v in items) or 1
    gap = 8
    bar_w = max(8, (plot_w - gap * (len(items) + 1)) / max(len(items), 1))
    bars = []
    for i, (label, value) in enumerate(items):
        x = margin_l + gap + i * (bar_w + gap)
        h = plot_h * value / max_value
        y = margin_t + plot_h - h
        bars.append('<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" rx="4" fill="{}"><title>{}: {:,}</title></rect>'.format(x, y, bar_w, h, color, html.escape(label), value))
        bars.append('<text x="{:.1f}" y="{}" text-anchor="middle" font-size="11" fill="#52645d">{}</text>'.format(x + bar_w / 2, height - 32, html.escape(label[:14])))
        bars.append('<text x="{:.1f}" y="{:.1f}" text-anchor="middle" font-size="10" fill="#18342b">{:,}</text>'.format(x + bar_w / 2, max(28, y - 5), value))
    return '<svg viewBox="0 0 {} {}" role="img" aria-label="{}"><text x="20" y="23" font-size="15" font-weight="700" fill="#18342b">{}</text><line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#b9c8c1"/>{}</svg>'.format(width, height, html.escape(title), html.escape(title), margin_l, margin_t + plot_h, width - margin_r, margin_t + plot_h, "".join(bars))


def _svg_pca(rows, width=760, height=360):
    points = [(x.get("PC1"), x.get("PC2"), str(x.get("sample_id") or "")) for x in rows]
    points = [(float(x), float(y), label) for x, y, label in points if x is not None and y is not None]
    if not points:
        return '<div class="empty">PCA没有足够的PC1/PC2结果</div>'
    xs, ys = [x[0] for x in points], [x[1] for x in points]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
    xspan, yspan = xmax - xmin or 1, ymax - ymin or 1
    circles = []
    for x, y, label in points:
        px = 55 + (x - xmin) / xspan * (width - 85)
        py = 25 + (ymax - y) / yspan * (height - 75)
        circles.append('<circle cx="{:.2f}" cy="{:.2f}" r="4" fill="#2f725f" opacity=".72"><title>{}: PC1={:.5g}, PC2={:.5g}</title></circle>'.format(px, py, html.escape(label), x, y))
    return '<svg viewBox="0 0 {} {}" role="img" aria-label="PCA PC1 PC2"><line x1="55" y1="{}" x2="{}" y2="{}" stroke="#9eb2aa"/><line x1="55" y1="25" x2="55" y2="{}" stroke="#9eb2aa"/><text x="{}" y="{}" text-anchor="middle">PC1</text><text x="16" y="{}" transform="rotate(-90 16 {})" text-anchor="middle">PC2</text>{}</svg>'.format(width, height, height - 50, width - 30, height - 50, height - 50, width / 2, height - 12, height / 2, height / 2, "".join(circles))


def _svg_ld_decay(rows, width=760, height=300):
    values = [
        (str(x.get("distance_bin_kb")), float(x.get("mean_r2")), int(x.get("pair_count") or 0))
        for x in rows if x.get("mean_r2") is not None
    ]
    if not values:
        return '<div class="empty">LD衰减没有可绘制的位点对</div>'
    maximum = max(value for _, value, _ in values)
    raw_limit = max(maximum * 1.16, 0.001)
    magnitude = 10 ** math.floor(math.log10(raw_limit))
    scaled = raw_limit / magnitude
    nice_scaled = next(step for step in (1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10) if scaled <= step)
    axis_max = nice_scaled * magnitude
    margin_left, margin_right, margin_top, margin_bottom = 72, 28, 32, 68
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    points = []
    for index, (label, value, pair_count) in enumerate(values):
        x = margin_left + index * plot_width / max(1, len(values) - 1)
        y = margin_top + (1 - max(0, value) / axis_max) * plot_height
        points.append((x, y, label, value, pair_count))
    polyline = " ".join("{:.2f},{:.2f}".format(x, y) for x, y, _, _, _ in points)
    grids = []
    for tick in range(5):
        value = axis_max * tick / 4
        y = margin_top + (1 - tick / 4) * plot_height
        grids.append('<line x1="{}" y1="{:.2f}" x2="{}" y2="{:.2f}" stroke="#dfe8e3"/><text x="{}" y="{:.2f}" text-anchor="end" font-size="10" fill="#52645d">{:.3g}</text>'.format(margin_left, y, width - margin_right, y, margin_left - 8, y + 3, value))
    marks = "".join(
        '<circle cx="{:.2f}" cy="{:.2f}" r="5" fill="#c76b27"><title>{} kb: mean r²={:.4f}; {:,} pairs</title></circle>'
        '<text x="{:.2f}" y="{:.2f}" text-anchor="middle" font-size="10" font-weight="700" fill="#8b4717">{:.4f}</text>'
        '<text x="{:.2f}" y="{}" text-anchor="middle" font-size="10">{}</text>'
        '<text x="{:.2f}" y="{}" text-anchor="middle" font-size="9" fill="#6d7e77">{:,}对</text>'.format(
            x, y, html.escape(label), value, pair_count, x, max(13, y - 10), value,
            x, height - 35, html.escape(label), x, height - 20, pair_count,
        ) for x, y, label, value, pair_count in points
    )
    return '<svg viewBox="0 0 {} {}" role="img" aria-label="LD decay"><text x="18" y="{}" transform="rotate(-90 18 {})" text-anchor="middle" font-size="11">平均 r²</text>{}<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#8fa39a"/><line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#8fa39a"/><polyline points="{}" fill="none" stroke="#c76b27" stroke-width="3"/>{}<text x="{}" y="15" text-anchor="end" font-size="10" fill="#6d7e77">纵轴自动缩放：0–{:.3g}</text></svg>'.format(
        width, height, height / 2, height / 2, "".join(grids), margin_left, margin_top + plot_height,
        width - margin_right, margin_top + plot_height, margin_left, margin_top, margin_left,
        margin_top + plot_height, polyline, marks, width - margin_right, axis_max,
    )


def _svg_sample_qc(rows, thresholds=None, width=760, height=340):
    points = [(x.get("missing_rate"), x.get("het_rate"), x.get("sample_id"), x.get("sample_status"), x.get("group") or x.get("batch") or "") for x in rows]
    points = [x for x in points if x[0] is not None and x[1] is not None]
    if not points:
        return '<div class="empty">没有可绘制的样本缺失率/杂合率</div>'
    thresholds = thresholds or {}
    missing_values = [x[0] for x in points]
    het_values = [x[1] for x in points]
    missing_warn = _safe_float(thresholds.get("sample_missing_warn"))
    het_warn = _safe_float(thresholds.get("het_rate_warn"))
    x_core = _quantile(missing_values, .98) if len(points) >= 20 else max(missing_values)
    y_core = _quantile(het_values, .98) if len(points) >= 20 else max(het_values)
    xmax = max((x_core or 0) * 1.18, (missing_warn or 0) * 1.15, .01)
    ymax = max((y_core or 0) * 1.18, (het_warn or 0) * 1.15, .05)
    left, right, top, bottom = 66, 24, 34, 54
    plot_w, plot_h = width - left - right, height - top - bottom
    colors = {"pass": "#2f725f", "warning": "#d18b00", "critical": "#a62d33"}
    grid = []
    for index in range(5):
        fraction = index / 4
        x = left + fraction * plot_w
        y = top + (1 - fraction) * plot_h
        grid.append('<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{base}" stroke="#dce6e1"/><text x="{x:.2f}" y="{label}" text-anchor="middle" font-size="11" fill="#65776f">{value}</text>'.format(
            x=x, top=top, base=top + plot_h, label=height - 29, value=_ratio(fraction * xmax)))
        grid.append('<line x1="{left}" y1="{y:.2f}" x2="{end}" y2="{y:.2f}" stroke="#dce6e1"/><text x="{label}" y="{y:.2f}" text-anchor="end" dominant-baseline="middle" font-size="11" fill="#65776f">{value}</text>'.format(
            left=left, end=left + plot_w, y=y, label=left - 8, value=_ratio(fraction * ymax)))
    threshold_marks = []
    if missing_warn is not None and 0 < missing_warn < xmax:
        x = left + missing_warn / xmax * plot_w
        threshold_marks.append('<line x1="{0:.2f}" y1="{1}" x2="{0:.2f}" y2="{2}" stroke="#c58a28" stroke-width="1.5" stroke-dasharray="5 4"/><text x="{0:.2f}" y="{3}" text-anchor="middle" font-size="10" fill="#8a611d">缺失提醒线</text>'.format(x, top, top + plot_h, top - 10))
    if het_warn is not None and 0 < het_warn < ymax:
        y = top + (1 - het_warn / ymax) * plot_h
        threshold_marks.append('<line x1="{0}" y1="{1:.2f}" x2="{2}" y2="{1:.2f}" stroke="#c58a28" stroke-width="1.5" stroke-dasharray="5 4"/><text x="{2}" y="{3:.2f}" text-anchor="end" font-size="10" fill="#8a611d">杂合提醒线</text>'.format(left, y, left + plot_w, y - 5))
    marks = []
    clipped_labels = []
    for missing, het, sample, status, group in points:
        clipped_x, clipped_y = missing > xmax, het > ymax
        px = left + min(missing, xmax) / xmax * plot_w
        py = top + (1 - min(het, ymax) / ymax) * plot_h
        tooltip = '{}{}：缺失率={}；杂合率={}'.format(html.escape(str(sample)), " · " + html.escape(group) if group else "", _ratio(missing), _ratio(het))
        color = colors.get(status, colors["pass"])
        if clipped_x or clipped_y:
            points_text = "{:.2f},{:.2f} {:.2f},{:.2f} {:.2f},{:.2f}".format(px, py - 6, px - 5.5, py + 4.5, px + 5.5, py + 4.5)
            marks.append('<polygon points="{}" fill="{}" stroke="#ffffff" stroke-width="1.2"><title>{}；超出主显示范围</title></polygon>'.format(points_text, color, tooltip))
            if len(clipped_labels) < 6:
                clipped_labels.append('<text x="{:.2f}" y="{:.2f}" text-anchor="end" font-size="10" fill="#53645d">{}</text>'.format(px - 7, max(top + 10, py - 8), html.escape(str(sample))))
        else:
            marks.append('<circle cx="{:.2f}" cy="{:.2f}" r="3.6" fill="{}" stroke="#ffffff" stroke-width=".8" opacity=".82"><title>{}</title></circle>'.format(px, py, color, tooltip))
    legend = '<g transform="translate({},{})" font-size="10" fill="#53645d"><circle cx="0" cy="0" r="3.5" fill="#2f725f"/><text x="8" y="3">通过</text><circle cx="50" cy="0" r="3.5" fill="#d18b00"/><text x="58" y="3">提醒</text><circle cx="100" cy="0" r="3.5" fill="#a62d33"/><text x="108" y="3">关键</text><polygon points="157,-5 152,4 162,4" fill="#65776f"/><text x="168" y="3">超出P98主区间</text></g>'.format(left + 8, 16)
    return '<svg viewBox="0 0 {width} {height}" role="img" aria-labelledby="sample-qc-title sample-qc-desc"><title id="sample-qc-title">样本缺失率与杂合率</title><desc id="sample-qc-desc">主坐标按百分之九十八分位缩放，超范围样本以三角形贴边并直接标注。</desc>{legend}{grid}{thresholds}<line x1="{left}" y1="{base}" x2="{end}" y2="{base}" stroke="#83978e"/><line x1="{left}" y1="{top}" x2="{left}" y2="{base}" stroke="#83978e"/><text x="{mid_x}" y="{axis_y}" text-anchor="middle" font-size="12">样本缺失率</text><text x="17" y="{mid_y}" transform="rotate(-90 17 {mid_y})" text-anchor="middle" font-size="12">样本杂合率</text>{marks}{labels}</svg>'.format(
        width=width, height=height, legend=legend, grid="".join(grid), thresholds="".join(threshold_marks),
        left=left, base=top + plot_h, end=left + plot_w, top=top, mid_x=left + plot_w / 2,
        axis_y=height - 7, mid_y=top + plot_h / 2, marks="".join(marks), labels="".join(clipped_labels))


def render_report(result):
    summary = result["summary"]
    profile = result["profile"]
    scan = result["scan"]
    site = result["site_metrics"]
    warnings = result["warnings"]
    sample_rows = []
    for item in result["samples"]:
        sample_rows.append("<tr data-sample='{}' data-status='{}'><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            html.escape(item["sample_id"].lower()), html.escape(item.get("sample_status") or "pass"),
            html.escape(item["sample_id"]), html.escape(item.get("group") or "—"), html.escape(item.get("batch") or "—"), _fmt(item.get("sample_score")), html.escape(item.get("sample_status") or "—"), _ratio(item["missing_rate"]), _ratio(item["het_rate"]),
            _fmt(item["median_dp"]), _fmt(item["median_gq"]), _ratio(item["ab_outlier_rate"]),
            _ratio(item.get("phase_rate"))))
    warning_rows = []
    for item in warnings:
        warning_rows.append("<tr data-level='{}'><td><span class='badge {}'>{}</span></td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            item["level"], item["level"], item["level"].upper(), html.escape(item["scope"]), html.escape(str(item["target"])),
            html.escape(item["message"]), html.escape(item["advice"])))
    threshold_rows = []
    for key, value in profile["thresholds"].items():
        threshold_rows.append("<tr><td>{}</td><td>{}</td><td>{}</td></tr>".format(html.escape(key), html.escape(str(value)), html.escape(profile["threshold_sources"].get(key, ""))))
    field_labels = {
        "species_name": "物种", "ploidy": "生物学倍性", "genotype_ploidy": "VCF GT编码倍性",
        "subgenomes": "亚基因组数量", "mating_system": "繁殖/材料类型",
    }
    field_rows = []
    for key, label in field_labels.items():
        field_rows.append("<tr><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            html.escape(label), html.escape(str(profile.get(key, "—"))),
            html.escape((profile.get("field_sources") or {}).get(key, ""))))
    status_label = {"pass": "通过", "warning": "需关注", "critical": "高风险"}[summary["status"]]
    population = result.get("population_analysis") or {}
    population_rows = []
    module_names = {"hwe": "HWE", "pca": "PCA", "kinship": "亲缘关系/IBD", "ld": "LD衰减", "roh": "ROH"}
    for key, module in (population.get("modules") or {}).items():
        details = module.get("summary") or {}
        compact = "；".join("{}={}".format(k, _fmt(v)) for k, v in list(details.items())[:6])
        files = " ".join("<a href='{}'>{}</a>".format(html.escape(name), html.escape(name)) for name in module.get("artifacts", []))
        population_rows.append("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            html.escape(module_names.get(key, key)), html.escape(module.get("status", "")),
            html.escape(compact or module.get("error") or "—"), files or "—"))
    population_section = ""
    if population:
        modules = population.get("modules") or {}
        pca_chart = _svg_pca((modules.get("pca") or {}).get("preview") or [])
        ld_chart = _svg_ld_decay((modules.get("ld") or {}).get("preview") or [])
        population_section = "<section class='card' id='section-population'><h2>群体遗传分析</h2><p class='note'>模块使用QC过滤后的二等位面板；PCA/亲缘使用LD剪枝标记，LD衰减使用未剪枝的均匀限量标记，避免人为压低r²。HWE是否参与质量解释由Profile决定，其余默认作为探索性证据。</p><div class='table'><table><thead><tr><th>模块</th><th>状态</th><th>摘要</th><th>结果文件</th></tr></thead><tbody>{}</tbody></table></div></section><div class='grid'><section class='card'><h2>PCA：PC1 × PC2</h2>{}</section><section class='card'><h2>LD衰减（平均r²）</h2>{}</section></div>".format("".join(population_rows) or "<tr><td colspan='4'>未启用</td></tr>", pca_chart, ld_chart)
    metric_rows = []
    for key, item in (site.get("quality_field_summaries") or {}).items():
        metric_rows.append("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            html.escape(key), _ratio(item.get("coverage")), _fmt(item.get("q05")),
            _fmt(item.get("median")), _fmt(item.get("q95")), _fmt(item.get("max"))))
    availability_rows = ["<tr><td>{}</td><td>{}</td><td>{}</td></tr>".format(
        html.escape(x["module"]), html.escape(x["status"]), html.escape(x["reason"]))
        for x in result.get("module_availability", [])]
    group_rows = ["<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
        html.escape(x["dimension"]), html.escape(x["label"]), _fmt(x["sample_count"]), _ratio(x["median_missing_rate"]),
        _ratio(x["median_het_rate"]), _fmt(x["median_dp"]), _fmt(x["median_gq"]))
        for x in result.get("cohort_metrics", {}).get("group_summaries", [])]
    subgenome_rows = ["<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
        html.escape(x["subgenome"]), _fmt(x["evaluated_sites"]), _fmt(x["called_genotypes"]), _ratio(x["heterozygosity_rate"]))
        for x in result.get("cohort_metrics", {}).get("subgenome_summary", [])]
    fake_window_rows = ["<tr><td>{}:{}-{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
        html.escape(str(x["chrom"])), _fmt(x["start"], 0), _fmt(x["end"], 0), _ratio(x["het_rate"]),
        _ratio(x["quarter_ab_rate"]), _fmt(x["mean_dp"]), "是" if x["fake_het_candidate"] else "")
        for x in result.get("fake_heterozygosity_windows", [])[:30]]
    qc_recommendation = result.get("qc_recommendation") or {}
    qc_parameter_rows = ["<tr><td>{}</td><td>{}</td></tr>".format(html.escape(str(key)), html.escape(str(value if value is not None else "关闭"))) for key, value in (qc_recommendation.get("parameters") or {}).items()]
    qc_reason_rows = "".join("<li>{}</li>".format(html.escape(str(item))) for item in qc_recommendation.get("reasons", []))
    readiness = result.get("analysis_readiness") or build_analysis_readiness(result)
    priority_rows = "".join(
        "<article class='priority {}'><div><span>{}</span><b>{}</b><small>{} · {}</small></div><a href='#{}'>查看证据</a></article>".format(
            html.escape(item.get("level") or "warning"), html.escape((item.get("level") or "warning").upper()),
            html.escape(item.get("message") or ""), html.escape(item.get("scope") or ""), html.escape(str(item.get("target") or "")),
            html.escape(item.get("anchor") or "section-warnings"),
        ) for item in readiness.get("top_priorities", [])
    ) or "<p class='empty'>没有需要优先处理的关键告警。</p>"
    dimension_rows = "".join(
        "<div class='dimension'><span>{}</span><b>{}</b></div>".format(html.escape(label), "NA" if value is None else "{}/100".format(value))
        for key, label in (("file_reference", "文件与参考"), ("sample_median", "样本中位"), ("site", "位点质量"), ("population_integrity", "群体完整性"), ("provenance", "来源证据"), ("requested_module_completion", "所选模块完成度"))
        for value in [readiness.get("score_dimensions", {}).get(key)]
    )
    download_names = [
        ("report_summary.json", "JSON摘要"), ("sample_metrics.tsv", "样本指标"),
        ("variant_metrics.tsv", "位点指标"), ("warnings.tsv", "告警表"),
        ("recommend_filters.tsv", "推荐过滤"),
    ]
    population_artifacts = set((result.get("population_analysis") or {}).get("artifacts") or [])
    for filename, label in (("pairwise_similarity.tsv", "样本对相似度"), ("roh_segments.tsv", "ROH区段"), ("pca_scores.tsv", "PCA坐标"), ("ld_decay.tsv", "LD衰减")):
        if filename in population_artifacts:
            download_names.append((filename, label))
    download_names.append(("GPA_Accelerator_VCF_QC_report.zip", "完整报告包"))
    download_links = "".join("<a href='{}'>{}</a>".format(html.escape(filename), html.escape(label)) for filename, label in download_names)
    embedded = json.dumps(result, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    template = """<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>GPA-Accelerator VCF质量评估报告</title><style>
    :root{--ink:#18342b;--muted:#667a72;--line:#d8e2dd;--brand:#2f725f;--soft:#f2f7f4;--warn:#a96200;--crit:#a62d33}*{box-sizing:border-box;scroll-behavior:smooth}body{margin:0;background:#edf3ef;color:var(--ink);font-family:"Microsoft YaHei",Arial,sans-serif}main{max-width:1180px;margin:auto;padding:30px}.hero,.card{background:white;border:1px solid var(--line);border-radius:18px;padding:24px;margin-bottom:18px}.hero{background:linear-gradient(135deg,#173c31,#347966);color:white}.hero h1{font-size:34px;margin:5px 0}.hero p{opacity:.82}.report-nav{position:sticky;top:0;z-index:20;display:flex;gap:6px;overflow:auto;margin:0 0 18px;padding:8px;border:1px solid var(--line);border-radius:14px;background:rgba(255,255,255,.94);box-shadow:0 8px 26px rgba(24,52,43,.08)}.report-nav a{padding:8px 11px;border-radius:9px;color:var(--brand);font-size:12px;font-weight:700;text-decoration:none;white-space:nowrap}.report-nav a:hover{background:var(--soft)}.kpis{display:grid;grid-template-columns:repeat(5,1fr);gap:12px}.kpi{background:var(--soft);border-radius:14px;padding:16px}.kpi b{display:block;font-size:25px;margin-top:6px}.score{font-size:64px;font-weight:800}.status-pass{color:#d6ffe8}.status-warning{color:#ffe09c}.status-critical{color:#ffb1b4}.readiness{border-left:6px solid var(--brand)}.readiness.caution{border-left-color:var(--warn)}.readiness.hold{border-left-color:var(--crit)}.readiness h2{margin-bottom:5px}.priority-list{display:grid;gap:8px;margin:14px 0}.priority{display:flex;align-items:center;justify-content:space-between;gap:12px;padding:11px 13px;border:1px solid var(--line);border-radius:11px;background:#fafcfb}.priority div{display:grid;grid-template-columns:auto 1fr;gap:3px 9px}.priority span{grid-row:1/3;padding:3px 6px;border-radius:7px;background:#fff0cc;color:#805000;font-size:10px;font-weight:800}.priority.critical span{background:#ffe0e1;color:#98242b}.priority small{color:var(--muted)}.priority a{color:var(--brand);font-size:12px;font-weight:700;white-space:nowrap}.dimensions{display:grid;grid-template-columns:repeat(6,1fr);gap:8px}.dimension{padding:10px;border-radius:10px;background:var(--soft)}.dimension span,.dimension b{display:block}.dimension span{color:var(--muted);font-size:10px}.dimension b{margin-top:4px;font-size:15px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:18px}h2{font-size:21px}table{width:100%;border-collapse:collapse;font-size:12px}th,td{text-align:left;padding:9px;border-bottom:1px solid #e8eeeb}th{position:sticky;top:0;background:#f6faf8}.table{max-height:520px;overflow:auto;border:1px solid var(--line);border-radius:10px}.badge{padding:3px 7px;border-radius:10px;font-weight:700}.badge.info{background:#e8eef3}.badge.warning{background:#fff0cc;color:#805000}.badge.critical{background:#ffe0e1;color:#98242b}.note{padding:12px;background:#fff7df;border-left:4px solid #d18b00}.filterbar{display:flex;gap:8px;flex-wrap:wrap;margin:10px 0}.filterbar button,.filterbar input{min-height:36px;padding:7px 10px;border:1px solid var(--line);border-radius:8px;background:white;color:var(--ink)}.filterbar button{cursor:pointer}.filterbar input{min-width:260px}.actions{display:flex;gap:8px;flex-wrap:wrap}.actions button{border:0;border-radius:9px;padding:10px 14px;background:#e8f1ed;color:var(--ink);cursor:pointer}.actions button:first-child{background:white}.download-grid{display:flex;gap:8px;flex-wrap:wrap}.download-grid a{padding:9px 11px;border:1px solid var(--line);border-radius:9px;background:var(--soft);color:var(--brand);font-size:12px;font-weight:700;text-decoration:none}.empty{padding:30px;color:var(--muted)}svg{width:100%;height:auto}@media(max-width:800px){.kpis,.grid{grid-template-columns:1fr 1fr}.dimensions{grid-template-columns:repeat(3,1fr)}}@media print{body{background:white}main{max-width:none;padding:0}.actions,.report-nav,.filterbar{display:none}.card,.hero{break-inside:avoid;border-color:#bbb}.table{max-height:none;overflow:visible}}
    </style></head><body><main><section class='hero' id='section-overview'><div class='actions'><button onclick='window.print()'>打印/另存为PDF</button><button onclick='downloadJson()'>下载嵌入JSON</button></div><p>GPA-Accelerator · VCF DEEP QUALITY REPORT</p><h1>VCF质量评估报告</h1><div class='score status-{status}'>{score} <small style='font-size:22px'>/ 100 · {status_label}</small></div><p>{name} · {generated}</p></section>
    <nav class='report-nav' aria-label='报告章节'><a href='#section-overview'>总览</a><a href='#section-input'>输入审计</a><a href='#section-site'>位点质量</a><a href='#section-samples'>样本质量</a><a href='#section-population'>群体与亲缘</a><a href='#section-warnings'>告警</a><a href='#section-recommendations'>质控建议</a><a href='#section-downloads'>下载</a></nav>
    <section class='kpis'><div class='kpi'>记录数<b>{records}</b></div><div class='kpi'>样本数<b>{samples}</b></div><div class='kpi'>染色体/Contig<b>{contigs}</b></div><div class='kpi'>严重告警<b>{critical}</b></div><div class='kpi'>一般告警<b>{warning}</b></div></section>
    <section class='card readiness {readiness_code}'><h2>{readiness_title}</h2><p>{readiness_description}</p><div class='priority-list'>{priority_rows}</div><h3>分层评分</h3><div class='dimensions'>{dimension_rows}</div></section>
    <section class='card' id='section-input'><h2>运行口径与输入审计</h2><p><b>Profile：</b>{profile_name}；<b>作物预设：</b>{crop_name}；<b>物种：</b>{species}；<b>生物学倍性：</b>{ploidy}；<b>VCF GT编码倍性：</b>{gt_ploidy}；<b>亚基因组：</b>{subgenomes}</p><p><b>建议基因组大小：</b>{genome_size} Mb；<b>参考版本提示：</b>{reference_hint}</p><p class='note'>{focus_note}</p><p><b>扫描：</b>{scan_mode}，评估 {evaluated}/{records} 条记录，耗时 {elapsed} 秒。</p><p class='note'>{method_note}</p></section>
    <div class='grid' id='section-site'><section class='card'><h2>变异类型</h2>{type_chart}</section><section class='card'><h2>位点缺失率分布</h2>{missing_chart}</section></div>
    <div class='grid'><section class='card'><h2>MAF分布</h2>{maf_chart}</section><section class='card'><h2>SV长度分布</h2>{sv_chart}</section></div>
    <section class='card'><h2>位点质量字段分布</h2><p>覆盖率表示抽样位点中该字段存在的比例；缺字段显示为NA，不按0分处理。</p><div class='table'><table><thead><tr><th>字段</th><th>覆盖率</th><th>P05</th><th>中位数</th><th>P95</th><th>最大值</th></tr></thead><tbody>{metric_rows}</tbody></table></div><p>Header质控证据分：<b>{evidence_score}/100</b>；可继续最简化INDEL：{nonminimal}；相邻重复记录：{duplicates}。</p></section>
    <section class='card'><h2>模块可用性</h2><p class='note'>需要外部输入的模块不会伪造结果，也不会因为用户未提供文件而扣质量分。</p><div class='table'><table><thead><tr><th>模块</th><th>状态</th><th>原因/所需输入</th></tr></thead><tbody>{availability_rows}</tbody></table></div></section>
    <div class='grid' id='section-samples'><section class='card'><h2>样本缺失率 × 杂合率</h2>{sample_qc_chart}</section><section class='card' id='section-annotation'><h2>参考、区域与注释审计</h2><p><b>REF核验：</b>{reference_status}；检查 {reference_checked} 条，错配 {reference_mismatch} 条。</p><p><b>区域BED：</b>{region_status}；重叠记录 {region_overlap}（{region_rate}）。</p><p><b>功能注释：</b>{annotation_tags}。</p></section></div>
    <section class='card'><h2>分组与批次质量</h2><div class='table'><table><thead><tr><th>维度</th><th>标签</th><th>样本数</th><th>中位缺失率</th><th>中位杂合率</th><th>中位DP</th><th>中位GQ</th></tr></thead><tbody>{group_rows}</tbody></table></div></section>
    <div class='grid'><section class='card'><h2>亚基因组概览</h2><div class='table'><table><thead><tr><th>亚基因组</th><th>评估位点</th><th>有效GT</th><th>杂合率</th></tr></thead><tbody>{subgenome_rows}</tbody></table></div></section><section class='card'><h2>假杂合候选窗口（前30）</h2><div class='table'><table><thead><tr><th>窗口</th><th>杂合率</th><th>0.25/0.75 AB比例</th><th>平均DP</th><th>候选</th></tr></thead><tbody>{fake_window_rows}</tbody></table></div></section></div>
    <section class='card' id='section-warnings'><h2>告警与建议</h2><div class='filterbar'><button onclick="filterWarnings('all')">全部</button><button onclick="filterWarnings('critical')">仅关键</button><button onclick="filterWarnings('warning')">仅提醒</button><button onclick="filterWarnings('info')">仅信息</button></div><div class='table'><table><thead><tr><th>级别</th><th>范围</th><th>对象</th><th>问题</th><th>建议</th></tr></thead><tbody id='warningTableBody'>{warning_rows}</tbody></table></div></section>
    <section class='card' id='section-recommendations'><h2>物种与VCF自适应质控建议</h2><p class='note'>这是保守起始方案，不会自动覆盖原VCF。预计触发位点过滤：{qc_estimated_removal}；高缺失候选样本：{qc_candidate_samples} 个，默认不自动删除。</p><div class='grid'><div class='table'><table><thead><tr><th>参数</th><th>推荐值</th></tr></thead><tbody>{qc_parameter_rows}</tbody></table></div><div><ul>{qc_reason_rows}</ul></div></div></section>
    <section class='card'><h2>样本质量指标</h2><div class='filterbar'><input id='sampleFilter' type='search' placeholder='搜索样本ID' oninput='filterSamples(this.value)'><button onclick="setSampleStatus('all')">全部状态</button><button onclick="setSampleStatus('critical')">关键</button><button onclick="setSampleStatus('warning')">提醒</button></div><div class='table'><table><thead><tr><th>样本</th><th>Group</th><th>Batch</th><th>样本分</th><th>状态</th><th>缺失率</th><th>杂合率</th><th>中位DP</th><th>中位GQ</th><th>AB异常</th><th>相位率</th></tr></thead><tbody id='sampleTableBody'>{sample_rows}</tbody></table></div></section>
    <section class='card'><h2>作物字段、有效阈值与来源</h2><p class='note'>作物预设是建议起点，不是锁定规则；“用户自定义（基于某作物）”表示该项已被手动覆盖。</p><div class='grid'><div class='table'><table><thead><tr><th>字段</th><th>有效值</th><th>来源</th></tr></thead><tbody>{field_rows}</tbody></table></div><div class='table'><table><thead><tr><th>阈值</th><th>有效值</th><th>来源</th></tr></thead><tbody>{threshold_rows}</tbody></table></div></div></section>
    {population_section}
    <section class='card'><h2>自动修复安全策略</h2><p>原始VCF永不被静默覆盖。覆盖源文件、改写REF/ALT/GT、坐标转换、染色体批量重命名和删除文件均被定义为危险操作，执行前必须再次确认。</p></section>
    <section class='card' id='section-downloads'><h2>机器可读结果与复核材料</h2><p>下列文件与HTML使用同一数据内核，可直接用于R、Excel或后续脚本。</p><div class='download-grid'>{download_links}</div></section>
    <script type='application/json' id='reportData'>{embedded}</script><script>let sampleStatus='all';function downloadJson(){const text=document.getElementById('reportData').textContent;const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([text],{type:'application/json;charset=utf-8'}));a.download='report_summary.json';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000)}function filterWarnings(level){document.querySelectorAll('#warningTableBody tr').forEach(row=>row.hidden=level!=='all'&&row.dataset.level!==level)}function setSampleStatus(status){sampleStatus=status;filterSamples(document.getElementById('sampleFilter').value)}function filterSamples(value){const query=(value||'').trim().toLowerCase();document.querySelectorAll('#sampleTableBody tr').forEach(row=>row.hidden=(sampleStatus!=='all'&&row.dataset.status!==sampleStatus)||(query&&!row.dataset.sample.includes(query)))}</script></main></body></html>"""
    values = {
        "status": summary["status"], "score": summary["score"], "status_label": status_label,
        "name": html.escape(result["input"]["name"]), "generated": html.escape(result["generated_at"]),
        "records": "{:,}".format(summary["record_count"]), "samples": summary["sample_count"], "contigs": summary["contig_count"],
        "critical": summary["critical_count"], "warning": summary["warning_count"], "profile_name": html.escape(profile["name"]),
        "species": html.escape(profile.get("species_name") or profile["kingdom"]), "ploidy": profile["ploidy"],
        "gt_ploidy": profile.get("effective_gt_ploidy") or "未识别", "subgenomes": profile["subgenomes"],
        "crop_name": html.escape(profile.get("crop_name") or "未选择（手动/Profile模式）"),
        "genome_size": html.escape(str(profile.get("genome_size_mb") or "—")),
        "reference_hint": html.escape(profile.get("reference_hint") or "未由作物预设提供"),
        "focus_note": html.escape(
            "非植物自定义模式：全部核心阈值由使用者提供并负责解释；GPA-Accelerator不提供动物默认参数。"
            if profile.get("analysis_scope") == "non_plant_custom"
            else ((profile.get("preset_disclaimer") or "") + " 所有作物参数均可在运行前自由修改，报告记录最终有效值与来源。" if profile.get("crop_id") else "植物分析模式：使用内置植物Profile，并记录全部用户覆盖参数及来源。")
        ),
        "scan_mode": "完整扫描" if scan["mode"] == "full" else "智能抽样", "evaluated": "{:,}".format(scan["evaluated_records"]),
        "elapsed": scan["elapsed_seconds"], "method_note": html.escape(scan["method_note"]),
        "type_chart": _svg_bars(sorted(site["variant_types"].items()), "SNP / INDEL / SV"),
        "missing_chart": _svg_bars([(x, site["missing_rate_bins"].get(x, 0)) for x in ["0–1%", "1–5%", "5–10%", "10–20%", ">20%"]], "位点缺失率"),
        "maf_chart": _svg_bars([(x, site["maf_bins"].get(x, 0)) for x in ["0", "(0,1%]", "(1%,5%]", "(5%,10%]", "(10%,30%]", "(30%,50%]"]], "次要等位基因频率"),
        "sv_chart": _svg_bars(list((site.get("sv_length_bins") or {}).items()), "SV长度"),
        "metric_rows": "".join(metric_rows) or "<tr><td colspan='6'>VCF未提供这些字段</td></tr>",
        "availability_rows": "".join(availability_rows),
        "sample_qc_chart": _svg_sample_qc(result.get("samples", []), result.get("profile", {}).get("thresholds", {})),
        "group_rows": "".join(group_rows) or "<tr><td colspan='7'>未提供样本元数据或没有可分组字段</td></tr>",
        "subgenome_rows": "".join(subgenome_rows) or "<tr><td colspan='4'>无可分配亚基因组</td></tr>",
        "fake_window_rows": "".join(fake_window_rows) or "<tr><td colspan='5'>没有候选窗口</td></tr>",
        "reference_status": html.escape(result.get("reference_audit", {}).get("status", "not_provided")),
        "reference_checked": _fmt(result.get("reference_audit", {}).get("checked_records", 0)),
        "reference_mismatch": _fmt(result.get("reference_audit", {}).get("mismatch_records", 0)),
        "region_status": html.escape(result.get("region_audit", {}).get("status", "not_provided")),
        "region_overlap": _fmt(result.get("region_audit", {}).get("overlap_records", 0)),
        "region_rate": _ratio(result.get("region_audit", {}).get("overlap_rate")),
        "annotation_tags": html.escape("；".join("{}={}".format(k, v) for k, v in result.get("annotation_audit", {}).get("tag_record_counts", {}).items()) or "未检测到ANN/CSQ/BCSQ"),
        "evidence_score": result.get("header_audit", {}).get("qc_evidence_score", "—"),
        "nonminimal": site.get("nonminimal_indel_count", 0), "duplicates": site.get("adjacent_duplicate_count", 0),
        "warning_rows": "".join(warning_rows) or "<tr><td colspan='5'>未触发告警</td></tr>", "sample_rows": "".join(sample_rows),
        "qc_estimated_removal": _ratio(qc_recommendation.get("estimated_site_removal_fraction")),
        "qc_candidate_samples": len(qc_recommendation.get("sample_exclusion_candidates") or []),
        "qc_parameter_rows": "".join(qc_parameter_rows) or "<tr><td colspan='2'>无可用参数</td></tr>", "qc_reason_rows": qc_reason_rows,
        "field_rows": "".join(field_rows), "threshold_rows": "".join(threshold_rows), "population_section": population_section, "embedded": embedded,
        "readiness_code": html.escape(readiness.get("code") or "caution"), "readiness_title": html.escape(readiness.get("title") or "需要复核"),
        "readiness_description": html.escape(readiness.get("description") or ""), "priority_rows": priority_rows, "dimension_rows": dimension_rows,
        "download_links": download_links,
    }
    for key, value in values.items():
        template = template.replace("{" + key + "}", str(value))
    return template


class QualityJobManager:
    def __init__(self, service):
        self.service = service
        self.evaluator = QualityEvaluator(service)
        self.population = PopulationAnalyzer()
        self._jobs = {}
        self._lock = threading.Lock()

    @staticmethod
    def default_report_root():
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / ".local" / "share")
        return Path(base) / "GPA-Accelerator" / "vcf-reports"

    def catalog(self):
        return {
            "profiles": profile_catalog(),
            "crops": crop_catalog(),
            "default_report_root": str(self.default_report_root()),
            "analysis_focus": "plant",
            "non_plant_mode": "custom_parameters_only",
        }

    def start(self, payload):
        path = self.service._validate_file(payload.get("path"))
        profile = resolve_profile(payload.get("profile"))
        output_text = str(payload.get("output_dir") or "").strip()
        root = Path(os.path.expandvars(os.path.expanduser(output_text))).resolve() if output_text else self.default_report_root().resolve()
        run_id = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:8]
        run_dir = root / run_id
        job = {
            "id": run_id, "status": "queued", "message": "等待运行", "progress": 0.0,
            "processed_records": 0, "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "run_dir": str(run_dir), "artifacts": [], "result": None, "error": None,
            "cancel": threading.Event(),
        }
        raw_population = payload.get("population_options") if isinstance(payload.get("population_options"), dict) else {}
        population_options = {key: bool(raw_population.get(key, False)) for key in ("hwe", "pca", "kinship", "ld", "roh")}
        for key in ("site_missing", "maf", "prune_window", "prune_step", "prune_r2", "pca_components", "ld_max_markers", "ld_window_kb"):
            if key in raw_population:
                population_options[key] = raw_population[key]
        config = {
            "profile": payload.get("profile") or {"profile_id": profile["id"]},
            "scan_mode": payload.get("scan_mode") or "smart",
            "target_records": payload.get("target_records") or 200000,
            "reference_path": payload.get("reference_path") or "",
            "sample_meta_path": payload.get("sample_meta_path") or "",
            "region_bed_path": payload.get("region_bed_path") or "",
            "population_options": population_options,
        }
        with self._lock:
            self._jobs[run_id] = job
        thread = threading.Thread(target=self._run, args=(job, str(path), config), daemon=True)
        job["thread"] = thread
        thread.start()
        return self._public(job)

    def _run(self, job, path, config):
        try:
            run_dir = Path(job["run_dir"])
            run_dir.mkdir(parents=True, exist_ok=False)
            job["status"] = "running"
            job["message"] = "正在读取VCF"

            def progress(records, percent, message):
                job["processed_records"] = records
                job["progress"] = round(percent * .55, 2) if percent is not None else None
                job["message"] = message

            result = self.evaluator.evaluate(path, config, progress, job["cancel"])
            population_options = config.get("population_options") or {}
            if any(population_options.get(key, False) for key in ("hwe", "pca", "kinship", "ld", "roh")):
                def population_progress(percent, message):
                    job["progress"] = round(percent, 2)
                    job["message"] = message
                result["population_analysis"] = self.population.run(
                    path, result["profile"], run_dir, population_options,
                    population_progress, job["cancel"],
                )
            else:
                result["population_analysis"] = self.population.run(
                    path, result["profile"], run_dir, population_options,
                    None, job["cancel"],
                )
            self._merge_population_evidence(result)
            variant_rows = result.pop("_variant_metrics_rows", [])
            report_path = run_dir / "report.html"
            json_path = run_dir / "report_summary.json"
            samples_path = run_dir / "sample_metrics.tsv"
            warnings_path = run_dir / "warnings.tsv"
            manifest_path = run_dir / "run_manifest.json"
            variants_path = run_dir / "variant_metrics.tsv"
            site_summary_path = run_dir / "site_metric_summary.tsv"
            density_path = run_dir / "density_windows.tsv"
            modules_path = run_dir / "module_availability.tsv"
            filters_path = run_dir / "recommend_filters.tsv"
            sv_path = run_dir / "sv_metrics.tsv"
            groups_path = run_dir / "group_batch_metrics.tsv"
            subgenome_path = run_dir / "subgenome_metrics.tsv"
            fake_het_path = run_dir / "fake_heterozygosity_windows.tsv"
            consequences_path = run_dir / "annotation_consequences.tsv"
            priorities_path = run_dir / "analysis_priorities.tsv"
            score_dimensions_path = run_dir / "score_dimensions.tsv"
            report_path.write_text(render_report(result), encoding="utf-8")
            json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
            with samples_path.open("w", encoding="utf-8-sig", newline="") as handle:
                fields = list(result["samples"][0].keys()) if result["samples"] else ["sample_id"]
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader(); writer.writerows(result["samples"])
            warning_fields = ["level", "scope", "target", "code", "message", "evidence", "advice"]
            with warnings_path.open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=warning_fields, delimiter="\t")
                writer.writeheader(); writer.writerows(result["warnings"])
            variant_fields = ["chrom", "pos", "id", "ref", "alt", "type", "qual", "filter", "missing_rate", "maf", "qd", "mq", "fs", "sor", "mq_rank_sum", "read_pos_rank_sum", "svtype", "svlen", "imprecise"]
            with variants_path.open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=variant_fields, delimiter="\t", extrasaction="ignore")
                writer.writeheader(); writer.writerows(variant_rows)
            with sv_path.open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=variant_fields, delimiter="\t", extrasaction="ignore")
                writer.writeheader(); writer.writerows(x for x in variant_rows if x.get("svtype") or x.get("type") in {"SV", "DEL", "DUP", "INV", "INS", "CNV", "BND"})
            summary_rows = []
            for metric, values in result["site_metrics"].get("quality_field_summaries", {}).items():
                summary_rows.append({"metric": metric, **values})
            with site_summary_path.open("w", encoding="utf-8-sig", newline="") as handle:
                fields = ["metric", "available", "n", "coverage", "quantile_sample_n", "min", "q05", "median", "q95", "max"]
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
                writer.writeheader(); writer.writerows(summary_rows)
            _write_rows = lambda path, rows, fields: self._write_tsv(path, rows, fields)
            _write_rows(density_path, result["site_metrics"].get("density_windows", []), ["chrom", "start", "end", "variant_count"])
            _write_rows(modules_path, result.get("module_availability", []), ["module", "status", "reason"])
            thresholds = result["profile"]["thresholds"]
            qc_recommendation = result.get("qc_recommendation") or {}
            qc_parameters = qc_recommendation.get("parameters") or {}
            estimated_removed = qc_recommendation.get("estimated_site_removal_fraction")
            recommended = [
                {"scope": "site", "rule": "F_MISSING", "threshold": "<={}".format(qc_parameters.get("max_site_missing", thresholds["site_missing_warn"])), "reason": "Profile与当前VCF生成的位点缺失起始线", "estimated_removed": estimated_removed, "execution": "建议生成新副本并比较前后指标"},
                {"scope": "sample", "rule": "missing_rate", "threshold": "<{}".format(qc_parameters.get("max_sample_missing", thresholds["sample_missing_warn"])), "reason": "Profile样本缺失提醒线", "estimated_removed": len(qc_recommendation.get("sample_exclusion_candidates") or []), "execution": "只列候选；先复核批次和深度，不自动删除"},
                {"scope": "population", "rule": "biallelic_only", "threshold": "analysis_panel_only", "reason": "HWE/IBS/PCA/LD统一使用QC双等位面板", "estimated_removed": "NA", "execution": "仅用于分析面板，不改写原VCF"},
            ]
            _write_rows(filters_path, recommended, ["scope", "rule", "threshold", "reason", "estimated_removed", "execution"])
            _write_rows(groups_path, result.get("cohort_metrics", {}).get("group_summaries", []), ["dimension", "label", "sample_count", "median_missing_rate", "median_het_rate", "median_dp", "median_gq"])
            _write_rows(subgenome_path, result.get("cohort_metrics", {}).get("subgenome_summary", []), ["subgenome", "evaluated_sites", "called_genotypes", "heterozygous_genotypes", "heterozygosity_rate"])
            _write_rows(fake_het_path, result.get("fake_heterozygosity_windows", []), ["chrom", "start", "end", "evaluated_sites", "het_rate", "quarter_ab_rate", "mean_dp", "fake_het_candidate"])
            consequence_rows = [{"consequence": key, "count": value} for key, value in sorted(result.get("annotation_audit", {}).get("consequence_counts", {}).items(), key=lambda x: x[1], reverse=True)]
            _write_rows(consequences_path, consequence_rows, ["consequence", "count"])
            readiness = result.get("analysis_readiness") or {}
            _write_rows(priorities_path, readiness.get("top_priorities", []), ["level", "code", "scope", "target", "message", "evidence", "advice", "anchor"])
            dimension_rows = [{"dimension": key, "score": value, "missing_is_not_zero": value is None} for key, value in (readiness.get("score_dimensions") or {}).items()]
            _write_rows(score_dimensions_path, dimension_rows, ["dimension", "score", "missing_is_not_zero"])
            stat = Path(path).stat()
            manifest = {
                "run_id": job["id"], "input_path": path, "input_size": stat.st_size,
                "input_mtime_ns": stat.st_mtime_ns, "input_fingerprint": hashlib.sha256((path + str(stat.st_size) + str(stat.st_mtime_ns)).encode()).hexdigest(),
                "schema_version": result["schema_version"], "profile": result["profile"], "scan": result["scan"],
                "ruleset_version": "callvcf-quality-rules-1.2", "analysis_readiness": readiness,
                "input_fingerprint_method": "sha256(path + size + mtime_ns); identity fingerprint, not a full-content checksum",
                "conditional_inputs": {
                    "reference": result.get("reference_audit"),
                    "sample_metadata": result.get("metadata_audit"),
                    "region_bed": result.get("region_audit"),
                },
                "repair_actions_executed": [],
            }
            manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
            zip_path = run_dir / "GPA_Accelerator_VCF_QC_report.zip"
            population_files = []
            for name in result.get("population_analysis", {}).get("artifacts", []):
                item = run_dir / name
                if item.is_file() and item.parent == run_dir:
                    population_files.append(item)
            core_files = [report_path, json_path, samples_path, warnings_path, variants_path, site_summary_path, density_path, modules_path, filters_path, sv_path, groups_path, subgenome_path, fake_het_path, consequences_path, priorities_path, score_dimensions_path, manifest_path]
            with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
                for item in core_files + population_files:
                    archive.write(item, item.name)
            artifacts = []
            for item in core_files + population_files + [zip_path]:
                artifacts.append({"name": item.name, "size": item.stat().st_size, "url": "/api/quality/artifact?run_id={}&name={}".format(job["id"], item.name)})
            job["result"] = result
            job["artifacts"] = artifacts
            job["status"] = "complete"
            job["progress"] = 100.0
            job["message"] = "报告已生成"
        except Exception as exc:
            job["status"] = "cancelled" if job["cancel"].is_set() else "failed"
            job["error"] = str(exc)
            job["message"] = str(exc)

    @staticmethod
    def _merge_population_evidence(result):
        """Promote population evidence into the layered report without over-scoring exploratory modules."""
        population = result.get("population_analysis") or {}
        modules = population.get("modules") or {}
        extra = []
        score_penalty = 0

        kinship = modules.get("kinship") or {}
        if kinship.get("status") == "complete":
            summary = kinship.get("summary") or {}
            critical_pairs = int(summary.get("critical_pairs") or 0)
            warning_pairs = int(summary.get("warning_pairs") or 0)
            if critical_pairs:
                extra.append(_warning(
                    "critical", "sample_pair", "cohort", "POSSIBLE_DUPLICATE_PAIRS",
                    "检测到疑似重复或极高相似样本对",
                    "{}对达到关键阈值；估计器为PLINK IBD/IBS，不等同于KING".format(critical_pairs),
                    "结合样本谱系、批次和测序来源复核；确认前不要自动删除样本",
                ))
                score_penalty += min(15, 4 + critical_pairs)
            if warning_pairs:
                extra.append(_warning(
                    "warning", "sample_pair", "cohort", "RELATED_SAMPLE_PAIRS",
                    "检测到需要复核的高相似样本对",
                    "{}对达到提醒阈值".format(warning_pairs),
                    "查看suspicious_sample_pairs.tsv与relatedness_components.tsv后结合设计解释",
                ))
                score_penalty += min(5, max(1, warning_pairs // 5 + 1))

        pca = modules.get("pca") or {}
        if pca.get("status") == "complete":
            outliers = int((pca.get("summary") or {}).get("robust_pc_outliers") or 0)
            if outliers:
                extra.append(_warning(
                    "warning", "cohort", "cohort", "PCA_ROBUST_OUTLIERS",
                    "PCA中存在稳健距离离群样本",
                    "{}个样本；该结果用于群体结构提示，不直接判定样本失败".format(outliers),
                    "结合分组元数据、批次与亲缘关系共同复核",
                ))

        hwe = modules.get("hwe") or {}
        hwe_summary = hwe.get("summary") or {}
        if hwe.get("status") == "complete" and hwe_summary.get("scoring_enabled"):
            tested = int(hwe_summary.get("tested_sites") or 0)
            failed = int(hwe_summary.get("p_below_profile_threshold") or 0)
            if tested and failed / tested >= .01:
                extra.append(_warning(
                    "warning", "cohort", "cohort", "HWE_EXCESS_DEVIATION",
                    "HWE显著偏离位点比例偏高",
                    "{}/{}（{:.2%}）低于Profile阈值".format(failed, tested, failed / tested),
                    "先排查分层、批次、基因分型错误和近交背景，再决定过滤",
                ))
                score_penalty += min(5, max(1, round(100 * failed / tested)))

        if not extra:
            result["analysis_readiness"] = build_analysis_readiness(result)
            return
        result.setdefault("warnings", []).extend(extra)
        order = {"critical": 0, "warning": 1, "info": 2}
        result["warnings"].sort(key=lambda x: (order.get(x.get("level"), 9), x.get("scope", ""), x.get("target", "")))
        summary = result.setdefault("summary", {})
        summary["critical_count"] = sum(x.get("level") == "critical" for x in result["warnings"])
        summary["warning_count"] = sum(x.get("level") == "warning" for x in result["warnings"])
        summary["score"] = max(0, int(summary.get("score") or 0) - score_penalty)
        result["analysis_readiness"] = build_analysis_readiness(result)
        summary["population_evidence_penalty"] = score_penalty
        summary["status"] = "critical" if summary["critical_count"] else "warning" if summary["warning_count"] else "pass"

    def status(self, run_id):
        with self._lock:
            job = self._jobs.get(str(run_id))
        if not job:
            raise VCFError("质量评估任务不存在或本地服务已重启")
        return self._public(job)

    def cancel(self, run_id):
        with self._lock:
            job = self._jobs.get(str(run_id))
        if not job:
            raise VCFError("质量评估任务不存在")
        job["cancel"].set()
        job["message"] = "正在取消"
        return self._public(job)

    def artifact(self, run_id, name):
        with self._lock:
            job = self._jobs.get(str(run_id))
        if not job:
            raise VCFError("报告任务不存在")
        allowed = {x["name"] for x in job.get("artifacts", [])}
        if name not in allowed or Path(name).name != name:
            raise VCFError("报告文件不存在")
        path = (Path(job["run_dir"]) / name).resolve()
        if path.parent != Path(job["run_dir"]).resolve() or not path.is_file():
            raise VCFError("报告文件不存在")
        return path

    @staticmethod
    def _public(job):
        return {key: value for key, value in job.items() if key not in {"thread", "cancel", "run_dir"}}

    @staticmethod
    def _write_tsv(path, rows, fields):
        with Path(path).open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader(); writer.writerows(rows)
