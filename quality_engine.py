"""Streaming VCF quality evaluation and self-contained report generation."""

import csv
import gzip
import hashlib
import html
import json
import math
import os
import shutil
import statistics
import subprocess
import threading
import time
import uuid
import zipfile
from collections import Counter
from pathlib import Path

from quality_profiles import profile_catalog, resolve_profile
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
                "ab_n": 0, "ab_out": 0,
            } for sample in samples
        }
        header = {
            "fileformat": None, "reference": None, "sources": [], "contigs": [],
            "filters": [], "info_ids": [], "format_ids": [], "has_chrom_header": False,
        }
        type_counts = Counter()
        filter_counts = Counter()
        svtype_counts = Counter()
        contig_counts = Counter()
        maf_bins = Counter()
        missing_bins = Counter()
        site_warn_count = 0
        site_critical_count = 0
        multiallelic = 0
        malformed = 0
        transitions = 0
        transversions = 0
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
                    header["contigs"].append(text.split("##contig=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
                elif text.startswith("##FILTER=<ID="):
                    header["filters"].append(text.split("##FILTER=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
                elif text.startswith("##INFO=<ID="):
                    header["info_ids"].append(text.split("##INFO=<ID=", 1)[1].split(",", 1)[0])
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
            ref, alt, filt, info_text = parts[3], parts[4], parts[6], parts[7]
            info = parse_info(info_text)
            variant_type = classify_variant(ref, alt, "", info_text)
            type_counts[variant_type] += 1
            contig_counts[chrom] += 1
            if "," in alt:
                multiallelic += 1
            for value in filt.split(";"):
                filter_counts[value or "."] += 1
            if variant_type == "SV":
                svtype_counts[info.get("SVTYPE") or "未标注"] += 1

            sampled = stride == 1 or ((sum(ord(x) for x in chrom) * 131 + pos) % stride == 0)
            if sampled:
                sampled_records += 1
                if variant_type == "SNP" and len(ref) == 1 and len(alt) == 1:
                    pair = (ref.upper(), alt.upper())
                    if pair in {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}:
                        transitions += 1
                    else:
                        transversions += 1
                sample_tail = parts[8].split("\t") if len(parts) > 8 else []
                format_keys = sample_tail[0].split(":") if sample_tail else []
                index = {key: i for i, key in enumerate(format_keys)}
                missing = 0
                allele_count = Counter()
                allele_number = 0
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
                    acc["ploidy_counts"][len(alleles)] += 1
                    gt_ploidy_counts[len(alleles)] += 1
                    for allele in alleles:
                        allele_count[allele] += 1
                        allele_number += 1
                    if "DP" in index and index["DP"] < len(values):
                        value = _safe_int(values[index["DP"]])
                        if value is not None and value >= 0:
                            acc["dp"][min(value, 10000)] += 1
                    if "GQ" in index and index["GQ"] < len(values):
                        value = _safe_int(values[index["GQ"]])
                        if value is not None and value >= 0:
                            acc["gq"][min(value, 999)] += 1
                    if category == "het" and "AD" in index and index["AD"] < len(values):
                        ads = [_safe_float(x) for x in values[index["AD"]].split(",")]
                        if len(ads) >= 2 and ads[0] is not None and ads[1] is not None and ads[0] + ads[1] > 0:
                            ab = ads[1] / (ads[0] + ads[1])
                            acc["ab_n"] += 1
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

            if progress and total_records % 5000 == 0:
                expected = metadata.get("record_count")
                pct = min(99.0, total_records * 100 / expected) if expected else None
                progress(total_records, pct, "正在流式扫描VCF")

        self.service._remember_record_count(path, total_records)
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
                "ploidy_mismatch_rate": (
                    sum(count for value, count in acc["ploidy_counts"].items() if value != effective_gt_ploidy) / called
                    if called and effective_gt_ploidy else None
                ),
            })

        warnings = []
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
        variant_classes = set(type_counts)
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
        result = {
            "schema_version": "callvcf-qc-1.0",
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
            },
            "header_audit": header,
            "site_metrics": {
                "variant_types": dict(type_counts), "filters": dict(filter_counts),
                "svtypes": dict(svtype_counts), "contigs": dict(contig_counts),
                "multiallelic_count": multiallelic, "malformed_count": malformed,
                "unsorted_count": unsorted_records, "site_warning_count": site_warn_count,
                "site_critical_count": site_critical_count,
                "missing_rate_bins": dict(missing_bins), "maf_bins": dict(maf_bins),
                "transitions": transitions, "transversions": transversions,
                "titv": transitions / transversions if transversions else None,
            },
            "cohort_metrics": {
                "median_het_rate": het_center, "het_mad": het_mad,
                "het_q01": het_q01, "het_q99": het_q99,
                "het_relative_upper": relative_upper,
                "het_relative_critical": relative_critical,
                "het_rule_mode": het_rule_mode,
            },
            "samples": report_samples,
            "warnings": warnings,
            "repair_policy": {
                "safe_actions": ["建立缺失索引", "生成新的排序副本", "生成新的标准化副本"],
                "dangerous_actions": ["覆盖原VCF", "改写REF/ALT或GT", "坐标转换", "染色体批量重命名", "删除原文件"],
                "dangerous_requires_second_confirmation": True,
                "original_is_never_silently_overwritten": True,
            },
        }
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


def render_report(result):
    summary = result["summary"]
    profile = result["profile"]
    scan = result["scan"]
    site = result["site_metrics"]
    warnings = result["warnings"]
    sample_rows = []
    for item in result["samples"]:
        sample_rows.append("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            html.escape(item["sample_id"]), _ratio(item["missing_rate"]), _ratio(item["het_rate"]),
            _fmt(item["median_dp"]), _fmt(item["median_gq"]), _ratio(item["ab_outlier_rate"]),
            _ratio(item["ploidy_mismatch_rate"])))
    warning_rows = []
    for item in warnings:
        warning_rows.append("<tr><td><span class='badge {}'>{}</span></td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            item["level"], item["level"].upper(), html.escape(item["scope"]), html.escape(str(item["target"])),
            html.escape(item["message"]), html.escape(item["advice"])))
    threshold_rows = []
    for key, value in profile["thresholds"].items():
        threshold_rows.append("<tr><td>{}</td><td>{}</td><td>{}</td></tr>".format(html.escape(key), html.escape(str(value)), html.escape(profile["threshold_sources"].get(key, ""))))
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
        population_section = "<section class='card'><h2>群体遗传分析</h2><p class='note'>这些模块使用过滤后的二等位标记面板；HWE是否参与质量解释由Profile决定，其余默认作为探索性证据，不自动给样本定性。</p><div class='table'><table><thead><tr><th>模块</th><th>状态</th><th>摘要</th><th>结果文件</th></tr></thead><tbody>{}</tbody></table></div></section><div class='grid'><section class='card'><h2>PCA：PC1 × PC2</h2>{}</section><section class='card'><h2>LD衰减（平均r²）</h2>{}</section></div>".format("".join(population_rows) or "<tr><td colspan='4'>未启用</td></tr>", pca_chart, ld_chart)
    embedded = json.dumps(result, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    template = """<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>CallVCF质量评估报告</title><style>
    :root{--ink:#18342b;--muted:#667a72;--line:#d8e2dd;--brand:#2f725f;--soft:#f2f7f4;--warn:#a96200;--crit:#a62d33}*{box-sizing:border-box}body{margin:0;background:#edf3ef;color:var(--ink);font-family:"Microsoft YaHei",Arial,sans-serif}main{max-width:1180px;margin:auto;padding:30px}.hero,.card{background:white;border:1px solid var(--line);border-radius:18px;padding:24px;margin-bottom:18px}.hero{background:linear-gradient(135deg,#173c31,#347966);color:white}.hero h1{font-size:34px;margin:5px 0}.hero p{opacity:.82}.kpis{display:grid;grid-template-columns:repeat(5,1fr);gap:12px}.kpi{background:var(--soft);border-radius:14px;padding:16px}.kpi b{display:block;font-size:25px;margin-top:6px}.score{font-size:64px;font-weight:800}.status-pass{color:#d6ffe8}.status-warning{color:#ffe09c}.status-critical{color:#ffb1b4}.grid{display:grid;grid-template-columns:1fr 1fr;gap:18px}h2{font-size:21px}table{width:100%;border-collapse:collapse;font-size:12px}th,td{text-align:left;padding:9px;border-bottom:1px solid #e8eeeb}th{position:sticky;top:0;background:#f6faf8}.table{max-height:520px;overflow:auto;border:1px solid var(--line);border-radius:10px}.badge{padding:3px 7px;border-radius:10px;font-weight:700}.badge.info{background:#e8eef3}.badge.warning{background:#fff0cc;color:#805000}.badge.critical{background:#ffe0e1;color:#98242b}.note{padding:12px;background:#fff7df;border-left:4px solid #d18b00}.actions{display:flex;gap:8px;flex-wrap:wrap}.actions button{border:0;border-radius:9px;padding:10px 14px;background:#e8f1ed;color:var(--ink);cursor:pointer}.actions button:first-child{background:white}.empty{padding:30px;color:var(--muted)}svg{width:100%;height:auto}@media(max-width:800px){.kpis,.grid{grid-template-columns:1fr 1fr}}@media print{body{background:white}main{max-width:none;padding:0}.actions{display:none}.card,.hero{break-inside:avoid;border-color:#bbb}.table{max-height:none;overflow:visible}}
    </style></head><body><main><section class='hero'><div class='actions'><button onclick='window.print()'>打印/另存为PDF</button><button onclick='downloadJson()'>下载嵌入JSON</button></div><p>CallVCF · VCF DEEP QUALITY REPORT</p><h1>VCF质量评估报告</h1><div class='score status-{status}'>{score} <small style='font-size:22px'>/ 100 · {status_label}</small></div><p>{name} · {generated}</p></section>
    <section class='kpis'><div class='kpi'>记录数<b>{records}</b></div><div class='kpi'>样本数<b>{samples}</b></div><div class='kpi'>染色体/Contig<b>{contigs}</b></div><div class='kpi'>严重告警<b>{critical}</b></div><div class='kpi'>一般告警<b>{warning}</b></div></section>
    <section class='card'><h2>运行口径</h2><p><b>Profile：</b>{profile_name}；<b>物种：</b>{species}；<b>生物学倍性：</b>{ploidy}；<b>VCF GT编码倍性：</b>{gt_ploidy}；<b>亚基因组：</b>{subgenomes}</p><p><b>扫描：</b>{scan_mode}，评估 {evaluated}/{records} 条记录，耗时 {elapsed} 秒。</p><p class='note'>{method_note}</p></section>
    <div class='grid'><section class='card'><h2>变异类型</h2>{type_chart}</section><section class='card'><h2>位点缺失率分布</h2>{missing_chart}</section></div>
    <section class='card'><h2>告警与建议</h2><div class='table'><table><thead><tr><th>级别</th><th>范围</th><th>对象</th><th>问题</th><th>建议</th></tr></thead><tbody>{warning_rows}</tbody></table></div></section>
    <section class='card'><h2>样本质量指标</h2><div class='table'><table><thead><tr><th>样本</th><th>缺失率</th><th>杂合率</th><th>中位DP</th><th>中位GQ</th><th>AB异常</th><th>倍性不符</th></tr></thead><tbody>{sample_rows}</tbody></table></div></section>
    <section class='card'><h2>有效阈值与来源</h2><div class='table'><table><thead><tr><th>阈值</th><th>有效值</th><th>来源</th></tr></thead><tbody>{threshold_rows}</tbody></table></div></section>
    {population_section}
    <section class='card'><h2>自动修复安全策略</h2><p>原始VCF永不被静默覆盖。覆盖源文件、改写REF/ALT/GT、坐标转换、染色体批量重命名和删除文件均被定义为危险操作，执行前必须再次确认。</p></section>
    <script type='application/json' id='reportData'>{embedded}</script><script>function downloadJson(){const text=document.getElementById('reportData').textContent;const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([text],{type:'application/json;charset=utf-8'}));a.download='report_summary.json';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000)}</script></main></body></html>"""
    values = {
        "status": summary["status"], "score": summary["score"], "status_label": status_label,
        "name": html.escape(result["input"]["name"]), "generated": html.escape(result["generated_at"]),
        "records": "{:,}".format(summary["record_count"]), "samples": summary["sample_count"], "contigs": summary["contig_count"],
        "critical": summary["critical_count"], "warning": summary["warning_count"], "profile_name": html.escape(profile["name"]),
        "species": html.escape(profile.get("species_name") or profile["kingdom"]), "ploidy": profile["ploidy"],
        "gt_ploidy": profile.get("effective_gt_ploidy") or "未识别", "subgenomes": profile["subgenomes"],
        "scan_mode": "完整扫描" if scan["mode"] == "full" else "智能抽样", "evaluated": "{:,}".format(scan["evaluated_records"]),
        "elapsed": scan["elapsed_seconds"], "method_note": html.escape(scan["method_note"]),
        "type_chart": _svg_bars(sorted(site["variant_types"].items()), "SNP / INDEL / SV"),
        "missing_chart": _svg_bars([(x, site["missing_rate_bins"].get(x, 0)) for x in ["0–1%", "1–5%", "5–10%", "10–20%", ">20%"]], "位点缺失率"),
        "warning_rows": "".join(warning_rows) or "<tr><td colspan='5'>未触发告警</td></tr>", "sample_rows": "".join(sample_rows),
        "threshold_rows": "".join(threshold_rows), "population_section": population_section, "embedded": embedded,
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
        return Path(base) / "CallVCF" / "reports"

    def catalog(self):
        return {"profiles": profile_catalog(), "default_report_root": str(self.default_report_root())}

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
            report_path = run_dir / "report.html"
            json_path = run_dir / "report_summary.json"
            samples_path = run_dir / "sample_metrics.tsv"
            warnings_path = run_dir / "warnings.tsv"
            manifest_path = run_dir / "run_manifest.json"
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
            stat = Path(path).stat()
            manifest = {
                "run_id": job["id"], "input_path": path, "input_size": stat.st_size,
                "input_mtime_ns": stat.st_mtime_ns, "input_fingerprint": hashlib.sha256((path + str(stat.st_size) + str(stat.st_mtime_ns)).encode()).hexdigest(),
                "schema_version": result["schema_version"], "profile": result["profile"], "scan": result["scan"],
                "repair_actions_executed": [],
            }
            manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
            zip_path = run_dir / "CallVCF_QC_report.zip"
            population_files = []
            for name in result.get("population_analysis", {}).get("artifacts", []):
                item = run_dir / name
                if item.is_file() and item.parent == run_dir:
                    population_files.append(item)
            with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
                for item in [report_path, json_path, samples_path, warnings_path, manifest_path] + population_files:
                    archive.write(item, item.name)
            artifacts = []
            for item in [report_path, json_path, samples_path, warnings_path, manifest_path] + population_files + [zip_path]:
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
