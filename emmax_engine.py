#!/usr/bin/env python3
"""Managed PLINK -> EMMAX mixed-model workflow for GPA-Accelerator."""

import base64
import csv
import gzip
import heapq
import html
import itertools
import json
import math
import os
import platform
import random
import re
import shutil
import statistics
import subprocess
import threading
import time
import uuid
import zipfile
from pathlib import Path

from association_engine import (
    LOCATION_ALIASES, SAMPLE_ALIASES, TRAIT_ALIASES, VALUE_ALIASES, YEAR_ALIASES,
    _composites, _number, _read_table, _resolve_column, _resolve_traits,
)
from tool_manager import (
    _wsl_operational, deploy_bundled_emmax, resolve_emmax, resolve_plink,
)
from vcf_service import VCFError


def _safe_float(value, default, low, high):
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = default
    return max(low, min(high, number))


def _safe_int(value, default, low, high):
    try:
        number = int(value)
    except (TypeError, ValueError):
        number = default
    return max(low, min(high, number))


def _slug(value):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "trait")).strip("_.")
    return text[:80] or "trait"


def _chrom_key(value):
    text = str(value or "")
    cleaned = re.sub(r"^(?:chr|chromosome)", "", text, flags=re.I)
    match = re.search(r"\d+", cleaned)
    return (0, int(match.group())) if match else (1, cleaned.casefold())


def _load_selected_phenotypes(payload):
    path = Path(str(payload.get("phenotype_path") or "")).expanduser().resolve()
    if not path.is_file():
        raise VCFError("表型文件不存在：{}".format(path))
    mode = str(payload.get("match_mode") or "exact").strip().casefold()
    if mode not in {"exact", "prefix"}:
        raise VCFError("表型匹配方式必须是 exact 或 prefix")
    query = str(payload.get("trait_query") or "").strip()
    include_replicates = bool(payload.get("include_replicates"))
    headers, rows, source = _read_table(path, str(payload.get("sheet") or "").strip() or None)
    columns = payload.get("columns") or {}
    sample_column = _resolve_column(headers, columns.get("sample"), SAMPLE_ALIASES, True, "材料/样本列")
    trait_column = _resolve_column(headers, columns.get("trait"), TRAIT_ALIASES)
    value_column = _resolve_column(headers, columns.get("value"), VALUE_ALIASES)
    year_column = _resolve_column(headers, columns.get("year"), YEAR_ALIASES)
    location_column = _resolve_column(headers, columns.get("location"), LOCATION_ALIASES)
    year_filter = str(payload.get("year_filter") or "").strip().casefold()
    location_filter = str(payload.get("location_filter") or "").strip().casefold()
    metadata_columns = {x for x in (sample_column, trait_column, value_column, year_column, location_column) if x}
    collected = {}

    if trait_column and value_column:
        available = sorted({str(row.get(trait_column) or "").strip() for row in rows if str(row.get(trait_column) or "").strip()})
        selected = _resolve_traits(available, query, mode, include_replicates)
        chosen = {name.casefold(): name for name in selected}
        collected = {name: {} for name in selected}
        for row in rows:
            sample = str(row.get(sample_column) or "").strip()
            trait = str(row.get(trait_column) or "").strip()
            if not sample or trait.casefold() not in chosen:
                continue
            if year_filter and str(row.get(year_column) or "").strip().casefold() != year_filter:
                continue
            if location_filter and str(row.get(location_column) or "").strip().casefold() != location_filter:
                continue
            number = _number(row.get(value_column))
            if number is not None:
                collected[chosen[trait.casefold()]].setdefault(sample, []).append(number)
    else:
        candidates = [header for header in headers if header not in metadata_columns and any(_number(row.get(header)) is not None for row in rows)]
        selected = _resolve_traits(candidates, query, mode, include_replicates)
        collected = {name: {} for name in selected}
        for row in rows:
            sample = str(row.get(sample_column) or "").strip()
            if not sample:
                continue
            if year_filter and year_column and str(row.get(year_column) or "").strip().casefold() != year_filter:
                continue
            if location_filter and location_column and str(row.get(location_column) or "").strip().casefold() != location_filter:
                continue
            for trait in selected:
                number = _number(row.get(trait))
                if number is not None:
                    collected[trait].setdefault(sample, []).append(number)

    values = {
        trait: {sample: statistics.fmean(items) for sample, items in by_sample.items() if items}
        for trait, by_sample in collected.items() if by_sample
    }
    selected = [trait for trait in selected if trait in values]
    if not selected:
        raise VCFError("所选表型没有有效数值；请检查列映射、NA和筛选条件")
    analyses = [(trait, "separate", values[trait]) for trait in selected]
    if bool(payload.get("combine", mode == "prefix")) and len(selected) > 1:
        coverage = _safe_float(payload.get("min_trait_coverage"), 0.5, 0.1, 1.0)
        for name, composite in _composites(selected, values, coverage).items():
            analyses.append(("{}__{}".format(query or "TRAITS", name), "combined", composite))
    return {
        "path": path, "source": source, "query": query, "mode": mode,
        "selected_traits": selected, "analyses": analyses,
        "columns": {"sample": sample_column, "trait": trait_column, "value": value_column,
                    "year": year_column, "location": location_column},
    }


def _parse_fam(path):
    rows = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.strip().split()
            if len(fields) >= 2:
                rows.append((fields[0], fields[1]))
    if not rows:
        raise VCFError("PLINK 未生成有效 TFAM/FAM 样本列表")
    return rows


def _write_phenotype(path, fam_rows, values):
    by_fold = {sample.casefold(): value for sample, value in values.items()}
    matched = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        for family, individual in fam_rows:
            value = by_fold.get(individual.casefold(), by_fold.get(family.casefold()))
            if value is None:
                rendered = "NA"
            else:
                rendered = "{:.12g}".format(value); matched += 1
            handle.write("{}\t{}\t{}\n".format(family, individual, rendered))
    return matched


def _read_covariate_rows(path):
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise VCFError("协变量文件为空")
    rows = [re.split(r"[\t, ]+", line.strip()) for line in lines]
    header = rows[0]
    has_header = len(rows) > 1 and any(_number(item) is None for item in header[2:]) and all(_number(item) is not None for item in rows[1][2:])
    data = rows[1:] if has_header else rows
    if any(len(row) < 3 for row in data):
        raise VCFError("协变量文件至少需要 FID、IID、截距/协变量 三列")
    return data


def _prepare_covariates(source, destination, fam_rows):
    rows = _read_covariate_rows(source)
    by_id = {row[1].casefold(): row[2:] for row in rows}
    width = max(len(values) for values in by_id.values())
    missing = []
    with destination.open("w", encoding="utf-8", newline="") as handle:
        for family, individual in fam_rows:
            values = by_id.get(individual.casefold()) or by_id.get(family.casefold())
            if not values or len(values) != width or any(_number(value) is None for value in values):
                missing.append(individual); continue
            rendered = list(values)
            if not rendered or _number(rendered[0]) != 1.0:
                rendered.insert(0, "1")
            handle.write("{}\t{}\t{}\n".format(family, individual, "\t".join(rendered)))
    if missing:
        destination.unlink(missing_ok=True)
        raise VCFError("协变量缺少或含非数值：{}{}".format(", ".join(missing[:10]), "…" if len(missing) > 10 else ""))


def _parse_reml(path):
    if not path.is_file():
        return {"pseudo_heritability": None}
    values = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.strip().split()
        for field in reversed(fields):
            try:
                values.append(float(field)); break
            except ValueError:
                continue
    return {
        "log_likelihood_with_variance": values[0] if len(values) > 0 else None,
        "log_likelihood_without_variance": values[1] if len(values) > 1 else None,
        "delta": values[2] if len(values) > 2 else None,
        "genetic_variance": values[3] if len(values) > 3 else None,
        "residual_variance": values[4] if len(values) > 4 else None,
        "pseudo_heritability": values[5] if len(values) > 5 else None,
    }


def _reservoir_add(items, value, seen, limit, generator):
    if len(items) < limit:
        items.append(value)
    else:
        index = generator.randrange(seen)
        if index < limit:
            items[index] = value


def _process_ps(ps_path, tped_path, output_path, trait):
    sample, top_heap, chromosome_max = [], [], {}
    tested, invalid = 0, 0
    generator = random.Random(83173)
    with ps_path.open("r", encoding="utf-8", errors="replace") as ps, tped_path.open("r", encoding="utf-8", errors="replace") as tped, gzip.open(output_path, "wt", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["trait", "chrom", "pos", "marker", "beta", "se", "pvalue", "effect_allele_code"])
        for ps_line, tped_line in itertools.zip_longest(ps, tped, fillvalue=""):
            ps_fields, tped_fields = ps_line.strip().split(), tped_line.strip().split()
            if len(ps_fields) < 4:
                invalid += 1; continue
            try:
                beta, se, pvalue = map(float, ps_fields[1:4])
            except ValueError:
                invalid += 1; continue
            if not all(math.isfinite(value) for value in (beta, se, pvalue)) or pvalue < 0 or pvalue > 1:
                invalid += 1; continue
            pvalue = max(pvalue, 1e-300)
            marker = ps_fields[0]
            chrom = tped_fields[0] if len(tped_fields) >= 4 else ""
            try:
                pos = int(float(tped_fields[3])) if len(tped_fields) >= 4 else 0
            except ValueError:
                pos = 0
            tested += 1
            chromosome_max[chrom] = max(chromosome_max.get(chrom, 0), pos)
            row = (pvalue, marker, chrom, pos, beta, se)
            _reservoir_add(sample, row, tested, 100000, generator)
            if len(top_heap) < 200:
                heapq.heappush(top_heap, (-pvalue, marker, chrom, pos, beta, se))
            elif pvalue < -top_heap[0][0]:
                heapq.heapreplace(top_heap, (-pvalue, marker, chrom, pos, beta, se))
            writer.writerow([trait, chrom, pos, marker, "{:.12g}".format(beta), "{:.12g}".format(se), "{:.12g}".format(pvalue), "2"])
    if not tested:
        raise VCFError("EMMAX .ps 未包含可解析的关联结果")
    top = sorted([(-item[0], item[1], item[2], item[3], item[4], item[5]) for item in top_heap])
    pvalues = [item[0] for item in sample]
    median_p = statistics.median(pvalues)
    try:
        chi_median = statistics.NormalDist().inv_cdf(1.0 - median_p / 2.0) ** 2
        lambda_gc = chi_median / 0.4549364231195727
    except (ValueError, statistics.StatisticsError):
        lambda_gc = None
    return {"tested_markers": tested, "invalid_rows": invalid, "sample": sample, "top": top,
            "chromosome_max": chromosome_max, "lambda_gc": lambda_gc}


def _manhattan_svg(trait, sample, chromosome_max, tested):
    width, height, left, right, top, bottom = 1280, 500, 75, 1245, 58, 430
    chromosomes = sorted(chromosome_max, key=_chrom_key)
    offsets, cursor = {}, 0
    for chrom in chromosomes:
        offsets[chrom] = cursor
        cursor += max(1, chromosome_max[chrom]) + max(1, int(chromosome_max[chrom] * .03))
    total = max(1, cursor)
    max_y = max(6.0, min(30.0, max((-math.log10(max(item[0], 1e-300)) for item in sample), default=6.0)))
    threshold = -math.log10(0.05 / max(1, tested))
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="#fbfdfb"/>',
             f'<text x="35" y="32" font-size="22" font-weight="700" fill="#173c34">{html.escape(trait)} · EMMAX Manhattan</text>',
             f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#687c75"/><line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#687c75"/>']
    for tick in range(0, int(max_y) + 1, max(1, int(max_y // 6))):
        y = bottom - (bottom - top) * tick / max_y
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" stroke="#e3ebe6"/><text x="65" y="{y+4:.1f}" text-anchor="end" font-size="12">{tick}</text>')
    if threshold <= max_y:
        y = bottom - (bottom - top) * threshold / max_y
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" stroke="#bd4d5b" stroke-dasharray="7 5"><title>Bonferroni 0.05/{tested}</title></line>')
    colors = ("#267064", "#d18b35")
    for pvalue, marker, chrom, pos, beta, se in sample:
        x = left + (right - left) * (offsets.get(chrom, 0) + pos) / total
        score = min(max_y, -math.log10(max(pvalue, 1e-300)))
        y = bottom - (bottom - top) * score / max_y
        color = colors[chromosomes.index(chrom) % 2] if chrom in chromosomes else colors[0]
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.2" fill="{color}" fill-opacity=".72"><title>{html.escape(marker)} p={pvalue:.4g}</title></circle>')
    for chrom in chromosomes:
        center = offsets[chrom] + chromosome_max[chrom] / 2
        x = left + (right - left) * center / total
        parts.append(f'<text x="{x:.1f}" y="454" text-anchor="middle" font-size="11" fill="#415850">{html.escape(chrom)}</text>')
    parts += ['<text x="660" y="486" text-anchor="middle" font-size="14">Chromosome</text>', '<text x="18" y="245" transform="rotate(-90 18 245)" text-anchor="middle" font-size="14">−log10(P)</text>', '</svg>']
    return "".join(parts)


def _qq_svg(trait, pvalues, lambda_gc):
    observed = sorted(-math.log10(max(p, 1e-300)) for p in pvalues)
    count = len(observed)
    expected = sorted(-math.log10((index + .5) / count) for index in range(count)) if count else []
    maximum = max(1.0, min(30.0, max(expected + observed, default=1.0)))
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="640" height="560" viewBox="0 0 640 560"><rect width="100%" height="100%" fill="#fbfdfb"/>',
             f'<text x="30" y="34" font-size="21" font-weight="700" fill="#173c34">{html.escape(trait)} · EMMAX Q-Q</text>',
             f'<text x="610" y="34" text-anchor="end" font-size="13" fill="#51665f">λGC={"—" if lambda_gc is None else f"{lambda_gc:.3f}"}</text>',
             '<line x1="70" y1="485" x2="590" y2="485" stroke="#667a73"/><line x1="70" y1="485" x2="70" y2="65" stroke="#667a73"/>',
             '<line x1="70" y1="485" x2="590" y2="65" stroke="#b46a62" stroke-dasharray="6 5"/>']
    for x_value, y_value in zip(expected, observed):
        x = 70 + 520 * x_value / maximum; y = 485 - 420 * min(maximum, y_value) / maximum
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.3" fill="#2e7568" fill-opacity=".65"/>')
    parts += ['<text x="330" y="530" text-anchor="middle" font-size="14">Expected −log10(P)</text>', '<text x="22" y="275" transform="rotate(-90 22 275)" text-anchor="middle" font-size="14">Observed −log10(P)</text>', '</svg>']
    return "".join(parts)


class EmmaxJobManager:
    def __init__(self, service):
        self.service = service
        self._jobs = {}
        self._lock = threading.Lock()

    def _runtime(self):
        plink = resolve_plink()
        if not plink:
            raise VCFError("EMMAX 需要 PLINK 1.9 准备 TPED/TFAM；请先在本地软件资源中一键安装 PLINK")
        system = platform.system()
        if system == "Windows" and not _wsl_operational():
            raise VCFError("EMMAX 已包含在安装包中，但 Ubuntu/WSL 尚未完成首次初始化；请先打开一次 Ubuntu 并创建 Linux 用户")
        emmax = resolve_emmax()
        if not emmax:
            try:
                deploy_bundled_emmax(); emmax = resolve_emmax()
            except OSError as exc:
                raise VCFError("EMMAX 部署到本地工具目录失败：{}".format(exc))
        if not emmax:
            raise VCFError("EMMAX 随包运行文件部署失败")
        return {"plink": plink, "emmax": emmax["emmax"], "kin": emmax["kin"], "system": system}

    def start(self, payload):
        vcf = Path(str(payload.get("vcf_path") or payload.get("path") or "")).expanduser().resolve()
        if not vcf.is_file():
            raise VCFError("VCF 文件不存在：{}".format(vcf))
        if vcf.suffix.casefold() == ".bcf":
            raise VCFError("当前 EMMAX 一键流程接受 VCF/VCF.GZ；BCF 请先用安全修复执行器导出为 VCF.GZ")
        phenotype = _load_selected_phenotypes(payload)
        output_root = Path(str(payload.get("output_dir") or (Path.home() / "GPA-Accelerator" / "emmax-reports"))).expanduser().resolve()
        run_id = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:8]
        run_dir = output_root / run_id
        run_dir.mkdir(parents=True, exist_ok=False)
        estimated = max(2 * 1024 ** 3, vcf.stat().st_size * 15)
        free = shutil.disk_usage(run_dir).free
        if free < estimated:
            run_dir.rmdir()
            raise VCFError("EMMAX TPED 中间文件预计需要约 {:.1f} GB；当前输出磁盘仅余 {:.1f} GB".format(estimated / 1024 ** 3, free / 1024 ** 3))
        job = {
            "id": run_id, "status": "queued", "progress": 0, "stage": "queued",
            "message": "等待启动", "run_dir": str(run_dir), "created_at": time.time(),
            "payload": dict(payload), "vcf": str(vcf), "phenotype": phenotype,
            "estimated_workspace_bytes": estimated, "cancel": threading.Event(),
            "artifacts": [], "result": None, "error": None, "warnings": [],
        }
        with self._lock:
            self._jobs[run_id] = job
        threading.Thread(target=self._execute, args=(job,), daemon=True, name="emmax-{}".format(run_id)).start()
        return self._public(job)

    def _set(self, job, progress, stage, message):
        with self._lock:
            job.update(progress=progress, stage=stage, message=message, updated_at=time.time())

    def _run_command(self, command, cwd, log_path, cancel):
        command = [str(item) for item in command]
        with log_path.open("a", encoding="utf-8", errors="replace") as log:
            log.write("\n$ {}\n".format(" ".join(command))); log.flush()
            process = subprocess.Popen(command, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       text=True, encoding="utf-8", errors="replace",
                                       creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            while process.poll() is None:
                if cancel.is_set():
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                    raise VCFError("用户取消了 EMMAX 任务")
                time.sleep(.25)
            output = process.stdout.read() if process.stdout else ""
            log.write(output or "")
            if process.returncode:
                raise VCFError("外部程序退出码 {}：{}".format(process.returncode, (output or "无输出")[-1200:]))

    @staticmethod
    def _wsl_path(path):
        wsl = shutil.which("wsl.exe") or shutil.which("wsl")
        proc = subprocess.run([str(wsl), "-e", "wslpath", "-a", str(Path(path).resolve())],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                              encoding="utf-8", errors="replace", timeout=20,
                              creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        if proc.returncode or not proc.stdout.strip():
            raise VCFError("无法把 Windows 路径映射到 WSL：{}".format(proc.stderr.strip()))
        return proc.stdout.strip().splitlines()[-1]

    def _linux_command(self, runtime, executable, arguments):
        normalized = [str(item)[5:] if str(item).startswith("PATH:") else str(item) for item in arguments]
        if runtime["system"] != "Windows":
            return [executable] + normalized
        wsl = shutil.which("wsl.exe") or shutil.which("wsl")
        converted = []
        for original, normalized_item in zip(arguments, normalized):
            text = self._wsl_path(normalized_item) if str(original).startswith("PATH:") else normalized_item
            converted.append(text)
        return [str(wsl), "-e", self._wsl_path(executable)] + converted

    def _execute(self, job):
        run_dir = Path(job["run_dir"]); log = run_dir / "emmax_run.log"
        payload, phenotype = job["payload"], job["phenotype"]
        try:
            runtime = self._runtime()
            self._set(job, 5, "preflight", "运行环境和磁盘空间检查完成")
            assoc_call = _safe_float(payload.get("association_call_rate"), .50, .1, 1.0)
            assoc_maf = _safe_float(payload.get("association_maf"), .001, 0.0, .5)
            kin_call = _safe_float(payload.get("kinship_call_rate"), .95, .5, 1.0)
            kin_maf = _safe_float(payload.get("kinship_maf"), .01, 0.0, .5)
            prune_window = _safe_int(payload.get("prune_window"), 50, 5, 10000)
            prune_step = _safe_int(payload.get("prune_step"), 5, 1, 1000)
            prune_r2 = _safe_float(payload.get("prune_r2"), .2, .01, .99)
            base, prune, kin_prefix, assoc_prefix = (run_dir / name for name in ("genotype", "kinship_prune", "kinship", "association"))
            vcf = Path(job["vcf"])
            base_command = [runtime["plink"], "--vcf", vcf, "--double-id", "--allow-extra-chr",
                            "--biallelic-only", "strict", "--geno", 1.0 - assoc_call, "--maf", assoc_maf,
                            "--set-missing-var-ids", "@:#:$1:$2", "--make-bed", "--out", base]
            self._set(job, 12, "plink_import", "PLINK 正在导入并过滤 VCF")
            self._run_command(base_command, run_dir, log, job["cancel"])
            fam_rows = _parse_fam(base.with_suffix(".fam"))
            self._set(job, 25, "ld_pruning", "构建 kinship 标记面板并进行 LD 剪枝")
            try:
                self._run_command([runtime["plink"], "--bfile", base, "--allow-extra-chr", "--geno", 1.0 - kin_call,
                                   "--maf", kin_maf, "--indep-pairwise", prune_window, prune_step, prune_r2, "--out", prune],
                                  run_dir, log, job["cancel"])
            except VCFError as pruning_error:
                if job["cancel"].is_set():
                    raise
                self._run_command([runtime["plink"], "--bfile", base, "--allow-extra-chr", "--geno", 1.0 - kin_call,
                                   "--maf", kin_maf, "--write-snplist", "--out", prune], run_dir, log, job["cancel"])
                snplist = Path(str(prune) + ".snplist")
                if not snplist.is_file() or not snplist.read_text(encoding="utf-8", errors="replace").strip():
                    raise pruning_error
                shutil.copy2(snplist, Path(str(prune) + ".prune.in"))
                job["warnings"].append("LD剪枝未能执行（通常因为过滤后标记过少）；kinship改用通过call rate/MAF的未剪枝标记")
            self._set(job, 35, "tped", "生成 kinship 与关联 TPED/TFAM")
            self._run_command([runtime["plink"], "--bfile", base, "--allow-extra-chr", "--extract", str(prune) + ".prune.in",
                               "--recode", "12", "transpose", "--output-missing-genotype", "0", "--out", kin_prefix],
                              run_dir, log, job["cancel"])
            self._run_command([runtime["plink"], "--bfile", base, "--allow-extra-chr", "--recode", "12", "transpose",
                               "--output-missing-genotype", "0", "--out", assoc_prefix],
                              run_dir, log, job["cancel"])
            existing_kin = str(payload.get("kinship_path") or "").strip()
            if existing_kin:
                kin_file = Path(existing_kin).expanduser().resolve()
                if not kin_file.is_file():
                    raise VCFError("指定的 kinship 文件不存在：{}".format(kin_file))
            else:
                mode = str(payload.get("kinship_mode") or "BN").upper()
                args = ["-v"] + (["-s"] if mode == "IBS" else []) + ["-d", "10", "PATH:" + str(kin_prefix)]
                self._set(job, 48, "kinship", "EMMAX 正在计算 {} 亲缘矩阵".format(mode))
                self._run_command(self._linux_command(runtime, runtime["kin"], args), run_dir, log, job["cancel"])
                kin_file = Path(str(kin_prefix) + (".aIBS.kinf" if mode == "IBS" else ".aBN.kinf"))
                if not kin_file.is_file():
                    raise VCFError("EMMAX-kin 未生成预期亲缘矩阵")

            covariate = None
            source_cov = str(payload.get("covariate_path") or "").strip()
            if source_cov:
                source_cov = Path(source_cov).expanduser().resolve()
                if not source_cov.is_file():
                    raise VCFError("协变量文件不存在：{}".format(source_cov))
                covariate = run_dir / "aligned_covariates.txt"
                _prepare_covariates(source_cov, covariate, fam_rows)

            trait_results, top_rows, result_files = [], [], []
            analyses = phenotype["analyses"]
            for index, (trait, kind, values) in enumerate(analyses):
                if job["cancel"].is_set():
                    raise VCFError("用户取消了 EMMAX 任务")
                trait_dir = run_dir / "trait_{:02d}_{}".format(index + 1, _slug(trait))
                trait_dir.mkdir()
                pheno_path = trait_dir / "phenotype.txt"
                matched = _write_phenotype(pheno_path, fam_rows, values)
                if matched < 3:
                    trait_results.append({"trait": trait, "analysis_kind": kind, "status": "skipped", "matched_samples": matched,
                                          "warning": "匹配材料少于3个"})
                    continue
                out_prefix = trait_dir / "emmax"
                progress = 52 + int(38 * index / max(1, len(analyses)))
                self._set(job, progress, "association", "EMMAX 正在分析 {}（{}/{}）".format(trait, index + 1, len(analyses)))
                args = ["-v", "-d", "10", "-t", "PATH:" + str(assoc_prefix), "-p", "PATH:" + str(pheno_path),
                        "-k", "PATH:" + str(kin_file), "-o", "PATH:" + str(out_prefix)]
                if covariate:
                    args += ["-c", "PATH:" + str(covariate)]
                self._run_command(self._linux_command(runtime, runtime["emmax"], args), run_dir, log, job["cancel"])
                ps_path, reml_path = Path(str(out_prefix) + ".ps"), Path(str(out_prefix) + ".reml")
                if not ps_path.is_file():
                    raise VCFError("EMMAX 未生成 {} 的 .ps 结果".format(trait))
                full_result = trait_dir / "emmax_results.tsv.gz"
                parsed = _process_ps(ps_path, Path(str(assoc_prefix) + ".tped"), full_result, trait)
                reml = _parse_reml(reml_path)
                manhattan = _manhattan_svg(trait, parsed["sample"], parsed["chromosome_max"], parsed["tested_markers"])
                qq = _qq_svg(trait, [item[0] for item in parsed["sample"]], parsed["lambda_gc"])
                (trait_dir / "manhattan.svg").write_text(manhattan, encoding="utf-8")
                (trait_dir / "qq.svg").write_text(qq, encoding="utf-8")
                top = []
                for rank, (pvalue, marker, chrom, pos, beta, se) in enumerate(parsed["top"], 1):
                    row = {"trait": trait, "rank": rank, "chrom": chrom, "pos": pos, "marker": marker,
                           "beta": beta, "se": se, "pvalue": pvalue}
                    top.append(row); top_rows.append(row)
                trait_results.append({"trait": trait, "analysis_kind": kind, "status": "complete",
                                      "matched_samples": matched, "tested_markers": parsed["tested_markers"],
                                      "invalid_rows": parsed["invalid_rows"], "lambda_gc": parsed["lambda_gc"],
                                      "bonferroni_threshold": .05 / parsed["tested_markers"], "top_hits": top[:20],
                                      "reml": reml, "result_file": full_result.relative_to(run_dir).as_posix(),
                                      "manhattan": (trait_dir / "manhattan.svg").relative_to(run_dir).as_posix(),
                                      "qq": (trait_dir / "qq.svg").relative_to(run_dir).as_posix()})
                result_files.append(full_result)
                ps_path.unlink(missing_ok=True)

            if not any(item["status"] == "complete" for item in trait_results):
                raise VCFError("所有表型都因匹配材料不足而跳过")
            self._set(job, 92, "report", "生成 EMMAX 汇总、图形与下载文件")
            parameters = {"association_call_rate": assoc_call, "association_maf": assoc_maf,
                          "kinship_call_rate": kin_call, "kinship_maf": kin_maf,
                          "prune_window": prune_window, "prune_step": prune_step, "prune_r2": prune_r2,
                          "kinship_mode": str(payload.get("kinship_mode") or "BN").upper(),
                          "covariates": bool(covariate), "keep_intermediates": bool(payload.get("keep_intermediates"))}
            result = {"schema_version": "gpa-emmax-1.0", "run_id": job["id"], "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                      "input": {"vcf": job["vcf"], "phenotype": str(phenotype["path"]), "source": phenotype["source"]},
                      "selected_traits": phenotype["selected_traits"], "parameters": parameters,
                      "warnings": list(job["warnings"]),
                      "summary": {"vcf_samples": len(fam_rows), "analyses_requested": len(analyses),
                                  "analyses_completed": sum(item["status"] == "complete" for item in trait_results)},
                      "traits": trait_results,
                      "method_notes": ["亲缘矩阵标记先按 call rate/MAF 过滤并 LD 剪枝；不使用 HWE 过滤，避免对自交或群体结构明显的植物材料造成系统性删除。",
                                       "表型按 TFAM 样本顺序写入；缺失值记为 NA。", "EMMAX .ps 的 beta 对应 TPED 中编码为2的等位基因。",
                                       "Manhattan 图在标记过多时使用确定性抽样，完整结果保存在每个表型的 emmax_results.tsv.gz。"]}
            (run_dir / "emmax_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
            with (run_dir / "emmax_top_hits.tsv").open("w", encoding="utf-8-sig", newline="") as handle:
                fields = ["trait", "rank", "chrom", "pos", "marker", "beta", "se", "pvalue"]
                writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader(); writer.writerows(top_rows)
            report = self._report_html(result, run_dir)
            (run_dir / "emmax_report.html").write_text(report, encoding="utf-8")
            report_names = ["emmax_summary.json", "emmax_top_hits.tsv", "emmax_report.html", "emmax_run.log"]
            for trait in trait_results:
                if trait["status"] == "complete":
                    report_names.extend([trait["manhattan"], trait["qq"]])
            zip_path = run_dir / "GPA_Accelerator_EMMAX_report.zip"
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
                for name in report_names:
                    path = run_dir / name
                    if path.is_file():
                        archive.write(path, name)
            artifacts = self._artifacts(job["id"], run_dir, report_names + [zip_path.name] + [str(path.relative_to(run_dir)) for path in result_files])
            if not bool(payload.get("keep_intermediates")):
                protected = {path.resolve() for path in result_files}
                protected.update((run_dir / name).resolve() for name in report_names + [zip_path.name])
                for path in run_dir.rglob("*"):
                    if path.is_file() and path.resolve() not in protected:
                        path.unlink(missing_ok=True)
                for path in sorted((item for item in run_dir.rglob("*") if item.is_dir()), reverse=True):
                    try:
                        path.rmdir()
                    except OSError:
                        pass
            with self._lock:
                job.update(status="complete", progress=100, stage="complete", message="EMMAX 分析完成",
                           result=result, artifacts=artifacts, updated_at=time.time())
        except Exception as exc:
            cancelled = job["cancel"].is_set()
            with self._lock:
                job.update(status="cancelled" if cancelled else "failed", stage="cancelled" if cancelled else "failed",
                           message="任务已取消" if cancelled else "EMMAX 分析失败", error=str(exc), updated_at=time.time())

    def _artifacts(self, run_id, run_dir, names):
        artifacts = []
        for name in dict.fromkeys(names):
            path = (run_dir / name).resolve()
            if not path.is_file() or run_dir not in path.parents:
                continue
            base = "/api/emmax/artifact?run_id={}&name={}".format(run_id, str(path.relative_to(run_dir)).replace("\\", "/"))
            artifacts.append({"name": str(path.relative_to(run_dir)).replace("\\", "/"), "size": path.stat().st_size,
                              "url": base, "view_url": base + "&view=1"})
        return artifacts

    @staticmethod
    def _report_html(result, run_dir):
        cards, sections = [], []
        for item in result["traits"]:
            if item["status"] != "complete":
                continue
            top = item["top_hits"][0] if item["top_hits"] else None
            cards.append("<div><b>{}</b><span>{} markers · λGC {}</span><strong>{}</strong></div>".format(
                html.escape(item["trait"]), item["tested_markers"], "—" if item["lambda_gc"] is None else "{:.3f}".format(item["lambda_gc"]),
                "—" if not top else "P={:.4g}".format(top["pvalue"])))
            manhattan = base64.b64encode((run_dir / item["manhattan"]).read_bytes()).decode("ascii")
            qq = base64.b64encode((run_dir / item["qq"]).read_bytes()).decode("ascii")
            rows = "".join("<tr><td>{}</td><td>{}:{}</td><td>{}</td><td>{:.5g}</td><td>{:.4g}</td></tr>".format(
                html.escape(hit["marker"]), html.escape(str(hit["chrom"])), hit["pos"], hit["rank"], hit["beta"], hit["pvalue"])
                for hit in item["top_hits"][:20])
            sections.append("<section><h2>{}</h2><div class='plots'><img src='data:image/svg+xml;base64,{}'><img src='data:image/svg+xml;base64,{}'></div><table><thead><tr><th>Marker</th><th>位置</th><th>Rank</th><th>Beta</th><th>P</th></tr></thead><tbody>{}</tbody></table></section>".format(
                html.escape(item["trait"]), manhattan, qq, rows))
        template = """<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width'><title>GPA-Accelerator EMMAX报告</title><style>body{margin:0;background:#edf4f0;color:#173c34;font-family:Arial,'Microsoft YaHei',sans-serif}.page{max-width:1320px;margin:auto;padding:30px}.summary{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:12px}.summary div,section{background:#fff;border:1px solid #d5e2da;border-radius:18px;padding:20px;margin:15px 0}.summary b,.summary span,.summary strong{display:block}.summary strong{font-size:23px;margin-top:8px}.plots{display:grid;grid-template-columns:2fr 1fr;gap:12px}.plots img{width:100%;min-width:0}table{border-collapse:collapse;width:100%}th,td{padding:9px;border-bottom:1px solid #e0e9e4;text-align:left}th{background:#edf5f1}@media(max-width:850px){.plots{grid-template-columns:1fr}}</style></head><body><main class='page'><h1>EMMAX 全基因组混合模型报告</h1><p>__GENERATED__ · __VCF__</p><div class='summary'>__CARDS__</div>__SECTIONS__<section><h2>方法说明</h2><ul>__NOTES__</ul></section></main></body></html>"""
        return (template.replace("__GENERATED__", html.escape(result["generated_at"]))
                .replace("__VCF__", html.escape(result["input"]["vcf"]))
                .replace("__CARDS__", "".join(cards)).replace("__SECTIONS__", "".join(sections))
                .replace("__NOTES__", "".join("<li>{}</li>".format(html.escape(note)) for note in result["method_notes"])))

    def _public(self, job):
        return {key: job.get(key) for key in ("id", "status", "progress", "stage", "message", "run_dir", "error", "warnings", "result", "artifacts", "estimated_workspace_bytes")}

    def status(self, run_id):
        with self._lock:
            job = self._jobs.get(str(run_id or ""))
            if not job:
                raise VCFError("EMMAX 任务不存在或本地服务已重启")
            return self._public(job)

    def cancel(self, run_id):
        with self._lock:
            job = self._jobs.get(str(run_id or ""))
            if not job:
                raise VCFError("EMMAX 任务不存在")
            if job["status"] in {"complete", "failed", "cancelled"}:
                return self._public(job)
            job["cancel"].set(); job["message"] = "正在取消外部程序…"
            return self._public(job)

    def artifact(self, run_id, name):
        with self._lock:
            job = self._jobs.get(str(run_id or ""))
            if not job:
                raise VCFError("EMMAX 任务不存在或本地服务已重启")
            allowed = {item["name"] for item in job.get("artifacts") or []}
            if name not in allowed:
                raise VCFError("EMMAX 结果文件不存在")
            run_dir = Path(job["run_dir"]).resolve()
        path = (run_dir / name).resolve()
        if run_dir not in path.parents or not path.is_file():
            raise VCFError("禁止访问 EMMAX 任务目录外的文件")
        return path
