"""Two-phase, no-overwrite VCF repair executor backed by bcftools."""

import os
import html
import json
import shutil
import subprocess
import threading
import time
import uuid
from pathlib import Path

from vcf_service import VCFError
from quality_engine import QualityEvaluator, render_report


ACTIONS = {
    "index": {"name": "建立缺失索引", "risk": "safe", "needs_output": False},
    "sort_copy": {"name": "生成排序后的新副本", "risk": "safe", "needs_output": True},
    "normalize_copy": {"name": "生成标准化的新副本", "risk": "dangerous", "needs_output": True, "needs_reference": True},
    "fill_tags_copy": {"name": "补全AC/AN/AF/MAF/NS/F_MISSING/HWE/ExcHet统计标签的新副本", "risk": "dangerous", "needs_output": True},
    "deduplicate_copy": {"name": "移除完全重复记录的新副本", "risk": "dangerous", "needs_output": True},
    "subset_samples_copy": {"name": "保留指定样本的新副本", "risk": "dangerous", "needs_output": True, "needs_samples": True},
    "mask_genotypes_copy": {"name": "按条件将低质量GT设为缺失的新副本", "risk": "dangerous", "needs_output": True, "needs_expression": True},
    "filter_copy": {"name": "生成过滤后的新副本", "risk": "dangerous", "needs_output": True, "needs_expression": True},
    "quality_control_copy": {"name": "按物种与VCF情况执行质控并智能复评", "risk": "dangerous", "needs_output": True, "needs_qc_parameters": True},
}


def _quality_comparison(before, after, parameters):
    def status_counts(result):
        counts = {"pass": 0, "warning": 0, "critical": 0}
        for item in result.get("samples", []):
            key = item.get("sample_status") or "pass"
            counts[key] = counts.get(key, 0) + 1
        return counts

    def warning_map(result):
        return {(item.get("code"), item.get("scope"), item.get("target")): item for item in result.get("warnings", [])}

    before_summary = before.get("summary") or {}
    after_summary = after.get("summary") or {}
    before_warnings, after_warnings = warning_map(before), warning_map(after)
    before_records = int(before_summary.get("record_count") or 0)
    after_records = int(after_summary.get("record_count") or 0)
    removed = max(0, before_records - after_records)
    removed_fraction = removed / before_records if before_records else None
    resolved = [before_warnings[key] for key in before_warnings.keys() - after_warnings.keys()]
    introduced = [after_warnings[key] for key in after_warnings.keys() - before_warnings.keys()]
    persistent = [after_warnings[key] for key in after_warnings.keys() & before_warnings.keys()]
    delta = int(after_summary.get("score") or 0) - int(before_summary.get("score") or 0)
    recommendations = []
    if removed_fraction is not None and removed_fraction > .50:
        recommendations.append("本次移除了超过50%的位点，参数可能过严；建议先核对缺失率、MAF和位点质量阈值。")
    elif removed_fraction is not None and removed_fraction < .001 and delta <= 0:
        recommendations.append("位点集合和评分几乎未变化；当前VCF可能已完成相似质控，或所选阈值未命中数据。")
    if delta > 0:
        recommendations.append("综合评分提高{}分；仍需结合保留下来的告警判断是否继续处理。".format(delta))
    elif delta < 0:
        recommendations.append("综合评分下降{}分；建议回退到原VCF并放宽参数，质控副本不会覆盖原文件。".format(abs(delta)))
    if after_summary.get("critical_count"):
        recommendations.append("质控后仍有{}个关键告警；优先处理参考版本、文件结构或关键样本问题。".format(after_summary["critical_count"]))
    if after_summary.get("warning_count"):
        recommendations.append("质控后仍有{}个提醒；请在复评报告中查看证据与建议。".format(after_summary["warning_count"]))
    for item in persistent:
        advice = str(item.get("advice") or "").strip()
        if advice and advice not in recommendations:
            recommendations.append(advice)
        if len(recommendations) >= 8:
            break
    if not recommendations:
        recommendations.append("未发现必须继续质控的证据；建议保留原VCF与本次参数清单以保证可追溯。")
    return {
        "parameters": parameters,
        "before": {"score": before_summary.get("score"), "status": before_summary.get("status"), "records": before_records, "samples": before_summary.get("sample_count"), "critical": before_summary.get("critical_count"), "warning": before_summary.get("warning_count"), "sample_status": status_counts(before)},
        "after": {"score": after_summary.get("score"), "status": after_summary.get("status"), "records": after_records, "samples": after_summary.get("sample_count"), "critical": after_summary.get("critical_count"), "warning": after_summary.get("warning_count"), "sample_status": status_counts(after)},
        "change": {"score_delta": delta, "records_removed": removed, "records_removed_fraction": removed_fraction, "resolved_warning_count": len(resolved), "new_warning_count": len(introduced), "persistent_warning_count": len(persistent)},
        "resolved_warnings": resolved[:100], "new_warnings": introduced[:100], "persistent_warnings": persistent[:100],
        "further_recommendations": recommendations,
    }


def _render_quality_comparison(comparison):
    before, after, change = comparison["before"], comparison["after"], comparison["change"]
    ratio = "—" if change["records_removed_fraction"] is None else "{:.2%}".format(change["records_removed_fraction"])
    recs = "".join("<li>{}</li>".format(html.escape(str(item))) for item in comparison["further_recommendations"])
    params = "".join("<tr><td>{}</td><td>{}</td></tr>".format(html.escape(str(key)), html.escape(str(value if value is not None else "关闭"))) for key, value in comparison["parameters"].items())
    return """<!doctype html><html lang='zh-CN'><meta charset='utf-8'><title>GPA-Accelerator质控前后比较</title><style>body{{font-family:system-ui,'Microsoft YaHei',sans-serif;color:#17362d;background:#f2f7f4;margin:0;padding:28px}}.page{{max-width:1080px;margin:auto}}.grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:16px}}.card{{background:#fff;border:1px solid #d7e4dd;border-radius:16px;padding:20px}}h1,h2{{margin-top:0}}.big{{font-size:32px;font-weight:700}}.delta{{color:#2f725f}}table{{width:100%;border-collapse:collapse}}td{{padding:8px;border-bottom:1px solid #e4ece8}}small{{color:#62746d}}@media(max-width:700px){{.grid{{grid-template-columns:1fr}}}}</style><body><main class='page'><h1>GPA-Accelerator 质控前后智能复评</h1><div class='grid'><section class='card'><small>质控前评分</small><div class='big'>{before_score}/100</div><p>{before_records:,} 条记录</p></section><section class='card'><small>质控后评分</small><div class='big'>{after_score}/100</div><p>{after_records:,} 条记录</p></section><section class='card'><small>变化</small><div class='big delta'>{delta:+d} 分</div><p>移除 {removed:,} 条（{ratio}）</p></section></div><section class='card'><h2>告警变化</h2><p>已解决 {resolved} 项；新增 {new} 项；仍存在 {persistent} 项。</p><ul>{recs}</ul></section><section class='card'><h2>本次质控参数</h2><table>{params}</table></section></main></body></html>""".format(
        before_score=before.get("score"), before_records=before.get("records", 0), after_score=after.get("score"), after_records=after.get("records", 0),
        delta=int(change.get("score_delta") or 0), removed=change.get("records_removed", 0), ratio=ratio,
        resolved=change.get("resolved_warning_count", 0), new=change.get("new_warning_count", 0), persistent=change.get("persistent_warning_count", 0), recs=recs, params=params,
    )


class RepairExecutor:
    def __init__(self, service):
        self.service = service
        self.quality_evaluator = QualityEvaluator(service)
        self.bcftools = getattr(service, "bcftools", None)
        self.backend_mode = "native" if self.bcftools else None
        self.wsl = None
        self.refresh_backend()
        self._plans = {}
        self._jobs = {}
        self._lock = threading.Lock()

    def refresh_backend(self):
        if self.bcftools and self.backend_mode == "native":
            return self.catalog() if hasattr(self, "_lock") else None
        self.bcftools = None
        self.backend_mode = None
        self.wsl = None
        if not self.bcftools and os.name == "nt":
            candidate = shutil.which("wsl.exe") or shutil.which("wsl")
            if candidate:
                try:
                    proc = subprocess.run(
                        [candidate, "-e", "sh", "-lc", "command -v bcftools"],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace", timeout=5,
                        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
                    )
                    resolved = proc.stdout.strip().splitlines()
                    if proc.returncode == 0 and resolved:
                        self.wsl = candidate
                        self.bcftools = resolved[-1].strip()
                        self.backend_mode = "wsl"
                except Exception:
                    pass
        return self.catalog() if hasattr(self, "_lock") else None

    def catalog(self):
        return {
            "available": bool(self.bcftools),
            "backend": ("WSL · {}".format(self.bcftools) if self.backend_mode == "wsl" else self.bcftools),
            "backend_mode": self.backend_mode,
            "actions": [{"id": key, **value} for key, value in ACTIONS.items()],
            "policy": {
                "source_overwrite_supported": False,
                "existing_target_overwrite_supported": False,
                "dangerous_requires_exact_phrase": True,
                "partial_output_is_removed_on_failure": True,
            },
        }

    @staticmethod
    def _resolve_existing(path, label):
        value = Path(os.path.expandvars(os.path.expanduser(str(path or "").strip()))).resolve()
        if not value.is_file():
            raise VCFError("{}不存在：{}".format(label, value))
        return value

    def plan(self, payload):
        action_id = str(payload.get("action") or "")
        action = ACTIONS.get(action_id)
        if not action:
            raise VCFError("未知修复操作")
        if not self.bcftools:
            raise VCFError("当前本地后端未找到bcftools；修复执行器不可用，但质量评估和PLINK群体分析仍可使用")
        source = self._resolve_existing(payload.get("path"), "输入VCF")
        output = None
        if action.get("needs_output"):
            raw = str(payload.get("output_path") or "").strip()
            if not raw:
                raise VCFError("请选择新输出文件路径")
            output = Path(os.path.expandvars(os.path.expanduser(raw))).resolve()
            if output == source:
                raise VCFError("安全策略禁止覆盖输入VCF")
            if output.exists():
                raise VCFError("目标文件已存在；为防止覆盖，请更换新文件名")
            if not str(output).lower().endswith((".vcf.gz", ".vcf.bgz")):
                raise VCFError("新副本必须使用.vcf.gz或.vcf.bgz后缀")
            if not output.parent.is_dir():
                raise VCFError("输出目录不存在：{}".format(output.parent))
        reference = None
        if action.get("needs_reference"):
            reference = self._resolve_existing(payload.get("reference_path"), "参考基因组FASTA")
        expression = str(payload.get("expression") or "").strip()
        if action.get("needs_expression"):
            if not expression:
                raise VCFError("过滤表达式不能为空")
            if len(expression) > 2000 or "\x00" in expression:
                raise VCFError("过滤表达式无效或过长")
        samples = []
        if action.get("needs_samples"):
            raw_samples = payload.get("samples") or ""
            if isinstance(raw_samples, list):
                raw_samples = "\n".join(str(x) for x in raw_samples)
            samples = [x.strip() for x in str(raw_samples).replace(",", "\n").replace(";", "\n").splitlines() if x.strip()]
            samples = list(dict.fromkeys(samples))
            if not samples:
                raise VCFError("请至少输入一个要保留的样本ID")
            if len(samples) > 10000 or any(len(x) > 300 or "\x00" in x for x in samples):
                raise VCFError("样本列表过大或包含无效ID")

        qc_parameters = None
        quality_profile = payload.get("quality_profile") or {}
        quality_inputs = payload.get("quality_inputs") if isinstance(payload.get("quality_inputs"), dict) else {}
        if action.get("needs_qc_parameters"):
            raw_qc = payload.get("qc_parameters") or {}
            if not isinstance(raw_qc, dict):
                raise VCFError("质控参数格式无效")

            def optional_number(name, minimum, maximum):
                value = raw_qc.get(name)
                if value in (None, ""):
                    return None
                try:
                    number = float(value)
                except (TypeError, ValueError):
                    raise VCFError("质控参数{}不是有效数字".format(name))
                if number < minimum or number > maximum:
                    raise VCFError("质控参数{}应在{}到{}之间".format(name, minimum, maximum))
                return number

            qc_parameters = {
                "mode": "custom" if str(raw_qc.get("mode")) == "custom" else "profile_adaptive",
                "max_site_missing": optional_number("max_site_missing", 0, 1),
                "max_sample_missing": optional_number("max_sample_missing", 0, 1),
                "min_dp": optional_number("min_dp", 0, 100000),
                "max_dp": optional_number("max_dp", 0, 100000),
                "min_gq": optional_number("min_gq", 0, 999),
                "min_qual": optional_number("min_qual", 0, 1000000),
                "min_maf": optional_number("min_maf", 0, .5),
                "min_qd": optional_number("min_qd", 0, 100000),
                "min_mq": optional_number("min_mq", 0, 100000),
                "max_fs": optional_number("max_fs", 0, 1000000),
                "max_sor": optional_number("max_sor", 0, 1000000),
                "pass_only": bool(raw_qc.get("pass_only")),
                "biallelic_only": bool(raw_qc.get("biallelic_only")),
                "normalize": bool(raw_qc.get("normalize")),
                "deduplicate": bool(raw_qc.get("deduplicate")),
                "remove_samples": bool(raw_qc.get("remove_samples")),
            }
            if qc_parameters["max_site_missing"] is None:
                raise VCFError("位点缺失率上限不能为空")
            if qc_parameters["min_dp"] is not None and qc_parameters["max_dp"] is not None and qc_parameters["min_dp"] >= qc_parameters["max_dp"]:
                raise VCFError("GT最小DP必须小于最大DP")
            if qc_parameters["normalize"]:
                reference = self._resolve_existing(payload.get("reference_path"), "参考基因组FASTA")
            if qc_parameters["remove_samples"]:
                raw_excluded = payload.get("exclude_samples") or []
                if not isinstance(raw_excluded, list):
                    raw_excluded = str(raw_excluded).replace(",", "\n").splitlines()
                samples = list(dict.fromkeys(str(x).strip() for x in raw_excluded if str(x).strip()))
                if not samples:
                    raise VCFError("已启用低质量样本排除，但当前没有候选样本；请关闭该选项或重新评估")

        plan_id = uuid.uuid4().hex
        phrase = "确认执行-{}".format(plan_id[:8]) if action["risk"] == "dangerous" else None
        command = self._command_preview(action_id, source, output, reference, expression, samples)
        plan = {
            "id": plan_id, "action": action_id, "action_name": action["name"], "risk": action["risk"],
            "source": str(source), "output": str(output) if output else None,
            "reference": str(reference) if reference else None, "expression": expression or None,
            "samples": samples,
            "qc_parameters": qc_parameters,
            "quality_profile": quality_profile,
            "quality_inputs": quality_inputs,
            "command_preview": command, "confirmation_phrase": phrase,
            "created_at": time.time(), "expires_at": time.time() + 900,
            "source_size": source.stat().st_size, "source_mtime_ns": source.stat().st_mtime_ns,
            "status": "planned",
        }
        with self._lock:
            self._plans[plan_id] = plan
        return self._public_plan(plan)

    def _command_preview(self, action, source, output, reference, expression, samples):
        bcftools = "wsl.exe -e {}".format(self.bcftools) if self.backend_mode == "wsl" else str(self.bcftools)
        if action == "index":
            return [bcftools, "index", "-t", str(source)]
        if action == "sort_copy":
            return [bcftools, "sort", "-Oz", "-o", "<temporary-output>", str(source)]
        if action == "normalize_copy":
            return [bcftools, "norm", "-f", str(reference), "-m", "-any", "-Oz", "-o", "<temporary-output>", str(source)]
        if action == "fill_tags_copy":
            return [bcftools, "+fill-tags", str(source), "-Oz", "-o", "<temporary-output>", "--", "-t", "AC,AN,AF,MAF,NS,F_MISSING,HWE,ExcHet"]
        if action == "deduplicate_copy":
            return [bcftools, "norm", "-d", "exact", "-Oz", "-o", "<temporary-output>", str(source)]
        if action == "subset_samples_copy":
            return [bcftools, "view", "-s", ",".join(samples), "-Oz", "-o", "<temporary-output>", str(source)]
        if action == "mask_genotypes_copy":
            return [bcftools, "+setGT", str(source), "-Oz", "-o", "<temporary-output>", "--", "-t", "q", "-n", ".", "-i", expression]
        if action == "quality_control_copy":
            return [bcftools, "受控多阶段质控", "→", "掩蔽低质量GT", "→", "补全统计标签", "→", "过滤位点", "→", "建立CSI", "→", "前后复评"]
        return [bcftools, "view", "-i", expression, "-Oz", "-o", "<temporary-output>", str(source)]

    def _translate_wsl_path(self, value):
        proc = subprocess.run(
            [self.wsl, "-e", "wslpath", "-a", str(value)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace", timeout=5,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        if proc.returncode or not proc.stdout.strip():
            raise VCFError("WSL无法转换路径：{}".format(value))
        return proc.stdout.strip().splitlines()[-1]

    def _bcftools_command(self, arguments, path_indexes=()):
        arguments = [str(x) for x in arguments]
        if self.backend_mode != "wsl":
            return [str(self.bcftools)] + arguments
        for index in path_indexes:
            arguments[index] = self._translate_wsl_path(arguments[index])
        return [self.wsl, "-e", self.bcftools] + arguments

    def execute(self, payload):
        plan_id = str(payload.get("plan_id") or "")
        with self._lock:
            plan = self._plans.get(plan_id)
        if not plan:
            raise VCFError("修复计划不存在或本地服务已重启")
        if time.time() > plan["expires_at"]:
            raise VCFError("修复计划已过期，请重新生成并核对")
        if plan["status"] != "planned":
            raise VCFError("该计划已执行或正在执行，不能重复提交")
        if plan["risk"] == "dangerous" and str(payload.get("confirmation") or "").strip() != plan["confirmation_phrase"]:
            raise VCFError("二次确认短语不匹配")
        source = Path(plan["source"])
        if not source.is_file():
            raise VCFError("输入文件在计划生成后已不存在")
        source_stat = source.stat()
        if source_stat.st_size != plan["source_size"] or source_stat.st_mtime_ns != plan["source_mtime_ns"]:
            raise VCFError("输入文件在计划生成后发生变化；请重新生成并核对计划")
        output = Path(plan["output"]) if plan.get("output") else None
        if output is not None and output.exists():
            raise VCFError("目标文件在计划生成后已出现；已阻止覆盖")
        plan["status"] = "queued"
        job = {
            "id": plan_id, "status": "queued", "message": "等待执行", "progress": 0,
            "plan": self._public_plan(plan), "error": None, "result": None,
            "cancel": threading.Event(),
        }
        with self._lock:
            self._jobs[plan_id] = job
        thread = threading.Thread(target=self._run_job, args=(job, plan), daemon=True)
        job["thread"] = thread
        thread.start()
        return self._public_job(job)

    def _run_process(self, command, log, cancel):
        creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
        proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, creationflags=creationflags)
        while proc.poll() is None:
            if cancel.is_set():
                proc.kill(); proc.wait()
                raise VCFError("修复任务已取消")
            time.sleep(.15)
        if proc.returncode:
            raise VCFError("bcftools执行失败（退出码{}）".format(proc.returncode))

    def _run_job(self, job, plan):
        output = Path(plan["output"]) if plan.get("output") else None
        partial = None
        partial_index = None
        qc_intermediates = []
        try:
            plan["status"] = "running"
            job["status"] = "running"; job["message"] = "bcftools正在执行"; job["progress"] = 15
            source = Path(plan["source"])
            audit_root = Path(os.environ.get("LOCALAPPDATA") or (Path.home() / ".local" / "share")) / "GPA-Accelerator" / "repairs" / plan["id"]
            audit_root.mkdir(parents=True, exist_ok=False)
            log_path = audit_root / "repair.log"
            with log_path.open("wb") as log:
                if plan["action"] == "index":
                    suffix = ".csi"
                    final_index = Path(str(source) + suffix)
                    alternate = Path(str(source) + ".tbi")
                    if final_index.exists() or alternate.exists():
                        raise VCFError("索引已存在；安全策略不会覆盖现有索引")
                    partial_index = audit_root / (source.name + suffix + ".partial")
                    self._run_process(self._bcftools_command(
                        ["index", "-c", "-o", str(partial_index), str(source)], path_indexes=(3, 4)
                    ), log, job["cancel"])
                    if final_index.exists() or alternate.exists():
                        raise VCFError("执行期间出现索引文件；已阻止覆盖")
                    os.replace(str(partial_index), str(final_index))
                    job["result"] = {"index_path": str(final_index), "audit_log": str(log_path)}
                else:
                    partial = output.parent / (".{}.{}.partial.vcf.gz".format(output.stem, plan["id"][:8]))
                    if partial.exists():
                        raise VCFError("临时输出路径已存在")
                    if plan["action"] == "quality_control_copy":
                        job["message"] = "正在执行多阶段质控"
                        qc_intermediates = self._run_quality_control(plan, partial, audit_root, log, job["cancel"])
                    else:
                        command = self._build_command(plan, partial)
                        self._run_process(command, log, job["cancel"])
                    job["progress"] = 75; job["message"] = "正在校验并建立新文件索引"
                    partial_index = Path(str(partial) + ".csi")
                    self._run_process(self._bcftools_command(
                        ["index", "-c", "-o", str(partial_index), str(partial)], path_indexes=(3, 4)
                    ), log, job["cancel"])
                    if not partial.is_file() or not partial_index.is_file() or not partial.stat().st_size:
                        raise VCFError("bcftools未生成完整输出和索引")
                    if output.exists() or Path(str(output) + ".tbi").exists() or Path(str(output) + ".csi").exists():
                        raise VCFError("执行期间目标路径已出现；已阻止覆盖")
                    before_stats = audit_root / "before.bcftools.stats.txt"
                    after_stats = audit_root / "after.bcftools.stats.txt"
                    self._run_to_file(self._bcftools_command(["stats", str(source)], path_indexes=(1,)), before_stats, log, job["cancel"])
                    self._run_to_file(self._bcftools_command(["stats", str(partial)], path_indexes=(1,)), after_stats, log, job["cancel"])
                    comparison = None
                    comparison_path = None
                    comparison_report = None
                    after_report = None
                    if plan["action"] == "quality_control_copy":
                        job["progress"] = 86; job["message"] = "正在重新评估质控前后变化"
                        quality_inputs = plan.get("quality_inputs") or {}
                        quality_config = {
                            "profile": plan.get("quality_profile") or {},
                            "scan_mode": "smart", "target_records": 20000,
                            "reference_path": plan.get("reference") or quality_inputs.get("reference_path") or "",
                            "sample_meta_path": quality_inputs.get("sample_meta_path") or "",
                            "region_bed_path": quality_inputs.get("region_bed_path") or "",
                        }
                        before_quality = self.quality_evaluator.evaluate(str(source), quality_config, cancel_event=job["cancel"])
                        after_quality = self.quality_evaluator.evaluate(str(partial), quality_config, cancel_event=job["cancel"])
                        after_quality["input"]["path"] = str(output)
                        after_quality["input"]["name"] = output.name
                        comparison = _quality_comparison(before_quality, after_quality, plan.get("qc_parameters") or {})
                        comparison_path = audit_root / "qc_before_after_comparison.json"
                        comparison_report = audit_root / "qc_before_after_comparison.html"
                        after_report = audit_root / "qc_after_report.html"
                        comparison_path.write_text(json.dumps(comparison, ensure_ascii=False, indent=2), encoding="utf-8")
                        comparison_report.write_text(_render_quality_comparison(comparison), encoding="utf-8")
                        after_report.write_text(render_report(after_quality), encoding="utf-8")
                    manifest_path = audit_root / "repair_manifest.json"
                    manifest = {
                        "plan_id": plan["id"], "action": plan["action"], "action_name": plan["action_name"],
                        "source": str(source), "source_size": source.stat().st_size,
                        "output": str(output), "output_size": partial.stat().st_size,
                        "backend": self.catalog().get("backend"), "completed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                        "before_stats": str(before_stats), "after_stats": str(after_stats),
                        "qc_parameters": plan.get("qc_parameters"),
                        "qc_comparison": str(comparison_path) if comparison_path else None,
                        "source_overwritten": False,
                    }
                    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
                    os.replace(str(partial), str(output))
                    os.replace(str(partial_index), str(output) + ".csi")
                    for item in qc_intermediates:
                        try:
                            Path(item).unlink()
                        except OSError:
                            pass
                    job["result"] = {"output_path": str(output), "index_path": str(output) + ".csi", "audit_log": str(log_path), "before_stats": str(before_stats), "after_stats": str(after_stats), "manifest": str(manifest_path), "qc_comparison": comparison, "qc_comparison_path": str(comparison_path) if comparison_path else None, "qc_comparison_report": str(comparison_report) if comparison_report else None, "qc_after_report": str(after_report) if after_report else None}
            job["status"] = "complete"; job["message"] = "安全修复已完成"; job["progress"] = 100
            plan["status"] = "complete"
        except Exception as exc:
            for item in (partial, partial_index, *qc_intermediates):
                if item and Path(item).is_file():
                    try:
                        Path(item).unlink()
                    except OSError:
                        pass
            if "audit_root" in locals() and Path(audit_root).is_dir():
                for item in Path(audit_root).glob(".0*_*.vcf.gz"):
                    try:
                        item.unlink()
                    except OSError:
                        pass
            job["status"] = "cancelled" if job["cancel"].is_set() else "failed"
            job["message"] = str(exc); job["error"] = str(exc)
            plan["status"] = job["status"]

    def _run_quality_control(self, plan, partial, audit_root, log, cancel):
        parameters = plan.get("qc_parameters") or {}
        current = Path(plan["source"])
        intermediates = []

        def stage(name):
            path = audit_root / (".{}.vcf.gz".format(name))
            intermediates.append(path)
            return path

        if parameters.get("normalize"):
            output = stage("01_normalized")
            args = ["norm", "-f", plan["reference"], "-m", "-any"]
            if parameters.get("deduplicate"):
                args.extend(["-d", "exact"])
            args.extend(["-Oz", "-o", str(output), str(current)])
            self._run_process(self._bcftools_command(args, path_indexes=(2, len(args) - 2, len(args) - 1)), log, cancel)
            current = output
        elif parameters.get("deduplicate"):
            output = stage("01_deduplicated")
            args = ["norm", "-d", "exact", "-Oz", "-o", str(output), str(current)]
            self._run_process(self._bcftools_command(args, path_indexes=(len(args) - 2, len(args) - 1)), log, cancel)
            current = output

        if parameters.get("remove_samples") and plan.get("samples"):
            output = stage("02_sample_subset")
            args = ["view", "-s", "^" + ",".join(plan["samples"]), "-Oz", "-o", str(output), str(current)]
            self._run_process(self._bcftools_command(args, path_indexes=(len(args) - 2, len(args) - 1)), log, cancel)
            current = output

        mask_clauses = []
        if parameters.get("min_dp") is not None:
            mask_clauses.append("FMT/DP<{}".format(parameters["min_dp"]))
        if parameters.get("max_dp") is not None:
            mask_clauses.append("FMT/DP>{}".format(parameters["max_dp"]))
        if parameters.get("min_gq") is not None:
            mask_clauses.append("FMT/GQ<{}".format(parameters["min_gq"]))
        if mask_clauses:
            output = stage("03_masked_gt")
            expression = " || ".join(mask_clauses)
            args = ["+setGT", str(current), "-Oz", "-o", str(output), "--", "-t", "q", "-n", ".", "-i", expression]
            self._run_process(self._bcftools_command(args, path_indexes=(1, 4)), log, cancel)
            current = output

        tagged = stage("04_filled_tags")
        args = ["+fill-tags", str(current), "-Oz", "-o", str(tagged), "--", "-t", "AC,AN,AF,MAF,NS,F_MISSING,HWE,ExcHet"]
        self._run_process(self._bcftools_command(args, path_indexes=(1, 4)), log, cancel)
        current = tagged

        clauses = ["F_MISSING<={}".format(parameters["max_site_missing"])]
        optional_rules = (("min_qual", "QUAL", ">="), ("min_qd", "QD", ">="), ("min_mq", "MQ", ">="), ("max_fs", "FS", "<="), ("max_sor", "SOR", "<="), ("min_maf", "MAF", ">="))
        for parameter, field, operator in optional_rules:
            value = parameters.get(parameter)
            if value is not None and not (parameter == "min_maf" and value <= 0):
                clauses.append("({field}=\".\" || {field}{operator}{value})".format(field=field, operator=operator, value=value))
        args = ["view", "-i", " && ".join(clauses)]
        if parameters.get("pass_only"):
            args.extend(["-f", "PASS"])
        if parameters.get("biallelic_only"):
            args.extend(["-m", "2", "-M", "2"])
        args.extend(["-Oz", "-o", str(partial), str(current)])
        self._run_process(self._bcftools_command(args, path_indexes=(len(args) - 2, len(args) - 1)), log, cancel)
        return intermediates

    def _build_command(self, plan, partial):
        source = plan["source"]
        if plan["action"] == "sort_copy":
            return self._bcftools_command(["sort", "-Oz", "-o", str(partial), source], path_indexes=(3, 4))
        if plan["action"] == "normalize_copy":
            return self._bcftools_command(
                ["norm", "-f", plan["reference"], "-m", "-any", "-Oz", "-o", str(partial), source],
                path_indexes=(2, 7, 8),
            )
        if plan["action"] == "fill_tags_copy":
            return self._bcftools_command(
                ["+fill-tags", source, "-Oz", "-o", str(partial), "--", "-t", "AC,AN,AF,MAF,NS,F_MISSING,HWE,ExcHet"],
                path_indexes=(1, 4),
            )
        if plan["action"] == "deduplicate_copy":
            return self._bcftools_command(["norm", "-d", "exact", "-Oz", "-o", str(partial), source], path_indexes=(5, 6))
        if plan["action"] == "subset_samples_copy":
            return self._bcftools_command(["view", "-s", ",".join(plan["samples"]), "-Oz", "-o", str(partial), source], path_indexes=(5, 6))
        if plan["action"] == "mask_genotypes_copy":
            return self._bcftools_command(["+setGT", source, "-Oz", "-o", str(partial), "--", "-t", "q", "-n", ".", "-i", plan["expression"]], path_indexes=(1, 4))
        return self._bcftools_command(
            ["view", "-i", plan["expression"], "-Oz", "-o", str(partial), source],
            path_indexes=(5, 6),
        )

    def _run_to_file(self, command, destination, log, cancel):
        creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
        with Path(destination).open("wb") as output:
            proc = subprocess.Popen(command, stdout=output, stderr=log, creationflags=creationflags)
            while proc.poll() is None:
                if cancel.is_set():
                    proc.kill(); proc.wait()
                    raise VCFError("修复任务已取消")
                time.sleep(.15)
        if proc.returncode:
            raise VCFError("bcftools stats审计失败（退出码{}）".format(proc.returncode))

    def status(self, plan_id):
        with self._lock:
            job = self._jobs.get(str(plan_id))
        if not job:
            raise VCFError("修复任务不存在")
        return self._public_job(job)

    def cancel(self, plan_id):
        with self._lock:
            job = self._jobs.get(str(plan_id))
        if not job:
            raise VCFError("修复任务不存在")
        job["cancel"].set(); job["message"] = "正在取消"
        return self._public_job(job)

    @staticmethod
    def _public_plan(plan):
        return {key: value for key, value in plan.items() if key not in {"created_at", "expires_at"}}

    @staticmethod
    def _public_job(job):
        return {key: value for key, value in job.items() if key not in {"thread", "cancel"}}
