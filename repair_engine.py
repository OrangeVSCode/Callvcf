"""Two-phase, no-overwrite VCF repair executor backed by bcftools."""

import os
import shutil
import subprocess
import threading
import time
import uuid
from pathlib import Path

from vcf_service import VCFError


ACTIONS = {
    "index": {"name": "建立缺失索引", "risk": "safe", "needs_output": False},
    "sort_copy": {"name": "生成排序后的新副本", "risk": "safe", "needs_output": True},
    "normalize_copy": {"name": "生成标准化的新副本", "risk": "dangerous", "needs_output": True, "needs_reference": True},
    "filter_copy": {"name": "生成过滤后的新副本", "risk": "dangerous", "needs_output": True, "needs_expression": True},
}


class RepairExecutor:
    def __init__(self, service):
        self.bcftools = getattr(service, "bcftools", None)
        self.backend_mode = "native" if self.bcftools else None
        self.wsl = None
        if not self.bcftools and os.name == "nt":
            candidate = shutil.which("wsl.exe") or shutil.which("wsl")
            if candidate:
                try:
                    proc = subprocess.run(
                        [candidate, "-e", "sh", "-lc", "command -v bcftools"],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=5,
                        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
                    )
                    resolved = proc.stdout.strip().splitlines()
                    if proc.returncode == 0 and resolved:
                        self.wsl = candidate
                        self.bcftools = resolved[-1].strip()
                        self.backend_mode = "wsl"
                except Exception:
                    pass
        self._plans = {}
        self._jobs = {}
        self._lock = threading.Lock()

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

        plan_id = uuid.uuid4().hex
        phrase = "确认执行-{}".format(plan_id[:8]) if action["risk"] == "dangerous" else None
        command = self._command_preview(action_id, source, output, reference, expression)
        plan = {
            "id": plan_id, "action": action_id, "action_name": action["name"], "risk": action["risk"],
            "source": str(source), "output": str(output) if output else None,
            "reference": str(reference) if reference else None, "expression": expression or None,
            "command_preview": command, "confirmation_phrase": phrase,
            "created_at": time.time(), "expires_at": time.time() + 900,
            "source_size": source.stat().st_size, "source_mtime_ns": source.stat().st_mtime_ns,
            "status": "planned",
        }
        with self._lock:
            self._plans[plan_id] = plan
        return self._public_plan(plan)

    def _command_preview(self, action, source, output, reference, expression):
        bcftools = "wsl.exe -e {}".format(self.bcftools) if self.backend_mode == "wsl" else str(self.bcftools)
        if action == "index":
            return [bcftools, "index", "-t", str(source)]
        if action == "sort_copy":
            return [bcftools, "sort", "-Oz", "-o", "<temporary-output>", str(source)]
        if action == "normalize_copy":
            return [bcftools, "norm", "-f", str(reference), "-m", "-any", "-Oz", "-o", "<temporary-output>", str(source)]
        return [bcftools, "view", "-i", expression, "-Oz", "-o", "<temporary-output>", str(source)]

    def _translate_wsl_path(self, value):
        proc = subprocess.run(
            [self.wsl, "-e", "wslpath", "-a", str(value)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=5,
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
        try:
            plan["status"] = "running"
            job["status"] = "running"; job["message"] = "bcftools正在执行"; job["progress"] = 15
            source = Path(plan["source"])
            audit_root = Path(os.environ.get("LOCALAPPDATA") or (Path.home() / ".local" / "share")) / "CallVCF" / "repairs" / plan["id"]
            audit_root.mkdir(parents=True, exist_ok=False)
            log_path = audit_root / "repair.log"
            with log_path.open("wb") as log:
                if plan["action"] == "index":
                    suffix = ".tbi"
                    final_index = Path(str(source) + suffix)
                    alternate = Path(str(source) + ".csi")
                    if final_index.exists() or alternate.exists():
                        raise VCFError("索引已存在；安全策略不会覆盖现有索引")
                    partial_index = audit_root / (source.name + suffix + ".partial")
                    self._run_process(self._bcftools_command(
                        ["index", "-t", "-o", str(partial_index), str(source)], path_indexes=(3, 4)
                    ), log, job["cancel"])
                    if final_index.exists() or alternate.exists():
                        raise VCFError("执行期间出现索引文件；已阻止覆盖")
                    os.replace(str(partial_index), str(final_index))
                    job["result"] = {"index_path": str(final_index), "audit_log": str(log_path)}
                else:
                    partial = output.parent / (".{}.{}.partial.vcf.gz".format(output.stem, plan["id"][:8]))
                    if partial.exists():
                        raise VCFError("临时输出路径已存在")
                    command = self._build_command(plan, partial)
                    self._run_process(command, log, job["cancel"])
                    job["progress"] = 75; job["message"] = "正在校验并建立新文件索引"
                    self._run_process(self._bcftools_command(
                        ["index", "-t", str(partial)], path_indexes=(2,)
                    ), log, job["cancel"])
                    partial_index = Path(str(partial) + ".tbi")
                    if not partial.is_file() or not partial_index.is_file() or not partial.stat().st_size:
                        raise VCFError("bcftools未生成完整输出和索引")
                    if output.exists() or Path(str(output) + ".tbi").exists():
                        raise VCFError("执行期间目标路径已出现；已阻止覆盖")
                    os.replace(str(partial), str(output))
                    os.replace(str(partial_index), str(output) + ".tbi")
                    job["result"] = {"output_path": str(output), "index_path": str(output) + ".tbi", "audit_log": str(log_path)}
            job["status"] = "complete"; job["message"] = "安全修复已完成"; job["progress"] = 100
            plan["status"] = "complete"
        except Exception as exc:
            for item in (partial, partial_index):
                if item and Path(item).is_file():
                    try:
                        Path(item).unlink()
                    except OSError:
                        pass
            job["status"] = "cancelled" if job["cancel"].is_set() else "failed"
            job["message"] = str(exc); job["error"] = str(exc)
            plan["status"] = job["status"]

    def _build_command(self, plan, partial):
        source = plan["source"]
        if plan["action"] == "sort_copy":
            return self._bcftools_command(["sort", "-Oz", "-o", str(partial), source], path_indexes=(3, 4))
        if plan["action"] == "normalize_copy":
            return self._bcftools_command(
                ["norm", "-f", plan["reference"], "-m", "-any", "-Oz", "-o", str(partial), source],
                path_indexes=(2, 7, 8),
            )
        return self._bcftools_command(
            ["view", "-i", plan["expression"], "-Oz", "-o", str(partial), source],
            path_indexes=(5, 6),
        )

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
