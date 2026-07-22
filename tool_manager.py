#!/usr/bin/env python3
"""Install and discover optional third-party command-line tools."""

import hashlib
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

from vcf_service import VCFError


PLINK_VERSION = "1.9 beta 7.11 (2025-08-19)"
PLINK_URLS = {
    "Windows": "https://s3.amazonaws.com/plink1-assets/plink_win64_20250819.zip",
    "Linux": "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip",
}
LDBLOCKSHOW_URL = "https://codeload.github.com/hewm2008/LDBlockShow/zip/refs/heads/main"
EMMAX_VERSION = "emmax-intel64-20120205 / binary distribution 20120210"
EMMAX_ARCHIVE_SHA256 = "E2A582851BA1BE908757D4EF436E98AD76664A0C55E00D13E55FA35FE2BA54DD"
RESOURCE_ROOT = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))


def tool_root():
    override = os.environ.get("GPA_ACCELERATOR_TOOL_DIR") or os.environ.get("CALLVCF_TOOL_DIR")
    if override:
        return Path(override).expanduser().resolve()
    if platform.system() == "Windows":
        base = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local")
    else:
        base = Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share")
    current = (base / "GPA-Accelerator" / "tools").resolve()
    legacy = (base / "CallVCF" / "tools").resolve()
    # Reuse existing managed tools after the product rename; new installations use the new root.
    return legacy if legacy.exists() and not current.exists() else current


def _download(url, destination, max_bytes=200 * 1024 * 1024):
    request = urllib.request.Request(url, headers={"User-Agent": "GPA-Accelerator/1.0"})
    digest = hashlib.sha256()
    total = 0
    with urllib.request.urlopen(request, timeout=90) as response, destination.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            if total > max_bytes:
                raise VCFError("工具下载超过 200 MB 安全限制")
            digest.update(chunk)
            handle.write(chunk)
    return total, digest.hexdigest()


def _safe_extract(archive, destination):
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as handle:
        for member in handle.infolist():
            target = (destination / member.filename).resolve()
            if destination.resolve() not in target.parents and target != destination.resolve():
                raise VCFError("工具压缩包包含不安全路径")
        handle.extractall(destination)


def _metadata_path(name):
    return tool_root() / name / "callvcf-tool.json"


def _write_metadata(name, payload):
    path = _metadata_path(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def resolve_plink(explicit=None):
    if explicit and Path(str(explicit)).expanduser().is_file():
        return str(Path(str(explicit)).expanduser().resolve())
    filename = "plink.exe" if platform.system() == "Windows" else "plink"
    installed = tool_root() / "plink" / filename
    if installed.is_file():
        return str(installed)
    return shutil.which("plink")


def resolve_ldblockshow(explicit=None):
    if explicit and Path(str(explicit)).expanduser().is_file():
        return str(Path(str(explicit)).expanduser().resolve())
    installed = tool_root() / "ldblockshow" / "bin" / "LDBlockShow"
    if installed.is_file():
        return str(installed)
    return shutil.which("LDBlockShow")


def bundled_emmax_dir():
    return (RESOURCE_ROOT / "vendor" / "emmax").resolve()


def deploy_bundled_emmax():
    source = bundled_emmax_dir()
    required = (source / "emmax-intel64", source / "emmax-kin-intel64", source / "LICENSE.txt")
    if not all(path.is_file() for path in required):
        raise VCFError("安装包中缺少 EMMAX 运行文件；请重新下载完整安装包")
    destination = tool_root() / "emmax"
    destination.mkdir(parents=True, exist_ok=True)
    for source_path in required + (source / "README.txt",):
        if source_path.is_file():
            target = destination / source_path.name
            if not target.is_file() or target.stat().st_size != source_path.stat().st_size:
                shutil.copy2(source_path, target)
            if not target.suffix:
                target.chmod(target.stat().st_mode | stat.S_IEXEC)
    _write_metadata("emmax", {
        "source": "https://csg.sph.umich.edu/kang/emmax/download/index.html",
        "distribution": "emmax-intel-binary-20120210.tar.gz",
        "version": EMMAX_VERSION, "archive_sha256": EMMAX_ARCHIVE_SHA256,
        "license": "MIT", "platform": "Ubuntu x86_64 via WSL on Windows",
    })
    return destination


def resolve_emmax():
    installed = tool_root() / "emmax" / "emmax-intel64"
    kin = tool_root() / "emmax" / "emmax-kin-intel64"
    if installed.is_file() and kin.is_file():
        return {"emmax": str(installed), "kin": str(kin)}
    if platform.system() != "Windows":
        command = shutil.which("emmax-intel64") or shutil.which("emmax")
        kin_command = shutil.which("emmax-kin-intel64") or shutil.which("emmax-kin")
        if command and kin_command:
            return {"emmax": command, "kin": kin_command}
    return None


def _probe(command, args):
    if not command:
        return None
    try:
        proc = subprocess.run([str(command)] + list(args), stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True, timeout=8)
        text = (proc.stdout + "\n" + proc.stderr).strip()
        return text.splitlines()[0][:300] if text else "可执行文件已找到"
    except Exception as exc:
        return "已找到，但版本检测失败：{}".format(exc)


def _wsl_operational():
    executable = shutil.which("wsl.exe") or shutil.which("wsl")
    if not executable:
        return False
    try:
        proc = subprocess.run([str(executable), "-e", "sh", "-lc", "true"],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=12)
        return proc.returncode == 0
    except Exception:
        return False


def _wsl_bcftools_path():
    if platform.system() != "Windows":
        return None
    executable = shutil.which("wsl.exe") or shutil.which("wsl")
    if not executable:
        return None
    try:
        proc = subprocess.run(
            [str(executable), "-e", "sh", "-lc", "command -v bcftools"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace", timeout=12,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        values = proc.stdout.strip().splitlines()
        return values[-1].strip() if proc.returncode == 0 and values else None
    except Exception:
        return None
def tools_status():
    system = platform.system()
    plink = resolve_plink()
    ldblockshow = resolve_ldblockshow()
    wsl_installed = _wsl_operational() if system == "Windows" else False
    native_bcftools = shutil.which("bcftools")
    wsl_bcftools = _wsl_bcftools_path()
    emmax = resolve_emmax()
    bundled_emmax = all((bundled_emmax_dir() / name).is_file() for name in ("emmax-intel64", "emmax-kin-intel64", "LICENSE.txt"))
    if bundled_emmax and not emmax and getattr(sys, "frozen", False):
        try:
            deploy_bundled_emmax()
            emmax = resolve_emmax()
        except (OSError, VCFError):
            pass
    return {
        "tool_root": str(tool_root()), "platform": system,
        "plink": {"installed": bool(plink), "path": plink, "version": _probe(plink, ["--version"]),
                  "bundled_version": PLINK_VERSION, "license": "GPL-3.0"},
        "ldblockshow": {"installed": bool(ldblockshow), "path": ldblockshow,
                        "version": "LDBlockShow (official hewm2008 build)" if ldblockshow else None,
                        "license": "MIT", "requires_wsl": system == "Windows",
                        "wsl_available": wsl_installed},
        "bcftools": {
            "installed": bool(native_bcftools or wsl_bcftools),
            "path": native_bcftools or ("WSL · {}".format(wsl_bcftools) if wsl_bcftools else None),
            "version": _probe(native_bcftools, ["--version"]) if native_bcftools else None,
            "license": "MIT/Expat（部分插件GPL）",
            "requires_wsl": system == "Windows" and not native_bcftools,
            "wsl_available": wsl_installed,
            "install_mode": "WSL apt managed runtime" if system == "Windows" else "system package",
        },
        "emmax": {
            "installed": bool(emmax), "bundled": bundled_emmax,
            "ready": bool(emmax) and (system != "Windows" or wsl_installed),
            "path": emmax["emmax"] if emmax else None,
            "kin_path": emmax["kin"] if emmax else None,
            "version": EMMAX_VERSION, "license": "MIT",
            "requires_wsl": system == "Windows", "wsl_available": wsl_installed,
            "install_mode": "随安装包部署；Windows 通过 WSL 调用" if system == "Windows" else "随安装包部署",
            "archive_sha256": EMMAX_ARCHIVE_SHA256,
        },
    }


def install_tool(name):
    name = str(name or "").lower()
    root = tool_root()
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="callvcf-tool-") as temp_name:
        temp = Path(temp_name)
        archive = temp / "download.zip"
        if name == "plink":
            url = PLINK_URLS.get(platform.system())
            if not url:
                raise VCFError("当前平台暂不支持自动安装 PLINK，请手动选择可执行文件")
            size, digest = _download(url, archive)
            extracted = temp / "extracted"
            _safe_extract(archive, extracted)
            filename = "plink.exe" if platform.system() == "Windows" else "plink"
            source = next((p for p in extracted.rglob(filename) if p.is_file()), None)
            if not source:
                raise VCFError("PLINK 官方压缩包中未找到 {}".format(filename))
            destination = root / "plink"
            destination.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination / filename)
            for candidate in extracted.rglob("LICENSE*"):
                if candidate.is_file():
                    shutil.copy2(candidate, destination / candidate.name)
            if platform.system() != "Windows":
                (destination / filename).chmod((destination / filename).stat().st_mode | stat.S_IEXEC)
            _write_metadata("plink", {"source": url, "version": PLINK_VERSION,
                                      "sha256": digest, "download_bytes": size, "license": "GPL-3.0"})
        elif name == "ldblockshow":
            size, digest = _download(LDBLOCKSHOW_URL, archive)
            extracted = temp / "extracted"
            _safe_extract(archive, extracted)
            source_root = next((p for p in extracted.iterdir() if p.is_dir()), None)
            binary = next((p for p in extracted.rglob("LDBlockShow")
                           if p.is_file() and p.parent.name == "bin"), None)
            if not source_root or not binary:
                raise VCFError("LDBlockShow 官方压缩包中未找到 bin/LDBlockShow")
            destination = root / "ldblockshow"
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(source_root, destination)
            installed = destination / "bin" / "LDBlockShow"
            installed.chmod(installed.stat().st_mode | stat.S_IEXEC)
            _write_metadata("ldblockshow", {"source": LDBLOCKSHOW_URL, "version": "main",
                                             "sha256": digest, "download_bytes": size, "license": "MIT"})
        elif name == "bcftools":
            if platform.system() != "Windows":
                raise VCFError("当前版本只在Windows中提供bcftools一键管理；Linux请使用系统包管理器安装")
            wsl = shutil.which("wsl.exe") or shutil.which("wsl")
            if not wsl or not _wsl_operational():
                raise VCFError("请先完成Ubuntu/WSL首次初始化（创建Linux用户名并进入一次终端），然后再点一键安装")
            command = [str(wsl), "-u", "root", "-e", "sh", "-lc",
                       "export DEBIAN_FRONTEND=noninteractive; apt-get update && apt-get install -y bcftools tabix"]
            try:
                proc = subprocess.run(
                    command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, timeout=1800,
                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
                )
            except subprocess.TimeoutExpired:
                raise VCFError("bcftools安装超过30分钟；请检查WSL网络后重试")
            if proc.returncode or not _wsl_bcftools_path():
                raise VCFError("bcftools安装失败：{}".format((proc.stdout or "无详细输出")[-1500:]))
            _write_metadata("bcftools", {
                "source": "Ubuntu apt repositories", "package": "bcftools + tabix",
                "license": "MIT/Expat（部分插件GPL）", "backend": "WSL",
            })
        elif name == "emmax":
            deploy_bundled_emmax()
        else:
            raise VCFError("无法识别的工具：{}".format(name))
    return tools_status()
