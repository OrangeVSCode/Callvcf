#!/usr/bin/env python3
import argparse
import json
import mimetypes
import os
import sys
import threading
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

from vcf_service import VCFError, create_service
from advanced_analysis import AdvancedAnalyzer


ROOT = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
STATIC_DIR = ROOT / "static"
SERVICE = None
ADVANCED = None


def select_local_file(initial_dir=None):
    try:
        import tkinter as tk
        from tkinter import filedialog
    except ImportError:
        raise VCFError("当前 Python 未包含 tkinter；请直接粘贴 VCF 文件路径")
    root = tk.Tk()
    root.withdraw()
    try:
        root.attributes("-topmost", True)
        path = filedialog.askopenfilename(
            title="选择 VCF / VCF.GZ / BCF 文件",
            initialdir=initial_dir if initial_dir and Path(initial_dir).is_dir() else None,
            filetypes=[
                ("Variant files", "*.vcf *.vcf.gz *.vcf.bgz *.bcf"),
                ("VCF", "*.vcf"),
                ("Compressed VCF", "*.vcf.gz *.vcf.bgz"),
                ("BCF", "*.bcf"),
                ("All files", "*.*"),
            ],
        )
    finally:
        root.destroy()
    return {"path": path or None}


def select_resource(kind="file", initial_dir=None):
    try:
        import tkinter as tk
        from tkinter import filedialog
    except ImportError:
        raise VCFError("当前 Python 未包含 tkinter，请直接粘贴资源路径")
    root = tk.Tk()
    root.withdraw()
    try:
        root.attributes("-topmost", True)
        kwargs = {"initialdir": initial_dir if initial_dir and Path(initial_dir).is_dir() else None}
        if kind in {"directory", "phenotype_dir", "output_dir"}:
            path = filedialog.askdirectory(title="选择目录", **kwargs)
        else:
            filters = {
                "gff": [("Gene annotation", "*.gff *.gff3 *.gtf *.gff.gz *.gff3.gz *.gtf.gz"), ("All files", "*.*")],
                "annotation": [("Annotation table", "*.tsv *.csv *.txt *.gz"), ("All files", "*.*")],
                "domain": [("Domain table", "*.tsv *.csv *.txt *.gz"), ("All files", "*.*")],
                "phenotype": [("Phenotype PS", "*.ps"), ("All files", "*.*")],
                "executable": [("Executable", "*.exe *"), ("All files", "*.*")],
            }
            path = filedialog.askopenfilename(title="选择资源文件", filetypes=filters.get(kind, [("All files", "*.*")]), **kwargs)
    finally:
        root.destroy()
    return {"path": path or None}


class AppServer(ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True


class Handler(BaseHTTPRequestHandler):
    server_version = "VCFQueryTool/1.0"

    def log_message(self, fmt, *args):
        sys.stdout.write("[%s] %s\n" % (self.log_date_time_string(), fmt % args))
        sys.stdout.flush()

    def _json(self, status, payload):
        body = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(body)

    def _read_json(self):
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            raise VCFError("无效的请求长度")
        if length <= 0 or length > 5 * 1024 * 1024:
            raise VCFError("请求为空或超过 5 MB")
        try:
            return json.loads(self.rfile.read(length).decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            raise VCFError("请求不是有效的 JSON")

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/api/health":
            return self._json(200, {
                "ok": True,
                "service": "VCF Query Tool",
                "backend": getattr(SERVICE, "backend", "bcftools"),
                "bcftools": SERVICE.bcftools,
            })
        relative = "index.html" if parsed.path in {"", "/"} else unquote(parsed.path.lstrip("/"))
        candidate = (STATIC_DIR / relative).resolve()
        if STATIC_DIR not in candidate.parents and candidate != STATIC_DIR:
            return self._json(403, {"ok": False, "error": "禁止访问"})
        if not candidate.is_file():
            candidate = STATIC_DIR / "index.html"
        body = candidate.read_bytes()
        mime = mimetypes.guess_type(str(candidate))[0] or "application/octet-stream"
        self.send_response(200)
        self.send_header("Content-Type", mime + ("; charset=utf-8" if mime.startswith("text/") else ""))
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-cache")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        try:
            payload = self._read_json()
            route = urlparse(self.path).path
            if route == "/api/select-file":
                result = select_local_file(payload.get("initial_dir"))
            elif route == "/api/select-resource":
                result = select_resource(payload.get("kind", "file"), payload.get("initial_dir"))
            elif route == "/api/shutdown":
                result = {"message": "CallVCF 正在关闭"}
                self._json(200, {"ok": True, "data": result})
                threading.Timer(0.2, self.server.shutdown).start()
                return
            elif route == "/api/files":
                result = SERVICE.discover(payload.get("root"), max_depth=int(payload.get("max_depth", 5)))
            elif route == "/api/inspect":
                result = SERVICE.inspect(payload.get("path"), force=bool(payload.get("force")))
            elif route == "/api/check":
                result = SERVICE.check_loci(payload.get("path"), payload.get("loci"))
            elif route == "/api/distribution":
                result = SERVICE.distributions(payload.get("path"), payload.get("loci"))
            elif route == "/api/sample-loci":
                result = SERVICE.sample_locus_matrix(
                    payload.get("path"), payload.get("loci"), payload.get("samples") or []
                )
            elif route == "/api/sample-stats":
                result = SERVICE.sample_stats(
                    payload.get("path"), payload.get("samples") or [], payload.get("loci")
                )
            elif route == "/api/lead-analysis":
                result = ADVANCED.analyze(payload)
            else:
                return self._json(404, {"ok": False, "error": "接口不存在"})
            return self._json(200, {"ok": True, "data": result})
        except VCFError as exc:
            return self._json(400, {"ok": False, "error": str(exc)})
        except Exception as exc:
            traceback.print_exc()
            return self._json(500, {"ok": False, "error": "内部错误：{}".format(exc)})


def main():
    parser = argparse.ArgumentParser(description="Interactive VCF query tool")
    parser.add_argument("--host", default=os.environ.get("VCF_TOOL_HOST", "127.0.0.1"))
    parser.add_argument("--port", type=int, default=int(os.environ.get("VCF_TOOL_PORT", "8765")))
    parser.add_argument("--bcftools", default=os.environ.get("BCFTOOLS"))
    args = parser.parse_args()

    global SERVICE, ADVANCED
    SERVICE = create_service(args.bcftools)
    ADVANCED = AdvancedAnalyzer(SERVICE)
    server = AppServer((args.host, args.port), Handler)
    print("VCF Query Tool: http://{}:{}".format(args.host, args.port), flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
