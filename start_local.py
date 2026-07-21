#!/usr/bin/env python3
import argparse
import json
import threading
import webbrowser
from urllib.request import urlopen

import app as app_module
from app import AppServer, Handler
from advanced_analysis import AdvancedAnalyzer
from quality_engine import QualityJobManager
from repair_engine import RepairExecutor
from vcf_service import create_service


def existing_callvcf(port):
    try:
        with urlopen("http://127.0.0.1:{}/api/health".format(port), timeout=1) as response:
            data = json.loads(response.read().decode("utf-8"))
        return bool(data.get("ok") and data.get("service") == "VCF Query Tool")
    except Exception:
        return False


def main():
    parser = argparse.ArgumentParser(description="Launch the local CallVCF desktop interface")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--bcftools", default=None)
    parser.add_argument("--no-browser", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if existing_callvcf(args.port):
        if not args.no_browser:
            webbrowser.open("http://127.0.0.1:{}".format(args.port))
        return

    app_module.SERVICE = create_service(args.bcftools)
    app_module.ADVANCED = AdvancedAnalyzer(app_module.SERVICE)
    app_module.QUALITY = QualityJobManager(app_module.SERVICE)
    app_module.REPAIR = RepairExecutor(app_module.SERVICE)
    try:
        server = AppServer(("127.0.0.1", args.port), Handler)
    except OSError:
        server = AppServer(("127.0.0.1", 0), Handler)
    port = server.server_address[1]
    url = "http://127.0.0.1:{}".format(port)
    print("CallVCF local interface: {}".format(url), flush=True)
    print("Backend: {}".format(getattr(app_module.SERVICE, "backend", "unknown")), flush=True)
    if not args.no_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
