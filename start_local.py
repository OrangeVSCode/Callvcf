#!/usr/bin/env python3
import argparse
import json
import threading
import webbrowser
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import app as app_module
from app import AppServer, Handler
from advanced_analysis import AdvancedAnalyzer
from quality_engine import QualityJobManager
from phenotype_engine import PhenotypeAnalyzer
from association_engine import VariantPhenotypeAnalyzer
from similarity_engine import SimilarityJobManager
from emmax_engine import EmmaxJobManager
from repair_engine import RepairExecutor
from vcf_service import create_service


def existing_service_status(port):
    """Return ``gpa``, ``other`` or ``free`` for the requested local port.

    Windows may allow two Python HTTP servers to share a port.  Detecting a
    legacy CallVCF/VCF Explorer service before binding prevents the new app
    from intermittently serving the old interface on the same address.
    """
    try:
        with urlopen("http://127.0.0.1:{}/api/health".format(port), timeout=1) as response:
            data = json.loads(response.read().decode("utf-8"))
        if data.get("ok") and data.get("service") == "GPA-Accelerator":
            return "gpa"
        return "other"
    except HTTPError:
        return "other"
    except URLError:
        return "free"
    except Exception:
        return "other"


def main():
    parser = argparse.ArgumentParser(description="Launch the local GPA-Accelerator desktop interface")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--bcftools", default=None)
    parser.add_argument("--no-browser", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    existing_status = existing_service_status(args.port)
    if existing_status == "gpa":
        if not args.no_browser:
            webbrowser.open("http://127.0.0.1:{}".format(args.port))
        return

    app_module.SERVICE = create_service(args.bcftools)
    app_module.ADVANCED = AdvancedAnalyzer(app_module.SERVICE)
    app_module.QUALITY = QualityJobManager(app_module.SERVICE)
    app_module.REPAIR = RepairExecutor(app_module.SERVICE)
    app_module.PHENOTYPE = PhenotypeAnalyzer()
    app_module.ASSOCIATION = VariantPhenotypeAnalyzer(app_module.SERVICE)
    app_module.EMMAX = EmmaxJobManager(app_module.SERVICE)
    app_module.SIMILARITY = SimilarityJobManager()
    preferred_port = 0 if existing_status == "other" else args.port
    try:
        server = AppServer(("127.0.0.1", preferred_port), Handler)
    except OSError:
        server = AppServer(("127.0.0.1", 0), Handler)
    port = server.server_address[1]
    url = "http://127.0.0.1:{}".format(port)
    print("GPA-Accelerator local interface: {}".format(url), flush=True)
    if existing_status == "other":
        print("Legacy local service detected on port {}; using an isolated port.".format(args.port), flush=True)
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
