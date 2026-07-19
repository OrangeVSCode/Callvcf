#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
exec python3 app.py --host "${VCF_TOOL_HOST:-127.0.0.1}" --port "${VCF_TOOL_PORT:-8765}"

