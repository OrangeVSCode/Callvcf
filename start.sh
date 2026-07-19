#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ -f callvcf.pid ]] && kill -0 "$(cat callvcf.pid)" 2>/dev/null; then
  echo "CallVCF is already running (PID $(cat callvcf.pid))"
  exit 0
fi

nohup python3 app.py \
  --host "${VCF_TOOL_HOST:-127.0.0.1}" \
  --port "${VCF_TOOL_PORT:-8765}" \
  > callvcf.log 2>&1 < /dev/null &
echo "$!" > callvcf.pid
echo "CallVCF started: http://${VCF_TOOL_HOST:-127.0.0.1}:${VCF_TOOL_PORT:-8765} (PID $!)"

