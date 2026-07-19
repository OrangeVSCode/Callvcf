#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ ! -f callvcf.pid ]]; then
  echo "CallVCF is not running (no PID file)"
  exit 0
fi

PID="$(cat callvcf.pid)"
if [[ ! "$PID" =~ ^[0-9]+$ ]]; then
  echo "Invalid PID file: $PID" >&2
  exit 1
fi
if kill -0 "$PID" 2>/dev/null; then
  kill "$PID"
  for _ in {1..20}; do
    kill -0 "$PID" 2>/dev/null || break
    sleep 0.1
  done
fi
rm -f callvcf.pid
echo "CallVCF stopped"

