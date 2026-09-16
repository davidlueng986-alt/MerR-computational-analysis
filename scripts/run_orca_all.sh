#!/usr/bin/env bash
set -euo pipefail
for d in jobs/*; do
  [ -d "$d" ] || continue
  echo "=== $d ==="
  (cd "$d" && ./run.sh)
done
