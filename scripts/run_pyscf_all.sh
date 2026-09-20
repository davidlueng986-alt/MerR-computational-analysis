#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
python3 -m mercury_merr.run_species --name all --engine pyscf --workdir .
