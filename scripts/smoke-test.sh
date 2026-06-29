#!/usr/bin/env bash
# Local smoke test: regenerate the tiny synthetic dataset, then run the
# pipeline end-to-end (no scheduler). Override cores with CORES=N.
set -eo pipefail

CORES="${CORES:-4}"

python test/generate-fixtures.py

snakemake --cores "${CORES}" \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --configfile test/config.yaml \
  "$@"
