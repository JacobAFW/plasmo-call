#!/usr/bin/env bash
# =============================================================================
# plasmo-call installer  —  clone -> ./install.sh -> run
#
# Strategy (LOCKED): build ON TOP of vvg-box. vvg-box provides the pixi
# instance and the auto-detecting PBS/SLURM snakemake profiles; this script
# bootstraps it, then installs plasmo-call's own pixi dependencies on top.
#
# NOTE: the PBS profile shipped by vvg-box does not work unmodified on every
# PBS HPC. plasmo-call overrides it with a known-working PBS setup (see
# profiles/pbs/). See README "Schedulers".
#
# Test ladder: (1) LOCAL, no scheduler  (2) PBS HPC  (3) SLURM HPC.
# Env vars honoured (passed through to vvg-box): VVG_BASEDIR, PIXI_ENVNAME, PYVER
# =============================================================================
set -euo pipefail

echo ">> plasmo-call install — SKELETON. Steps below are stubbed, not yet active."

# 1. Bootstrap vvg-box (environment + scheduler profiles)
#    "${SHELL}" <(curl -sSL https://raw.githubusercontent.com/vivaxgen/vvg-box/main/install.sh)
echo "[TODO] bootstrap vvg-box into \${VVG_BASEDIR:-./box}"

# 2. Install plasmo-call's pixi dependencies on top (pixi.toml -> pixi.lock)
echo "[TODO] pixi install (pins gatk4/bcftools/samtools/bwa/picard/snakemake)"

# 3. Select scheduler profile: vvg-box auto-detects qsub/sbatch; for PBS we
#    prefer plasmo-call's profiles/pbs override.
echo "[TODO] wire scheduler profile (local | pbs | slurm)"

echo ">> Skeleton install complete (no-op). Real steps land in the build phase."
