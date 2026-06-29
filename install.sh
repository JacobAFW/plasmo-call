#!/usr/bin/env bash
# =============================================================================
# plasmo-call installer  —  clone -> ./install.sh -> activate -> run
#
# Strategy (LOCKED): build ON TOP of vvg-box.
#   1. Bootstrap vvg-box into a self-contained dir inside this repo (./box by
#      default; the path is gitignored). vvg-box provides pixi + the auto-
#      detecting PBS/SLURM Snakemake profiles.
#   2. Install plasmo-call's own pinned pixi workspace at the repo root
#      (./pixi.toml + ./pixi.lock). This is where gatk4/bcftools/samtools/
#      bwa/bwa-mem2/picard/snakemake-minimal live.
#   3. Verify tool versions + Snakefile parses.
#
# vvg-box ships a PBS profile that doesn't run unmodified on every PBS HPC.
# plasmo-call overrides it via profiles/pbs/. See profiles/README.md.
#
# Test ladder: (1) LOCAL, no scheduler  (2) PBS HPC  (3) SLURM HPC.
#
# Env vars honoured:
#   VVG_BASEDIR   default: ./box   (where vvg-box installs itself)
#   PIXI_ENVNAME  default: vvg-box (vvg-box's internal pixi env name)
#   PYVER         passed through to vvg-box (default 3.12)
# =============================================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo
echo ">> plasmo-call install — bootstrapping vvg-box, then layering tool pins."
echo

# Refuse to run inside an active conda/pixi env (mirrors vvg-box's own guard).
if [[ ${CONDA_SHLVL:-0} -ge 1 ]]; then
  echo "ERROR: deactivate the active conda/mamba environment first." >&2
  exit 1
fi
if [[ -n "${PIXI_ENVIRONMENT_NAME:-}" ]]; then
  echo "ERROR: deactivate the active pixi environment first." >&2
  exit 1
fi

# 1. Bootstrap vvg-box into ${VVG_BASEDIR} (default ./box, kept local to repo).
export VVG_BASEDIR="${VVG_BASEDIR:-${REPO_ROOT}/box}"
export PIXI_ENVNAME="${PIXI_ENVNAME:-vvg-box}"
export PYVER="${PYVER:-3.12}"

if [[ -x "${VVG_BASEDIR}/bin/activate" ]]; then
  echo ">> vvg-box already installed at ${VVG_BASEDIR} (skipping bootstrap)."
else
  # vvg-box's installer calls `realpath "${VVG_BASEDIR}"` before creating the
  # directory itself, so it must already exist.
  mkdir -p "${VVG_BASEDIR}"

  # vvg-box's installer uses GNU `ln -sr` (relative symlinks). BSD ln on macOS
  # lacks `-r`, which breaks the bootstrap there. If we detect BSD ln, drop a
  # tiny shim in front of PATH that translates `-sr` to a python-computed
  # relative target and forwards everything else to /bin/ln untouched.
  if ! ln --version 2>/dev/null | grep -q GNU; then
    SHIMS_DIR="${VVG_BASEDIR}/.shims"
    mkdir -p "${SHIMS_DIR}"
    cat > "${SHIMS_DIR}/ln" <<'SHIM'
#!/usr/bin/env bash
# plasmo-call ln(1) shim: emulate GNU `ln -s -r [-f]` on BSD ln (macOS).
# Splits merged short-opt clusters (-sr, -srf, -rs, etc) so we can detect -r.
# Anything that doesn't end up as a -s -r combo passes straight to /bin/ln.
want_rel=0; want_sym=0; want_force=0; rest=()
for a in "$@"; do
  case "$a" in
    --relative)  want_rel=1 ;;
    --symbolic)  want_sym=1 ;;
    --force)     want_force=1 ;;
    -[srfRSF]*)
      # short-opt cluster: split into individual letters.
      sub="${a#-}"
      pass_through=""
      for (( i=0; i<${#sub}; i++ )); do
        c="${sub:i:1}"
        case "$c" in
          s) want_sym=1 ;;
          r) want_rel=1 ;;
          f) want_force=1 ;;
          *) pass_through+="$c" ;;
        esac
      done
      [[ -n "$pass_through" ]] && rest+=("-$pass_through")
      ;;
    *)           rest+=("$a") ;;
  esac
done
if [[ $want_rel -eq 1 && $want_sym -eq 1 && ${#rest[@]} -ge 2 ]]; then
  src="${rest[${#rest[@]}-2]}"
  tgt="${rest[${#rest[@]}-1]}"
  rel="$(python3 -c 'import os,sys; print(os.path.relpath(sys.argv[1], os.path.dirname(sys.argv[2]) or "."))' "$src" "$tgt")"
  flags="-s"
  [[ $want_force -eq 1 ]] && flags="-sf"
  exec /bin/ln "$flags" "$rel" "$tgt"
fi
# Not a -s -r call; rebuild from any leftover and pass through.
[[ $want_sym -eq 1 ]] && rest=("-s" "${rest[@]}")
[[ $want_force -eq 1 ]] && rest=("-f" "${rest[@]}")
exec /bin/ln "${rest[@]}"
SHIM
    chmod +x "${SHIMS_DIR}/ln"
    export PATH="${SHIMS_DIR}:${PATH}"
    echo ">> macOS detected — installed BSD-ln shim at ${SHIMS_DIR}"
  fi

  echo ">> Bootstrapping vvg-box into ${VVG_BASEDIR}"
  # Pipe vvg-box's upstream installer through bash. Env vars above flow in.
  if command -v curl >/dev/null 2>&1; then
    FETCH=(curl -fsSL)
  elif command -v wget >/dev/null 2>&1; then
    FETCH=(wget -qO-)
  else
    echo "ERROR: need curl or wget to fetch vvg-box installer." >&2
    exit 1
  fi
  # Run in a subshell so the installer's `set -e` doesn't kill us; bash reads
  # from /dev/stdin so it stays non-interactive.
  #
  # cwd must be OUTSIDE plasmo-call: pixi walks up looking for a pixi.toml, and
  # would otherwise mutate our repo-root manifest instead of vvg-box's own
  # workspace at ${VVG_BASEDIR}/opt/pixi/${PIXI_ENVNAME}/pixi.toml.
  ( cd /tmp && "${FETCH[@]}" https://raw.githubusercontent.com/vivaxgen/vvg-box/main/install.sh | bash )
fi

# 2. Activate vvg-box to get pixi onto PATH, then install plasmo-call's pixi
#    workspace at the repo root.
# shellcheck disable=SC1091
source "${VVG_BASEDIR}/bin/activate"

if ! command -v pixi >/dev/null 2>&1; then
  echo "ERROR: pixi not on PATH after sourcing vvg-box activate." >&2
  exit 1
fi

echo ">> Installing plasmo-call pixi workspace at ${REPO_ROOT}"
# Override PIXI_PROJECT_MANIFEST so pixi uses our manifest cleanly without
# warning about overriding vvg-box's activation-time value.
# conda-forge compiler activation scripts reference unset vars; soften `-u`
# while pixi enters the env so they don't blow up the install.
set +u
export PIXI_PROJECT_MANIFEST="${REPO_ROOT}/pixi.toml"
( cd "${REPO_ROOT}" && pixi install )

# 3. Verify: print tool versions, confirm local scheduler resolution.
echo
echo ">> Tool versions (from plasmo-call pixi workspace):"
( cd "${REPO_ROOT}" && pixi run versions ) || true

echo
echo ">> Scheduler profile resolution (vvg-box auto-detect):"
"${VVG_BASEDIR}/envs/vvg-box/bin/set-snakemake-profile.py" || true
set -u

echo
echo ">> plasmo-call install complete."
echo "   Activate with:  source ${VVG_BASEDIR}/bin/activate"
echo "   Then run:       pixi run run-local CORES=4"
echo
# EOF
