#!/usr/bin/env bash
# BQSR acceptance: runs bootstrap end-to-end, copies its filtered VCFs into
# test/known_sites/, then re-runs with mode=known_sites against them. Also
# dry-runs auto in both empty and populated states so the resolver's loud
# banner is visible in both directions.
#
# CORES overrides the snakemake --cores value (default 4).
set -eo pipefail

CORES="${CORES:-4}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

snakemake_run() {
  local overlay="$1"
  snakemake --cores "${CORES}" \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml test/config.yaml "${overlay}"
}

snakemake_dry() {
  # Capture first so `head` closing the pipe doesn't SIGPIPE snakemake under
  # `set -eo pipefail`.
  local overlay="$1" out
  out=$(snakemake -n \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml test/config.yaml "${overlay}" 2>&1)
  printf '%s\n' "${out}" | head -12
}

echo
echo "########################################################################"
echo "# 1. mode=bootstrap end-to-end                                          "
echo "########################################################################"
python test/generate-fixtures.py
rm -rf output .snakemake
snakemake_run test/config-bootstrap.yaml

echo
echo "########################################################################"
echo "# 2. capture bootstrap VCFs into test/known_sites/                      "
echo "########################################################################"
mkdir -p test/known_sites
cp output/bqsr/bootstrap/snvs_pass.vcf.gz       test/known_sites/snvs.vcf.gz
cp output/bqsr/bootstrap/snvs_pass.vcf.gz.tbi   test/known_sites/snvs.vcf.gz.tbi
cp output/bqsr/bootstrap/indels_pass.vcf.gz     test/known_sites/indels.vcf.gz
cp output/bqsr/bootstrap/indels_pass.vcf.gz.tbi test/known_sites/indels.vcf.gz.tbi
ls -la test/known_sites/

echo
echo "########################################################################"
echo "# 3. mode=known_sites end-to-end (uses the VCFs captured in step 2)     "
echo "########################################################################"
rm -rf output .snakemake
snakemake_run test/config-known-sites.yaml

echo
echo "########################################################################"
echo "# 4. mode=auto banner: known_variants empty (should pick bootstrap)     "
echo "########################################################################"
snakemake_dry test/config-auto-empty.yaml

echo
echo "########################################################################"
echo "# 5. mode=auto banner: known_variants populated (should pick known_sites)"
echo "########################################################################"
snakemake_dry test/config-auto-populated.yaml

echo
echo ">> BQSR acceptance: all 3 modes + auto resolution complete."
