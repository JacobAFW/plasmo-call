#!/usr/bin/env bash
# Print resolved versions of every plasmo-call calling tool. Sourced by
# install.sh's verification step and by the `versions` pixi task.
set -eo pipefail

printf '%-10s : %s\n' "snakemake" "$(snakemake --version)"
printf '%-10s : %s\n' "gatk"      "$(gatk --version 2>&1 | tr '\n' ' ' | sed 's/  */ /g')"
printf '%-10s : %s\n' "bcftools"  "$(bcftools --version | head -n1)"
printf '%-10s : %s\n' "samtools"  "$(samtools --version | head -n1)"
printf '%-10s : %s\n' "bwa"       "$(bwa 2>&1 | awk '/^Version:/{print $2; exit}')"
if command -v bwa-mem2 >/dev/null 2>&1; then
  printf '%-10s : %s\n' "bwa-mem2" "$(bwa-mem2 version 2>&1 | tail -n1)"
else
  printf '%-10s : %s\n' "bwa-mem2" "(not available on this platform)"
fi
printf '%-10s : %s\n' "picard"    "$(picard MarkDuplicates --version 2>&1 | head -n1)"
printf '%-10s : %s\n' "python"    "$(python --version)"
