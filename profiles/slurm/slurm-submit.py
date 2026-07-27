#!/usr/bin/env python3
# =============================================================================
# profiles/slurm/slurm-submit.py
# plasmo-call SLURM submit script for Snakemake 8's `cluster-generic` executor.
#
# Reads site-specific values from `profiles/slurm/site.local.yaml`
# (gitignored; copy from `site.example.yaml`) and issues:
#
#   sbatch -J <name> \
#          --cpus-per-task=<t> \
#          --mem=<m>M \
#          --time=<HH:MM:SS> \
#          -p <partition> -A <account> [--qos=<qos>] \
#          [--mail-user=<email> --mail-type=<mailtype>] \
#          -o slurm-%j.out -e slurm-%j.err \
#          <jobscript>
#
# Symmetric with pbs-submit.py — same site-yaml → command mechanism.
# =============================================================================

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import yaml
from snakemake.utils import read_job_properties


PROFILE_DIR = Path(__file__).resolve().parent
SITE_FILE     = PROFILE_DIR / "site.local.yaml"
SITE_EXAMPLE  = PROFILE_DIR / "site.example.yaml"


def _load_site() -> dict:
    if SITE_FILE.exists():
        return yaml.safe_load(SITE_FILE.read_text()) or {}
    if SITE_EXAMPLE.exists():
        sys.stderr.write(
            f"WARN: {SITE_FILE.name} not found; falling back to "
            f"{SITE_EXAMPLE.name} placeholders. Copy the example and fill in "
            "your site values before running for real.\n"
        )
        return yaml.safe_load(SITE_EXAMPLE.read_text()) or {}
    sys.stderr.write(f"ERROR: neither {SITE_FILE} nor {SITE_EXAMPLE} exists.\n")
    sys.exit(2)


def _fmt_walltime(runtime_min: int) -> str:
    total = int(runtime_min)
    h, rem = divmod(total, 60)
    m, s = rem, 0
    return f"{h:02d}:{m:02d}:{s:02d}"


def _job_name(props: dict) -> str:
    base = f"smk.{props.get('rule', 'job')}"
    wc = props.get("wildcards", {}) or {}
    if wc:
        parts = ".".join(str(v).replace(":", "_") for v in wc.values())
        base = f"{base}.{parts}"
    return base[:120]


def main() -> None:
    site = _load_site()
    jobscript = sys.argv[-1]
    props = read_job_properties(jobscript)

    threads = int(props.get("threads", 1))
    resources = props.get("resources", {}) or {}
    mem_mb   = int(resources.get("mem_mb", 8192))
    runtime  = int(resources.get("runtime", 720))    # minutes

    cmd: list[str] = ["sbatch", "--parsable"]

    cmd += ["-J", _job_name(props)]
    cmd += [f"--cpus-per-task={threads}"]
    cmd += [f"--mem={mem_mb}M"]
    cmd += [f"--time={_fmt_walltime(runtime)}"]

    if site.get("partition"):
        cmd += ["-p", str(site["partition"])]
    if site.get("account"):
        cmd += ["-A", str(site["account"])]
    if site.get("qos"):
        cmd += [f"--qos={site['qos']}"]
    if site.get("email"):
        cmd += [f"--mail-user={site['email']}"]
        cmd += [f"--mail-type={site.get('mailtype', 'FAIL')}"]

    cmd += ["-o", "slurm-%j.out", "-e", "slurm-%j.err"]

    extras = os.environ.get("SNAKEMAKE_CLUSTER_EXTRA_FLAGS", "").strip()
    if extras:
        cmd += extras.split()

    cmd.append(jobscript)

    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(
            f"ERROR: sbatch failed (exit {e.returncode})\n"
            f"       cmd: {' '.join(cmd)}\n"
            f"       stderr: {e.stderr}\n"
        )
        sys.exit(e.returncode)

    # `--parsable` prints just the job id, which is what Snakemake expects.
    print(res.stdout.strip())


if __name__ == "__main__":
    main()
