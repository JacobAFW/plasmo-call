#!/usr/bin/env python3
# =============================================================================
# profiles/pbs/pbs-submit.py
# plasmo-call PBS submit script for Snakemake 8's `cluster-generic` executor.
#
# Reads site-specific values from `profiles/pbs/site.local.yaml` (gitignored;
# copy from `site.example.yaml`) and issues the known-working qsub
# invocation from JacobAFW/Variant_Calling_Pipeline:
#
#   qsub -N <name> \
#        -l ncpus=<t>,mem=<m>MB,walltime=<hh:mm:ss> \
#        -A <account> -q <queue> \
#        -l storage=<storage> \
#        -M <email> -m <mailon> \
#        -j oe -V -S /bin/sh \
#        <jobscript>
#
# Nothing site-specific is baked into this file — account/queue/storage/email
# all come from the gitignored site YAML. Missing site values are simply
# omitted (so blank email → no `-M/-m` flags).
#
# Packaged so it can be offered back to the vvg-box authors as a drop-in
# replacement for `etc/snakemake-profiles/pbspro/pbs-submit.py`.
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
    """Site params from site.local.yaml. If missing, fall back to the example
    (with a loud stderr warning) — this keeps `snakemake -n` dry-runs working
    on a fresh clone that hasn't yet copied the template, without silently
    letting a real submission use placeholder values."""
    if SITE_FILE.exists():
        return yaml.safe_load(SITE_FILE.read_text()) or {}
    if SITE_EXAMPLE.exists():
        sys.stderr.write(
            f"WARN: {SITE_FILE.name} not found; falling back to "
            f"{SITE_EXAMPLE.name} placeholders. Copy the example and fill in "
            "your site values before running for real.\n"
        )
        return yaml.safe_load(SITE_EXAMPLE.read_text()) or {}
    sys.stderr.write(
        f"ERROR: neither {SITE_FILE} nor {SITE_EXAMPLE} exists.\n"
    )
    sys.exit(2)


def _fmt_walltime(runtime_min: int) -> str:
    """Convert Snakemake `runtime` (minutes) into PBS `HH:MM:SS`."""
    total = int(runtime_min)
    h, rem = divmod(total, 60)
    m, s = rem, 0
    return f"{h:02d}:{m:02d}:{s:02d}"


def _job_name(props: dict) -> str:
    """Compose a PBS jobname from rule + wildcards. PBS limits are ~236 chars
    but many schedulers reject much shorter names; be conservative."""
    base = f"smk.{props.get('rule', 'job')}"
    wc = props.get("wildcards", {}) or {}
    if wc:
        # Sanitise: PBS jobnames must not contain whitespace or `:` (colons
        # confuse some qsub versions). Our intervals look like
        # `LT727648:1-89247`; replace `:` with `_`.
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

    cmd: list[str] = ["qsub"]

    cmd += ["-N", _job_name(props)]
    cmd += ["-l", f"ncpus={threads},mem={mem_mb}MB,walltime={_fmt_walltime(runtime)}"]

    # Site-specific flags — each omitted if the site value is missing/blank.
    if site.get("account"):
        cmd += ["-A", str(site["account"])]
    if site.get("queue"):
        cmd += ["-q", str(site["queue"])]
    if site.get("storage"):
        cmd += ["-l", f"storage={site['storage']}"]
    if site.get("email"):
        cmd += ["-M", str(site["email"])]
        cmd += ["-m", str(site.get("mailon", "a"))]

    # Always-on flags matching the known-working invocation.
    cmd += ["-o", "/dev/null", "-e", "/dev/null", "-V", "-S", "/bin/sh"]

    # Extra flags from env for one-off overrides (matches vvg-box's convention).
    extras = os.environ.get("SNAKEMAKE_CLUSTER_EXTRA_FLAGS", "").strip()
    if extras:
        cmd += extras.split()

    cmd.append(jobscript)

    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(
            f"ERROR: qsub failed (exit {e.returncode})\n"
            f"       cmd: {' '.join(cmd)}\n"
            f"       stderr: {e.stderr}\n"
        )
        sys.exit(e.returncode)

    # Snakemake reads the job id from stdout.
    print(res.stdout.strip())


if __name__ == "__main__":
    main()
