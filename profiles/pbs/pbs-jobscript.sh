#!/bin/sh
# plasmo-call PBS jobscript.
# Snakemake substitutes {properties} and {exec_job} at submit time.
# properties = {properties}
#
# Logging policy — keep a log ONLY when the job fails. A clean per-interval
# scatter run would otherwise drop hundreds of PBS .o files; here PBS's own
# stdout/stderr go to /dev/null (see pbs-submit.py), and this wrapper captures
# the job's output to logs/plasmo.<jobid>.log, deleting it on success. So the
# logs/ dir contains exactly the jobs that failed. Resource usage for a killed
# job is still available via `qstat -xf <jobid>`.
cd "${PBS_O_WORKDIR:-.}" 2>/dev/null || true
mkdir -p logs
_log="logs/plasmo.${PBS_JOBID:-$$}.log"
{ {exec_job} ; } > "$_log" 2>&1
_rc=$?
if [ "$_rc" -eq 0 ]; then
    rm -f "$_log"
else
    echo "plasmo-call: job ${PBS_JOBID:-?} FAILED (exit $_rc) — log kept at $_log" >&2
fi
exit $_rc
