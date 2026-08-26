#!/bin/sh
# plasmo-call PBS jobscript -- keep a per-job log ONLY on failure.
# PBS stdout/stderr go to /dev/null (see pbs-submit.py); this captures the job
# to logs/plasmo.<jobid>.log and deletes it on success.
# IMPORTANT: Snakemake runs str.format on this file, so the only curly-brace
# fields allowed are properties and exec_job below. Everything else uses plain
# $VAR (no braces) and a subshell in parentheses -- adding any other curly
# brace will raise a KeyError at submit time.
# properties = {properties}
cd "$PBS_O_WORKDIR" || exit 1
mkdir -p logs
_log="logs/plasmo.$PBS_JOBID.log"
( {exec_job} ) > "$_log" 2>&1
_rc=$?
if [ "$_rc" -eq 0 ]; then
    rm -f "$_log"
else
    echo "plasmo-call: FAILED (exit $_rc) -- kept $_log" >&2
fi
exit $_rc
