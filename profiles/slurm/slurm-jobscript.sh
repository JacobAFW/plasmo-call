#!/bin/bash
# plasmo-call SLURM jobscript — minimal wrapper.
# Snakemake substitutes {properties} and {exec_job} at submit time.
# properties = {properties}
{exec_job}
