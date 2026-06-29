# Scheduler profiles

plasmo-call runs in three modes. **Local is the first-class default and the
first rung of the test ladder.**

1. **local** — no scheduler: `snakemake --cores N`. Build + smoke-test here.
2. **pbs** — `profiles/pbs/` (this repo's known-working PBS setup).
3. **slurm** — `profiles/slurm/` (adapted from vvg-box).

## Why PBS lives here and not just in vvg-box
vvg-box ships a `pbspro` profile, but it does **not** run unmodified on every
PBS HPC (the vvg-box authors could not get it working on Jacob's PBS cluster).
The original Variant_Calling_Pipeline, however, submits successfully on that
same cluster using the classic Snakemake cluster string:

    qsub -N {jobname} -l ncpus={ncpus},mem={mem},walltime={walltime},storage={storage} \
         -A {account} -q {queue} -M {email} -m {mailon} -j oe -V -S /bin/sh

`profiles/pbs/` is seeded from that working invocation. All site-specific values
(account, queue, storage, email) are parameters — never hard-coded here.
Candidate to upstream back to the vvg-box authors once validated.
