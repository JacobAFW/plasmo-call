# PBS profile — known-working reference (sanitised from the predecessor)

This is the **mechanism that actually submits successfully on Jacob's PBS HPC**,
captured from JacobAFW/Variant_Calling_Pipeline (`config/pbs-submission/`). It is
the ground truth to port/reconcile with vvg-box's cluster-generic PBS profile.

Site-specific values are shown as `<placeholders>` — they are PARAMETERS, never
committed. (In the predecessor they were: account `jw1542`, queue `normalbw`,
storage `gdata/<proj>+scratch/<proj>`, a Menzies email — all stripped here.)

## Submission (classic Snakemake cluster string — Snakemake <8)
```
cluster-config: "cluster.yaml"
cluster: "qsub -N {cluster.jobname} -l ncpus={cluster.ncpus},mem={cluster.mem},walltime={cluster.walltime},storage={cluster.storage} -A {cluster.account} -q {cluster.queue} -M {cluster.email} -m {cluster.mailon} -j {cluster.jobout} -V -S /bin/sh"
jobs: 100
notemp: true
verbose: true
```

## Per-rule resources (from the working cluster.yaml)
| scope                | ncpus | mem   | walltime  |
|----------------------|-------|-------|-----------|
| __default__          | 1     | 8GB   | 12:00:00  |
| bwa_map (mapping)    | 5     | 50GB  | 24:00:00  |
| haplotype_caller     | 3     | 24GB  | 48:00:00  |
| combine_gvcfs        | 3     | 24GB  | 48:00:00  |
| joint_genotyping     | 3     | 24GB  | 48:00:00  |
| bcftools_caller      | 2     | 16GB  | 48:00:00  |
| concat_bcftools      | 2     | 16GB  | 48:00:00  |
| consensus_of_vcfs    | 2     | 16GB  | 48:00:00  |
| concat_vcfs          | 2     | 16GB  | 12:00:00  |

Default flags also set: `-j oe`, `-m a` (mail on abort), `-V`, `-S /bin/sh`,
`storage=gdata/<proj>+scratch/<proj>`.

## Note (GATK4 migration)
`realigner_target_creator` existed in the GATK3 predecessor and is **dropped**.
New rule names (mapping/markdup/bqsr/etc.) inherit `__default__` unless given an
explicit entry. Reconcile this submission string with vvg-box's Snakemake-8
`cluster-generic` executor (pbs-submit.py/pbs-status.py); if the generic path
keeps failing on this HPC, fall back to this proven cluster-string mechanism.
