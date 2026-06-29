# PBS profile (seed = Jacob's known-working setup)

Source of truth: the classic `qsub ... -l ncpus=,mem=,walltime=,storage=` cluster
string from JacobAFW/Variant_Calling_Pipeline, which is confirmed working on the
target PBS HPC. To be reconciled with vvg-box's Snakemake-8 cluster-generic
profile during the build.

Site values are PARAMETERS (set per site, not committed):
  account, queue, storage, email, walltime/mem/ncpus per rule.

TODO (build phase):
  - config.yaml  (cluster-generic submit/status wiring OR legacy cluster string)
  - per-rule resources (port cluster.yaml: bwa_map, haplotype_caller, etc.)
