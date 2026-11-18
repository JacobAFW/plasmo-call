# Running on an HPC (PBS / Gadi) — run-book & gotchas

Local runs need none of this. These are the things that bite on a scheduler,
learned running on NCI Gadi. Update `<proj>` etc. for your site.

## Install
- Run `bash install.sh` on a **login node** (compute nodes have no internet for
  the vvg-box / pixi downloads).
- Run it in a **clean shell with no conda/pixi module loaded** — vvg-box's
  installer refuses to run inside an active conda/pixi environment.
- Clone into `/g/data/<proj>` (persistent), not home (small quota) or scratch
  (purged).

## Selecting the plasmo-call profile (not vvg-box's)
`source box/bin/activate` makes vvg-box **auto-detect the scheduler and set its
own `SNAKEMAKE_PROFILE`** (its stock `pbspro` profile). That profile is NOT
Gadi-tuned (it omits the storage flag — see below), so for plasmo-call you must
override it:

```bash
unset SNAKEMAKE_PROFILE
pixi run -- snakemake --profile profiles/pbs \
  --configfile config/config.yaml \
  --configfile config/config.run-local.yaml
```

**Run through `pixi run`, not bare `snakemake`.** plasmo-call's tools
(gatk4/samtools/bwa/picard) live in plasmo-call's pixi env
(`.pixi/envs/default/bin`), NOT in the vvg-box env. `source box/bin/activate`
gives you vvg-box's snakemake but not the tools — so bare `snakemake` submits
jobs whose PATH (propagated via qsub `-V`) has no tools, and every job dies with
`picard: command not found` (exit 127). `pixi run -- snakemake ...` runs in the
tool env, and `-V` then carries that PATH to the compute nodes.

The same `unset` is needed to run the **local smoke test** without submitting
jobs: `unset SNAKEMAKE_PROFILE; pixi run smoke-test`.

## Run via a CONTROLLER job — never on the login node
Do NOT run the Snakemake orchestrator interactively (login node / tmux) — the
login-node reaper or a dropped session will kill it, orphaning jobs (you'll see
`touch: .snakemake/tmp.XXX/NN.jobfinished: No such file or directory`). Instead
submit a lightweight **controller job** that runs Snakemake from a compute node
and submits the workers via nested qsub (Gadi allows this). Example
`run-controller.pbs` (fill in project/queue/storage):

```bash
#!/bin/bash
#PBS -P <proj>
#PBS -q <queue>
#PBS -N plasmo-call-ctl
#PBS -l ncpus=1,mem=8GB,walltime=48:00:00
#PBS -l storage=gdata/<proj>+scratch/<proj>
#PBS -j oe
set -euo pipefail
cd /path/to/plasmo-call
./box/bin/pixi run -- snakemake --unlock \
  --profile profiles/pbs --configfile config/config.yaml \
  --configfile config/config.run-local.yaml || true
./box/bin/pixi run -- snakemake \
  --profile profiles/pbs --configfile config/config.yaml \
  --configfile config/config.run-local.yaml \
  --keep-going --rerun-incomplete
```

Submit with `qsub run-controller.pbs`; the controller's own `plasmo-call-ctl.o<id>`
log holds the Snakemake summary. (Use `./box/bin/pixi run` — not
`source box/bin/activate` — so tools are on PATH without `SNAKEMAKE_PROFILE`
contamination.)

## Worker job logs — only failures are kept
Worker jobs send PBS stdout/stderr to `/dev/null`; the jobscript captures each
job to `logs/plasmo.<jobid>.log` and deletes it on success. So `logs/` holds
only the jobs that FAILED — no hundreds of `.o` files from a clean scatter run.
Resource usage for a killed job: `qstat -xf <jobid>`.

## Site config — `storage` is mandatory on Gadi
Fill `profiles/pbs/site.local.yaml` (gitignored). **If `storage` is missing, jobs
submit but die instantly** with `cd: /g/data/<proj>: No such file or directory` —
Gadi compute nodes can't see `/g/data` or `/scratch` unless the job requests them.

```yaml
account:  "<proj>"                    # e.g. pq84
queue:    "normalbw"                  # or your queue
storage:  "gdata/<proj>+scratch/<proj>"
email:    ""
mailon:   "a"
```

## Reference: a pre-indexed READ-ONLY reference works directly
**No copy needed.** Point `config.reference.fasta` straight at the curated,
read-only, already-indexed fasta on `/g/data`.

At parse time the pipeline probes for each index next to the fasta — `.fai`,
`.dict`, and the aligner's sidecar (`.bwt` for `bwa`, `.bwt.2bit.64` for
`bwa-mem2`) — and **defines an index-build rule only for the ones that are
missing**. An index that already exists has no rule producing it, so Snakemake
consumes it as a plain static input and never tries to write it. Nothing is
written to the reference directory, so there is no `ProtectedOutputException`.

The startup banner tells you exactly what it decided:

```
========================================================================
[plasmo-call] reference = /g/data/<proj>/ref/PvP01.fasta
              faidx            found    -> reusing as-is
              dict             found    -> reusing as-is
              bwa index        found    -> reusing as-is
              reference is fully pre-indexed; no index rules scheduled
              (a READ-ONLY reference directory is fine — nothing is written)
========================================================================
```

Three cases, all supported:

| Reference state | What happens |
|---|---|
| Fully indexed, read-only | Nothing built, nothing written — runs as-is |
| Fully indexed, writable  | Nothing built (indices reused, not regenerated) |
| Fresh / unindexed        | Missing indices built next to the fasta, as before |

Only the aligner selected by `config.aligner` gets an index — the unused
aligner's sidecars are never required or built.

**Caveats.**
- A *partially* indexed read-only reference still fails: the missing index has a
  build rule, and the rule can't write to the directory. Either ask the data
  custodian to complete the index, or copy the fasta somewhere writable.
- An existing index is used **as-is and never refreshed**, even if it predates
  the fasta — that is precisely what makes the read-only case work. If the index
  is older than the fasta the banner prints a `WARNING: index file(s) OLDER than
  the FASTA — reused anyway` line. To force a rebuild, delete the stale index
  (on a writable reference) — mtime alone will not trigger one.

## Disk: throwaway intermediates ARE cleaned during the run
The scheduler profiles no longer set `notemp: true`, so Snakemake's `temp()`
markers take effect on cluster runs the same way they do locally. At ~580
samples that is the difference between keeping and discarding hundreds of GB of
scratch mid-run.

Deleted as soon as nothing needs them: `output/mapped_reads/*.bam`,
`*_dupmarked.bam`, recal tables, per-interval GATK/bcftools VCFs, the combined
GVCF, and the bootstrap pass's intermediates.

**Kept** (the deliverables): `output/bam/{sample}.bam` + `.bam.bai`, the
recalibrated BAMs when BQSR is on, per-sample GVCFs, the per-chromosome joint
VCFs, the GenomicsDB workspaces, all three Consensus tiers + the filter stats,
and `output/bam/{sample}_picard_metrics.txt` — the per-sample MarkDuplicates
duplicate-rate metrics, deliberately NOT `temp()` so the cohort's dup rates
survive the run as a QC artefact (no rule consumes the file).

## No jobfs: GenomicsDBImport tmp dir is `/tmp`
**No rule requests node-local scratch (jobfs).** `profiles/pbs/config.yaml`
pins `disk_mb=0` in `default-resources` and no rule overrides it, and
GenomicsDBImport's `--tmp-dir` is the committed default `/tmp`:
```yaml
# config/config.yaml (already the default — nothing to override)
joint_calling:
  genomicsdb:
    tmp_dir: "/tmp"
```

Why the explicit `disk_mb=0` matters: Snakemake's built-in default is
`disk_mb = max(2 * input.size_mb, 1000)`, and `pbs-submit.py` turns any non-zero
`disk_mb` into `,jobfs=<n>MB` on the `-l` string — so *every* job would request
jobfs sized to 2× its input. For `bcftools_caller`, whose input is ALL BAMs,
that is ~2 × 900 × ~2 GB ≈ 3.6 TB at ~900 samples: unschedulable on Gadi
(nodes carry ~400 GB jobfs), and it only surfaces *after* mapping and calling
have finished. The 2-sample test hid this because 2× two small BAMs fits
anywhere. Pinning `disk_mb=0` matches the legacy GATK 3.8 pipeline, which ran
~900 samples requesting no jobfs with Java temp on `/tmp`.

Why `/tmp` is enough: on the real 2-sample run GenomicsDBImport used ~0 bytes
of `--tmp-dir` scratch. Its TileDB workspace is written to the **output path**
on `/g/data` (measured 95 MB for 2 samples, ~1.3 MB per interval), not to
`--tmp-dir`.

**Opting back in** (only if a run actually shows node-local scratch pressure —
the jobfs mechanism in `pbs-submit.py` is kept dormant for exactly this):
set `disk_mb` on that one rule *and* point the tmp dir at the scratch var, since
`$PBS_JOBFS` only exists when the job requested jobfs:
```yaml
# profiles/pbs/config.yaml → set-resources
  genomicsdb_import:
    disk_mb: 100000                   # ~100 GB of $PBS_JOBFS
# config/config.run-local.yaml
joint_calling:
  genomicsdb:
    tmp_dir: "$PBS_JOBFS"
```
Then check real use with `qstat -xf <jobid> | grep jobfs` and tune. On
**SLURM** the same `disk_mb` becomes `--tmp=<n>` (MB) and the scratch is
exposed as `$TMPDIR` / `$SLURM_TMPDIR`, **not** `$PBS_JOBFS` — a `$PBS_JOBFS`
value silently expands to an empty string there.

## JVM heap is derived from `mem_mb` — not left to the JVM
Every GATK and Picard rule passes an **explicit** heap sized from its own
`resources.mem_mb` (`xmx_mb()` in `workflow/rules/common.smk`, default 80%):

* `gatk --java-options "-Xmx<n>m" <Tool> …`
* `_JAVA_OPTIONS="-Xmx<n>m" picard <Tool> …`

Why this matters on a cluster:

* Left alone, the JVM sizes its heap off **physical node RAM** (~1/4). On a
  big shared node that is far more than the job's `mem=` request, so PBS/SLURM
  kills the job for exceeding its allocation instead of throwing a clean Java
  OOM you could diagnose.
* bioconda's `picard` wrapper hard-codes `-Xms512m -Xmx2g` unless an explicit
  `-Xm*` arg or `$_JAVA_OPTIONS` is present. Without the `_JAVA_OPTIONS`
  prefix, MarkDuplicates sits at a 2 GB ceiling no matter what you requested.
  (`_JAVA_OPTIONS` is used rather than `-Xmx` as an argument because the
  wrapper both *checks for* it and the JVM *itself* picks it up — it survives
  wrapper quirks. You'll see `Picked up _JAVA_OPTIONS: …` in the job log;
  that line is the confirmation it took.)

Because the heap tracks `mem_mb`, bumping a rule's `set-resources` entry in the
profile is all it takes — the heap follows. Local (no-profile) runs declare no
`mem_mb`, so they fall back to `LOCAL_FALLBACK_MEM_MB` (4096 → `-Xmx3276m`).

**One deliberate exception: `genomicsdb_import` gets 25% capped at 8192m.**
GenomicsDBImport's heavy memory is **native / off-heap** (TileDB write
buffers), not JVM heap; a large heap crowds the native allocations out of the
`mem=` envelope and OOMs the job. Broad recommends a modest heap regardless of
node size — footprint is controlled by `--batch-size` instead. At
`mem_mb: 32768` that renders as `-Xmx8192m` (the cap), not `-Xmx26214m`.

## Species priors — `knowlesi` is now a committed, data-derived preset
`config/species/knowlesi.yaml` is no longer a TBD stub: it carries
`heterozygosity: 0.005` / `indel_heterozygosity: 0.000625`, derived from the
808-sample monoclonal core-genome Pk cohort (see the file's `source` field).
**Production runs should use `species: knowlesi`** — switch
`config.run-local.yaml` off `knowlesi-provisional`.

Both priors are now passed to **GenotypeGVCFs as well as HaplotypeCaller**. In a
GVCF workflow the final genotype QUAL is computed by GenotypeGVCFs, which has
its own `--heterozygosity` / `--indel-heterozygosity` defaulting to the human
`0.001` — passing the species priors only to HaplotypeCaller left them largely
inert in joint calling.

## Test / provisional species priors
The `config/species/*-provisional.yaml` files are **gitignored** (test-only), so
they are NOT in a fresh clone — recreate on the cluster if you use one, e.g.:
```bash
cat > config/species/knowlesi-provisional.yaml <<'YAML'
species_name: "Plasmodium knowlesi (PROVISIONAL)"
heterozygosity: 0.003
indel_heterozygosity: 0.0017
status: "PROVISIONAL — test only"
YAML
```

## Updating
Use `git pull` in the existing clone — do NOT re-clone. Your env (`box/`,
`.pixi/`), data, `site.local.yaml`, `config.run-local.yaml`, provisional species
files, and `reference/` are all gitignored and survive a pull.
