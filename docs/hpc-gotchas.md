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

## GenomicsDBImport tmp dir
Point it at node-local scratch — shared filesystems make it crawl:
```yaml
# config.run-local.yaml
joint_calling:
  genomicsdb:
    tmp_dir: "$PBS_JOBFS"
```

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
