# plasmo-call — PBS profile

Snakemake 8 profile for running plasmo-call on PBS Pro / OpenPBS clusters
(tested lineage: the classic `qsub -N -l ncpus,mem,walltime,storage -A -q -M
-m -j oe -V -S /bin/sh` submission from `JacobAFW/Variant_Calling_Pipeline`,
known-working on Gadi (NCI) and adapted here for the CURRENT rule set).

Packaged self-contained so it can be offered upstream to the vvg-box
authors as a drop-in replacement for `etc/snakemake-profiles/pbspro/`.

## Layout

| File | Committed | Purpose |
| --- | --- | --- |
| `config.yaml`         | yes | Snakemake profile config. Uses `cluster-generic` + per-rule `set-threads` / `set-resources` (memory + walltime). |
| `pbs-submit.py`       | yes | Self-contained qsub wrapper. Reads `site.local.yaml`; nothing site-specific baked in. |
| `pbs-jobscript.sh`    | yes | Minimal jobscript. |
| `site.example.yaml`   | yes | TEMPLATE with placeholders. Copy to `site.local.yaml`. |
| `site.local.yaml`     | **no** (gitignored via `profiles/**/*.local.yaml`) | Your cluster's actual values. |
| `legacy-reference.md` | yes | Original known-working invocation. |

## One-time setup on the HPC

```sh
cp profiles/pbs/site.example.yaml profiles/pbs/site.local.yaml
${EDITOR:-nano} profiles/pbs/site.local.yaml     # fill in account, queue, storage, email
```

Also override the GenomicsDBImport tmp-dir for Gadi (node-local scratch;
GenomicsDBImport writes a lot and shared filesystems make it crawl). Put
this in a run-local config overlay so nothing site-specific goes in the
committed `config/config.yaml`:

```sh
cat > config/config.run-local.yaml <<'YAML'
joint_calling:
  genomicsdb:
    tmp_dir: "$PBS_JOBFS"
YAML
```

(`config/config.run-local.yaml` is gitignored — see `.gitignore`.)

## Submitting

```sh
# Activate vvg-box env, then:
pixi run snakemake \
    --profile profiles/pbs \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --configfile config/config.run-local.yaml
```

The profile injects site.local.yaml values into each `qsub`. What actually
gets submitted for e.g. `haplotype_caller` looks like:

```
qsub -N smk.haplotype_caller.<sample> \
     -l ncpus=3,mem=24576MB,walltime=48:00:00 \
     -A <account> -q <queue> \
     -l storage=<storage> \
     -M <email> -m a \
     -j oe -V -S /bin/sh \
     <jobscript>
```

Nothing site-specific reaches the committed tree — verify with:

```sh
grep -rE 'gdata/|scratch/|@.*\.edu' profiles/pbs/config.yaml \
    profiles/pbs/pbs-submit.py profiles/pbs/site.example.yaml
# ← expect no hits beyond the placeholders inside site.example.yaml
```

## Per-rule resources — where they come from

`config.yaml` in this directory carries the `set-threads` / `set-resources`
map, ported from the legacy `cluster.yaml` (`legacy-reference.md`) and
adapted for the CURRENT rules:

- `bwa_map` — 5 cpu / 50 GB / 24h (unchanged)
- `haplotype_caller` — 3 cpu / 24 GB / 48h (unchanged)
- `genomicsdb_import` — **new**; memory-heavy (writes many TileDB fragments) — 4 cpu / 32 GB / 24h
- `genotype_gvcfs` — **now per-INTERVAL**; modest 2 cpu / 16 GB / 12h
- `bcftools_caller` — 2 cpu / 16 GB / 48h (unchanged)
- `combine_gvcfs` — 3 cpu / 24 GB / 48h (legacy backend, only for small cohorts)
- `filter_hard` — 2 cpu / 16 GB / 6h (GATK VariantFiltration on the cohort)
- All concat/consensus/light-filter/stats rules — 2 cpu / 8 GB / 6h

If your cluster differs, override via `--set-threads` / `--set-resources`
on the CLI, or edit `config.yaml`.

## Not tested here

This environment has no PBS. Everything ships parse-verified via
`snakemake -n --profile profiles/pbs …` and the actual submission is
Jacob's to run on the HPC. Report back and we'll adjust.
