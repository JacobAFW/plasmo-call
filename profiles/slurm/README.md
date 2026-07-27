# plasmo-call — SLURM profile

Snakemake 8 profile for running plasmo-call on SLURM clusters. Mirrors the
PBS profile in shape (same submit-script contract, same per-rule
threads / mem / walltime map) so the pipeline costs the same resources
whichever scheduler you're on.

## Layout

| File | Committed | Purpose |
| --- | --- | --- |
| `config.yaml`         | yes | Snakemake profile config. Uses `cluster-generic` + per-rule `set-threads` / `set-resources`. |
| `slurm-submit.py`     | yes | Self-contained sbatch wrapper. Reads `site.local.yaml`. |
| `slurm-jobscript.sh`  | yes | Minimal jobscript. |
| `site.example.yaml`   | yes | TEMPLATE with placeholders. Copy to `site.local.yaml`. |
| `site.local.yaml`     | **no** (gitignored via `profiles/**/*.local.yaml`) | Your cluster's actual values. |

## One-time setup on the HPC

```sh
cp profiles/slurm/site.example.yaml profiles/slurm/site.local.yaml
${EDITOR:-nano} profiles/slurm/site.local.yaml   # partition, account, qos, email
```

If your cluster exposes per-job node-local scratch (usually `$SLURM_TMPDIR`
or `$TMPDIR`), point GenomicsDBImport at it — same pattern as PBS:

```sh
cat > config/config.run-local.yaml <<'YAML'
joint_calling:
  genomicsdb:
    tmp_dir: "$SLURM_TMPDIR"
YAML
```

(`config/config.run-local.yaml` is gitignored — see `.gitignore`.)

## Submitting

```sh
pixi run snakemake \
    --profile profiles/slurm \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --configfile config/config.run-local.yaml
```

What actually gets submitted for e.g. `haplotype_caller`:

```
sbatch --parsable -J smk.haplotype_caller.<sample> \
       --cpus-per-task=3 --mem=24576M --time=48:00:00 \
       -p <partition> -A <account> [--qos=<qos>] \
       [--mail-user=<email> --mail-type=FAIL] \
       -o slurm-%j.out -e slurm-%j.err \
       <jobscript>
```

## Not tested here

No SLURM in this environment; the profile parses via `snakemake -n` and the
actual submission is Jacob's to run on the cluster.
