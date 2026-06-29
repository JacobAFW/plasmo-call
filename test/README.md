# test/

Small bundled smoke-test dataset for the LOCAL (no-scheduler) run — rung 1 of
the test ladder. Tiny synthetic FASTQs + a small reference so the pipeline
completes fast end-to-end.

## What's tracked

| File                     | Purpose                                              |
|--------------------------|------------------------------------------------------|
| `generate-fixtures.py`   | Stdlib-only generator: 2 synthetic chroms + 2 samples |
| `config.yaml`            | Overlay config (source_dir, reference, aligner=bwa)  |
| `README.md`              | This file                                            |

Everything else (`reference/`, `fastq/`, `output/`) is generated and gitignored.

## Run the smoke test

```bash
source box/bin/activate                       # vvg-box env (gets pixi onto PATH)
pixi run smoke-test                           # generate fixtures + end-to-end pipeline
```

Or by hand:

```bash
pixi run -- python test/generate-fixtures.py  # produces test/reference/ + test/fastq/
pixi run -- snakemake --cores 4 \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --configfile test/config.yaml
```

## What the smoke test exercises

- Reference prep: `samtools faidx`, `picard CreateSequenceDictionary`, `bwa index`
- Mapping: `bwa mem` → `samtools sort` (coord) → Picard `MarkDuplicates` → RG reheader → index
- Calling: GATK4 `HaplotypeCaller -ERC GVCF` (Anto params + vivax priors)
  → `CombineGVCFs` → `GenotypeGVCFs` per chromosome

BQSR, bcftools arm, and consensus + concat are deferred to later milestones.

## Dataset shape (deterministic, seed 20260629)

- 2 chromosomes: `chr1` (5 kb), `chr2` (4 kb)
- 2 samples (`sample01`, `sample02`), ~15× coverage each
- 100 bp paired-end reads, ~300 bp insert, 0.2% per-base error rate

The error rate is intentional — pure-identity reads would produce an empty
final VCF and obscure regressions in HaplotypeCaller params.
