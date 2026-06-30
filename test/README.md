# test/

Small bundled smoke-test dataset for the LOCAL (no-scheduler) run — rung 1 of
the test ladder. Tiny synthetic FASTQs + a small reference so the pipeline
completes fast end-to-end.

## What's tracked

| File                          | Purpose                                              |
|-------------------------------|------------------------------------------------------|
| `generate-fixtures.py`        | Stdlib-only generator: 2 synthetic chroms + 2 samples |
| `config.yaml`                 | Base overlay: source_dir, reference, aligner=bwa, mode=off |
| `config-bootstrap.yaml`       | BQSR mode = bootstrap (no external VCF)              |
| `config-known-sites.yaml`     | BQSR mode = known_sites (points at `test/known_sites/`) |
| `config-auto-empty.yaml`      | BQSR mode = auto, no known_variants → expects bootstrap |
| `config-auto-populated.yaml`  | BQSR mode = auto, known_variants set → expects known_sites |
| `README.md`                   | This file                                            |

Everything else (`reference/`, `fastq/`, `output/`, `known_sites/`) is generated
and gitignored.

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

bcftools arm + consensus + concat are deferred to a later milestone.

## BQSR acceptance (runs all four modes)

```bash
pixi run bqsr-acceptance     # bootstrap → known_sites → auto×2 banners
```

That script:
1. Runs `bqsr.mode=bootstrap` end-to-end on the synthetic data
   (27 jobs: bootstrap HC + filter → BR/ApplyBQSR → main HC chain).
2. Copies the bootstrap-produced filtered VCFs into `test/known_sites/`.
3. Re-runs with `bqsr.mode=known_sites` against those VCFs
   (21 jobs — no bootstrap stage).
4. Dry-runs `bqsr.mode=auto` with both empty and populated `known_variants`
   to surface the resolver's banner (`>> BQSR mode = bootstrap (from auto ...)`
   vs `>> BQSR mode = known_sites (from auto ...)`).

You can also run each mode by hand:

```bash
pixi run -- snakemake --cores 4 \
  --snakefile workflow/Snakefile \
  --configfile config/config.yaml test/config.yaml test/config-bootstrap.yaml
```

## Dataset shape (deterministic, seed 20260629)

- 2 chromosomes: `chr1` (5 kb), `chr2` (4 kb)
- 2 samples (`sample01`, `sample02`), ~15× coverage each
- 100 bp paired-end reads, ~300 bp insert, 0.2% per-base error rate

The error rate is intentional — pure-identity reads would produce an empty
final VCF and obscure regressions in HaplotypeCaller params.
