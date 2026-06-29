# plasmo-call

A reproducible, portable **GATK4 consensus variant-calling pipeline** for
malaria (*Plasmodium*) short-read data: `fastq → Consensus.vcf.gz`.
Clone → run one install script → run, on **PBS or SLURM**.

> Status: **scaffold / work in progress.** The repository layout, config, and
> scheduler/BQSR design are in place; pipeline rule logic is being ported from
> the predecessor (see below). Steps marked `TODO` are not yet active.

## What it does

Per sample: `bwa`/`bwa-mem2 mem` → sort → Picard `MarkDuplicates` → read-group
reheader → **BQSR** (`BaseRecalibrator` → `ApplyBQSR`) → GATK4 `HaplotypeCaller`
(GVCF) → `CombineGVCFs`/`GenomicsDBImport` → `GenotypeGVCFs`. In parallel,
`bcftools mpileup | call`. The two are combined into a **consensus** call set —
GATK genotypes restricted to bcftools-called positions — for a conservative,
high-confidence VCF.

Built on [vvg-box](https://github.com/vivaxgen/vvg-box) (pixi environment +
PBS/SLURM Snakemake profiles). Predecessor (reference, read-only):
[JacobAFW/Variant_Calling_Pipeline](https://github.com/JacobAFW/Variant_Calling_Pipeline)
— note that was **GATK 3.8**; plasmo-call is a full GATK4 migration and drops the
GATK3 indel-realignment steps (replaced by GATK4 HaplotypeCaller reassembly).

## Install

```bash
git clone https://github.com/JacobAFW/plasmo-call.git
cd plasmo-call
./install.sh        # bootstraps vvg-box, then pixi-installs the tools  (TODO)
```

## Run

Build and validate **locally first (no scheduler)**, then test on each cluster:

| Rung | Mode  | Command                                   |
|------|-------|-------------------------------------------|
| 1    | local | `snakemake --cores N`  (on `test/` data)  |
| 2    | PBS   | use `profiles/pbs/`                        |
| 3    | SLURM | use `profiles/slurm/`                      |

See [`docs/schedulers.md`](docs/schedulers.md) and
[`profiles/README.md`](profiles/README.md).

## Configure

Everything site-specific is config — **no institute paths, accounts, storage, or
emails live in the code.** Edit `config/config.yaml`:

- **Reference** — `reference.fasta` + `reference.bed` (malaria default expected,
  fully configurable). Genome data is **not** committed.
- **BQSR** — `bqsr.mode`: `auto` (blank `known_variants` → bootstrap, else
  known-sites; the choice is logged loudly), `known_sites`, `bootstrap`, or
  `off`. GATK4's `BaseRecalibrator` still requires known sites; for *Plasmodium*
  (no dbSNP) the pipeline can **bootstrap** them (call → hard-filter → reuse).
- **Species priors** — `species:` selects a file in `config/species/`. The
  HaplotypeCaller `--heterozygosity` / `--indel-heterozygosity` priors are
  per-species config, **not** hard-coded.

### ⚠ Species presets

| Species  | Priors            | Status                              |
|----------|-------------------|-------------------------------------|
| vivax    | 0.0029 / 0.0017   | **Documented** (only validated one) |
| knowlesi | unset             | **TBD** — set before use            |
| malariae | unset             | **TBD** — set before use            |
| ovale    | unset             | **TBD** — set before use            |

Non-vivax priors are deliberately left blank — supply values from literature or
estimate them before running those species. The remaining HaplotypeCaller flags
(`--kmer-size 10/25/40`, `--dont-use-soft-clipped-bases`,
`--min-assembly-region-size 100`, `--do-not-run-physical-phasing`,
`--base-quality-score-threshold 12`, `-mbq 5`, `-DF MappingQualityReadFilter`)
are species-agnostic shared defaults in `config/params.yaml`.

## What's included / deliberately excluded

**Included:** workflow code, config templates, species presets, scheduler
profiles, install script, docs.
**Excluded** (see `.gitignore`): all sequencing/variant data (FASTQ/BAM/CRAM/VCF),
reference genomes & indices, run outputs/logs, the pixi env dir, and any
filled-in site profile or secret. `pixi.lock` **is** committed for reproducibility.

## License

[MIT](LICENSE) © 2026 Jacob A. F. Westaway.
