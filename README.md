# plasmo-call

A reproducible, portable **GATK4 consensus variant-calling pipeline** for malaria
(*Plasmodium*) short-read data. Clone → run one install script → run, locally or
on **PBS / SLURM**. Validated end-to-end on real *P. knowlesi* whole-genome data.

## What it does

`fastq → consensus VCFs (+ filtered tiers)`, calling variants with **two callers
and taking their consensus** for a conservative, high-confidence set.

Per sample: `bwa` / `bwa-mem2 mem` → sort → Picard `MarkDuplicates` → read-group
reheader → **BQSR** (`BaseRecalibrator` → `ApplyBQSR`) → GATK4 `HaplotypeCaller`
(GVCF). GVCFs are consolidated with **GenomicsDBImport** (default; scales to
thousands of samples) or `CombineGVCFs`, then genotyped **per interval** with
`GenotypeGVCFs` and gathered per chromosome. In parallel, `bcftools mpileup | call`
runs per interval. The two arms are combined into a **consensus** — GATK genotypes
restricted to bcftools-called positions.

Built on [vvg-box](https://github.com/vivaxgen/vvg-box) (pixi environment +
auto-detecting PBS/SLURM Snakemake profiles). Predecessor (reference, read-only):
[JacobAFW/Variant_Calling_Pipeline](https://github.com/JacobAFW/Variant_Calling_Pipeline)
— that was **GATK 3.8**; plasmo-call is a full GATK4 migration and drops the GATK3
indel-realignment steps (subsumed by GATK4 HaplotypeCaller local reassembly).

## Install

```bash
git clone https://github.com/JacobAFW/plasmo-call.git
cd plasmo-call
./install.sh        # bootstraps vvg-box into ./box, then pixi-installs the tools
```

## Run

Validate **locally first (no scheduler)**, then run on a cluster:

| Rung | Mode  | How                                                        |
|------|-------|------------------------------------------------------------|
| 1    | local | `pixi run smoke-test` (bundled data) · `pixi run run-local` |
| 2    | PBS   | fill `profiles/pbs/site.local.yaml`, run with `profiles/pbs/`   |
| 3    | SLURM | fill `profiles/slurm/site.local.yaml`, run with `profiles/slurm/` |

The scheduler profiles ship with the submission logic and per-rule resources; you
supply your site's account/queue/storage/scratch in a gitignored `site.local.yaml`
(copy the `site.example.yaml` template). See [`docs/schedulers.md`](docs/schedulers.md)
and [`profiles/README.md`](profiles/README.md).

## Configure

Everything site-specific is config — **no institute paths, accounts, storage, or
emails live in the code.** Edit `config/config.yaml`:

- **Reference** — `reference.fasta` + `reference.bed` (malaria default expected,
  fully configurable). Missing `.fai`/`.dict`/aligner indices are built
  automatically next to the fasta; any that **already exist are reused as-is**,
  and no rule is scheduled to rebuild them. A curated, **read-only, pre-indexed**
  reference therefore works directly — nothing is written to it, so no copy to a
  writable location is needed (see [`docs/hpc-gotchas.md`](docs/hpc-gotchas.md)).
  Genome data is not committed.
- **Species** — `species:` selects `config/species/<name>.yaml`, which supplies the
  per-species HaplotypeCaller `--heterozygosity` / `--indel-heterozygosity` priors.
  It is **required** — the run stops if the selected species' priors are unset, so
  data is never silently called with another species' priors. Shared, species-agnostic
  HaplotypeCaller flags live in `config/params.yaml`.
- **BQSR** — `bqsr.mode`: `auto` (blank `known_variants` → bootstrap, else
  known-sites; the choice is logged), `known_sites`, `bootstrap`, or `off`. GATK4's
  `BaseRecalibrator` still requires known sites; for *Plasmodium* (no dbSNP) the
  pipeline can **bootstrap** them (call → hard-filter → reuse).
- **Joint calling** — `joint_calling.consolidation`: `genomicsdb` (default,
  scale-first) or `combine_gvcfs` (small cohorts).
- **Filtering** — `filter.light.min_qual` and `filter.hard.{qd,fs,mq}` (see Outputs).

## Outputs

The consensus is emitted **raw plus two filter tiers**, side by side — the raw set
is never replaced, so any tier can be re-derived cheaply without re-calling:

| Output | Filter | Use |
|--------|--------|-----|
| `Consensus.vcf.gz` | none (GATK ∩ bcftools) | reproducible base |
| `Consensus.light.vcf.gz` | `QUAL ≥ 30` | permissive / discovery |
| `Consensus.hard.vcf.gz` | `QD < 15`, `FS > 1`, `MQ < 40` | stringent, popgen/GWAS-grade |

The **hard** tier follows the Sanger/MalariaGEN, non-model-organism recipe (QD/FS/MQ
only). GATK's rank-sum filters (`MQRankSum`, `ReadPosRankSum`, `SOR`) are deliberately
*not* applied: they are diploid-tuned tests defined only at heterozygous sites, and
*Plasmodium* is haploid (apparent hets are mixed/polyclonal infections, handled
downstream via Fws/COI). They are shown in the report for transparency but not filtered.

`output/filter/stats.{txt,tsv}` is a **funnel report** — record counts at
raw → light → hard by SNP/indel, plus per-site QUAL/QD/FS/MQ (and the rank-sum/SOR)
distributions with both tiers' thresholds marked, so thresholds can be reviewed and
retuned per dataset.

> Downstream population-genetics QC — sample/site missingness, MAF, Fws/clonality,
> and core-genome / hypervariable-region masking — is intentionally left to the
> analysis workflow, not this caller.

## What's included / deliberately excluded

**Included:** workflow code, config templates, scheduler profiles, install script,
docs. **Excluded** (see `.gitignore`): all sequencing/variant data
(FASTQ/BAM/CRAM/VCF), reference genomes & indices, run outputs/logs, the pixi env
and vvg-box `box/`, and any filled-in `site.local.yaml` or secret. `pixi.lock` **is**
committed for reproducibility.

## License

[MIT](LICENSE) © 2026 Jacob A. F. Westaway.
