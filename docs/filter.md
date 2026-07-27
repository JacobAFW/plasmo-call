# plasmo-call — post-calling filter tiers

Prompt H redesigns the filter layer as **three tiers** consumed by anything
downstream. All three VCFs are emitted side by side; the raw is never
replaced or hidden.

## Tiers

| Tier   | Path                                              | Filter                                                | Source                                                                 |
| ------ | ------------------------------------------------- | ----------------------------------------------------- | ---------------------------------------------------------------------- |
| RAW    | `output/calling/consensus/Consensus.vcf.gz`       | none — GATK ∩ bcftools intersection only              | plasmo-call consensus arm (Prompt D)                                   |
| LIGHT  | `output/calling/consensus/Consensus.light.vcf.gz` | `QUAL >= light.min_qual`                              | vivaxgen-equivalent floor                                              |
| HARD   | `output/calling/consensus/Consensus.hard.vcf.gz`  | `QD >= hard.qd && FS <= hard.fs && MQ >= hard.mq`     | JacobAFW/Pk_Malaysian_Population_Genetics `05_Variant_filtering`       |

Plus:

- `output/filter/stats.txt` — human-readable **funnel report**: raw → light →
  hard counts split by SNP / indel, with per-site distributions of QUAL, QD,
  FS, MQ, MQRankSum, ReadPosRankSum, SOR — LIGHT and HARD thresholds
  annotated so you can see where each cut lands on the real distribution.
- `output/filter/stats.tsv` — per-site machine-readable table (one row per
  raw site; `IN_LIGHT` / `IN_HARD` flags).

The filter/stats rules consume ONLY the raw `Consensus.vcf.gz`. Amending a
threshold in `config.filter` and re-running regenerates the tiers + report
in seconds — no re-calling.

## Config

```yaml
filter:
  light:
    min_qual: 30           # LIGHT: QUAL floor
  hard:
    qd: 15.0               # HARD: QD < qd → fail
    fs:  1.0               #       FS > fs → fail
    mq: 40.0               #       MQ < mq → fail
```

All thresholds are overridable. No `enabled: true|false` toggle — both
tiers are always produced, because the funnel report needs to show what
each filter drops.

## Why HARD is QD / FS / MQ only

Deliberately omitted from the HARD filter:

| Annotation      | Why NOT filtered on                                                                                          |
| --------------- | ------------------------------------------------------------------------------------------------------------ |
| MQRankSum       | Rank-sum test defined only at **HET sites** — compares MQ between REF-supporting and ALT-supporting reads.   |
| ReadPosRankSum  | Rank-sum test defined only at **HET sites** — compares read position between REF and ALT.                    |
| SOR             | Strand-odds-ratio calibrated on mammalian diploid data.                                                      |

Plasmodium is **haploid**. plasmo-call calls at diploid ploidy so that
mixed / polyclonal infections show as heterozygous — the mixture level
(Fws / COI) is then quantified downstream. But the "heterozygotes" here are
not biological heterozygotes: they're mixed-lineage reads at a single
locus, and the rank-sum thresholds calibrated on diploid heterozygotes
don't transfer. Per GATK's non-model-organism guidance and the
Sanger / MalariaGEN malaria variant-filter recipes those annotations are
useful for review but not appropriate for hard filtering.

The **funnel report shows all three** so you can inspect them for
diagnostic patterns, but the HARD tier does not gate on them.

## Why LIGHT has no DP filter

The vivaxgen `--mindepth` cut is a **per-sample FORMAT/DP** filter applied
in the single-sample report step; it is not an INFO/DP site cut on the
joint cohort VCF. Applying an INFO/DP floor on the joint VCF would drop
sites purely because a few samples in the cohort have low local coverage,
which is a per-sample QC decision, not a cohort-level variant-quality
decision. Per-sample depth belongs in the downstream popgen tier, out of
scope here.

## Sources

- **Jacob's published filter (HARD tier):**
  `JacobAFW/Pk_Malaysian_Population_Genetics`, script
  `05_Variant_filtering` — QD < 15, FS > 1, MQ < 40 via GATK
  VariantFiltration then PASS-only extraction. Same threshold for SNPs and
  indels.
- **vivaxgen (LIGHT tier floor):** the vivaxgen malaria calling pipeline
  applies `QUAL >= 30` as its joint-VCF site floor before per-sample cuts.
- **Rank-sum / SOR rationale:**
  - GATK's "Hard-filtering germline short variants" best-practice doc —
    explicitly documents MQRankSum / ReadPosRankSum / SOR as defined at
    het sites and provides the mammalian-tuned thresholds.
  - GATK's "How to filter variants for a non-model organism" — advises
    dropping the rank-sum tests for non-diploid callsets.
  - Sanger / MalariaGEN malaria variant-filter recipes — QD / FS / MQ
    only for the hard tier on Plasmodium cohorts.

## Amending thresholds

`config.filter` is the single source of truth. Edit the threshold, then:

```
pixi run snakemake --cores <N> --rerun-triggers mtime \
  --snakefile workflow/Snakefile --configfile config/config.yaml
```

Snakemake will re-run only `filter_light`, `filter_hard`, and
`consensus_stats` — the calling side of the DAG is not touched.
