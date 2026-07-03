# GATK 3.8 → GATK 4 argument changes

Captured during the port from `JacobAFW/Variant_Calling_Pipeline` (GATK 3.8)
to `plasmo-call` (GATK 4.6.2.x). Every entry here was actually hit during the
port — nothing speculative.

## Tool invocation

| GATK 3.8                                  | GATK 4                            |
|-------------------------------------------|-----------------------------------|
| `java -jar GenomeAnalysisTK.jar -T <Tool>`| `gatk <Tool>` (tool is positional, no `-T`) |
| `-o output`                               | `-O output` (capital O)           |
| `--variant` (long) or implicit            | `-V` / `--variant` (still works)  |

`gatk` is a wrapper that picks the bundled jar and a sensible JVM config; no
need to hand-craft `java -Djava.iodir=… -Xms… -Xmx…` lines. Pass JVM flags via
`--java-options "…"` if you really need them.

## Argument naming convention

GATK 3.8 used `camelCase` / `snake_case` mixed; GATK 4 is strict `kebab-case`.
Direct renames seen in the port:

| GATK 3.8                          | GATK 4                                       |
|-----------------------------------|----------------------------------------------|
| `--emitRefConfidence`             | `--emit-ref-confidence`                      |
| `--minPruning`                    | `--min-pruning`                              |
| `--maxNumHaplotypesInPopulation`  | `--max-num-haplotypes-in-population`         |
| `--max_alternate_alleles`         | `--max-alternate-alleles`                    |
| `--heterozygosity`                | `--heterozygosity` (unchanged)               |
| `--indel_heterozygosity`          | `--indel-heterozygosity`                     |
| `-knownSites`                     | `--known-sites`                              |
| `--consensusDeterminationModel`   | (gone — see "Removed" below)                 |
| `-mbq`                            | `-mbq` (short form stable, long is `--min-base-quality-score`) |
| `-DF MappingQualityReadFilter`    | `-DF MappingQualityReadFilter` (stable)      |

## Removed / no replacement (just drop them)

- `--variant_index_type LINEAR` and `--variant_index_parameter 128000`
  GATK 4 auto-indexes GVCF output; the manual index hints are gone.
- `-G Standard` (annotation group)
  GATK 4 has a different default annotation set; `StandardAnnotation` is
  applied implicitly when needed. Drop the flag.
- `--genotyping_mode DISCOVERY`
  Default. Just drop it.
- `-stand_emit_conf`
  Removed — emission threshold is no longer separate from calling threshold.
- `-stand_call_conf`
  Renamed to `--standard-min-confidence-threshold-for-calling`; default is
  usually fine for malaria data.
- `-nt` and `-nct` (multi-threading)
  GATK 4 HaplotypeCaller is **single-threaded** by design — parallelism comes
  from per-sample / per-interval scatter (which we already do via Snakemake).
  CombineGVCFs / GenotypeGVCFs likewise have no `-nt`.
- `-contamination 0.0`
  Long form is now `--contamination-fraction-to-filter`; default is 0.0, so
  dropping the flag entirely is the cleanest move.

## Tools removed entirely (intentional drop in plasmo-call)

- `RealignerTargetCreator` + `IndelRealigner`
  Gone from GATK 4. HaplotypeCaller's local de novo reassembly subsumes
  indel realignment. Just remove the rules; no replacement command is needed.
- `PrintReads` for the "apply BQSR" step
  Replaced by **`ApplyBQSR`**. `BaseRecalibrator` is still the same name.
  Switched in Prompt C.

## BQSR specifics (Prompt C)

Args hit during the BQSR port:

| GATK 3.8                            | GATK 4                                       |
|-------------------------------------|----------------------------------------------|
| `-T BaseRecalibrator`               | `gatk BaseRecalibrator` (positional, no -T)  |
| `-knownSites known.vcf`             | `--known-sites known.vcf` (repeatable)       |
| `-T PrintReads … -BQSR table`       | `gatk ApplyBQSR --bqsr-recal-file table`     |
| `-T SelectVariants --selectType SNP`| `gatk SelectVariants --select-type-to-include SNP` |
| `-T VariantFiltration --filterExpression "EXPR" --filterName "name"` | `gatk VariantFiltration --filter-expression "EXPR" --filter-name name` |

### Quirks worth knowing

- **ApplyBQSR writes `<name>.bai`, not `<name>.bam.bai`.**
  GATK 4 follows the Picard naming convention (drop `.bam`, append `.bai`)
  rather than the samtools convention. Tools that expect the samtools form
  (htslib readers, some VCF callers) silently miss the index. plasmo-call
  adds a tiny `recal_bam_index` rule that drops a `.bam.bai` symlink next
  to it so both lookups work.

- **BaseRecalibrator on tiny data.** With very small inputs (the bundled
  smoke fixture is 9 kb × 2 samples × 15× cov ≈ 270 kb of read data) and a
  single bootstrap SNV in the known-sites mask, BaseRecalibrator still
  produces a usable recal table — but the recalibrated quality scores can
  collapse far enough that HaplotypeCaller stops calling the bootstrap
  variant in the main pass. That's correct GATK behaviour, not a pipeline
  bug; it surfaces because the synthetic reads have a flat Q40 quality
  string that bears no relation to the injected 0.2 % error. Real-world
  WGS depths and quality strings make this a non-issue.

- **`--filter-expression` is JEXL, not bash.** The original
  Plasmodium-tuned expressions (`QD<2.0 || FS>60.0 || ...`) are valid JEXL
  and need to be passed wrapped in single quotes through the Snakemake
  shell directive — Snakemake's default `bash -c` will otherwise try to
  interpret `||` as a shell operator.

- **`SelectVariants --exclude-filtered` is the "keep PASS" gate.** GATK 4
  did not change this from 3.x; just flagging because plasmo-call's
  bootstrap-filter chain uses the three-step `select → filter → select`
  pattern, and the third call's `--exclude-filtered` is the one that
  actually drops the FAILs.

- **`-L reference.bed` is honoured by both BaseRecalibrator and
  ApplyBQSR.** The predecessor passed it to both; plasmo-call does the
  same. Without it, BR scans every contig in the BAM header (slow on real
  references with decoys/alts).

## bcftools arm + consensus (Prompt D)

The bcftools arm and the consensus/concat step aren't GATK renames, but
they came with their own quirks worth capturing.

### bcftools mpileup|call

- `bcftools call -m -Oz -a FORMAT/GQ,FORMAT/GP,INFO/PV4 -v` matches the
  predecessor exactly. `-v` is variants-only output; on synthetic /
  low-coverage data this can yield an empty VCF for a whole chromosome
  and force downstream handling (see below).
- The predecessor scattered per-`CHROMOSOME_INTERVALS` entry (10-way
  split for long chroms, whole for short); plasmo-call preserves the
  same scatter using Snakemake wildcards. The wildcard value contains
  `:` and `-` (e.g. `chr1:1-5000`), so an explicit
  `wildcard_constraints: chromosome = r"[^:/]+"` is needed to keep
  Snakemake from ambiguously matching `bcftools_genotyped_chr1:1-5000.vcf.gz`
  under the per-chromosome concat rule.
- Predecessor's `bam_input_list` used `glob.glob('output/bam_recal/*_recalibrated.bam')`
  at rule-runtime to build the `-b` list. That's fragile — any stray
  BAM in the dir gets picked up. plasmo-call builds the list
  deterministically from `all_sample_bams()` (mode-aware helper in
  common.smk) so both calling arms are guaranteed to consume the same
  set of BAMs.

### Position-based consensus

- The predecessor's method — `bcftools query -f '%CHROM\t%POS\n' bcf.vcf > pos.txt`
  then `bcftools filter -R pos.txt gatk.vcf` — is preserved exactly.
  Two-column (chrom, pos) is a valid `-R` input to `bcftools filter`;
  `-R` treats each pair as a 1-bp region.
- **`bcftools filter -R` fails hard on an empty regions file** with
  "Failed to read the regions" and a non-zero exit — which snakemake's
  `set -euo pipefail` correctly propagates. When the bcftools arm calls
  nothing on a chromosome (common on synthetic flat-Q fixtures, rare on
  real WGS), the pos.txt is empty and the rule dies. plasmo-call guards
  with `if [ -s pos.txt ]` and falls back to `bcftools view --header-only`
  for the consensus VCF, preserving the intersection semantics
  (`empty ∩ anything = empty`) without introducing an allele-aware
  join. The guard fires transparently on real data because pos.txt
  will not be empty there.
- The `-R` intersection is **position-only**, not allele-aware. A GATK
  record at a bcftools position passes even if its ALT allele differs.
  Predecessor behaved the same; do not change without an explicit
  decision.

## Output indexing quirk worth knowing

`CombineGVCFs` in GATK 4.6.x on `osx-arm64` doesn't always drop the `.tbi`
sidecar — the bundled Intel native compression lib is x86-only, GATK falls
back to Java zip, and the index occasionally gets skipped. plasmo-call works
around this by chaining `gatk IndexFeatureFile -I {out}` immediately after
`CombineGVCFs` so `GenotypeGVCFs` always finds an index. Same pattern is
safe to use for `HaplotypeCaller` output if you ever hit it.

## Sanity check at the call site

Every HC flag plasmo-call uses today (rendered from `params.yaml` +
`species/vivax.yaml`):

```
--emit-ref-confidence GVCF
--kmer-size 10 --kmer-size 25 --kmer-size 40
--dont-use-soft-clipped-bases
--min-assembly-region-size 100
--do-not-run-physical-phasing
--base-quality-score-threshold 12
-mbq 5
-DF MappingQualityReadFilter
--heterozygosity 0.0029
--indel-heterozygosity 0.0017
```

All accepted by `gatk-4.6.2.0` without warnings (other than the macOS
IntelGKL native-lib fallback noise above, which is harmless).
