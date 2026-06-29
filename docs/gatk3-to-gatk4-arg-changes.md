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
  (We're not wiring BQSR until the next prompt, but flagging this here.)

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
