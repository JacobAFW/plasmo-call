# =============================================================================
# workflow/rules/calling_gatk.smk
# GATK4 HaplotypeCaller (-ERC GVCF) → CombineGVCFs → GenotypeGVCFs (per chrom).
#
# Shared params from config["params"]["haplotypecaller"] (Anto's tuned set).
# Per-species priors (--heterozygosity / --indel-heterozygosity) from
# config["species_priors"]. common.smk has already failed loudly if either is
# unset on the selected species.
#
# Differences vs GATK3 predecessor (worth knowing if you compare):
#   * tool name is positional in GATK4 (HaplotypeCaller, not -T HaplotypeCaller)
#   * all flags are kebab-case (--emit-ref-confidence not --emitRefConfidence)
#   * -nt is gone — GATK4 HC is single-threaded by design
#   * --variant_index_type / --variant_index_parameter are gone (auto)
#   * -G Standard is gone (default annotations differ; no longer needed)
#   * -o is now -O (capital), and --variant is -V
#   * RealignerTargetCreator / IndelRealigner are gone entirely — HC's local
#     reassembly subsumes them. Hence no rules for them in this file.
# =============================================================================

# ---- Helpers ---------------------------------------------------------------

def _hc_input_bam(wildcards):
    """Pick the BAM that feeds HaplotypeCaller, gated on the resolved BQSR_MODE.

    off  -> dedup'd / reheadered BAM straight from mapping.smk.
    else -> recalibrated BAM produced by rules/bqsr.smk's ApplyBQSR.
    """
    if BQSR_MODE == "off":
        return f"output/bam/{wildcards.sample}.bam"
    return f"output/bam_recal/{wildcards.sample}_recalibrated.bam"

def _hc_input_bai(wildcards):
    """Companion to _hc_input_bam: list the .bai alongside so snakemake tracks it."""
    if BQSR_MODE == "off":
        return f"output/bam/{wildcards.sample}.bam.bai"
    return f"output/bam_recal/{wildcards.sample}_recalibrated.bam.bai"

# HC_PARAMS is defined in common.smk so bqsr.smk (bootstrap HC) can share it.

# ---- HaplotypeCaller (per sample → GVCF) -----------------------------------

rule haplotype_caller:
    input:
        bam   = _hc_input_bam,
        bai   = _hc_input_bai,
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
    output:
        gvcf = "output/calling/gatk/gvcf/{sample}.g.vcf.gz",
        tbi  = "output/calling/gatk/gvcf/{sample}.g.vcf.gz.tbi",
    params:
        hc_args = HC_PARAMS,
    shell:
        "gatk HaplotypeCaller "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-O {output.gvcf} "
        "{params.hc_args}"

# ---- CombineGVCFs (all samples → one combined GVCF) ------------------------

rule combine_gvcfs:
    input:
        gvcfs = expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz", sample=SAMPLES),
        tbis  = expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz.tbi", sample=SAMPLES),
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
    output:
        gvcf = temp("output/calling/gatk/gvcf/GATK_combined.g.vcf.gz"),
        tbi  = temp("output/calling/gatk/gvcf/GATK_combined.g.vcf.gz.tbi"),
    params:
        v_args = lambda w, input: " ".join(f"-V {g}" for g in input.gvcfs),
    shell:
        # GATK4 CombineGVCFs is supposed to write the .tbi alongside the output,
        # but on osx-arm64 the bundled IntelGKL native lib falls back to Java
        # zip and the index sidecar is occasionally not produced. Re-index
        # explicitly so GenotypeGVCFs always finds a usable .tbi.
        "gatk CombineGVCFs "
        "-R {input.fasta} "
        "{params.v_args} "
        "-O {output.gvcf} && "
        "gatk IndexFeatureFile -I {output.gvcf}"

# ---- GenotypeGVCFs (per chromosome) ----------------------------------------

rule genotype_gvcfs:
    input:
        gvcf  = "output/calling/gatk/gvcf/GATK_combined.g.vcf.gz",
        tbi   = "output/calling/gatk/gvcf/GATK_combined.g.vcf.gz.tbi",
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
    output:
        vcf = "output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz",
        tbi = "output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz.tbi",
    shell:
        "gatk GenotypeGVCFs "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-L {wildcards.chromosome} "
        "-O {output.vcf}"
