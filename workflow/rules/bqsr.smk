# =============================================================================
# workflow/rules/bqsr.smk
# BaseRecalibrator + ApplyBQSR, plus the bootstrap pass that derives a
# known-sites VCF from the data itself when none is supplied.
#
# Mode wiring (resolved in common.smk → BQSR_MODE):
#   off          -> bqsr.smk is included but no rule fires; main HC reads
#                   output/bam/<sample>.bam directly.
#   known_sites  -> BaseRecalibrator/ApplyBQSR with the user's VCF(s).
#   bootstrap    -> bootstrap_* rules first produce output/bqsr/bootstrap/
#                   snvs_pass.vcf.gz + indels_pass.vcf.gz; those then feed
#                   BaseRecalibrator/ApplyBQSR.
#
# GATK3 → GATK4 renames hit in this file:
#   PrintReads -BQSR              -> ApplyBQSR --bqsr-recal-file
#   -knownSites                   -> --known-sites
#   --filterExpression / Name     -> --filter-expression / --filter-name
#   SelectVariants --selectType   -> SelectVariants --select-type-to-include
# =============================================================================

# ---- BaseRecalibrator (per sample → recal table) ---------------------------

rule base_recalibrator:
    input:
        bam         = "output/bam/{sample}.bam",
        bai         = "output/bam/{sample}.bam.bai",
        fasta       = REF_FASTA,
        fai         = REF_FAI,
        dict_       = REF_DICT,
        bed         = REF_BED,
        known_sites = known_sites_vcfs(),
        known_tbis  = known_sites_tbis(),
    output:
        table = temp("output/bam_recal/{sample}.recal.table"),
    params:
        ks_args = lambda w, input: " ".join(f"--known-sites {v}" for v in input.known_sites),
    shell:
        "gatk BaseRecalibrator "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-L {input.bed} "
        "{params.ks_args} "
        "-O {output.table}"

# ---- ApplyBQSR (per sample → recalibrated BAM) -----------------------------

rule apply_bqsr:
    input:
        bam   = "output/bam/{sample}.bam",
        bai   = "output/bam/{sample}.bam.bai",
        table = "output/bam_recal/{sample}.recal.table",
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
        bed   = REF_BED,
    output:
        bam = "output/bam_recal/{sample}_recalibrated.bam",
        # ApplyBQSR writes the BAI alongside (Picard-style .bai, not .bam.bai).
        bai = "output/bam_recal/{sample}_recalibrated.bai",
    shell:
        "gatk ApplyBQSR "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-L {input.bed} "
        "--bqsr-recal-file {input.table} "
        "-O {output.bam}"

# A samtools-style .bam.bai sidecar makes downstream tools happy regardless
# of whether they expect the GATK or Picard naming convention.
rule recal_bam_index:
    input:
        bam = "output/bam_recal/{sample}_recalibrated.bam",
        bai = "output/bam_recal/{sample}_recalibrated.bai",
    output:
        bai = "output/bam_recal/{sample}_recalibrated.bam.bai",
    shell:
        "ln -sf $(basename {input.bai}) {output.bai}"

# ============================================================================
# Bootstrap pass: only fires when BQSR_MODE == 'bootstrap'.
#
# Why a separate haplotype_caller rule:
#   * the main `haplotype_caller` in calling_gatk.smk reads
#     output/bam_recal/<sample>_recalibrated.bam (when mode != off);
#   * a single rule that switched its own input file based on mode would
#     loop back on itself (recal BAM ← recal table ← known-sites ← HC ← ...
#     ← recal BAM). So this bootstrap HC explicitly targets
#     output/bam/<sample>.bam — the dedup'd, reheadered, non-recal BAM.
# ============================================================================

rule bootstrap_haplotype_caller:
    input:
        bam   = "output/bam/{sample}.bam",
        bai   = "output/bam/{sample}.bam.bai",
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
    output:
        gvcf = temp("output/bqsr/bootstrap/gvcf/{sample}.g.vcf.gz"),
        tbi  = temp("output/bqsr/bootstrap/gvcf/{sample}.g.vcf.gz.tbi"),
    params:
        hc_args = HC_PARAMS,
    shell:
        "gatk HaplotypeCaller "
        "-R {input.fasta} "
        "-I {input.bam} "
        "-O {output.gvcf} "
        "{params.hc_args}"

rule bootstrap_combine_gvcfs:
    input:
        gvcfs = expand("output/bqsr/bootstrap/gvcf/{sample}.g.vcf.gz", sample=SAMPLES),
        tbis  = expand("output/bqsr/bootstrap/gvcf/{sample}.g.vcf.gz.tbi", sample=SAMPLES),
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
    output:
        gvcf = temp("output/bqsr/bootstrap/combined.g.vcf.gz"),
        tbi  = temp("output/bqsr/bootstrap/combined.g.vcf.gz.tbi"),
    params:
        v_args = lambda w, input: " ".join(f"-V {g}" for g in input.gvcfs),
    shell:
        # Index explicitly — same osx-arm64 IntelGKL quirk as the main pipeline.
        "gatk CombineGVCFs "
        "-R {input.fasta} "
        "{params.v_args} "
        "-O {output.gvcf} && "
        "gatk IndexFeatureFile -I {output.gvcf}"

rule bootstrap_genotype:
    input:
        gvcf  = "output/bqsr/bootstrap/combined.g.vcf.gz",
        tbi   = "output/bqsr/bootstrap/combined.g.vcf.gz.tbi",
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
    output:
        vcf = temp("output/bqsr/bootstrap/genotyped.vcf.gz"),
        tbi = temp("output/bqsr/bootstrap/genotyped.vcf.gz.tbi"),
    shell:
        # No -L: bootstrap operates on the whole genome (one cohort VCF).
        "gatk GenotypeGVCFs "
        "-R {input.fasta} "
        "-V {input.gvcf} "
        "-O {output.vcf}"

# ---- Hard-filter to high-confidence known sites (SNVs + indels) ------------
# Pipeline per type: SelectVariants (type) → VariantFiltration (mark FAILs)
# → SelectVariants --exclude-filtered (drop them). Output is bgzipped + tbi.

def _bootstrap_filter_rule(var_type: str, filter_expr_key: str, out_basename: str):
    """Factory for the two near-identical SNV/indel filter rules."""
    return var_type, filter_expr_key, out_basename

rule bootstrap_filter_snvs:
    input:
        vcf   = "output/bqsr/bootstrap/genotyped.vcf.gz",
        tbi   = "output/bqsr/bootstrap/genotyped.vcf.gz.tbi",
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
    output:
        vcf = "output/bqsr/bootstrap/snvs_pass.vcf.gz",
        tbi = "output/bqsr/bootstrap/snvs_pass.vcf.gz.tbi",
    params:
        expr     = config["bqsr"]["bootstrap"]["snv_filter"],
        var_type = "SNP",
    shell:
        "tmp=$(mktemp -d); "
        "gatk SelectVariants -R {input.fasta} -V {input.vcf} "
        "  --select-type-to-include {params.var_type} -O $tmp/typed.vcf.gz && "
        "gatk VariantFiltration -R {input.fasta} -V $tmp/typed.vcf.gz "
        "  --filter-expression '{params.expr}' --filter-name BootstrapFail "
        "  -O $tmp/marked.vcf.gz && "
        "gatk SelectVariants -R {input.fasta} -V $tmp/marked.vcf.gz "
        "  --exclude-filtered -O {output.vcf} && "
        "gatk IndexFeatureFile -I {output.vcf} && "
        "rm -rf $tmp"

rule bootstrap_filter_indels:
    input:
        vcf   = "output/bqsr/bootstrap/genotyped.vcf.gz",
        tbi   = "output/bqsr/bootstrap/genotyped.vcf.gz.tbi",
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
    output:
        vcf = "output/bqsr/bootstrap/indels_pass.vcf.gz",
        tbi = "output/bqsr/bootstrap/indels_pass.vcf.gz.tbi",
    params:
        expr     = config["bqsr"]["bootstrap"]["indel_filter"],
        var_type = "INDEL",
    shell:
        "tmp=$(mktemp -d); "
        "gatk SelectVariants -R {input.fasta} -V {input.vcf} "
        "  --select-type-to-include {params.var_type} -O $tmp/typed.vcf.gz && "
        "gatk VariantFiltration -R {input.fasta} -V $tmp/typed.vcf.gz "
        "  --filter-expression '{params.expr}' --filter-name BootstrapFail "
        "  -O $tmp/marked.vcf.gz && "
        "gatk SelectVariants -R {input.fasta} -V $tmp/marked.vcf.gz "
        "  --exclude-filtered -O {output.vcf} && "
        "gatk IndexFeatureFile -I {output.vcf} && "
        "rm -rf $tmp"
