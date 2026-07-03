# =============================================================================
# workflow/rules/consensus.smk
# Position-based intersection of the two arms, per chromosome, then concat.
#
# EXACT semantics preserved from the GATK3 predecessor:
#   1. `bcftools query -f '%CHROM\t%POS\n' bcftools.vcf.gz > pos.txt`
#      pulls a 2-column (chrom, pos) file from the bcftools calls.
#   2. `bcftools filter -R pos.txt gatk.vcf.gz` keeps only GATK records whose
#      (chrom, pos) matches — a conservative POSITION intersection.
#   3. Concat the per-chrom consensus VCFs → Consensus.vcf.gz.
#
# This is NOT allele-aware — a GATK record at a bcftools position passes even
# if the alleles differ. Do NOT change without an explicit decision, per the
# prompt's guidance.
#
# Notes on `bcftools filter -R`:
#   * accepts either 3-col BED or 2-col (chrom, pos); we use the 2-col form
#     the predecessor produced.
#   * requires the input VCF to be indexed (.tbi) — which GenotypeGVCFs
#     already provides via calling_gatk.smk.
# =============================================================================

localrules: concat_consensus

# ---- Per-chromosome position intersection ----------------------------------

rule consensus_of_vcfs:
    input:
        bcftools     = "output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz",
        bcftools_tbi = "output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz.tbi",
        gatk         = "output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz",
        gatk_tbi     = "output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz.tbi",
    output:
        pos = temp("output/calling/consensus/{chromosome}.pos.txt"),
        vcf = temp("output/calling/consensus/{chromosome}_consensus.vcf.gz"),
        tbi = temp("output/calling/consensus/{chromosome}_consensus.vcf.gz.tbi"),
    shell:
        # If bcftools called nothing on this chromosome, an empty positions file
        # would make `bcftools filter -R` bail ("Failed to read the regions"),
        # so short-circuit to a header-only consensus VCF. Semantics preserved:
        # empty ∩ anything = empty.
        r"bcftools query -f '%CHROM\t%POS\n' {input.bcftools} > {output.pos} && "
        "if [ -s {output.pos} ]; then "
        "  bcftools filter -R {output.pos} -Oz -o {output.vcf} {input.gatk}; "
        "else "
        "  bcftools view --header-only -Oz -o {output.vcf} {input.gatk}; "
        "fi && "
        "bcftools index -t -o {output.tbi} {output.vcf}"

# ---- Final concat -> Consensus.vcf.gz --------------------------------------

rule concat_consensus:
    input:
        vcfs = expand("output/calling/consensus/{chromosome}_consensus.vcf.gz",
                      chromosome=CHROMOSOME),
        tbis = expand("output/calling/consensus/{chromosome}_consensus.vcf.gz.tbi",
                      chromosome=CHROMOSOME),
    output:
        vcf = "output/calling/consensus/Consensus.vcf.gz",
        tbi = "output/calling/consensus/Consensus.vcf.gz.tbi",
    shell:
        "bcftools concat -Oz -o {output.vcf} {input.vcfs} && "
        "bcftools index -t -o {output.tbi} {output.vcf}"
