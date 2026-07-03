# =============================================================================
# workflow/rules/calling_bcftools.smk
# bcftools arm of the consensus caller: per-interval mpileup|call → per-chrom
# concat. Uses the same mode-aware BAM list as the GATK arm
# (all_sample_bams() in common.smk).
#
# Ported verbatim from the GATK3 predecessor apart from:
#   * bam list is built from all_sample_bams() (deterministic — the predecessor
#     used glob.glob on disk, which is fragile if any stray *_recalibrated.bam
#     is present),
#   * threads and call annotations flow through config.params.bcftools,
#   * `bcftools index --threads N -t` produces the .tbi in one shell block.
# =============================================================================

# ---- BAM list for mpileup -b ------------------------------------------------
# localrule: just a text file with one BAM path per line.
localrules: bam_input_list

rule bam_input_list:
    input:
        bams = all_sample_bams(),
        bais = all_sample_bais(),
    output:
        lst = temp("output/calling/bcftools/input_bam_files.list"),
    run:
        with open(output.lst, "w") as fh:
            for b in input.bams:
                fh.write(f"{b}\n")

# ---- mpileup | call per genomic interval ------------------------------------

rule bcftools_caller:
    input:
        bams  = all_sample_bams(),
        bais  = all_sample_bais(),
        lst   = "output/calling/bcftools/input_bam_files.list",
        fasta = REF_FASTA,
        fai   = REF_FAI,
    output:
        vcf = temp("output/calling/bcftools/intervals/bcftools_genotyped_intervals_{intervals}.vcf.gz"),
        tbi = temp("output/calling/bcftools/intervals/bcftools_genotyped_intervals_{intervals}.vcf.gz.tbi"),
    threads: config["params"]["bcftools"]["mpileup_threads"]
    params:
        annots = config["params"]["bcftools"]["call_annotations"],
    shell:
        "bcftools mpileup --threads {threads} "
        "  -f {input.fasta} "
        "  -b {input.lst} "
        "  -r {wildcards.intervals} "
        "| bcftools call --threads {threads} "
        "    -m -Oz "
        "    -a {params.annots} "
        "    -v -o {output.vcf} && "
        "bcftools index --threads {threads} -t -o {output.tbi} {output.vcf}"

# ---- Concat a chromosome's intervals into one per-chrom VCF -----------------

def _bcftools_intervals_for_chrom(wildcards):
    """Return every interval VCF whose CHROMOSOME_INTERVALS entry starts with 'chr:'."""
    prefix = f"{wildcards.chromosome}:"
    return [f"output/calling/bcftools/intervals/bcftools_genotyped_intervals_{iv}.vcf.gz"
            for iv in CHROMOSOME_INTERVALS if iv.startswith(prefix)]

def _bcftools_interval_tbis_for_chrom(wildcards):
    return [f"{v}.tbi" for v in _bcftools_intervals_for_chrom(wildcards)]

rule concat_bcftools:
    input:
        vcfs = _bcftools_intervals_for_chrom,
        tbis = _bcftools_interval_tbis_for_chrom,
    output:
        vcf = temp("output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz"),
        tbi = temp("output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz.tbi"),
    shell:
        "bcftools concat -Oz -o {output.vcf} {input.vcfs} && "
        "bcftools index -t -o {output.tbi} {output.vcf}"
