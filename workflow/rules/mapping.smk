# =============================================================================
# workflow/rules/mapping.smk
# Reference prep + read mapping + dedup + RG reheader + index.
#
# Ported from the GATK3 predecessor. Two deliberate simplifications:
#   * single sort (coord) out of the aligner — the original sorted by name
#     then re-sorted by coord; Picard MarkDuplicates works on coord-sorted
#     directly, so the name-sort hop is dropped.
#   * aligner is a config toggle (bwa | bwa-mem2). Both produce identical
#     SAM/BAM at the command-line level, so the shell line is shared.
# =============================================================================

# ---- Reference indices ------------------------------------------------------

rule samtools_faidx:
    input:
        REF_FASTA
    output:
        REF_FAI
    shell:
        "samtools faidx {input}"

rule sequence_dict:
    input:
        REF_FASTA
    output:
        REF_DICT
    shell:
        "picard CreateSequenceDictionary R={input} O={output}"

rule bwa_index:
    input:
        REF_FASTA
    output:
        # bwa writes 5 sidecars; .bwt is the cheapest to use as a sentinel.
        sentinel = BWA_INDEX_SENTINEL,
    shell:
        "bwa index {input}"

rule bwa_mem2_index:
    input:
        REF_FASTA
    output:
        sentinel = BWA_MEM2_INDEX_SENTINEL,
    shell:
        "bwa-mem2 index {input}"

def _aligner_index_input(_wildcards):
    """Choose which index sidecar gates bwa_map based on config.aligner."""
    return BWA_INDEX_SENTINEL if ALIGNER == "bwa" else BWA_MEM2_INDEX_SENTINEL

# ---- Map + sort -------------------------------------------------------------

rule bwa_map:
    input:
        fasta = REF_FASTA,
        fai   = REF_FAI,
        index = _aligner_index_input,
        r1    = f"{SOURCE_DIR}/{{sample}}_1.fastq.gz",
        r2    = f"{SOURCE_DIR}/{{sample}}_2.fastq.gz",
    output:
        bam = temp("output/mapped_reads/{sample}.bam"),
    threads: config.get("threads", {}).get("bwa", 5)
    params:
        aligner = ALIGNER,
        # @RG written at map time. ID + PL only; SM/LB are normalised by the
        # reheader rule downstream (matches the predecessor's flow).
        rg = lambda w: rf"@RG\tID:{w.sample}\tPL:ILLUMINA",
    shell:
        # `-M` marks split alignments as secondary (Picard-compatible).
        "{params.aligner} mem -t {threads} -M -R '{params.rg}' "
        "{input.fasta} {input.r1} {input.r2} "
        "| samtools sort -@ {threads} -o {output.bam} -"

# ---- MarkDuplicates ---------------------------------------------------------

rule mark_duplicates:
    input:
        bam = "output/mapped_reads/{sample}.bam",
    output:
        bam     = temp("output/bam/{sample}_dupmarked.bam"),
        metrics = temp("output/bam/{sample}_picard_metrics.txt"),
    shell:
        # AS=TRUE: assume coord-sorted input. VALIDATION_STRINGENCY=LENIENT
        # matches the predecessor.
        "picard MarkDuplicates AS=TRUE VALIDATION_STRINGENCY=LENIENT "
        "I={input.bam} O={output.bam} M={output.metrics}"

# ---- Reheader RG + index ----------------------------------------------------
# Overwrites whatever @RG the aligner wrote with a clean ID/SM/LB/PL block.
# Matches the predecessor exactly (sed in-place on the header text).

rule reheader_and_index:
    input:
        bam = "output/bam/{sample}_dupmarked.bam",
    output:
        bam = "output/bam/{sample}.bam",
        bai = "output/bam/{sample}.bam.bai",
    params:
        sed_expr = lambda w: (
            rf"'s,^@RG.*,@RG\tID:{w.sample}\tSM:{w.sample}\tLB:None\tPL:Illumina,g'"
        ),
    shell:
        "samtools view -H {input.bam} | sed {params.sed_expr} "
        "| samtools reheader - {input.bam} > {output.bam} && "
        "samtools index {output.bam}"
