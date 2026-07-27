# =============================================================================
# workflow/rules/calling_gatk.smk
# GATK4 arm: HaplotypeCaller per-sample → consolidation → GenotypeGVCFs per
# INTERVAL → concat per chromosome. Consolidation backend is selected by
# config.joint_calling.consolidation (see common.smk / JOINT_CALLING):
#
#   * genomicsdb    — GenomicsDBImport per interval, GenotypeGVCFs on gendb://
#   * combine_gvcfs — CombineGVCFs whole-genome, GenotypeGVCFs -L per interval
#
# Both backends produce identical per-interval outputs
# (`output/calling/gatk/joint/intervals/gatk_genotyped_intervals_{interval}.vcf.gz`),
# which are then concat'd to per-chrom VCFs consumed unchanged by consensus.smk.
#
# Shared HC params + per-species priors come from common.smk (HC_PARAMS).
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

localrules: sample_name_map

# ---- Helpers ---------------------------------------------------------------

# Mode-aware BAM path lives in common.smk (sample_bam_path); the wrappers below
# just adapt it for Snakemake's wildcards-callable-input contract, and share
# the choice with calling_bcftools.smk so the two arms can never diverge.
def _hc_input_bam(wildcards):
    return sample_bam_path(wildcards.sample)

def _hc_input_bai(wildcards):
    return f"{sample_bam_path(wildcards.sample)}.bai"

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

# ---- Sample-name-map TSV (GenomicsDBImport input at scale) ------------------
# GenomicsDBImport accepts --sample-name-map <TSV> instead of a long chain of
# `-V`s. At ~980 samples the CLI would otherwise blow past ARG_MAX; the map
# also lets us swap GVCF paths without touching the rule.

rule sample_name_map:
    input:
        gvcfs = expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz", sample=SAMPLES),
        tbis  = expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz.tbi", sample=SAMPLES),
    output:
        tsv = "output/calling/gatk/genomicsdb/sample_name_map.tsv",
    run:
        from pathlib import Path
        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output.tsv, "w") as fh:
            for s, g in zip(SAMPLES, input.gvcfs):
                fh.write(f"{s}\t{g}\n")

# ---- CombineGVCFs (legacy backend — whole genome → one combined GVCF) -------
# Only fires when JOINT_CALLING.consolidation == "combine_gvcfs". Semantics
# unchanged from the previous milestone; kept for small cohorts + as the
# behaviour-preserving baseline for validating the genomicsdb swap.

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

# ---- GenomicsDBImport (scale backend — per interval workspace) -------------
# GenomicsDBImport quirks handled here:
#   * workspace MUST NOT pre-exist — Snakemake's normal directory creation
#     would collide; we `rm -rf` first so re-runs (and interrupted runs) work.
#   * one contig per workspace — our intervals are single-contig by construction
#     (see _load_intervals in common.smk), so this is naturally satisfied.
#   * sample-name-map TSV instead of --variant/-V — see rule sample_name_map.
# The workspace directory is the tracked output (`directory(...)`). Colons in
# `{interval}` (e.g. "LT727648:1-89247") are legal on APFS/ext4 and Snakemake
# passes them through unchanged — same pattern the bcftools arm already uses.

rule genomicsdb_import:
    input:
        gvcfs           = expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz", sample=SAMPLES),
        tbis            = expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz.tbi", sample=SAMPLES),
        sample_name_map = "output/calling/gatk/genomicsdb/sample_name_map.tsv",
        fasta           = REF_FASTA,
        fai             = REF_FAI,
        dict_           = REF_DICT,
    output:
        workspace = directory("output/calling/gatk/genomicsdb/workspaces/{interval}"),
    params:
        batch_size     = JOINT_CALLING["batch_size"],
        reader_threads = JOINT_CALLING["reader_threads"],
        # tmp_dir: shell-expanded at runtime, so "$PBS_JOBFS" / "$TMPDIR"
        # resolve inside the scheduler job. See config.joint_calling.genomicsdb.tmp_dir.
        tmp_dir        = JOINT_CALLING["tmp_dir"],
    shell:
        # GenomicsDBImport refuses to write into an existing workspace; wipe
        # then let it create fresh so re-runs on the same interval succeed.
        "rm -rf {output.workspace} && "
        "gatk GenomicsDBImport "
        "--genomicsdb-workspace-path {output.workspace} "
        "--sample-name-map {input.sample_name_map} "
        "--batch-size {params.batch_size} "
        "--reader-threads {params.reader_threads} "
        "-L {wildcards.interval} "
        "--tmp-dir {params.tmp_dir}"

# ---- GenotypeGVCFs (per interval, backend-aware input) ---------------------
# One rule, two input shapes:
#   * genomicsdb:    -V gendb://<workspace>       (workspace is interval-scoped)
#   * combine_gvcfs: -V <combined.g.vcf.gz> -L    (subset the combined GVCF)
# Selection is done at DAG-build time via a callable input (Snakemake evaluates
# the callable once per (rule, wildcards) — the wrong backend's inputs are
# never requested, so no phantom jobs get scheduled).

def _genotype_gvcfs_input(wildcards):
    common = {
        "fasta": REF_FASTA,
        "fai":   REF_FAI,
        "dict_": REF_DICT,
    }
    if JOINT_CALLING["consolidation"] == "genomicsdb":
        return {
            **common,
            "workspace": f"output/calling/gatk/genomicsdb/workspaces/{wildcards.interval}",
        }
    return {
        **common,
        "gvcf": "output/calling/gatk/gvcf/GATK_combined.g.vcf.gz",
        "tbi":  "output/calling/gatk/gvcf/GATK_combined.g.vcf.gz.tbi",
    }

def _genotype_gvcfs_variant_arg(wildcards):
    if JOINT_CALLING["consolidation"] == "genomicsdb":
        return f"gendb://output/calling/gatk/genomicsdb/workspaces/{wildcards.interval}"
    return "output/calling/gatk/gvcf/GATK_combined.g.vcf.gz"

rule genotype_gvcfs:
    input:
        unpack(_genotype_gvcfs_input),
    output:
        vcf = temp("output/calling/gatk/joint/intervals/gatk_genotyped_intervals_{interval}.vcf.gz"),
        tbi = temp("output/calling/gatk/joint/intervals/gatk_genotyped_intervals_{interval}.vcf.gz.tbi"),
    params:
        variant_arg = _genotype_gvcfs_variant_arg,
    shell:
        # -L is redundant for the gendb backend (the workspace is already
        # interval-scoped) but harmless — GATK filters to the specified
        # interval either way, so keeping -L for BOTH backends makes the two
        # code paths symmetric and easier to reason about.
        "gatk GenotypeGVCFs "
        "-R {input.fasta} "
        "-V {params.variant_arg} "
        "-L {wildcards.interval} "
        "-O {output.vcf}"

# ---- Concat a chromosome's interval VCFs → per-chrom joint VCF -------------
# Same shape as calling_bcftools.smk's concat_bcftools — the two arms stay
# aligned on the per-chrom output, which is what consensus.smk consumes.

def _gatk_intervals_for_chrom(wildcards):
    prefix = f"{wildcards.chromosome}:"
    return [f"output/calling/gatk/joint/intervals/gatk_genotyped_intervals_{iv}.vcf.gz"
            for iv in CHROMOSOME_INTERVALS if iv.startswith(prefix)]

def _gatk_interval_tbis_for_chrom(wildcards):
    return [f"{v}.tbi" for v in _gatk_intervals_for_chrom(wildcards)]

rule concat_gatk_chrom:
    input:
        vcfs = _gatk_intervals_for_chrom,
        tbis = _gatk_interval_tbis_for_chrom,
    output:
        vcf = "output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz",
        tbi = "output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz.tbi",
    shell:
        "bcftools concat -Oz -o {output.vcf} {input.vcfs} && "
        "bcftools index -t -o {output.tbi} {output.vcf}"
