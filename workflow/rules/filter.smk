# =============================================================================
# workflow/rules/filter.smk
# Post-calling filter tiers + funnel stats report (Prompt H).
#
# DECOUPLED from calling: all rules here consume ONLY the finished
# `output/calling/consensus/Consensus.vcf.gz` (plus its .tbi). No BAMs, no
# per-chrom VCFs, no calling inputs. Amending a threshold in
# config.filter and re-running regenerates the tiers + funnel report in
# seconds WITHOUT re-calling.
#
# THREE VCFs are ALWAYS emitted side by side (no enabled/on-off toggle —
# that would hide the reproducible base each tier was built from):
#   * output/calling/consensus/Consensus.vcf.gz        — RAW, unchanged.
#   * output/calling/consensus/Consensus.light.vcf.gz  — LIGHT tier.
#   * output/calling/consensus/Consensus.hard.vcf.gz   — HARD tier.
# Plus the funnel stats report:
#   * output/filter/stats.txt   — RAW → LIGHT → HARD counts by SNP/indel,
#                                  per-site QUAL/QD/FS/MQ/MQRankSum/
#                                  ReadPosRankSum/SOR distributions with
#                                  tier thresholds annotated.
#   * output/filter/stats.tsv   — per-site machine-readable table.
#
# Tier recipes (thresholds in config.filter; defaults in common.smk):
#   * LIGHT = QUAL >= light.min_qual (bcftools view -e "QUAL<T"). vivaxgen-
#            equivalent floor. NO per-site DP cut — per-sample depth is
#            tier-2 downstream, out of scope here.
#   * HARD  = Sanger-style GATK hard filter applied to SNPs AND indels
#            alike, via GATK VariantFiltration then `bcftools view -f PASS,.`:
#              QD < hard.qd  → filter "QD<qd"
#              FS > hard.fs  → filter "FS>fs"
#              MQ < hard.mq  → filter "MQ<mq"
#            The three expressions run independently — a record only fails a
#            given expression if that expression evaluates true; a record
#            with QD missing is left unfiltered by that expression (GATK
#            leaves records with undefined JEXL vars UNTAGGED).
#
# Why HARD is QD/FS/MQ ONLY (see docs/filter.md for the long form):
#   MQRankSum, ReadPosRankSum, SOR are diploid-tuned rank-sum / strand-bias
#   tests defined only at HET sites. Plasmodium is HAPLOID (called at
#   diploid ploidy to capture mixed/polyclonal infections, quantified
#   downstream by Fws / COI); those thresholds don't transfer. Per GATK's
#   non-model-organism + Sanger/MalariaGEN malaria guidance they are SHOWN
#   in the funnel report but NOT filtered.
# =============================================================================

localrules: consensus_stats

# ---- LIGHT tier -------------------------------------------------------------
# QUAL floor only. bcftools view -e drops records failing the expression.
# We use -e (exclude) rather than -i (include) so a site with QUAL="."
# (very rare on the joint VCF, but possible on odd sites) is KEPT — the
# safer default. Users who want strict-numeric semantics can override to -i.

rule filter_light:
    input:
        vcf = "output/calling/consensus/Consensus.vcf.gz",
        tbi = "output/calling/consensus/Consensus.vcf.gz.tbi",
    output:
        vcf = "output/calling/consensus/Consensus.light.vcf.gz",
        tbi = "output/calling/consensus/Consensus.light.vcf.gz.tbi",
    params:
        min_qual = FILTER["light"]["min_qual"],
    shell:
        "bcftools view -e 'QUAL<{params.min_qual}' -Oz -o {output.vcf} {input.vcf} && "
        "bcftools index -t -o {output.tbi} {output.vcf}"

# ---- HARD tier --------------------------------------------------------------
# Sanger-style GATK hard filter: QD<qd | FS>fs | MQ<mq. Three independent
# --filter-expression flags so each contributes its own FILTER name and
# missing-annotation semantics (undefined JEXL var → record NOT filtered by
# that expression, which is what we want — conservative).
#
# We keep the tagged intermediate as temp() so someone can inspect FILTER
# tags on borderline sites during a review, without leaving it lying around.
# The final .hard.vcf.gz is the PASS-only cohort.

rule filter_hard:
    input:
        vcf   = "output/calling/consensus/Consensus.vcf.gz",
        tbi   = "output/calling/consensus/Consensus.vcf.gz.tbi",
        fasta = REF_FASTA,
        fai   = REF_FAI,
        dict_ = REF_DICT,
    output:
        tagged     = temp("output/filter/Consensus.hard.tagged.vcf.gz"),
        tagged_tbi = temp("output/filter/Consensus.hard.tagged.vcf.gz.tbi"),
        vcf = "output/calling/consensus/Consensus.hard.vcf.gz",
        tbi = "output/calling/consensus/Consensus.hard.vcf.gz.tbi",
    params:
        qd = FILTER["hard"]["qd"],
        fs = FILTER["hard"]["fs"],
        mq = FILTER["hard"]["mq"],
    shell:
        "gatk VariantFiltration "
        "-R {input.fasta} "
        "-V {input.vcf} "
        "--filter-expression 'QD < {params.qd}' --filter-name 'QD{params.qd}' "
        "--filter-expression 'FS > {params.fs}' --filter-name 'FS{params.fs}' "
        "--filter-expression 'MQ < {params.mq}' --filter-name 'MQ{params.mq}' "
        "-O {output.tagged} && "
        "bcftools view -f 'PASS,.' -Oz -o {output.vcf} {output.tagged} && "
        "bcftools index -t -o {output.tbi} {output.vcf}"

# ---- Funnel stats report ----------------------------------------------------
# Consumes ALL THREE VCFs so it can print the RAW → LIGHT → HARD funnel
# with SNP/indel splits, and pull QUAL/QD/FS/MQ/MQRankSum/ReadPosRankSum/SOR
# distributions from the raw callset (with LIGHT + HARD thresholds annotated
# so the user can see where each cut lands on the real distribution).

rule consensus_stats:
    input:
        raw       = "output/calling/consensus/Consensus.vcf.gz",
        raw_tbi   = "output/calling/consensus/Consensus.vcf.gz.tbi",
        light     = "output/calling/consensus/Consensus.light.vcf.gz",
        light_tbi = "output/calling/consensus/Consensus.light.vcf.gz.tbi",
        hard      = "output/calling/consensus/Consensus.hard.vcf.gz",
        hard_tbi  = "output/calling/consensus/Consensus.hard.vcf.gz.tbi",
    output:
        txt = "output/filter/stats.txt",
        tsv = "output/filter/stats.tsv",
    params:
        light_min_qual = FILTER["light"]["min_qual"],
        hard_qd        = FILTER["hard"]["qd"],
        hard_fs        = FILTER["hard"]["fs"],
        hard_mq        = FILTER["hard"]["mq"],
    shell:
        "python scripts/consensus_stats.py "
        "--raw {input.raw} --light {input.light} --hard {input.hard} "
        "--txt {output.txt} --tsv {output.tsv} "
        "--light-min-qual {params.light_min_qual} "
        "--hard-qd {params.hard_qd} "
        "--hard-fs {params.hard_fs} "
        "--hard-mq {params.hard_mq}"
