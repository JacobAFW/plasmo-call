# =============================================================================
# workflow/rules/common.smk
# Config load + validation, sample discovery, chromosome/interval wrangling.
# Imported first by the Snakefile; every other rule file relies on what's
# defined here.
# =============================================================================
import math
import os
import sys
import yaml
from pathlib import Path

# ---- Load params.yaml + the selected species file ---------------------------
# config (config.yaml) is already populated by Snakemake's `configfile:` line.
# params.yaml and config/species/<species>.yaml are loaded here and merged onto
# `config` under nested keys so the rest of the workflow has a single namespace.

_repo_root = Path(workflow.basedir).parent     # .../plasmo-call
_cfg_dir   = _repo_root / "config"

with open(_cfg_dir / "params.yaml") as _fh:
    config["params"] = yaml.safe_load(_fh)

_species = config.get("species")
if not _species:
    sys.exit("ERROR: config.species is not set. Pick one of: vivax | knowlesi | malariae | ovale")
_species_file = _cfg_dir / "species" / f"{_species}.yaml"
if not _species_file.exists():
    sys.exit(f"ERROR: species file not found: {_species_file}")
with open(_species_file) as _fh:
    config["species_priors"] = yaml.safe_load(_fh)

# ---- Validate species priors (fail fast on TBD stubs) -----------------------
# GATK4 HaplotypeCaller needs concrete numeric values for both priors. The
# stub species files (knowlesi/malariae/ovale) ship with `heterozygosity:`
# left blank (yaml-null), so we bail before any rule fires.
_het  = config["species_priors"].get("heterozygosity")
_ihet = config["species_priors"].get("indel_heterozygosity")
if _het is None or _ihet is None:
    sys.exit(
        f"ERROR: species '{_species}' has unset priors in {_species_file}\n"
        f"       heterozygosity={_het!r}, indel_heterozygosity={_ihet!r}\n"
        f"       Set both before running. See config/species/vivax.yaml for the\n"
        f"       only documented preset; do NOT copy vivax values to other species."
    )

# ---- Sample discovery -------------------------------------------------------
SOURCE_DIR = config.get("source_dir", "")
if not SOURCE_DIR:
    sys.exit("ERROR: config.source_dir is empty. Set it to the directory of paired FASTQs.")

SAMPLES, = glob_wildcards(f"{SOURCE_DIR}/{{sample}}_1.fastq.gz")
if not SAMPLES:
    sys.exit(f"ERROR: no samples found in {SOURCE_DIR} matching {{sample}}_1.fastq.gz")

# ---- Reference paths --------------------------------------------------------
REF_FASTA = config["reference"]["fasta"]
REF_BED   = config["reference"]["bed"]
if not REF_FASTA or not REF_BED:
    sys.exit("ERROR: config.reference.fasta and config.reference.bed must both be set.")
REF_FAI   = f"{REF_FASTA}.fai"
REF_DICT  = f"{os.path.splitext(REF_FASTA)[0]}.dict"
# bwa index sidecars share the FASTA prefix; bwa-mem2 adds .0123 / .bwt.2bit.64
BWA_INDEX_SENTINEL      = f"{REF_FASTA}.bwt"
BWA_MEM2_INDEX_SENTINEL = f"{REF_FASTA}.bwt.2bit.64"

ALIGNER = config.get("aligner", "bwa")
if ALIGNER not in ("bwa", "bwa-mem2"):
    sys.exit(f"ERROR: config.aligner must be 'bwa' or 'bwa-mem2' (got {ALIGNER!r})")

# ---- Chromosome / interval wrangling ---------------------------------------
# Ported from the GATK3 predecessor (pandas wide_to_long). Same semantics,
# stdlib-only:
#   * read 3-column BED (chrom, start, end)
#   * chromosomes with end <= SPLIT_THRESHOLD stay as a single interval
#   * chromosomes with end >  SPLIT_THRESHOLD are split into N_SEGMENTS equal
#     pieces; the last piece absorbs the remainder so the final coord is `end`.
#   * intervals are emitted in GATK 1-based inclusive form ("chrom:start-end").
SPLIT_THRESHOLD = 100_000   # bp; matches the original threshold
N_SEGMENTS      = 10        # matches the original 10-way split

def _load_intervals(bed_path: str):
    chromosomes: list[str] = []
    intervals:   list[str] = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                sys.exit(f"ERROR: malformed BED line in {bed_path}: {line!r}")
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            chromosomes.append(chrom)
            if end > SPLIT_THRESHOLD:
                # Predecessor used floor(end/10*i) for cut points; replicate.
                step = math.floor(end / N_SEGMENTS)
                cuts = [math.floor(end / N_SEGMENTS * i) for i in range(1, N_SEGMENTS + 1)]
                seg_start = start + 1                        # 1-based inclusive
                for i, cut in enumerate(cuts):
                    seg_end = end if i == N_SEGMENTS - 1 else cut
                    intervals.append(f"{chrom}:{seg_start}-{seg_end}")
                    seg_start = seg_end + 1
            else:
                intervals.append(f"{chrom}:{start + 1}-{end}")
    return chromosomes, intervals

CHROMOSOME, CHROMOSOME_INTERVALS = _load_intervals(REF_BED)
if not CHROMOSOME:
    sys.exit(f"ERROR: no chromosomes parsed from {REF_BED}")

# ---- BQSR mode resolution + LOUD logging -----------------------------------
# bqsr.mode: off | known_sites | bootstrap | auto
#   auto  -> bootstrap if known_variants is empty, else known_sites
# The chosen path MUST be visible; bootstrap adds a full extra calling pass,
# so we never want it to fire silently.
_VALID_BQSR_MODES = ("off", "known_sites", "bootstrap", "auto")

def _resolve_bqsr_mode() -> str:
    bqsr_cfg = config.get("bqsr", {}) or {}
    requested = bqsr_cfg.get("mode", "auto")
    if requested not in _VALID_BQSR_MODES:
        sys.exit(
            f"ERROR: config.bqsr.mode must be one of {list(_VALID_BQSR_MODES)}; "
            f"got {requested!r}"
        )

    known_variants = bqsr_cfg.get("known_variants") or []

    if requested == "auto":
        resolved = "known_sites" if known_variants else "bootstrap"
        origin = f"auto (resolved → {resolved})"
    else:
        resolved = requested
        origin = requested

    if resolved == "known_sites" and not known_variants:
        sys.exit(
            "ERROR: bqsr.mode resolved to 'known_sites' but bqsr.known_variants is empty.\n"
            "       Either point known_variants at one or more VCF(s), or set bqsr.mode\n"
            "       to 'bootstrap' / 'auto' / 'off'."
        )

    if resolved == "bootstrap":
        iterations = (bqsr_cfg.get("bootstrap") or {}).get("iterations", 1)
        if iterations != 1:
            sys.exit(
                f"ERROR: bqsr.bootstrap.iterations={iterations} is not yet implemented.\n"
                f"       Only single-pass bootstrap is wired today. Set iterations: 1."
            )

    # ---- LOUD banner (stderr so it shows even when stdout is captured) ------
    bar = "=" * 72
    lines = [bar, f"[plasmo-call] BQSR mode = {resolved}   (from {origin})"]
    if resolved == "off":
        lines.append("              skipping recalibration; HC reads output/bam/<sample>.bam")
    elif resolved == "known_sites":
        lines.append(f"              using {len(known_variants)} known-sites VCF(s):")
        for v in known_variants:
            lines.append(f"                - {v}")
    elif resolved == "bootstrap":
        lines.append("              no external known-sites VCF; will call+hard-filter on")
        lines.append("              non-recal BAMs to derive one. THIS ADDS A FULL EXTRA")
        lines.append("              CALLING PASS BEFORE THE MAIN HC RUN.")
    lines.append(bar)
    print("\n" + "\n".join(lines) + "\n", file=sys.stderr)
    return resolved

BQSR_MODE = _resolve_bqsr_mode()

# Helper used by both BaseRecalibrator (in bqsr.smk) and any rule that needs
# to know which VCFs back the recal table. Returns absolute or relative paths
# straight from the user's config in known_sites mode; in bootstrap mode the
# paths are the filtered VCFs the bootstrap pass will produce.
BOOTSTRAP_KNOWN_SITES = [
    "output/bqsr/bootstrap/snvs_pass.vcf.gz",
    "output/bqsr/bootstrap/indels_pass.vcf.gz",
]

def known_sites_vcfs() -> list[str]:
    if BQSR_MODE == "known_sites":
        return list(config["bqsr"]["known_variants"])
    if BQSR_MODE == "bootstrap":
        return list(BOOTSTRAP_KNOWN_SITES)
    return []                                # off — no recalibration

def known_sites_tbis() -> list[str]:
    return [f"{v}.tbi" for v in known_sites_vcfs()]

# ---- Shared HaplotypeCaller flag string -------------------------------------
# Used by both the main HC rule (calling_gatk.smk) and the bootstrap HC rule
# (bqsr.smk). Lives here because bqsr.smk loads before calling_gatk.smk.

def _hc_params_str() -> str:
    p  = config["params"]["haplotypecaller"]
    sp = config["species_priors"]
    parts = [
        f"--emit-ref-confidence {p['emit_ref_confidence']}",
        *[f"--kmer-size {k}" for k in p["kmer_size"]],
        "--dont-use-soft-clipped-bases" if p.get("dont_use_soft_clipped_bases") else "",
        f"--min-assembly-region-size {p['min_assembly_region_size']}",
        "--do-not-run-physical-phasing" if p.get("do_not_run_physical_phasing") else "",
        f"--base-quality-score-threshold {p['base_quality_score_threshold']}",
        f"-mbq {p['min_base_quality_score']}",
        *[f"-DF {f}" for f in p.get("disable_read_filter", [])],
        f"--heterozygosity {sp['heterozygosity']}",
        f"--indel-heterozygosity {sp['indel_heterozygosity']}",
    ]
    return " ".join(x for x in parts if x)

HC_PARAMS = _hc_params_str()
