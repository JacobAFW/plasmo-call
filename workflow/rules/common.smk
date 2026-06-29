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
