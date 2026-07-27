#!/usr/bin/env python3
# =============================================================================
# scripts/consensus_stats.py
# Funnel stats report for the raw + LIGHT + HARD Consensus tiers (Prompt H).
#
# What this produces:
#   * stats.txt — human-readable
#       - RAW → LIGHT → HARD funnel: total + SNP + indel counts at each
#         tier, absolute drop and percent-remaining of the previous tier
#         (that's the "funnel").
#       - Per-site distributions of QUAL, QD, FS, MQ, MQRankSum,
#         ReadPosRankSum, SOR — pulled from the RAW consensus. The LIGHT
#         and HARD thresholds are printed next to each distribution so
#         the reader can see where each cut lands on the real numbers.
#         MQRankSum / ReadPosRankSum / SOR are SHOWN but NOT filtered on
#         (see docs/filter.md).
#   * stats.tsv — per-site machine-readable table for downstream plotting.
#         Columns: CHROM POS TYPE QUAL QD FS MQ MQRankSum ReadPosRankSum
#                  SOR IN_LIGHT IN_HARD
#         where TYPE is "snp" or "indel" and IN_LIGHT / IN_HARD are 0/1
#         flags derived by joining on (CHROM, POS) against the tier VCFs.
#
# Args (argparse; the smk rule wires them):
#   --raw / --light / --hard   paths to the three Consensus VCFs
#   --txt / --tsv              output paths
#   --light-min-qual           LIGHT tier's QUAL cutoff (for annotation)
#   --hard-qd / --hard-fs / --hard-mq   HARD tier's thresholds (annotation)
#
# We call `bcftools query` for per-site fields (streaming) and
# `bcftools stats` for the summary block. Pandas is already pinned in
# pixi.toml — no new deps.
# =============================================================================
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

# ---- Per-site query columns pulled from the RAW consensus ------------------
# Kept aligned with the -f format string below.
SITE_COLS = [
    "CHROM", "POS", "QUAL",
    "TYPE",                      # "snp" | "indel" | "mnp" | "other"
    "DP", "QD", "FS", "MQ",
    "MQRankSum", "ReadPosRankSum", "SOR",
]

QUANTILE_LABELS = ["min", "p05", "p25", "med", "p75", "p95", "max"]
QUANTILES       = [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0]


def run(cmd: list[str]) -> str:
    """Run a subprocess, capture stdout, die loudly on non-zero exit."""
    res = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if res.returncode != 0:
        sys.stderr.write(f"ERROR running: {' '.join(cmd)}\n{res.stderr}\n")
        sys.exit(res.returncode)
    return res.stdout


def _to_num(s: pd.Series) -> pd.Series:
    """Coerce bcftools query output ('.' for missing) to numeric NaN."""
    return pd.to_numeric(s.replace(".", pd.NA), errors="coerce")


def _classify_type(alt: str, ref_len: int) -> str:
    """SNP if all ALTs are single nucleotides and REF is length 1; indel if
    ANY ALT differs in length from REF; else 'other'. Matches bcftools stats
    'number of SNPs'/'number of indels' bucketing to first order."""
    if alt in (".", "", None):
        return "other"
    alts = alt.split(",")
    if ref_len == 1 and all(len(a) == 1 for a in alts):
        return "snp"
    if any(len(a) != ref_len for a in alts):
        return "indel"
    return "other"


def build_raw_site_table(raw_vcf: str) -> pd.DataFrame:
    """Per-site DataFrame from the RAW consensus.

    Uses bcftools -u so undefined tags print '.' instead of erroring
    (MQRankSum / ReadPosRankSum / SOR are missing on most homozygous
    or single-allele sites).
    """
    fmt = (
        "%CHROM\t%POS\t%QUAL\t"
        "%REF\t%ALT\t"
        "%INFO/DP\t%INFO/QD\t%INFO/FS\t%INFO/MQ\t"
        "%INFO/MQRankSum\t%INFO/ReadPosRankSum\t%INFO/SOR\n"
    )
    stdout = run(["bcftools", "query", "-u", "-f", fmt, raw_vcf])

    rows = []
    for line in stdout.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        chrom, pos, qual, ref, alt = parts[:5]
        dp, qd, fs, mq, mqrs, rprs, sor = parts[5:12]
        rows.append({
            "CHROM": chrom, "POS": pos, "QUAL": qual,
            "TYPE": _classify_type(alt, len(ref)),
            "DP": dp, "QD": qd, "FS": fs, "MQ": mq,
            "MQRankSum": mqrs, "ReadPosRankSum": rprs, "SOR": sor,
        })
    df = pd.DataFrame(rows, columns=SITE_COLS)
    for c in ("QUAL", "DP", "QD", "FS", "MQ",
              "MQRankSum", "ReadPosRankSum", "SOR"):
        df[c] = _to_num(df[c])
    df["POS"] = _to_num(df["POS"]).astype("Int64")
    return df


def positions_set(vcf: str) -> set[tuple[str, int]]:
    """Set of (CHROM, POS) tuples in a VCF — used to flag IN_LIGHT / IN_HARD."""
    stdout = run(["bcftools", "query", "-f", "%CHROM\t%POS\n", vcf])
    out: set[tuple[str, int]] = set()
    for line in stdout.splitlines():
        if not line:
            continue
        chrom, pos = line.split("\t", 1)
        try:
            out.add((chrom, int(pos)))
        except ValueError:
            continue
    return out


def tier_counts(df: pd.DataFrame, in_col: str) -> dict[str, int]:
    """Return {'total','snp','indel'} for records where df[in_col] is truthy.

    Guards for empty inputs — on the synthetic flat-Q fixture (Prompt D)
    the consensus is legitimately 0 records, and empty-DF boolean masking
    can drop columns in pandas."""
    if df.empty:
        return {"total": 0, "snp": 0, "indel": 0}
    sub = df[df[in_col]] if in_col else df
    if sub.empty or "TYPE" not in sub.columns:
        return {"total": int(len(sub)), "snp": 0, "indel": 0}
    return {
        "total": int(len(sub)),
        "snp":   int((sub["TYPE"] == "snp").sum()),
        "indel": int((sub["TYPE"] == "indel").sum()),
    }


def _quantile_line(name: str, s: pd.Series, threshold_note: str = "") -> str:
    """One-line ASCII quantile summary with an optional threshold annotation."""
    s = s.dropna()
    if s.empty:
        return f"  {name:<16}  (all missing){('  ' + threshold_note) if threshold_note else ''}"
    qs = s.quantile(QUANTILES)
    parts = [f"{q:.3g}" for q in qs.values]
    body = "  ".join(f"{lbl}={val}" for lbl, val in zip(QUANTILE_LABELS, parts))
    tail = f"  {threshold_note}" if threshold_note else ""
    return f"  {name:<16}  {body}{tail}"


def write_report(raw: pd.DataFrame,
                 funnel: dict,
                 bcf_stats_txt: str,
                 raw_vcf: str,
                 light_vcf: str,
                 hard_vcf: str,
                 light_min_qual: float,
                 hard_qd: float,
                 hard_fs: float,
                 hard_mq: float,
                 out_txt: Path) -> None:
    lines: list[str] = []
    add = lines.append
    bar = "=" * 72

    add(bar)
    add("plasmo-call — consensus filter funnel report")
    add(f"raw   VCF : {raw_vcf}")
    add(f"light VCF : {light_vcf}   (QUAL >= {light_min_qual})")
    add(f"hard  VCF : {hard_vcf}   (QD >= {hard_qd}, FS <= {hard_fs}, MQ >= {hard_mq})")
    add(bar)
    add("")
    add("FUNNEL — record counts by tier, split by type")
    add(f"  {'tier':<8}{'total':>12}{'SNPs':>12}{'indels':>12}"
        f"{'%prev':>10}{'%raw':>10}")
    prev_total = None
    raw_total  = funnel["raw"]["total"]
    for tier in ("raw", "light", "hard"):
        c = funnel[tier]
        pct_prev = "" if prev_total is None else (
            "100.0%" if prev_total == 0
            else f"{c['total']/prev_total*100:.1f}%"
        )
        pct_raw  = "" if raw_total == 0 else f"{c['total']/raw_total*100:.1f}%"
        add(f"  {tier:<8}{c['total']:>12}{c['snp']:>12}{c['indel']:>12}"
            f"{pct_prev:>10}{pct_raw:>10}")
        prev_total = c["total"]
    add("")
    add(bar)
    add("PER-SITE DISTRIBUTIONS  (from RAW; quantiles: min p05 p25 med p75 p95 max)")
    add(f"  cuts annotated per row — LIGHT: QUAL>={light_min_qual};  "
        f"HARD: QD>={hard_qd}, FS<={hard_fs}, MQ>={hard_mq}.")
    add("  MQRankSum / ReadPosRankSum / SOR: shown for review; NOT filtered")
    add("  (see docs/filter.md — rank-sum tests are diploid-tuned and don't")
    add("   transfer to haploid Plasmodium called at diploid ploidy).")
    add("")
    add(_quantile_line("QUAL",           raw["QUAL"],
                       f"[LIGHT cut: QUAL<{light_min_qual}]"))
    add(_quantile_line("QD",             raw["QD"],
                       f"[HARD cut: QD<{hard_qd}]"))
    add(_quantile_line("FS",             raw["FS"],
                       f"[HARD cut: FS>{hard_fs}]"))
    add(_quantile_line("MQ",             raw["MQ"],
                       f"[HARD cut: MQ<{hard_mq}]"))
    add(_quantile_line("MQRankSum",      raw["MQRankSum"],
                       "[review only, not filtered]"))
    add(_quantile_line("ReadPosRankSum", raw["ReadPosRankSum"],
                       "[review only, not filtered]"))
    add(_quantile_line("SOR",            raw["SOR"],
                       "[review only, not filtered]"))
    add("")
    add(bar)
    add("bcftools stats (raw consensus — see 'SN' section for totals)")
    add(bar)
    add(bcf_stats_txt)
    out_txt.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Funnel stats over the RAW + LIGHT + HARD consensus tiers.")
    p.add_argument("--raw",   required=True, help="RAW consensus VCF (indexed)")
    p.add_argument("--light", required=True, help="LIGHT-tier VCF (indexed)")
    p.add_argument("--hard",  required=True, help="HARD-tier VCF (indexed)")
    p.add_argument("--txt",   required=True, help="Output human-readable report")
    p.add_argument("--tsv",   required=True, help="Output per-site TSV")
    p.add_argument("--light-min-qual", type=float, required=True,
                   help="LIGHT tier's QUAL floor (for annotation)")
    p.add_argument("--hard-qd", type=float, required=True,
                   help="HARD tier's QD threshold (for annotation)")
    p.add_argument("--hard-fs", type=float, required=True,
                   help="HARD tier's FS threshold (for annotation)")
    p.add_argument("--hard-mq", type=float, required=True,
                   help="HARD tier's MQ threshold (for annotation)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_txt = Path(args.txt)
    out_tsv = Path(args.tsv)
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Site table from the RAW consensus — all downstream tier tests join on
    # (CHROM, POS) from this frame.
    raw = build_raw_site_table(args.raw)

    # Flag each raw site with whether it survives LIGHT / HARD.
    light_positions = positions_set(args.light)
    hard_positions  = positions_set(args.hard)
    raw["IN_LIGHT"] = [
        (c, int(p)) in light_positions if pd.notna(p) else False
        for c, p in zip(raw["CHROM"], raw["POS"])
    ]
    raw["IN_HARD"] = [
        (c, int(p)) in hard_positions if pd.notna(p) else False
        for c, p in zip(raw["CHROM"], raw["POS"])
    ]

    funnel = {
        "raw":   tier_counts(raw, in_col=""),
        "light": tier_counts(raw, in_col="IN_LIGHT"),
        "hard":  tier_counts(raw, in_col="IN_HARD"),
    }

    bcf_stats_txt = run(["bcftools", "stats", args.raw])

    # Machine-readable per-site table.
    tsv_cols = [
        "CHROM", "POS", "TYPE", "QUAL", "QD", "FS", "MQ",
        "MQRankSum", "ReadPosRankSum", "SOR", "IN_LIGHT", "IN_HARD",
    ]
    raw[tsv_cols].to_csv(out_tsv, sep="\t", index=False, na_rep="NA")

    write_report(
        raw=raw,
        funnel=funnel,
        bcf_stats_txt=bcf_stats_txt,
        raw_vcf=args.raw,
        light_vcf=args.light,
        hard_vcf=args.hard,
        light_min_qual=args.light_min_qual,
        hard_qd=args.hard_qd,
        hard_fs=args.hard_fs,
        hard_mq=args.hard_mq,
        out_txt=out_txt,
    )


if __name__ == "__main__":
    main()
