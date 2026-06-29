#!/usr/bin/env python3
"""Generate a tiny synthetic Plasmodium-like test dataset for plasmo-call.

Produces:
  test/reference/synth.fasta        (2 chromosomes, ~5 kb each)
  test/reference/synth.bed          (3-col BED matching the FASTA)
  test/fastq/<sample>_1.fastq.gz    paired-end reads, gzipped
  test/fastq/<sample>_2.fastq.gz

Stdlib-only, deterministic (seeded). Reads are exact substrings of the
reference with a small random per-base error rate, so HaplotypeCaller has
a few candidate variants to chew on without the test taking forever.

Run:
    python test/generate-fixtures.py
"""
from __future__ import annotations

import gzip
import random
from pathlib import Path

SEED         = 20260629
CHROM_SIZES  = {"chr1": 5_000, "chr2": 4_000}    # bp
N_SAMPLES    = 2
COVERAGE     = 15                                # × per sample
READ_LEN     = 100                               # bp
INSERT_MEAN  = 300                               # bp (fragment length)
INSERT_SD    = 30                                # bp
ERROR_RATE   = 0.002                             # per-base substitution prob

BASES = "ACGT"
QUAL  = "I" * READ_LEN                            # ASCII 73 = Q40 (Phred+33)

ROOT      = Path(__file__).resolve().parent
REF_DIR   = ROOT / "reference"
FASTQ_DIR = ROOT / "fastq"
FASTA     = REF_DIR / "synth.fasta"
BED       = REF_DIR / "synth.bed"


def _rand_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choices(BASES, k=length))


def _mutate(rng: random.Random, seq: str) -> str:
    if ERROR_RATE <= 0:
        return seq
    out = []
    for base in seq:
        if rng.random() < ERROR_RATE:
            alt = rng.choice([b for b in BASES if b != base])
            out.append(alt)
        else:
            out.append(base)
    return "".join(out)


def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def _write_fasta(path: Path, chroms: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for name, seq in chroms.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def _write_bed(path: Path, chroms: dict[str, str]) -> None:
    with path.open("w") as fh:
        for name, seq in chroms.items():
            fh.write(f"{name}\t0\t{len(seq)}\n")


def _simulate_pairs(
    rng: random.Random,
    chroms: dict[str, str],
    sample: str,
):
    """Yield FASTQ records (id, seq, qual) for R1 then R2 alternately."""
    total_bp = sum(len(s) for s in chroms.values())
    n_pairs = max(1, (total_bp * COVERAGE) // (2 * READ_LEN))

    chrom_names = list(chroms)
    weights = [len(chroms[c]) for c in chrom_names]
    pair_idx = 0

    for _ in range(n_pairs):
        chrom = rng.choices(chrom_names, weights=weights, k=1)[0]
        seq = chroms[chrom]
        insert = max(
            READ_LEN * 2 + 10,
            int(rng.gauss(INSERT_MEAN, INSERT_SD)),
        )
        if insert > len(seq):
            continue
        start = rng.randint(0, len(seq) - insert)
        fragment = seq[start : start + insert]

        r1 = _mutate(rng, fragment[:READ_LEN])
        r2 = _mutate(rng, _revcomp(fragment[-READ_LEN:]))

        pair_idx += 1
        rid = f"{sample}:{chrom}:{start}:{pair_idx}"
        yield (f"{rid}/1", r1, QUAL), (f"{rid}/2", r2, QUAL)


def _write_fastq_gz(path: Path, records) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", compresslevel=4) as fh:
        for rid, seq, qual in records:
            fh.write(f"@{rid}\n{seq}\n+\n{qual}\n")


def main() -> None:
    rng = random.Random(SEED)

    # Build reference.
    chroms = {name: _rand_seq(rng, size) for name, size in CHROM_SIZES.items()}
    _write_fasta(FASTA, chroms)
    _write_bed(BED, chroms)
    print(f"reference  : {FASTA} ({sum(len(s) for s in chroms.values())} bp)")
    print(f"bed        : {BED}")

    # Simulate paired-end reads per sample.
    for i in range(1, N_SAMPLES + 1):
        sample = f"sample{i:02d}"
        rng_s = random.Random(SEED + i)
        r1_records, r2_records = [], []
        for r1, r2 in _simulate_pairs(rng_s, chroms, sample):
            r1_records.append(r1)
            r2_records.append(r2)
        r1_path = FASTQ_DIR / f"{sample}_1.fastq.gz"
        r2_path = FASTQ_DIR / f"{sample}_2.fastq.gz"
        _write_fastq_gz(r1_path, r1_records)
        _write_fastq_gz(r2_path, r2_records)
        print(f"sample {sample} : {len(r1_records)} pairs -> {r1_path.name} / {r2_path.name}")


if __name__ == "__main__":
    main()
