# Schedulers & the test ladder

**Build and validate locally first (no PBS/SLURM), then test on each cluster.**

| Rung | Mode  | How                                   |
|------|-------|---------------------------------------|
| 1    | local | `snakemake --cores N` on test/ data   |
| 2    | PBS   | `profiles/pbs/` (known-working seed)  |
| 3    | SLURM | `profiles/slurm/` (from vvg-box)      |

See `../profiles/README.md` for why the PBS profile is overridden here rather
than inherited straight from vvg-box.
