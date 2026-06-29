# Resolved tool versions

The exact versions resolved by `./install.sh` on a clean clone. Captured for
reproducibility; the authoritative source is `pixi.lock` at the repo root.

Regenerate this table with `pixi run versions`.

## Pins (pixi.toml)

| Tool                | Pin              | Channel              |
|---------------------|------------------|----------------------|
| snakemake-minimal   | `>=8.27,<9`      | conda-forge/bioconda |
| gatk4               | `4.6.2.*`        | bioconda             |
| bcftools            | `1.21.*`         | bioconda             |
| samtools            | `1.21.*`         | bioconda             |
| bwa                 | `0.7.18.*`       | bioconda             |
| bwa-mem2            | `2.2.1.*` (linux-64, osx-64 only — no osx-arm64 build in bioconda) | bioconda |
| picard              | `3.3.*`          | bioconda             |
| python              | `3.12.*`         | conda-forge          |
| pandas              | `2.2.*`          | conda-forge          |

## Resolved on `osx-arm64` (2026-06-29, dev box)

| Tool                | Resolved version                                    |
|---------------------|-----------------------------------------------------|
| snakemake           | 8.30.0                                              |
| gatk4               | 4.6.2.0 (HTSJDK 4.2.0, Picard 3.4.0 — embedded)     |
| bcftools            | 1.21                                                |
| samtools            | 1.21                                                |
| bwa                 | 0.7.18-r1243-dirty                                  |
| bwa-mem2            | (not available — osx-arm64 has no bioconda build)   |
| picard              | 3.3.0                                               |
| python              | 3.12.13                                             |
| Java (gatk4 runtime)| OpenJDK Zulu 17.0.18+8-LTS                          |

## Notes

- **GATK4 only.** GATK3.8 tools (RealignerTargetCreator, IndelRealigner) are
  intentionally absent — see `MEMORY.md`.
- **bwa-mem2 on Apple Silicon.** bioconda does not ship an `osx-arm64` build.
  The pixi manifest scopes `bwa-mem2` to `linux-64`/`osx-64` only, so
  `pixi install` succeeds on `osx-arm64` without it. Set
  `aligner: "bwa"` in `config/config.yaml` for local dev on M-series Macs;
  HPC runs (linux-64) can use either.
- **Snakemake major version.** vvg-box installs snakemake 9.x into its own
  pixi workspace (used only for the scheduler-profile autodetect). plasmo-call
  pins snakemake-minimal `>=8.27,<9` in its own workspace because the
  cluster-generic executor and the existing PBS profile are validated against
  the 8.x line. Both envs coexist; `pixi run` always uses plasmo-call's.
- **`WARN Using local manifest …`** during `pixi run` is harmless: pixi
  prefers the cwd's `pixi.toml` over `PIXI_PROJECT_MANIFEST` from vvg-box's
  activation. That preference is exactly what we want.
