## Common rules (applies to every BT-lane)
- Work folder. Write ONLY under results/<LANE>/. Read-only everywhere else (other results folders, references/, external_mount, ~/projects/...). No git commit/push, no emails, no publishing, no downloading unless it is stated in the assignment.
- Swedish in README/RESULTS. Each digit in results/<LANE>/*.json.
- PREREG.md is written BEFORE the first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged last in PREREG and in RESULTS.
- Cargo: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, at most 2 processes. Before heavy driving: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, maximum 15 min). Top ≥ 4 GB RAM or GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Intermediate results > 50 MB → external_media<LANE>/. The system disk (/) is almost full.
- Import ALDRIG other agents' modules with `from X import *` and never write to their folders: set own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill other processes, never `nvidia-smi -q`.
- Report honestly: cases are as valuable as wins. No claims without figure and source file.
- End with results/<LANE>/RESULTS.md which starts with the line `# <LANE>` and contains: what was done, outcome per criterion (pass/fail), deviations, files.
