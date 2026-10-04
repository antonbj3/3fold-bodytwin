## Shared rules (apply to every BT lane)
- Working folder . Write ONLY under results/<LANE>/. Read-only everywhere else (other results folders, references/, external_mount, ~/projects/...). No git commit/push, no email, no publishing, no downloading unless stated in the task.
- Swedish in README/RESULTS. Every number in results/<LANE>/*.json.
- PREREG.md is written BEFORE the first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged at the end of PREREG and in RESULTS.
- Load: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, at most 2 processes. Before a heavy run: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, at most 15 min). Peak ≥ 4 GB RAM or GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Intermediate results > 50 MB → external_media<LANE>/. The system disk (/) is almost full.
- NEVER import other agents' modules with `from X import *` and never write to their folders: set your own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill others' processes, never `nvidia-smi -q`.
- Report honestly: failures are as valuable as gains. No claims without a number and source file.
- End with results/<LANE>/RESULTS.md starting with the line `# <LANE>` and containing: what was done, outcome per criterion (holds/failed), deviations, files.
