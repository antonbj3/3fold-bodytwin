## Common rules (applies to every BT-lane)
- Arbetsmapp .. Skriv ENDAST under results/<LANE>/. Read-only everywhere else (andra results-mappar, references/, /mnt/..., ~/projects/...). Ingen git commit/push, no e-mails, no publishing, no download without it being stated in the assignment.
- Swedish in README/RESULTS. Each digit in results/<LANE>/*.json.
- PREREG.md is written BEFORE the first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged last in PREREG and in RESULTS.
- Last: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, not more than 2 processes. Before heavy driving: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, not more than 15 min). Topp ≥ 4 GB RAM eller GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Mellanresultat > 50 MB → /media/anton/sdc1-tmp/bodytwin/<LANE>/. Systemdisken (/) I'm almost drunk.
- Import ALDRIG other agents' modules with `from X import *` and never write to their folders: set own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill other processes, never `nvidia-smi -q`.
- Report honestly: cases are as valuable as wins. No claims without figure and source file.
- End with results/<LANE>/RESULTS.md which starts with the line `# <LANE>` and contains: what was done, outcome per criterion (pass/fail), deviations, files.
