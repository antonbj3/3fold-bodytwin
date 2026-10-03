# Load Rule (2026-09-23 00:30, after the machine froze and the session died twice)
DISK: BodyTwin roof on sdc1-tmp 12 GB totalt (agreed with dental 05:45; hold ≥15 GB free where). /mnt/shared_data almost drunk (7 GB, 01:05) — new dataset and large intermediate results to /media/anton/sdc1-tmp/bodytwin/ (45 GB).
Every calculation from BodyTwin agents (addition 00:50: memory counts too — swap was almost finished; anything with top ≥ 4 GB via heavy_run.sh, measure pilot first):
1. Before start: `awk '{print $1/20}' /proc/loadavg` — over 1.0: wait (sleep 60 in the foreground, max 15 min), then run a minor case.
2. `nice -n 19`, OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=2, at most 2 processes (no pools > 2).
3. ≥ 8 GB RAM eller GPU: tasks/heavy_run.sh (bigmem.lock).
4. Keep runs short and repeatable: write partial results continuously.
