## Common rules (applies to every BT-lane)
- Arbetsmapp .. Skriv ENDAST under results/<LANE>/. Read-only everywhere else (andra results-mappar, references/, /mnt/..., ~/projects/...). Ingen git commit/push, no e-mails, no publishing, no download without it being stated in the assignment.
- Swedish in README/RESULTS. Each digit in results/<LANE>/*.json.
- PREREG.md is written BEFORE the first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged last in PREREG and in RESULTS.
- Last: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, not more than 2 processes. Before heavy driving: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, not more than 15 min). Topp ≥ 4 GB RAM eller GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Mellanresultat > 50 MB → /media/anton/sdc1-tmp/bodytwin/<LANE>/. Systemdisken (/) I'm almost drunk.
- Import ALDRIG other agents' modules with `from X import *` and never write to their folders: set own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill other processes, never `nvidia-smi -q`.
- Report honestly: cases are as valuable as wins. No claims without figure and source file.
- End with results/<LANE>/RESULTS.md which starts with the line `# <LANE>` and contains: what was done, outcome per criterion (pass/fail), deviations, files.

# BT_IM2_CONTACT_CONFOUND_SPLIT_20260923

User renewed reserve_worker authorization 2026-09-23 ~16:08; breakthrough hunt across domains.

OVERRIDING RUN SCOPE: write only ./results/BT_IM2_CONTACT_CONFOUND_SPLIT_20260923/; source modules and other lanes read-only. Read applicable AGENTS and project START before work. Hash PREREG before first new numerical evaluation; exploratory changes separately documented. CPU2 (all BLAS/OMP/MKL/NUMEXPR), nice19, <=2GB unless the task explicitly authorizes one gated solve. Maximum new root-disk artifacts50MB; no bulk downloads, package/toolchain builds, GPU, cloud compute, source edits, public output or email. No child agents. No process killing. Never read credentials or billing settings. Do not wait on other lanes; snapshot completed referenced outputs and record missing dependencies. All numbers in machine-readable files. RESULTS must begin '# BT_IM2_CONTACT_CONFOUND_SPLIT_20260923' in first5lines; report negative or incomplete results honestly.

Every innovation lane must deliver CONNECT.md: source observation -> mechanism -> equation -> operator -> representation -> assumptions -> discriminating test, and one concrete cross-domain mapping with units and failure condition. Separate audited theorem, empirical evidence, hypothesis and open obligation. Audits must remain independent; do not turn a new candidate into approval of its parent result.

TASK:
Join the saved per-pair residual, force-consistency and contact-audit rows by exact subject/motion ID. Freeze a DJ1 versus STS1 contrast and quantify whether residual failure tracks suspected missing contact rather than joint-geometry update. Compare J-versus-best-baseline residual changes within each motion and inspect pelvic moment-only diagnostic as a negative control; do not rerun IK/ID. If saved rows cannot identify contact independently, report that as a structural nonidentifiability result instead of imputing chair force.
EVIDENCE TO VERIFY: IM2 joint fit missed its held-out 20% pelvic-residual gate on 0/13 and 0/6 pairs; fully held-out baseline vertical-force inconsistency was 304 N (45.7% BW), with possible unmeasured chair contact during STS.
STRONG BASELINE: Best per-pair B0/Breg/B1/B2 residual from IM2 and within-subject DJ1-versus-STS1 contrast; retain the original 20%/10% joint gate.
FALSIFIER: Contact-stratified effect vanishes, reverses, or cannot be estimated from saved observations, so contact omission cannot explain the failed geometry fit.
LIMITS: CPU1, RSS <=1 GB, <=20 MB new files; saved JSON only, no OpenSim, FE, GPU, download or patient claim.
Distinguish a finite saved-grid identification failure from structural or population nonidentifiability; report source gaps rather than proving unobserved states from missing data. Read completed upstream full report before using any quoted number.
