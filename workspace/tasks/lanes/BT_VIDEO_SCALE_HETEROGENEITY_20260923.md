## Common rules (applies to every BT-lane)
- Work folder. Write ONLY under results/<LANE>/. Read-only everywhere else (other results folders, references/, external_mount, ~/projects/...). No git commit/push, no emails, no publishing, no downloading unless it is stated in the assignment.
- Swedish in README/RESULTS. Each digit in results/<LANE>/*.json.
- PREREG.md is written BEFORE the first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged last in PREREG and in RESULTS.
- Cargo: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, at most 2 processes. Before heavy driving: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, maximum 15 min). Top ≥ 4 GB RAM or GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Intermediate results > 50 MB → external_media<LANE>/. The system disk (/) is almost full.
- Import ALDRIG other agents' modules with `from X import *` and never write to their folders: set own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill other processes, never `nvidia-smi -q`.
- Report honestly: cases are as valuable as wins. No claims without figure and source file.
- End with results/<LANE>/RESULTS.md which starts with the line `# <LANE>` and contains: what was done, outcome per criterion (pass/fail), deviations, files.

# BT_VIDEO_SCALE_HETEROGENEITY_20260923

User renewed reserve_worker authorization 2026-09-23 ~16:08; breakthrough hunt across domains.

OVERRIDING RUN SCOPE: write only results/BT_VIDEO_SCALE_HETEROGENEITY_20260923/; source modules and other lanes read-only. Read applicable AGENTS and project START before work. Hash PREREG before first new numerical evaluation; exploratory changes separately documented. CPU2 (all BLAS/OMP/MKL/NUMEXPR), nice19, <=2GB unless the task explicitly authorizes one gated solve. Maximum new root-disk artifacts50MB; no bulk downloads, package/toolchain builds, GPU, cloud compute, source edits, public output or email. No child agents. No process killing. Never read credentials or billing settings. Do not wait on other lanes; snapshot completed referenced outputs and record missing dependencies. All numbers in machine-readable files. RESULTS must begin '# BT_VIDEO_SCALE_HETEROGENEITY_20260923' in first5lines; report negative or incomplete results honestly.

Every innovation lane must deliver CONNECT.md: source observation -> mechanism -> equation -> operator -> representation -> assumptions -> discriminating test, and one concrete cross-domain mapping with units and failure condition. Separate audited theorem, empirical evidence, hypothesis and open obligation. Audits must remain independent; do not turn a new candidate into approval of its parent result.

TASK:
Using only BT-V7's saved per-subject/per-trial rows, separate event-timing error from the failed gravity-ruler scale criterion (8/10) and V3 DJ-MAE criterion (19.51 %BW). Freeze a subject-held-out calibration/control comparison: one shared bias, chain-specific bias, and no correction at equal fitted parameter count. Identify whether the two failed subjects are explained by a common optical-scale mechanism or remain irreducible. Preserve the subject6 failure and compare to the 18.3 %BW gate; no pose reprocessing or new population claim.
FALSIFIER: A chain-specific calibration does not improve held-out scale and force error simultaneously over shared/no-bias controls, or improvement depends on the two held-out subjects' labels.
WHY NOW: Tests a broad BodyTwin observation failure outside the collaborator's geometry path; saved rows permit a cheap independent mechanism check rather than another video-model run.
LIMITS: CPU 2, RSS <=2 GB, new artifacts <=25 MB; no video decode, GPU, FE, download, source edit or clinical inference.
EXPLICIT COMPLETED SOURCE FILES: ["results/BT-V7/RESULTS.md"]
