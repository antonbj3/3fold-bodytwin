## Common rules (applies to each BT-lane)
- Work Folder . Type ONLY under results/<LANE>/. Read-only everywhere else (other results folders, references/, external_mount, ~/projects/...). No git commit/push, no emails, no publishing, no downloading unless it says so in the assignment.
- Swedish in README/RESULTS. Every digit in results/<LANE>/*.json.
- PREREG.md is written BEFORE first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged last in PREREG and in RESULTS.
- Load: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, at most 2 processes. Before heavy driving: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, maximum 15 min). Top ≥ 4 GB RAM or GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Intermediate results > 50 MB → external_media<LANE>/. The system disk (/) is almost full.
- Import ALDRIG other agents' modules with `from X import *` and never write to their folders: set own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill other processes, never `nvidia-smi -q`.
- Report honestly: cases are as valuable as wins. No claims without figure and source file.
- End with results/<LANE>/RESULTS.md starting with the line `# <LANE>` and containing: what was done, outcome per criterion (passed/failed), deviations, files.

# BT_PLASTIC_HISTORY_ORDER_20260923

User renewed swarm_worker authorization 2026-09-23 ~16:08; breakthrough hunt across domains.

OVERRIDING RUN SCOPE: write only results/BT_PLASTIC_HISTORY_ORDER_20260923/; source modules and other lanes read-only. Read applicable AGENTS and project START before work. Hash PREREG before first new numerical evaluation; exploratory changes separately documented. CPU2 (all BLAS/OMP/MKL/NUMEXPR), nice19, <=2GB unless the task explicitly authorizes one gated solve. Maximum new root-disk artifacts50MB; no bulk downloads, package/toolchain builds, GPU, cloud compute, source edits, public output or email. No child agents. No process killing. Never read credentials or billing settings. Do not wait on other lanes; snapshot completed referenced outputs and record missing dependencies. All numbers in machine-readable files. RESULTS must begin '# BT_PLASTIC_HISTORY_ORDER_20260923' in first5lines; report negative or incomplete results honestly.

Every innovation lane must deliver CONNECT.md: source observation -> mechanism -> equation -> operator -> representation -> assumptions -> discriminating test, and one concrete cross-domain mapping with units and failure condition. Separate audited theorem, empirical evidence, hypothesis and open obligation. Audits must remain independent; do not turn a new candidate into approval of its parent result.

TASK:
Use an independent scalar elastoplastic return-map oracle to compare two equal-total load sequences A→B and B→A, with perfect plasticity, isotropic hardening and a kinematic-hardening control. Freeze loads that cause at least one branch switch and a no-yield control; report stress, plastic strain, backstress, dissipated work and subsequent probe response. Test whether a representation retaining only final load/geometry is sufficient versus minimal history state. Keep X3c’s damage-order observation as a separate motivating analogue, not evidence that bone has this plastic law.
FALSIFIER: The proposed history state cannot predict the held-out probe, or apparent order dependence vanishes under an exact return map and matched final state.
WHY NOW: A cheap broad-BodyTwin constitutive experiment on path dependence and noncommutativity, distinct from the collaborator geometry and from current thin-layer/operator certificates.
LIMITS: CPU 2, RSS <=2 GB, new artifacts <=25 MB; scalar analytic model only, no biomechanical/clinical plasticity claim, FE, GPU or download.
EXPLICIT COMPLETED SOURCE FILES: ["../3fold-motion-engine/_private/romi_collab/build/U24/RESULTS.md", "results/X3c/RESULTS.md"]
