## Shared rules (apply to every BT lane)
- Working folder . Write ONLY under results/<LANE>/. Read-only everywhere else (other results folders, references/, external_mount, ~/projects/...). No git commit/push, no emails, no publication, no download unless stated in the assignment.
- Swedish in README/RESULTS. Every number in results/<LANE>/*.json.
- PREREG.md is written BEFORE the first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterward; deviations are logged at the end of PREREG and in RESULTS.
- Load: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, at most 2 processes. Before a heavy run: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, at most 15 min). Peak ≥ 4 GB RAM or GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Intermediate results > 50 MB → external_media<LANE>/. The system disk (/) is almost full.
- NEVER import other agents' modules with `from X import *` and never write to their folders: set your own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill others' processes, never `nvidia-smi -q`.
- Report honestly: failures are as valuable as gains. No claims without a number and source file.
- Finish with results/<LANE>/RESULTS.md starting with the line `# <LANE>` and containing: what was done, outcome per criterion (holds/failed), deviations, files.

# BT_HISTORY_MINIMAL_STATE_AUDIT_20260923

User renewed swarm_worker authorization 2026-09-23 ~16:08; breakthrough hunt across domains.

OVERRIDING RUN SCOPE: write only results/BT_HISTORY_MINIMAL_STATE_AUDIT_20260923/; source modules and other lanes read-only. Read applicable AGENTS and project START before work. Hash PREREG before first new numerical evaluation; exploratory changes separately documented. CPU2 (all BLAS/OMP/MKL/NUMEXPR), nice19, <=2GB unless the task explicitly authorizes one gated solve. Maximum new root-disk artifacts50MB; no bulk downloads, package/toolchain builds, GPU, cloud compute, source edits, public output or email. No child agents. No process killing. Never read credentials or billing settings. Do not wait on other lanes; snapshot completed referenced outputs and record missing dependencies. All numbers in machine-readable files. RESULTS must begin '# BT_HISTORY_MINIMAL_STATE_AUDIT_20260923' in first5lines; report negative or incomplete results honestly.

Every innovation lane must deliver CONNECT.md: source observation -> mechanism -> equation -> operator -> representation -> assumptions -> discriminating test, and one concrete cross-domain mapping with units and failure condition. Separate audited theorem, empirical evidence, hypothesis and open obligation. Audits must remain independent; do not turn a new candidate into approval of its parent result.

TASK:
Independent audit of BT_PLASTIC_HISTORY_ORDER_20260923. Concrete theoretical concern: for linear Prager kinematic hardening with fixed known initial conditions, eps_p=eps-sigma/E and alpha=Hkin*eps_p+constant, so alpha may be redundant given eps/sigma. Test this independently and distinguish arbitrary internally inconsistent initialstates from states actually reachable under this model. Audit claimed minimal G/S/SP/SA representations, strain increment vs load vocabulary, samefinalstress vs samefinalstrain, reversal dissipation and plasticwork/hardeningenergysplit. Exact Fraction independent update for <=24 rational histories and independent energy integrals. Strong baseline reconstructs all algebraically implied state at same information; twohistories are not a universal minimality proof. No biological transferclaim, no source edits, CPU1<=1GB/20MB. Output corrected scope, audit.json and one constructive witness or redundancy proof.
