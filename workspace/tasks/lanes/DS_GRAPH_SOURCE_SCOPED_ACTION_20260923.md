## Common rules (applies to each BT-lane)
- Work Folder . Type ONLY under results/<LANE>/. Read-only everywhere else (other results folders, references/, external_mount, ~/projects/...). No git commit/push, no emails, no publishing, no downloading unless it says so in the assignment.
- Swedish in README/RESULTS. Every digit in results/<LANE>/*.json.
- PREREG.md is written BEFORE first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged last in PREREG and in RESULTS.
- Load: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, at most 2 processes. Before heavy driving: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, maximum 15 min). Top ≥ 4 GB RAM or GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Intermediate results > 50 MB → external_media<LANE>/. The system disk (/) is almost full.
- Import ALDRIG other agents' modules with `from X import *` and never write to their folders: set own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill other processes, never `nvidia-smi -q`.
- Report honestly: cases are as valuable as wins. No claims without figure and source file.
- End with results/<LANE>/RESULTS.md starting with the line `# <LANE>` and containing: what was done, outcome per criterion (passed/failed), deviations, files.

# DS_GRAPH_SOURCE_SCOPED_ACTION_20260923

User renewed swarm_worker authorization 2026-09-23 ~16:08; breakthrough hunt across domains.

OVERRIDING RUN SCOPE: write only results/DS_GRAPH_SOURCE_SCOPED_ACTION_20260923/; source modules and other lanes read-only. Read applicable AGENTS and project START before work. Hash PREREG before first new numerical evaluation; exploratory changes separately documented. CPU2 (all BLAS/OMP/MKL/NUMEXPR), nice19, <=2GB unless the task explicitly authorizes one gated solve. Maximum new root-disk artifacts50MB; no bulk downloads, package/toolchain builds, GPU, cloud compute, source edits, public output or email. No child agents. No process killing. Never read credentials or billing settings. Do not wait on other lanes; snapshot completed referenced outputs and record missing dependencies. All numbers in machine-readable files. RESULTS must begin '# DS_GRAPH_SOURCE_SCOPED_ACTION_20260923' in first5lines; report negative or incomplete results honestly.

Every innovation lane must deliver CONNECT.md: source observation -> mechanism -> equation -> operator -> representation -> assumptions -> discriminating test, and one concrete cross-domain mapping with units and failure condition. Separate audited theorem, empirical evidence, hypothesis and open obligation. Audits must remain independent; do not turn a new candidate into approval of its parent result.

TASK:
Freeze the same 24 cases and source hashes. Add a source-scoped query-aware raw index that selects workspace/domain before BM25, then ask whether it retains raw-index fact recall while reducing cross-domain leakage and improving evidence-backed next-experiment choice. Score action choice from a preregistered exact fact/negative-result/gate checklist without new model calls; compare with original raw index, packet routing and a same-byte domain-only ablation. Keep human-authored summary text separate from typed links.
EVIDENCE TO VERIFY: Across 24 questions, query-aware raw index recalled 0.780 source facts versus packet routing 0.583, but raw index included 3.29 cross-domain sources per question; typed-link ablation recalled only 0.292.
STRONG BASELINE: Original query-aware raw index and packet arm at the same <=12 KB context budget; domain-only index ablation distinguishes routing from graph information.
FALSIFIER: Domain scoping loses relevant evidence or does not improve supported action selection after the strong raw-index baseline is allowed equal bytes.
LIMITS: CPU1, RSS <=1 GB, <=20 MB new files; deterministic retrieval/procedural proxy only, no model calls, graph refresh or status promotion.
Distinguish a finite saved-grid identification failure from structural or population nonidentifiability; report source gaps rather than proving unobserved states from missing data. Read completed upstream full report before using any quoted number.
