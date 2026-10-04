## Common rules (applies to every BT-lane)
- Work folder. Write ONLY under results/<LANE>/. Read-only everywhere else (other results folders, references/, external_mount, ~/projects/...). No git commit/push, no emails, no publishing, no downloading unless it is stated in the assignment.
- Swedish in README/RESULTS. Each digit in results/<LANE>/*.json.
- PREREG.md is written BEFORE the first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged last in PREREG and in RESULTS.
- Cargo: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, at most 2 processes. Before heavy driving: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, maximum 15 min). Top ≥ 4 GB RAM or GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Intermediate results > 50 MB → external_media<LANE>/. The system disk (/) is almost full.
- Import ALDRIG other agents' modules with `from X import *` and never write to their folders: set own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill other processes, never `nvidia-smi -q`.
- Report honestly: cases are as valuable as wins. No claims without figure and source file.
- End with results/<LANE>/RESULTS.md which starts with the line `# <LANE>` and contains: what was done, outcome per criterion (pass/fail), deviations, files.

# DS_GRAPH_CITATION_TRACE_20260923

User renewed reserve_worker authorization 2026-09-23 ~16:08; breakthrough hunt across domains.

OVERRIDING RUN SCOPE: write only results/DS_GRAPH_CITATION_TRACE_20260923/; source modules and other lanes read-only. Read applicable AGENTS and project START before work. Hash PREREG before first new numerical evaluation; exploratory changes separately documented. CPU2 (all BLAS/OMP/MKL/NUMEXPR), nice19, <=2GB unless the task explicitly authorizes one gated solve. Maximum new root-disk artifacts50MB; no bulk downloads, package/toolchain builds, GPU, cloud compute, source edits, public output or email. No child agents. No process killing. Never read credentials or billing settings. Do not wait on other lanes; snapshot completed referenced outputs and record missing dependencies. All numbers in machine-readable files. RESULTS must begin '# DS_GRAPH_CITATION_TRACE_20260923' in first5lines; report negative or incomplete results honestly.

Every innovation lane must deliver CONNECT.md: source observation -> mechanism -> equation -> operator -> representation -> assumptions -> discriminating test, and one concrete cross-domain mapping with units and failure condition. Separate audited theorem, empirical evidence, hypothesis and open obligation. Audits must remain independent; do not turn a new candidate into approval of its parent result.

TASK:
Real agent performance repair from frozen paired trial, not more genericpreamble. Read external_research_path and each workspace GRAPH_WORKING_VIEW_20260923 packet/schema; frozen trial files read-only. Failure mode: packet annotation became citationdestination, original file/lines lost. Build an isolated shared evidence-card renderer for BodyTwin and dental: original exactpath/hash/line span/briefquote separate from derived summary, status/regime and negative findings. Validate each card againstsourcebytes; changingaquote, wrongline/hash/scope mustfail. Use4actualpackets perproject wherepossible includingunlinkedmaterial and failedIM2/GEOM. Show explicitunknown for unsupportedannotation ratherthaninventline. Compare source-trace coverage with currentpacket and equal12KB cap, but no claimagentgain from selftest; prepare frozen newpairedinput only if useful, do NOT callmodels. Preserve all originalfacts/sourceadmission; no workingview/globalCLI edits. Output reusable renderer,cards,actualnegativefixturecounts and proposedintegrationpatch in ownresultdir. <=2GB/25MB.
