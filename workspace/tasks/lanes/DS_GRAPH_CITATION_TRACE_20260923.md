## Common rules (apply to every BT lane)
- Working directory . Write ONLY under results/<LANE>/. Read-only everywhere else (other results directories, references/, external_mount, ~/projects/...). No git commit/push, no emails, no publishing, no downloads unless specified in the assignment.
- Swedish in README/RESULTS. Every number in results/<LANE>/*.json.
- PREREG.md is written BEFORE the first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged at the end of PREREG and in RESULTS.
- Load: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, at most 2 processes. Before a heavy run: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, at most 15 min). Peak ≥ 4 GB RAM or GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Intermediate results > 50 MB → external_media<LANE>/. The system disk (/) is almost full.
- NEVER import other agents’ modules with `from X import *` and never write to their directories: set your own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill other people’s processes, never `nvidia-smi -q`.
- Report honestly: failures are as valuable as gains. No claims without a number and source file.
- Finish with results/<LANE>/RESULTS.md starting with the line `# <LANE>` and containing: what was done, outcome per criterion (holds/failed), deviations, files.

# DS_GRAPH_CITATION_TRACE_20260923

User renewed swarm_worker authorization 2026-09-23 ~16:08; breakthrough hunt across domains.

OVERRIDING RUN SCOPE: write only results/DS_GRAPH_CITATION_TRACE_20260923/; source modules and other lanes read-only. Read applicable AGENTS and project START before work. Hash PREREG before first new numerical evaluation; exploratory changes separately documented. CPU2 (all BLAS/OMP/MKL/NUMEXPR), nice19, <=2GB unless the task explicitly authorizes one gated solve. Maximum new root-disk artifacts50MB; no bulk downloads, package/toolchain builds, GPU, cloud compute, source edits, public output or email. No child agents. No process killing. Never read credentials or billing settings. Do not wait on other lanes; snapshot completed referenced outputs and record missing dependencies. All numbers in machine-readable files. RESULTS must begin '# DS_GRAPH_CITATION_TRACE_20260923' in first5lines; report negative or incomplete results honestly.

Every innovation lane must deliver CONNECT.md: source observation -> mechanism -> equation -> operator -> representation -> assumptions -> discriminating test, and one concrete cross-domain mapping with units and failure condition. Separate audited theorem, empirical evidence, hypothesis and open obligation. Audits must remain independent; do not turn a new candidate into approval of its parent result.

TASK:
Real agent performance repair from frozen paired trial, not more genericpreamble. Read external_research_path and each workspace GRAPH_WORKING_VIEW_20260923 packet/schema; frozen trial files read-only. Failure mode: packet annotation became citationdestination, original file/lines lost. Build an isolated shared evidence-card renderer for BodyTwin and dental: original exactpath/hash/line span/briefquote separate from derived summary, status/regime and negative findings. Validate each card againstsourcebytes; changingaquote, wrongline/hash/scope mustfail. Use4actualpackets perproject wherepossible includingunlinkedmaterial and failedIM2/GEOM. Show explicitunknown for unsupportedannotation ratherthaninventline. Compare source-trace coverage with currentpacket and equal12KB cap, but no claimagentgain from selftest; prepare frozen newpairedinput only if useful, do NOT callmodels. Preserve all originalfacts/sourceadmission; no workingview/globalCLI edits. Output reusable renderer,cards,actualnegativefixturecounts and proposedintegrationpatch in ownresultdir. <=2GB/25MB.
