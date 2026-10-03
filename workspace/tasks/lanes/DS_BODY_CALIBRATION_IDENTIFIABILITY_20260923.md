## Common rules (applies to every BT-lane)
- Arbetsmapp .. Skriv ENDAST under results/<LANE>/. Read-only everywhere else (andra results-mappar, references/, /mnt/..., ~/projects/...). Ingen git commit/push, no e-mails, no publishing, no download without it being stated in the assignment.
- Swedish in README/RESULTS. Each digit in results/<LANE>/*.json.
- PREREG.md is written BEFORE the first run and hashed: `sha256sum PREREG.md > PREREG.sha256`. Never change criteria afterwards; deviations are logged last in PREREG and in RESULTS.
- Last: `nice -n 19`, OMP_NUM_THREADS=2, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, not more than 2 processes. Before heavy driving: `awk '{print $1/20}' /proc/loadavg` > 1,0 → wait (sleep 60, not more than 15 min). Topp ≥ 4 GB RAM eller GPU → `tasks/heavy_run.sh <ram_gb> <vram_gb> <kommando>`. Mellanresultat > 50 MB → /media/anton/sdc1-tmp/bodytwin/<LANE>/. Systemdisken (/) I'm almost drunk.
- Import ALDRIG other agents' modules with `from X import *` and never write to their folders: set own OUTDIR explicitly and check that no path variables are shadowed (incident H2b 2026-09-23 overwrote results/H2).
- Never `pgrep -f`, never kill other processes, never `nvidia-smi -q`.
- Report honestly: cases are as valuable as wins. No claims without figure and source file.
- End with results/<LANE>/RESULTS.md which starts with the line `# <LANE>` and contains: what was done, outcome per criterion (pass/fail), deviations, files.

# DS_BODY_CALIBRATION_IDENTIFIABILITY_20260923

User renewed reserve_worker authorization 2026-09-23 ~16:08; breakthrough hunt across domains.

OVERRIDING RUN SCOPE: write only ./results/DS_BODY_CALIBRATION_IDENTIFIABILITY_20260923/; source modules and other lanes read-only. Read applicable AGENTS and project START before work. Hash PREREG before first new numerical evaluation; exploratory changes separately documented. CPU2 (all BLAS/OMP/MKL/NUMEXPR), nice19, <=2GB unless the task explicitly authorizes one gated solve. Maximum new root-disk artifacts50MB; no bulk downloads, package/toolchain builds, GPU, cloud compute, source edits, public output or email. No child agents. No process killing. Never read credentials or billing settings. Do not wait on other lanes; snapshot completed referenced outputs and record missing dependencies. All numbers in machine-readable files. RESULTS must begin '# DS_BODY_CALIBRATION_IDENTIFIABILITY_20260923' in first5lines; report negative or incomplete results honestly.

Every innovation lane must deliver CONNECT.md: source observation -> mechanism -> equation -> operator -> representation -> assumptions -> discriminating test, and one concrete cross-domain mapping with units and failure condition. Separate audited theorem, empirical evidence, hypothesis and open obligation. Audits must remain independent; do not turn a new candidate into approval of its parent result.

TASK:
Expand measured force-channel idea from /mnt/shared_data/research_compute/AB_HUNT_WAVE2_20260923/jobs/W2B02/W2B04: one added channel'sgain mayvanish with unknown sensor scale,offset orsharedcalibration. Build exact sensitivity/nuisance model with explicit actualobservationtiming. Compare candidatechannels underboundedcalibration andunknowncalibration; constructnullspace collisions and one minimalcalibration intervention that restores a targetdirection. No counterfactual sameperson features, no assumed independenterrors becausemodalitynamesdiffer. Use a strict held-paramfamily baseline and countmeasurementscost. Readgraphworking rank/packet onexisting material/geometryobservability node ifavailable; emitpendingfeedback inownlane only. CPU1<=60s noFE/downloads.
