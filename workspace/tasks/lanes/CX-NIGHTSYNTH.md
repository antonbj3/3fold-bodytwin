# CX-NIGHTSYNTH — the night's summary, register entries and the next The swarm wave (runs on its own ~03:00; the coordinator is saving its coordinator quota)

## 1. Compile (read everything, recompute critical numbers yourself)
- `results/BT-HB-*` (113; ε-band independent of CX-EPSBAND + hamstring EMG vs meas−lo; NOTE the addendum at 23:15 in the later packets: the partial correlation and the loading phase are PRIMARY). Coverage/width per person and ε. Compare against `results/CX-EPSBAND/results.json` and Field's certified `the public staging tree/3fold-motion-engine/_private/romi_collab/build/U460/raw/coverage.json` (read-only). Where The swarm differs from those: say which is wrong and why (The swarm is weaker; check its code).
- `results/BT-KC-*` (8; SCoRE/SARA per person ×2 runs) against `results/CX-KNEECENTER/geometry.json`.
- `results/CX-L1ARMS`, `CX-EMGLEARN`, `CX-GAUGE2D` if they are finished (otherwise report their status).
- Field's night results in romi_collab/build/U47x–U48x (read-only), where they touch the knee.

## 2. Register
Append ONE row per new result to `notes/RESULTS_INDEX.md`. The next A-number = max(A-number in the file)+1; compute it immediately before writing. Use the same column format as the last rows, and mark "ej granskad" or "OGRANSKAD (The swarm)" as appropriate.

## 3. The morning summary
Write `~/research/ANYBODY_PERFORMANCE_20260924/session_input/bodytwin_night_20260925.md`, max 60 lines, in Swedish, for Anton. Cover:
- the state of the knee line (N1g 0.396 is the reference to beat; the model error = frames below Fmin; what L1ARMS showed);
- what holds and what does not;
- 3 proposed next nodes. Rules: no process narration,
  no verdict words ("kill"/"dead end").

## 4. The next The swarm wave
Build 40–100 packets in `results/BT-N1-*` that follow from the night's results. Each packet MUST compute something new against measured implant force or against an independent reference, with a frozen criterion in its BRIEF.md. No reanalysis/meta questions.
- Model on `tasks/hamband_packets.py`: real files only, NO symlinks, < 20 MB, `inputs/NIGHT_PREAMBLE.md` copied in, explicit metadata (implant side JW R/DM R/SC L/PS L; BW; the force column), and no jw_lungef1.
- Write QUEUE_PROPOSAL.txt with a value justification per packet group.
- Then append the lines to `tasks/lanes/bt_queue.txt` (format `A swarm <ID>`, round-robin over A/B/C). The coordinator has pre-approved this under the rule above.

lane runner has full permissions in the workspace; `~/projects/bodytwin` and romi_collab are read-only. No emails, pushes or credentials. Internal data stays local. Commit via `bash tasks/lanegit.sh add ... && bash tasks/lanegit.sh commit -m ...`.

## Supplement 23:55
Read also results/CX-L1ARMS2 and results/CX-DESIGNLOOP (Astras nattprioriteringar §1/§4 i ~/research/JOHN_PRESENTATION_20260924/NIGHT_PRIORITIES.md — Background for Anton's presentation). The morning summary shall have its own section per priority: what was completed, callable function, paired figures, measured time.

## Supplement 00:45
A1449: quadriceps arm alone → 9,71 % under lo. Raise this first in the summary (§1). Read results/CX-QUADARM if done. The swarm wave is welcome to contain independent reproductions of the solo change (per trial: just change the quad arm in A[3] according to results/CX-L1ARMS2/public_curves.json and count lo with HiGHS) — copy the npz + curve file as real files.

## Supplement 00:55
A1450 CX-QUADARM (Im curve confirms underway; JW's patellar release ρ shear → 4,4 mm arm at 50°) and A1451 CX-DESIGNLOOP (callable design loop works). Both in §1/§4 in the summary.
