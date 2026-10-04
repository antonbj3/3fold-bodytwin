# CX-NIGHTSYNTH2 — the night's summary, register entries and the next The swarm wave (runs on its own ~03:00; the coordinator is saving its coordinator quota)

## 1. Compile (read everything, recompute critical numbers yourself)
- `results/BT-HB-*` (113; ε-band independent of CX-EPSBAND + hamstring EMG vs meas−lo; NOTE the addendum at 23:15 in the later packets: the partial correlation and the loading phase are PRIMARY). Coverage/width per person and ε. Compare against `results/CX-EPSBAND/results.json` and Field's certified `../3fold-motion-engine/_private/romi_collab/build/U460/raw/coverage.json` (read-only). Where The swarm differs from those: say which is wrong and why (The swarm is weaker; check its code).
- `results/BT-KC-*` (8; SCoRE/SARA per person ×2 runs) against `results/CX-KNEECENTER/geometry.json`.
- `results/CX-L1ARMS`, `CX-EMGLEARN`, `CX-GAUGE2D` if they are finished (otherwise report their status).
- Field's night results in romi_collab/build/U47x–U48x (read-only), where they touch the knee.

## 2. Register
Append ONE row per new result to `notes/RESULTS_INDEX.md`. The next A-number = max(A-number in the file)+1; compute it immediately before writing. Use the same column format as the last rows, and mark "ej granskad" or "OGRANSKAD (The swarm)" as appropriate.

## 3. The morning summary
Write `external_research_path`, max 60 lines, in Swedish, for Anton. Cover:
- the state of the knee line (N1g 0.396 is the reference to beat; the model error = frames below Fmin; what L1ARMS showed);
- what holds and what does not;
- 3 proposed next nodes. Rules: no process narration,
  no verdict words ("kill"/"dead end").

## 4. The next The swarm wave
Build 40–100 packets in `results/BT-N2-*` that follow from the night's results. Each packet MUST compute something new against measured implant force or against an independent reference, with a frozen criterion in its BRIEF.md. No reanalysis/meta questions.
- Model on `tasks/hamband_packets.py`: real files only, NO symlinks, < 20 MB, `inputs/NIGHT_PREAMBLE.md` copied in, explicit metadata (implant side JW R/DM R/SC L/PS L; BW; the force column), and no jw_lungef1.
- Write QUEUE_PROPOSAL.txt with a value justification per packet group.
- Then append the lines to `tasks/lanes/bt_queue.txt` (format `A swarm <ID>`, round-robin over A/B/C). The coordinator has pre-approved this under the rule above.

lane runner has full permissions in the workspace; `~/projects/bodytwin` and romi_collab are read-only. No emails, pushes or credentials. Internal data stays local. Commit via `bash tasks/lanegit.sh add ... && bash tasks/lanegit.sh commit -m ...`.

## Supplement 23:55
Also read results/CX-L1ARMS2 and results/CX-DESIGNLOOP (the proof lane's night priorities §1/§4 in external_research_path — basis for Anton's presentation). The morning summary should have its own section per priority: what was completed, callable function, paired numbers, measured time.

## Second Run (05:00)
First read results/CX-NIGHTSYNTH/RESULTS.md and the morning summary; UPDATE it (do not write a new file) with everything added since then (BT-N1-results, L1ARMS2, DESIGNLOOP, Fields U47x–U48x). The next The_swarm wave is called BT-N2-*.

## Addendum 01:55
CX-NONGAITOPS (A1574) has built L1 operators for 34 DM/SC/PS non-walk trials in the external_media The_swarm wave BT-N2 should include the quad-arm test (CX-QUADARM geometry, correct Mkx characters: medial = 0,5·cj + Mkx/d) on these trials, one packet per trial, copied npz as real files; note the team cross check deviation in each BRIEF.
