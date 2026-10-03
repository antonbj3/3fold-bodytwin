# CX-DMCOMPUTE — compute the data matrix DETERMINISTICALLY (no agents), then design analysis packets across all trials

Coordinator's review of CX-DATAMATRIX: 605 packets that each compute 5 simple statistics from ONE raw CSV (peaks, impulses, RMS, active fraction, widths) are pure feature extraction. A script does it for all trials in seconds, so it must not be spent on 605 agent runs. The swarm's value is in ANALYSING the complete tables across many trials.

## 1. Compute (script, locally, ≤ a few min)
- Implement the 8 protocols in `results/BT-DM-*/inputs/PROTOCOL.json` exactly: P-GRF-PEAK, P-GRF-IMPULSE, P-GRF-BALANCE, P-EMG-TIMING, P-EMG-AMPLITUDE, P-MARKER-GEOM, P-IMU-VIRTUAL, P-CONTACT-COMP.
- Use `inputs/raw.csv` + reader.py per packet. Implement in `tasks/datamatrix_compute.py`.
- Compute all 605. Write `results/DATAMATRIX_TABLES/<protocol>.csv`, one row per trial with provenance (raw_sha256), and the odd/even check value.
- Add a stronger check where possible:
  - repeated trials of the same activity/person → repeatability (CV);
  - left/right for bilateral quantities;
  - implant force vs |GRF| for P-CONTACT-COMP.
- Sanity: units, missing channels (NaN, not 0), outliers flagged (not removed).
- Also mark the 605 BT-DM packets as already computed (write TABLE.csv in each and a RESULTS.md stating "computed deterministically by CX-DMCOMPUTE"), so they are never queued.
- The F-8 exception: `jw_lungef1` implant/knee_forces must not be computed or copied (Field's sealed facit). Check.

## 2. Analysis packets for The swarm (the "several flies" packets) — `results/BT-DA-001…` (up to 40)
Each packet has:
- the complete relevant table(s) + a list of the trials' activity/person labels in inputs;
- a question that gives new knowledge or model input across MANY trials;
- ≥ 2 consumers;
- a check that is not circular (hold out a person/activity; permute labels).
Examples:
- (a) N1g: fit k per activity against the MEASURED implant force (P-CONTACT-COMP) over ALL GC activities and persons, LOPO → a new, extended N1g for non-gait activities (lunge excluded: F-8);
- (b) which GRF features beyond |GRF| explain the implant force residual (the peak ratio, impulse, balance), held-out person;
- (c) EMG timing/amplitude vs implant force: does co-activation carry information beyond GRF (compare A324 H3/A325);
- (d) the virtual IMU (pelvis) → implant force estimator across activities (A321 sensor track), held out;
- (e) marker geometry (knee/ankle/pelvis width) per person vs the CT/anthropometric priors of the geometry manager;
- (f) the variation between repeated trials as the measurement-noise floor for all our models (what error is irreducible?).
Write `results/CX-DMCOMPUTE/QUEUE_PROPOSAL.txt` (`<A|B|C> swarm BT-DA-xxx | question | consumers | check`). **Do NOT queue** — the coordinator reviews.

`results/CX-DMCOMPUTE/RESULTS.md` starting with `# CX-DMCOMPUTE`: rows per protocol, the check statistics, flagged outliers, and the analysis-packet list. lane runner has full permissions in the workspace; the raw zips are read-only.
