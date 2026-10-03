# CX-KNEENET — knee contact force as a CONSTRAINT NET: variables in quantified relations, where each measurement tightens an edge

Idea (Anton: "constraint nets are hot"): stop choosing between a null model and a muscle model. Build the knee's quantities as a net of variables with quantified relations. The contact force is the variable that all edges tighten together.
- Variables: |GRF|, GRF direction/CoP, knee extension/flexion moment (ID), patellar tendon lever arm, knee flexion angle, speed, activity, EMG co-contraction index, body mass, contact force F.
- Edge status is TIGHT / OPEN / UNKNOWN, with evidence as {file, key, value}.
- The null model k·|GRF| is ONE edge. N1g beats everything because that edge is tight. The question is which other edges are tight enough to tighten F further.

## Format and mechanisms (read first)
- Format for constraint nets: `~/projects/cad-to-simulation-I/data/CONSTRAINT_NET_VEHICLE.json` (variables, edges, status, evidence), and `~/projects/bodytwin/data/CONSTRAINT_NET_PROJECTOR.json`. Map: `~/projects/graph_workspace/GRAPH_MAP.md` ("variables in relation" = the core).
- Graph engine: `~/projects/3fold_staging/3fold-graph-engine/src/graph_engine/` (identifiability_oed, observability_sigma_min_cert, fleet_oed, coverage_collapse_alarm, admission_bands, decorrelated_*, graph_hole_engine). Borrow mechanisms and name them in PREREG.
- Data and code: `results/L1` (prep/*.npz, code/, EMG in /mnt/shared_data/bodytwin_work/L1/bundle/emg/), BT-N12/N12b/N39, BT-B122/B128. Do not redo `results/CX-BEATN1G` (H1–H5 there, running in parallel) or `results/CX-EARLYSTANCE`; read them if they exist.

## Hypotheses (PREREG.md + sha256 BEFORE the first run; LOPO over the Grand Challenge persons, N39/L1 windows)
- K1 edge census: estimate each pairwise relation to F (and between inputs) on the training persons. Tightness = LOPO residual SD per stance-phase decile. Which edges are TIGHT in which phase? Hypothesis: |GRF| is tightest overall, while the knee moment and co-contraction are tighter in early stance.
- K2 factor-graph fusion: F = MAP over all edges, weighted by measured tightness (Gaussian factor graph; phase-dependent variances). Compare against N1g.
- K3 interval intersection: each edge gives an interval for F, from a calibrated band on the training persons. Prediction = the midpoint of the intersection. Report coverage and width vs N1g's band.
- K4 conflict = alarm: where edges disagree (the intersection is empty or narrow), does the disagreement predict N1g's error (AUC, correlation)? Borrow coverage_collapse_alarm/decision_cert.
- K5 densification/OED: which single additional measurement (EMG, ID moment, knee angle, speed, geometry) tightens F the most per unit of cost? Compute with identifiability/observability from the graph engine. This is "which measurement to buy".
- K6 net structure across activities: does activity enter as a variable that switches edge strength (k(activity)), or can a continuous variable (speed, flexion) replace the activity label?

## Criteria and counter-tests
- K2/K3: median RMSE ≥ 5 % below N1g in ≥ 3/4 LOPO rotations, and not worse in early stance. Placebo: permuted/time-shifted secondary variables must lose the gain.
- K4: AUC ≥ 0.7 for N1g error > median, against a placebo.
- Also report N1, N1g and B24. Symmetric skepticism.

## Deliverables: `results/CX-KNEENET/` only
- `KNEE_CONSTRAINT_NET.json` in the same format as the vehicle net: variables, edges with status and evidence, updated from the measurements.
- `RESULTS.md` starting with `# CX-KNEENET`, plus `results.json` and code with pytest.
- Local: nice, 2 threads, ≤ 60 s per test. Heavier runs via Modal or RUN_ON_CLOUD.md. Everything outside your directory is read-only.

## From the graph session (anton-df, 24/9 ~13:00) — mechanisms to use (read the code: ~/projects/graph_workspace/pub/src/graph_engine/)
Principle: N1g is the default. Correct only when a decision certificate says with confidence 1−α that N1g lies outside the physically feasible set. Otherwise N1g stands (at most an α-fraction of frames can get worse, which can be tested).
1. margin_net.py (H1/H3/H5, K2/K3): each physics bound is a requirement with a margin, m_lo=(N1g−F_min)/scale, m_hi=(F_max−N1g)/scale. Inputs (moments, capacity, EMG) are evidence sources with lineage. Declare shared model error with `shares: {group: ρ}`, otherwise the feasible set becomes too narrow. Report N_eff: are several physics terms really one piece of information?
2. decision_cert.py (the correction rule): the decision "N1g violates the margin in phase φ" → certify() gives p_flip. Correct only if p_flip ≤ α. Report the fraction of frames corrected and the error on exactly those frames vs N1g.
3. alarm.py (H2/H4, K4): sequential test with false-alarm guarantee on N1g's residuals per phase/person. It decides WHERE phase weighting or a personal k is justified.
4. Phases: no hard labels. Use posterior marginals from the HMM chain `~/projects/3fold_public/kernel-lane_runner-hmm-laplacian-20260923/src/kernel_engine/inference/hmm_laplacian.py` (XD2, 60/60 exact paths vs BT-Viterbi). k-best spread is NOT calibrated uncertainty; use the marginals.
Lesson (graph G34/G36/LC1): signals with a shared error source counted as independent give false certainty. Measure the design effect between moments, capacity and EMG before trusting a narrow feasible set.
