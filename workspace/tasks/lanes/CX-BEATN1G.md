# CX-BEATN1G — beat the null model N1g (k(activity)·|GRF|) with physics as a CONSTRAINT/SUPPLEMENT, not as a replacement

State: L1 (A312, `results/L1/`) shows that no recruitment criterion beats N1g (0.378 BW median LOPO). Earlier tries: N52 rule hybrid (A281), N56 physics+data correction (A294, the model term carried no information), N57 calibration (A299), B141 ranking. The error sits in early stance (+0.67–1.02 BW). All of these tried to make the muscle model better than N1g. Here N1g is the base and physics adds only what N1g cannot know.

## Build on (read first; do not redo)
- `results/L1` (code/, prep/*.npz, scores.json, FOR_CX_MUSCLE.md; EMG in /mnt/shared_data/bodytwin_work/L1/bundle/emg/).
- BT-N12/N12b/N39 (N1/N1g scoring, COMAK frames), BT-B122/B128 (N1g), BT-N52, BT-N56, BT-N57.
- `results/CX-EARLYSTANCE` (running in parallel: early-stance attribution — read it if it exists, do not duplicate).
- Inspiration, read before designing:
  - The graph engine: `~/projects/3fold_staging` and the graph workspace `~/projects/3fold-workspaces/` (grep decision_cert, margin_net, alarm, constraint net, set value/exact set value, "densification").
  - The field engine: certified error band/decision layer, shape derivatives (grep FIELD_RESULTS, U370/U374, shape derivative).
  - Write in PREREG which mechanism you borrow and from where.

## Hypotheses (PREREG.md + sha256 BEFORE the first run; LOPO over the Grand Challenge persons, same windows and normalisation as N39/L1)
- H1 feasible set: from the knee moments (ID) plus the muscle model's capacity, compute the smallest/largest physically possible contact force per frame, F_min/F_max (LP with muscle forces 0…Fmax). Prediction = clip(N1g, F_min, F_max). Also report how often N1g falls outside the set, and where.
- H2 phase-wise uncertainty weighting: F = w(phase)·N1g + (1−w)·model, where w is estimated on the training persons from the error variance per stance-phase decile.
- H3 co-contraction from EMG: F = k·|GRF| + c·CCI(t), where CCI is quadriceps/hamstrings/gastrocnemius co-contraction from normalised EMG (MaxEMG trials exist). k and c are fitted on the training persons.
- H4 personal k: k_p from one calibration trial per person (standing/slow gait) instead of the population k(activity).
- H5 moment term: F = k·|GRF| + b·|M_knee_ext|/r_pt, with the moment from GRF × lever arm and a patellar lever arm from geometry, but no muscle distribution.
- Your own variants are allowed but must be preregistered before being run.

## Criteria and counter-tests
- Median RMSE over persons at least 5 % lower than N1g in ≥ 3/4 LOPO rotations, and not worse in early stance.
- For each hypothesis, a placebo: permuted/time-shifted EMG/moment/limits. It must lose the gain. Otherwise the gain comes from the degree of freedom, not the physics.
- Also report B24 and N1. Symmetric skepticism: if nothing wins, say so and give the best explanation of why N1g is so hard to beat, with a number (e.g. how much of the variance |GRF| explains).

## Resources
- Local: nice, 2 threads, ≤ 60 s per test. Heavier runs via Modal (`tasks/modal_run.py`) or RUN_ON_CLOUD.md.
- Write only in `results/CX-BEATN1G/`. `~/projects/bodytwin` and other lanes' directories are read-only.
- `RESULTS.md` starting with `# CX-BEATN1G`, plus `results.json` and code with pytest.

## From the graph session (anton-df, 24/9 ~13:00) — mechanisms to use (read the code: ~/projects/graph_workspace/pub/src/graph_engine/)
Principle: N1g is the default. Correct only when a decision certificate says with confidence 1−α that N1g lies outside the physically feasible set. Otherwise N1g stands (at most an α-fraction of frames can get worse, which can be tested).
1. margin_net.py (H1/H3/H5, K2/K3): each physics bound is a requirement with a margin, m_lo=(N1g−F_min)/scale, m_hi=(F_max−N1g)/scale. Inputs (moments, capacity, EMG) are evidence sources with lineage. Declare shared model error with `shares: {group: ρ}`, otherwise the feasible set becomes too narrow. Report N_eff: are several physics terms really one piece of information?
2. decision_cert.py (the correction rule): the decision "N1g violates the margin in phase φ" → certify() gives p_flip. Correct only if p_flip ≤ α. Report the fraction of frames corrected and the error on exactly those frames vs N1g.
3. alarm.py (H2/H4, K4): sequential test with false-alarm guarantee on N1g's residuals per phase/person. It decides WHERE phase weighting or a personal k is justified.
4. Phases: no hard labels. Use posterior marginals from the HMM chain `~/projects/3fold_public/kernel-lane_runner-hmm-laplacian-20260923/src/kernel_engine/inference/hmm_laplacian.py` (XD2, 60/60 exact paths vs BT-Viterbi). k-best spread is NOT calibrated uncertainty; use the marginals.
Lesson (graph G34/G36/LC1): signals with a shared error source counted as independent give false certainty. Measure the design effect between moments, capacity and EMG before trusting a narrow feasible set.
