# MECHANISM CMC AS A SECOND, DECORRELATED MUSCLE-FORCE SOLVE (2026-07-21)

Executes the capstone's #1 recommendation (`docs/MECHANISM_LAYER_COMPLETENESS.md`): nearly every
downstream twin layer reuses the SAME single subject2/`walking1` Static-Optimization (SO) solve, so
their mutual "agreement" is internal consistency, not independent validation. `docs/
MECHANISM_ANKLE_RESERVE_FIX.md` made Computed Muscle Control (CMC) — a FORWARD-DYNAMIC,
tracking-controller-driven muscle solution, mechanistically independent of SO's per-frame
quasi-static optimization — complete the full 1.57s trial for the first time
(`data/msk_smoketest/subject2_walking1/ankle_reserve_fix/boost150/cmc/full_1p57s/`). This session
uses that completed run as a genuinely second, decorrelated solve and re-derives the knee-r/hip-r
contact-force headline — the SAME "in-vivo-anchored" number the capstone traces the common-mode
chain through — plus the metabolic-cost headline, FROM CMC INSTEAD OF SO, holding everything else
(model, IK kinematics, GRF, free-body-cut geometry, OrthoLoad anchor) fixed. Every number below is
machine-measured this session (`scripts/msk/cmc_second_solve.py`), reusing
`static_opt_knee.py`'s/`validate_hip_force.py`'s already-committed, purely-geometric crossing-muscle
machinery UNCHANGED, and `metabolic_cost.py`'s probe-builder UNCHANGED. Isolation respected:
bodytwin only, all upstream SO/CMC outputs read in place, all new outputs under
`data/msk_smoketest/subject2_walking1/cmc_second_solve/`, no git operations.

**Headline: CMC AGREES with SO on all 6 pre-registered comparisons (knee force, hip force,
Umberger2010 COT, Bhargava2004 COT — each via 2 extraction methods where applicable), within a
pre-registered 25%-relative tolerance, with real margin (largest gap 10.9%, smallest 1.2%). CMC
consistently reads 8-11% HIGHER than SO on 5/6 metrics (further over the OrthoLoad anchor, not
closer to it, not reversed) — a genuine, disclosed, non-random-looking pattern, not scatter. This
is a real strengthening of the SO-based knee/hip headline: a mechanistically decorrelated
forward-dynamics solve does NOT reveal the number to have been an SO-optimization artifact. It does
NOT fully break the common-mode: both solves still share the same subject-scaled model, the same
single IK/GRF trial, and (necessarily) a different reserve-actuator configuration — disclosed below,
argued not to matter for this specific headline.**

## Headline result

| comparison | SO (existing, fresh-reparsed) | CMC (this session) | rel. diff. | verdict (pre-registered <25%) |
|---|---:|---:|---:|---|
| knee-r, self-computed (BFS crossing-muscle) | 391.10 %BW | 427.31 %BW | 9.3% | **AGREE** |
| knee-r, official `opensim.JointReaction` | 391.11 %BW | 427.54 %BW | 9.3% | **AGREE** |
| hip-r, self-computed (BFS crossing-muscle) | 386.77 %BW | 428.53 %BW | 10.8% | **AGREE** |
| hip-r, official `opensim.JointReaction` | 387.04 %BW | 429.26 %BW | 10.9% | **AGREE** |
| Umberger2010, net COT | 7.6834 J/kg/m | 8.3198 J/kg/m | 8.3% | **AGREE** |
| Bhargava2004, net COT | 6.1977 J/kg/m | 6.1239 J/kg/m | 1.2% | **AGREE** |

Ratio to the OrthoLoad in-vivo anchor (re-verified live this session, not trusted from prior JSON:
knee median 258.22 %BW/n=72, hip median 273.93 %BW/n=162): SO sits at 1.41-1.52×; CMC sits at
1.56-1.66× — CMC clears the SAME `PASS-substantial` bar (ratio > 0.6, this cert family's own
convention) SO clears, on both joints and both extraction methods. Ratio to the Koelewijn 2019
in-vivo indirect-calorimetry anchor (3.95 J/kg/m interpolated level-walking): SO sits at 1.57-1.95×;
CMC sits at 1.55-2.11×.

## 1. Pre-registered claims and the forced adversaries (both stated before computing the table above)

**C1/C2 (force agreement):** `|CMC_X - SO_X| / SO_X < 0.25` for X in {self-computed, official-JR},
per joint. **C3 (metabolic agreement):** same 25%-relative rule for {Umberger2010, Bhargava2004} net
COT. All six held. The exploratory work that preceded fixing this threshold (confirming CMC's
states.sto double-header quirk, §4, and timing a probe-evaluation loop) never touched the actual
force or metabolic comparison numbers, so this is a genuine pre-registration, not a post-hoc-loosened
one — and the 25% band was cleared with real margin (worst case 10.9%, more than 2× inside the band),
not a knife-edge pass.

**Forced adversary, agreement direction (this result leans positive, so this is the adversary that
must be forced, not skipped):** could SO and CMC "agree" for an uninteresting reason — e.g. CMC's own
`use_fast_optimization_target` tracking-optimizer converging on an activation pattern similar to SO's
minimum-effort solution *by construction*, not because two independently-mechanised solves found the
same physical answer? Forced check: per-crossing-muscle activation correlation, SO's
`activation.sto` vs CMC's `states.sto` activation trace, interpolated onto the shared IK grid, for
all 31 muscles detected as knee-r-or-hip-r-crossing:

- Mean r = **0.857** across 31 muscles; **not** a uniform near-1.0 rubber stamp — range **0.276
  (`grac_r`) to 0.995** (several muscles, e.g. `glmin1_r`, `glmed3_r`).
- The low-correlation muscles are the ones with the LEAST anatomical/mechanical reason to converge:
  `grac_r` (0.276), `addmagDist_r` (0.321), `addmagMid_r` (0.336) are all secondary/redundant
  hip muscles with substantial null-space in any force-allocation problem (many combinations of
  their activation levels satisfy the same net joint moment). The HIGH-correlation muscles are the
  major, load-bearing, low-redundancy prime movers of stance-phase gait: `glmax1_r` (0.984),
  `glmed1/2/3_r` (0.987-0995), `iliacus_r` (0.989), `vaslat_r`/`vasmed_r` (0.989/0.975). This is the
  geometrically-expected signature of two DIFFERENT solvers converging on the same answer for the
  well-determined (load-bearing) part of the problem while disagreeing on the underdetermined
  (redundant-muscle) part — not the signature of one method silently copying the other (which would
  show uniformly high correlation everywhere, muscle-by-muscle, including the redundant ones).
  **Conclusion: real, non-trivial allocation differences exist per-muscle even though the aggregate
  net crossing-force converges to within ~9-11% — the "different mechanisms, similar net outcome"
  signature that strengthens, not undermines, the corroboration.**

**Forced adversary, divergence-guard direction (kept active even though the result agreed, per the
task's own instruction to be maximally symmetric):**
- Crossing-muscle SET identical SO-fed vs CMC-fed at both joints (**PASS** — pure model+pose
  geometry, independent of which force series is supplied; a mismatch would have been a bug, not a
  finding).
- Anatomical anchor (all 12 expected knee-r muscles, all 25 expected hip-r muscles) present in BOTH
  the SO-fed and CMC-fed detector runs (**PASS**).
- Hip negative control (`bfsh_r`, biceps-femoris-short-head, anatomically must not cross the hip):
  excluded in both SO-fed and CMC-fed runs (**PASS**).
- CMC force-reconstruction self-consistency (§4): replaying CMC's own recorded
  activation+fiber_length states through the model reproduces CMC's own reported
  `Actuation_force.sto` to **median 0.0003% / p90 0.0017% relative error** across 221,908
  muscle-frames (gate <2%, **PASS** with enormous margin) — the metabolic replay pipeline is not
  silently mis-extracting CMC's own data.
- Peak times fall strictly inside CMC's actually-recorded window (t=[0.03,1.56]s): knee peak at
  t=0.510s (identical to SO's own peak time), hip self-computed peak at t=0.550s (identical to SO),
  hip official-JR peak at t=0.540s (SO: 0.550s, a 1-sample/10ms shift on a broad plateau — not a red
  flag). **No peak lands in the unrecorded, would-be-extrapolated edge region** (**PASS**).

## 2. Disclosed methodological asymmetry, argued (not just asserted) not to drive the result

CMC's only-ever-successful run required BOOSTED ankle/subtalar/mtp reserves
(`SO_reserve_boosted.xml`, 150 N·m) to converge past 56%GC at all (`docs/
MECHANISM_ANKLE_RESERVE_FIX.md`); the SO-based "no-RRA baseline" this compares against used the
ORIGINAL weak reserve set (`walking1_reserveActuators.xml`) because SO never needed the boost. This
is a real, disclosed difference this script did not introduce and cannot remove (CMC categorically
fails without the boost). Two lines of evidence argue it does not explain the observed ~9-11% gap:

1. **Reserve actuators never enter the self-computed method at all** —
   `knee_crossing_muscles_and_forces` iterates `model.getMuscles()` only; reserves are
   generalized-coordinate actuators with no anatomical path, structurally invisible to the
   crossing-muscle detector. Yet the self-computed method shows the SAME ~9-11% gap as the
   official-JointReaction method (which DOES sum every actuator including reserves). If the
   reserve-boost asymmetry were the driver, self-computed and official-JR should have diverged from
   EACH OTHER; instead both stay within 0.1-0.2 %BW of each other on the CMC side (427.31 vs 427.54
   knee; 428.53 vs 429.26 hip) — mirroring SO's own near-perfect self/JR internal consistency
   (391.10 vs 391.11; 386.77 vs 387.04).
2. **Measured reserve magnitude is small next to crossing-muscle tension at these joints**: peak
   reserve usage during the CMC run tops out at 23.8% of a 150 N·m budget (`ankle_angle_l_reserve`,
   re-verified live this session), i.e. tens of N·m, next to peak crossing-muscle tension in the
   hundreds of N at the knee/hip specifically (this session's own
   `muscle_crossing_force_peak_pct_bw` ~300 %BW ≈ 2300 N).

## 3. Consistent direction, not scatter — disclosed, not mechanistically decomposed

5 of 6 metrics move the SAME direction under CMC (higher, i.e. further over both external anchors):
knee self +9.3%, knee JR +9.3%, hip self +10.8%, hip JR +10.9%, Umberger +8.3%. Only Bhargava2004 is
flat/slightly lower (−1.2%). A consistent-direction shift across independently-extracted quantities
is more informative than symmetric scatter would have been — it is the signature of a systematic
difference between the two solve methods, not measurement noise. **This session did not isolate the
mechanism** (open gap, §6) but notes it sits in the SAME direction as this repo's OTHER decorrelated
finding: `docs/MECHANISM_EMG_DRIVEN.md` (real surface EMG driving 7/80 muscles) found SO under-reads
co-contraction, with **+48.4% at the knee** vs SO. Two independently-built decorrelation efforts
(real EMG on a 7-muscle subset; full-80-muscle CMC forward dynamics) now agree on SIGN (SO reads
LOW) even though they disagree sharply on MAGNITUDE (48% vs 9%) — expected, since EMG-driven's 7/80
coverage and CMC's full-80-muscle forward-dynamics-with-tracking are very different partial/full
decorrelations of the same question, not repeated measurements of the same thing. This is reported
as a real, mutually-reinforcing triangulation on direction only, not overclaimed as magnitude
agreement.

## 4. Metabolic secondary analysis — feasibility, method, and a parsing bug found+fixed

**Feasibility discovery (forced, not assumed):** CMC's `states.sto` carries BOTH `activation` AND
`fiber_length` as genuine, time-varying integrated ODE states for all 80 muscles (verified live:
e.g. `vasmed_r` fiber_length std=0.0129m over the trial, not frozen) — richer than what SO provides
(SO has no fiber-length state at all; `metabolic_cost.py` must approximate it via
`equilibrateMuscles` under a forced rigid-tendon assumption, itself measured there at 8.27% median
force-reconstruction error). For CMC, replaying BOTH recorded states directly (no
`equilibrateMuscles` needed) and calling `realizeVelocity`/`realizeDynamics` reproduces CMC's own
reported tendon force almost exactly (§1) — a strictly MORE self-consistent reconstruction than SO's
own metabolic pipeline required.

**A real parsing bug found and fixed before trusting anything downstream:** CMC's full-trial
`states.sto` contains a duplicate, degenerate LEADING header block (`nRows=0 nColumns=1`, itself
terminated by its own `endheader`) immediately followed by the real header. Both this repo's
existing `vjf.parse_mot` AND `opensim.Storage`/`osim.StatesTrajectory.createFromStatesStorage`
silently stop at the FIRST `endheader` (`Storage`'s own log line reads `(nr=0 nc=1)`, and
`StatesTrajectory.createFromStatesStorage` then raises `"Table has no column labels"` — caught live,
not guessed). Fix: `parse_sto_last_endheader()` anchors on the LAST `endheader` occurrence — a
strict generalization that behaves identically on every OTHER `.sto` in this pipeline (confirmed via
`grep -c '^endheader$'`: `Actuation_force.sto`=1, `states.sto`=2) and correctly skips the bogus stub
here. Flagged for any future script touching a full-trial CMC `states.sto` in this repo.

**Result (all 2786 recorded frames replayed, no subsampling needed — 3.8s wall time, 1.35ms/frame):**
Umberger2010 net COT 8.32 J/kg/m (SO: 7.68, +8.3%), Bhargava2004 net COT 6.12 J/kg/m (SO: 6.20,
−1.2%). CMC's own COM-displacement-derived distance/speed (decorrelated from SO's own number):
1.6225m / 1.530s / 1.0605 m/s, vs SO's 1.6716m / 1.57s / 1.0647 m/s — consistent to within ~3%
despite being derived from an entirely separate (forward-dynamics-realized, not IK-prescribed)
kinematic trajectory; a useful, if coarse, partial answer to `docs/MECHANISM_ANKLE_RESERVE_FIX.md`'s
own disclosed gap #10 ("CMC's own forward-dynamics tracking accuracy was not separately measured").

## 5. What this DOES and does NOT resolve about the common-mode

**Resolved (narrowed):** the specific worry that the knee/hip "twin over-predicts OrthoLoad by
~1.4-1.6×" headline was an artifact of SO's particular per-frame minimum-effort optimization is NOT
supported — a mechanistically independent forward-dynamics solve, with its own separate tracking
controller, activation-contraction dynamics, and tendon compliance, lands in the same regime (both
methods clear `PASS-substantial` against OrthoLoad; CMC if anything overshoots by MORE, not less).

**NOT resolved (still common-mode):** both SO and CMC share the identical subject-scaled model
(same Fmax/geometry/scaling defects the Muscles-tier row of `docs/
MECHANISM_LAYER_COMPLETENESS.md` already flags as method-only-no-external-anchor for
subject-specificity), the identical single IK/GRF trial (subject2/`walking1`, no claim of
generality across subjects/trials/speeds), and CMC carries its OWN un-cleared defect
(`docs/MECHANISM_RRA_TASK_GAINS.md`: the Hicks moment-residual gate still fails for this exact
model/trial family, even after an exhaustive 9-configuration task-gain sweep) — so CMC is a
genuinely DIFFERENT solve, not a validated ground truth. The honest scope statement: this session
breaks the muscle-force-SOLVE-METHOD axis of the common-mode (the one the capstone's headline
sentence was about) and finds it holds up; it does not and cannot break the shared-model/shared-trial
axis, which remains the twin's deepest, still-unaddressed common-mode risk.

## 6. Honest gaps (full list)

1. **Self-computed and official-JointReaction are not mutually independent tests** — both consume
   the identical underlying force.sto (SO's or CMC's). The genuinely distinct quantities tested here
   are 4 (knee force, hip force, Umberger COT, Bhargava COT), not 6; the self/JR pairing is an
   internal-consistency check on THIS script's extraction, not a second independent physical test.
   Reported as 6 rows above for completeness, but not to be read as 6 independent confirmations.
2. **The consistent +8-11% direction (§3) is disclosed, not mechanistically isolated** — which
   specific difference between SO's minimum-effort optimization and CMC's forward-dynamics tracking
   (activation dynamics lag, real tendon compliance, the tracking controller's own error-correction
   demand, or something else) drives it was not decomposed this session.
3. **Umberger-vs-Bhargava asymmetry (+8.3% vs −1.2%) is reported, not explained** — which specific
   term in each probe's cost equation is more/less sensitive to CMC's real fiber-velocity-from-
   tendon-compliance dynamics vs SO's reconstructed-under-rigid-tendon fiber velocity was not
   isolated.
4. **Single CMC run available** (`boost150`, the only configuration that ever completed the full
   trial) — no replicate/seed variability or alternate boost-level (`boost300`, per
   `docs/MECHANISM_ANKLE_RESERVE_FIX.md`'s own disclosed gap, was never re-run for CMC) check was
   possible; this comparison uses the ONE existing successful CMC trial, not an ensemble.
5. **Single trial, right-leg-primary** (subject2 `walking1`) — same scope caveat as every upstream
   cert in this family; no claim of generality across subjects/trials/gait speeds/models.
6. **CMC's own kinematic-tracking fidelity was only coarsely checked** (net COM displacement, §4),
   not per-coordinate (the way `docs/MECHANISM_ANKLE_RESERVE_FIX.md` §3 characterized RRA's tracking
   error) — a full per-coordinate CMC-vs-IK tracking-error table remains unmeasured.
7. **The activation-correlation adversary check (§1) covers only the 31 muscles already flagged as
   knee-r-or-hip-r-crossing**, not all 80 muscles in the model.
8. **CMC's own Hicks-moment-gate-failing pelvis residual** (`docs/MECHANISM_RRA_TASK_GAINS.md`,
   inherited into the CMC pipeline via the shared `RRA_tasks.xml`/reserve-actuator convention) means
   CMC is a genuinely different, but not bias-free, solve — it is NOT treated as ground truth
   anywhere in this document, only as a decorrelated second opinion.
9. **The reserve-boost asymmetry argument (§2) is a measured, disclosed argument, not a
   counterfactual experiment** — SO was never re-run with the boosted reserve set (nor CMC with the
   weak set, since CMC cannot converge with it), so the claim rests on the reserve-magnitude/
   crossing-invisibility argument, not a direct ablation.

## Files

- `scripts/msk/cmc_second_solve.py` — the full, re-runnable pipeline. Reuses
  `validate_joint_force.py`'s (`vjf`) `parse_mot`/`get_descendant_bodies`/OrthoLoad-knee-parser/
  `MODEL_FILE`/`IK_MOT`/`GRF_MOT`/`G`; `static_opt_knee.py`'s (`sok`) crossing-muscle geometry,
  `sub_tag`/`sub_tag_allow_selfclosing`/`patch_external_loads`/`assert_no_template_paths_survived`/
  `JR_SETUP_TEMPLATE`/`check_so_convergence_and_sanity`/`extract_jr_knee_r_force_pct_bw`/anatomical
  and absurdity-ceiling constants, all UNCHANGED; `validate_hip_force.py`'s (`vhf`) generic
  `compute_self_cross_check_generic`/OrthoLoad-hip-parser/`extract_jr_hip_r_force_pct_bw`/hip
  anatomical+negative-control constants, all UNCHANGED; `metabolic_cost.py`'s (`mc`) `build_model`/
  `trapz_mean`/Koelewijn-anchor constants, UNCHANGED. New code added: `parse_sto_last_endheader`
  (the states.sto double-header fix, §4), `run_joint_reaction_generic` (a parametrized sibling of
  `sok.run_joint_reaction` — wrap, don't edit in place, per `COORDINATOR.md` Sec.2 — needed because the
  CMC run requires a different reserve-actuator file and a narrower, non-extrapolated time window),
  `classify_agreement`, and `run_metabolic_from_cmc_states` (the CMC-states replay loop).
- `data/msk_smoketest/subject2_walking1/cmc_second_solve/` — all NEW outputs: `jr_cmc/` (the one new
  `opensim.JointReaction` AnalyzeTool pass, fed CMC's force.sto), `cmc_second_solve_results.json`
  (every number in this document, traceable back to this one file).
- Inputs read in place, never modified: `data/msk_smoketest/subject2_walking1/static_optimization/`
  (SO's existing force.sto + JointReaction output), `data/msk_smoketest/subject2_walking1/
  ankle_reserve_fix/boost150/cmc/full_1p57s/` (CMC's existing force.sto + states.sto), `data/
  msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json` (SO's existing
  metabolic headline).

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` Sec.1 and the task).
