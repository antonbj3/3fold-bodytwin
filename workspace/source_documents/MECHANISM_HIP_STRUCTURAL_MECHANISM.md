# MECHANISM HIP STRUCTURAL MECHANISM — digging the SO-then-subtract over-count (2026-07-21)

Digs the Moissenet et al. 2014 (PMID 24210475) "SO-then-subtract... usually overestimates" critique
concretely on subject2/`walking1`'s hip peak (t=0.55s, 386.77 %BW self-computed vs OrthoLoad hip
median **273.93 %BW**, ratio **1.412x**) — the one joint the batch could not explain by strength
(`docs/MECHANISM_FMAX_PCSA_VALIDATION.md`: correction nets to ~0), by any single global knob
(`docs/MECHANISM_OVERPREDICTION_DECOMP.md`: nothing clears the 0.10-ratio bar), and already shown to be
DIFFUSE (`docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`: top-2 groups = 47.5%, glmed1_r over-recruited 5.9x
vs real EMG).

**Headline (joint-differentiated, not a single verdict — matching this doc family's own established
discipline): a specific, nameable, externally-EMG-anchored over-count IS localized — gluteus medius
(`glmed1_r`) is over-recruited by Static Optimization relative to this subject's own real surface EMG
(585.3N SO vs 99.76N EMG-implied, ratio 0.17) — and a MOMENT-BUDGET-PRESERVING correction of it closes
a MATERIAL, ROBUST fraction of the gap: 53-79% of the distance from 1.412 toward 1.0, depending on
whether the remaining 24 muscles are reallocated via a smooth (physiologically graded, ratio→1.194)
or bang-bang (LP-optimal, ratio→1.086) redistribution — both clear the pre-registered 0.10-ratio
materiality bar with a comfortable margin. This is DIFFERENT from, and stronger than, the diffuse
sub-threshold-recruitment trim: the 15 near-floor-activation muscles' own DIRECT contribution to the
contact force is small (only 3.3% of the total, naive-trim ratio 1.365, sub-material) and their
"valid" LP-rebalanced improvement turns out to be a RE-DISCOVERY of the same glmed1_r/recfem_r-zeroing
lever, not an independent mechanism — and its own smooth-objective cross-check is knife-edge (Delta
0.102, barely above the 0.10 bar). The bare geometric floor of the vector-sum (ignoring which muscle
is at fault) is even lower (ratio 1.065) but requires a non-physiological bang-bang pattern that a
smooth alternative objective does NOT reach (Delta 0.087, sub-material) — so "vector geometry alone"
is a real but WEAK, largely glmed1_r-attributable effect, not a broad diffuse-summing artifact. Even
after the most generous, valid correction, a 9-19% residual over-prediction remains — a genuine,
disclosed diagnosed-gap on top of the localized mechanism.**

## 0. Pre-registration (stated before any LP was solved)

- **MATERIALITY GATE** (reused verbatim, same convention as `docs/MECHANISM_OVERPREDICTION_DECOMP.md`):
  a mechanism is material iff it moves `hip_ratio` (vs OrthoLoad 273.93 %BW) by **>=0.10 absolute**,
  in the CLOSING direction (toward 1.0).
- **VALIDITY RULE** (this task's own symmetric-QC clause, verbatim): a trimmed/corrected muscle set
  that does NOT reproduce the SAME 3-DOF hip moment budget is **INVALID** — reported as a diagnostic
  only, never as evidence of "closing the gap."
- **ACTIVATION-FLOOR THRESHOLD**: read from the MEASURED distribution, never tuned to the answer —
  `contact_muscle_decomp_results.json`'s hip activation column has a clean bimodal gap: 15 muscles in
  [0.0100,0.0153] (the SO numerical floor) vs 10 in [0.1126,0.5740] — a 7.3x gap with nothing between.
  Primary threshold = 0.05 (centered in the gap); swept at {0.02,0.05,0.10} to confirm the trimmed SET
  is identical everywhere in the gap (non-knife-edge on the THRESHOLD — separate from the knife-edge
  found on the smooth-vs-bang-bang CHOICE, Sec.3).
- **FALSIFIER**: if trimming diffuse recruitment OR correcting glmed1_r to EMG moves the ratio
  materially, that localizes a specific, nameable over-count; if nothing structural moves it, the
  over-prediction is intrinsic to the free-body-cut method — a diagnosed-gap pointing at the method.
- **Symmetric-QC commitment**: no "fix" is motivated by hitting 273.93 %BW — every correction is
  anchored to an independent principle (real EMG, or a moment-balance identity) and its own validity
  is machine-checked (does the correction still satisfy the hip's actual moment budget?).

## 1. Geometric method (exact bounded-variable linear program, not a heuristic)

Every hip-crossing muscle *m*'s contribution to the compressive contact-force axis is `c_m = k_m *
tension_m`, `k_m = -alignment_with_contact_axis_m` — reused EXACTLY from
`contact_muscle_decomp_results.json` (an already-verified exact identity: `R_pure + sum(c_m)`
reconstitutes the published 386.77 %BW to 1e-11 relative). Every muscle also contributes `tension_m *
r_m^(k)` to the net moment about each of the hip's **3 generalized coordinates** — `hip_flexion_r`,
`hip_adduction_r`, `hip_rotation_r` (a true 3-DOF ball joint in this model, verified live via
`computeMomentArm` at the same peak pose; `hip_flexion_r` cross-checked exactly against the JSON's own
value, 0.00e+00 relative diff). The ORIGINAL 25-muscle SO solution's own moment output,
`M_k = sum_m tension_m^SO * r_m^(k)`, is the **moment budget** already self-consistently satisfied by
the model's own dynamics — no re-derivation needed (reserve-actuator check, Sec.2, confirms this is a
good proxy for flexion/adduction; the rotation axis is separately caveated, Sec.5). ANY alternative
recruitment of the SAME 25 muscles reproducing `M_k` for all 3 axes is an equally physically-admissible
way to carry the same joint torque. Given per-muscle bounds `[passive_m, Fmax_m]`, the **global
minimum achievable compressive contact force** is:

```
minimize   sum_m k_m * x_m
subject to sum_m x_m * r_m^(k) = M_k     for k in {flexion, adduction, rotation}
           lo_m <= x_m <= hi_m
```

— an exact LP (`scipy.optimize.linprog`, method=`highs`), the SAME machinery reused for every test
below with different bounds/pins. Two independent adversaries are forced on every result:
1. **Feasibility gating**: SO's own actual tension vector must itself satisfy every bound used
   (Sec.2) — a self-consistency requirement caught and fixed one real issue (Sec.2).
2. **Bang-bang-vs-smooth**: a plain LP's optimum is generically a *vertex* (at most 3 muscles
   interior, the rest pinned at a bound) — not obviously a plausible recruitment. Every material LP
   finding is cross-checked against a **smooth, SO-style quadratic objective**
   (`minimize sum((x_m/Fmax_m)^2)`, the SAME family SO itself uses, scipy SLSQP) on the identical
   moment-budget constraint, to test whether the improvement survives a non-pathological adversary.

## 2. Validity gates and one real fix forced along the way (OODA, not swept under the rug)

| gate | result |
|---|---|
| Reconstitution from JSON's own `k_coef`/`tension_so` reproduces published 386.7679 %BW | PASS (rel diff 7.91e-12) |
| Fresh `hip_flexion_r` moment-arm recompute vs JSON's stored value, same peak pose | PASS (0.00e+00) |
| `opd.passive_active_decomposition`'s nearest-frame tension vs `cmd`'s exact-interp tension | PASS (rel diff 0.0000) |
| SO's own tension vector feasible under the LP's bounds | **FAIL on first pass** → fixed (below) → **PASS** |
| LP-min/LP-max envelope contains SO's actual contact force (both bound variants) | PASS |

**The forced fix**: the first LP run failed its own feasibility gate — `equilibrateMuscles`'s passive-
force estimate for 4 low-tension muscles (`sart_r`, `addlong_r`, `addbrev_r`, `grac_r`) came out
**above** those muscles' own SO-reported total tension (e.g. `sart_r`: passive-est 26.01N > SO-actual
25.13N), which would make SO's own reference solution infeasible under its own LP. This is a symptom
of the **already-disclosed ~7-11% active/passive reconstruction error** at this joint
(`docs/MECHANISM_OVERPREDICTION_DECOMP.md` Sec.6, inherited from `metabolic_cost.py`'s rigid-tendon-
equilibration approximation) — not a new defect, but left uncorrected it would silently distort every
downstream LP. **Fix**: clamp the lower bound to `min(passive_estimate, tension_so)` — the lower bound
must never exclude SO's own reference point, mirroring the same self-consistency clamp already applied
to the upper bound (activation-extrapolated ceiling never allowed to fall below SO's own choice
either). After the fix: **all machine validity gates PASS** (`n_lower_bound_clamped=4/25`,
`n_tight_bound_clamped_to_so_value=0`).

**Moment budget** (self-consistent, from SO's own 25-muscle solution): `hip_flexion_r`=+27.9 N·m,
`hip_adduction_r`=**-69.6 N·m** (dominant — physiologically expected: this is mid/late single-limb
stance, 47.1%GC, where the hip ABDUCTOR moment countering pelvis drop is the largest hip moment
component, a textbook cross-check, not assumed), `hip_rotation_r`=-2.5 N·m (**caveat**: the reserve
actuator itself contributes 3.0 N·m here — 122% of the muscle-only budget — so this axis is noise-
dominated, not physiologically meaningful; a 2-DOF-only robustness re-run of the baseline LP, dropping
rotation, moves LP_min by only 4.6 %BW (287.25 vs 291.86), confirming this caveat does not change the
headline finding).

## 3. TEST (c) — does vector geometry alone force the excess? (baseline, no muscle-specific anomaly assumed)

| | contact force | ratio vs OrthoLoad | Delta toward Ortho | material (>=0.10)? |
|---|---:|---:|---:|---|
| SO actual | 386.77 %BW | 1.412 | — | (anchor) |
| **LP-min (generous Fmax bound)** | **291.86 %BW** | **1.065** | **0.346** | **YES** |
| LP-min (tight, activation-extrapolated bound) | 292.78 %BW | 1.069 | 0.343 | YES |
| LP-max (generous) | 1562.15 %BW | 5.70 | — | (envelope ceiling) |
| **Smooth QP cross-check** (same budget, quadratic objective) | 362.81 %BW | 1.324 | **0.087** | **NO — falls short** |

**Envelope check** (SO's actual contact force lies inside [LP-min, LP-max] for both bound variants):
PASS — confirms the LP is a valid, non-buggy characterization of the SAME feasible set SO's own
solution lives in.

**Forced adversary result**: the LP-min IS real and mathematically decisive — it proves the moment
budget does NOT intrinsically force 386.77 %BW; a moment-preserving reallocation exists at 291.86 %BW
(84.1% of the way from 1.412 to 1.0). But its solution is a **vertex** (2 muscles at their own Fmax
ceiling — `tfl_r`, `glmin2_r` — 3 interior — `iliacus_r`, `glmin1_r`, `glmin3_r` — and **20 of 25 at
floor, including `glmed1_r` and `recfem_r` driven to exactly 0**). A smooth quadratic-objective
reoptimization on the identical constraint (correlates 0.863 with SO's own actual tension pattern —
a fair, non-degenerate proxy) reaches only ratio 1.324 (Delta 0.087) — **short of materiality**. So:
**geometry alone permits a much lower number, but no smooth/plausible alternative objective, absent a
specific muscle-level anomaly, actually reaches it.**

## 4. TEST (a) — diffuse-recruitment (activation-floor) trim

| | contact force | ratio | Delta | gap closed | valid? |
|---|---:|---:|---:|---:|---|
| Naive trim (zero the 15 floor muscles, no rebalance) | 373.96 %BW | 1.365 | 0.047 | 11.4% | **INVALID** (breaks `hip_rotation_r` moment by 54.5%, `hip_flexion_r` by 9.0%) |
| Valid LP rebalance (10 "real" muscles absorb the full budget) | 288.89 %BW | 1.055 | 0.357 | 86.7% | material, but... |
| **Smooth QP cross-check** (same trim, quadratic objective) | 358.77 %BW | 1.310 | **0.102** | 24.8% | **knife-edge** (0.102 vs the 0.10 bar — a 2% margin) |

Threshold-choice robustness: the floor-muscle SET is **identical** at {0.02, 0.05, 0.10} — not
knife-edge on the threshold itself. But the **direct** answer to the task's literal question ("zero
out sub-floor muscles and re-sum, does force drop toward OrthoLoad?") is: only 3.3% of the total
contact force (12.8 of 386.8 %BW) comes from the 15 diffuse muscles directly, and removing them alone
is sub-material (Delta 0.047). The larger "valid LP" improvement (86.7% of the gap) is **not
independent evidence** for a diffuse-recruitment mechanism — inspecting its own solution shows it
re-defunds `glmed1_r`, `glmed2_r`, `psoas_r`, `recfem_r`, `sart_r` to zero among the surviving 10
muscles, i.e. it re-discovers the SAME lever Test (b) finds directly, rather than crediting the
15 floor muscles' removal per se. **Verdict for (a), taken on its own terms: NOT a materially
independent mechanism** — a real, clean negative for that specific sub-question.

## 5. TEST (b) — gluteus-medius (`glmed1_r`) EMG correction — the localized mechanism

This subject's own real surface EMG (reused verbatim, `emg_hybrid_force.sto`, exact interpolation at
t=0.55s): `glmed1_r` = **99.76 N**, vs SO's **585.31 N** (ratio 0.170, matching
`docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md` Sec.9's independently-established ~0.17 finding).

| | contact force | ratio | Delta | gap closed | valid? |
|---|---:|---:|---:|---:|---|
| Naive (glmed1_r only, no rebalance) | 323.81 %BW | 1.182 | 0.230 | 55.8% | **INVALID** (moment hole 20-43% of budget, dominated by `hip_adduction_r` at 43.4%) |
| **LP-optimal rebalance (bang-bang, generous bound)** | **297.58 %BW** | **1.086** | **0.326** | **79.0%** | **VALID, material** |
| LP-optimal rebalance (tight bound) | 298.49 %BW | 1.090 | 0.322 | 78.2% | VALID, material |
| **Smooth QP rebalance (graded, generous bound)** | **327.04 %BW** | **1.194** | **0.218** | **52.9%** | **VALID, material** |
| Combined with floor-trim (a)+(b), LP | 294.60 %BW | 1.075 | 0.336 | 81.7% | VALID, material (only +2.6pp beyond (b) alone — confirms (a)'s redundancy, Sec.4) |

**This is the one mechanism that is simultaneously**: (i) independently motivated (this subject's own
real EMG, not tuned to OrthoLoad, established in a prior session), (ii) moment-balance valid (machine-
verified, `moment_eq_max_abs_resid` ~1e-15 N·m), and (iii) **robust to the smooth-vs-bang-bang
adversary** — both realizations clear the materiality bar by a comfortable margin (0.218 and 0.326,
both well clear of 0.10), unlike test (a) and (c).

**Where does the smooth reallocation put the freed-up moment?** (physiologically coherent, not
arbitrary): `glmin1_r` +205N, `tfl_r` +150N, `glmed2_r` +145N, `glmin2_r` +63N — **the SAME functional
synergist group** (the other hip abductors) picks up gluteus medius's slack — while `recfem_r` -203N,
`psoas_r` -210N, `iliacus_r` -38N relax (the flexion-moment budget doesn't need the extra help the
adduction/abduction budget does). No muscle in this smooth solution hits its own Fmax ceiling
(`glmin1_r` reaches 85% of its own Fmax, `tfl_r` 69%) — a graded, anatomically sensible pattern, unlike
test (c)'s bang-bang vertex.

**Symmetric-QC forced check — does a "nearby," independently-motivated single-muscle refill work?**
`docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md` Sec.7/10 already flagged gluteus MAXIMUS as anomalously silent
(SO floor activation, 0.7% of hip contact force) despite broad literature expecting tone — an obvious
independent candidate to refill the glmed1_r-correction's moment hole. **Forced and it fails,
geometrically, not arbitrarily**: an exact 3-muscle solve (glmax1/2/3_r moment arms are non-singular,
rank 3, condition number 18.3) demands **negative tensions** for 2 of the 3 compartments
(`glmax1_r` → -341N, `glmax3_r` → -477N) — physically impossible. A properly-bounded LP version
(glmax free in `[passive,Fmax]`, everything else pinned) confirms this is **INFEASIBLE**, not merely
numerically awkward. **Why, geometrically**: the moment hole is adduction-dominated (43% of it), and
glmax's compartments have small-to-mixed-sign adduction moment arms (-0.013, -0.005, **+0.062** m) —
NOT the consistently-negative -0.039 to -0.063 m range the true abductor group
(`glmed1/glmed2/tfl/glmin1/2/3_r`) shares — gluteus
maximus is anatomically the wrong muscle group for a frontal-plane-dominated deficit. This infeasible
result is **excluded from the materiality tally** (an infeasible reallocation is not evidence of
anything closing — it is itself part of the diagnosed-gap finding, precisely locating which muscle
group the moment geometry says CAN vs CANNOT plausibly absorb this specific correction).

## 6. Pre-registered falsifier — resolved

**"Trimming diffuse recruitment OR correcting glmed1_r to EMG moves the ratio materially → localizes
a nameable over-count"**: **CONFIRMED specifically for the glmed1_r correction** (robust across both a
bang-bang and a smooth redistribution, Sec.5) — **NOT independently confirmed for the diffuse-
recruitment trim** (Sec.4: its direct contribution is sub-material, and its "valid" improvement is a
re-discovery of the same glmed1_r-type lever, not new evidence). **"If nothing structural moves it →
intrinsic to the method"**: **partially true even after the localized correction** — the fully-
corrected ratio remains 1.086-1.194 (a residual 9-19% over-prediction), and the bare geometric floor
(1.065, Sec.3) is only reachable via an implausible bang-bang pattern, not by any smooth alternative
tested. **Net reading**: this is BOTH a localized-mechanism finding (glute medius over-recruitment, a
specific, nameable, EMG-verified over-count) AND a diagnosed-gap finding (a residual over-prediction
that neither this correction nor the bare vector-geometry floor fully explains) — a joint-
differentiated, two-part answer, not a single up/down verdict, consistent with this doc family's own
established discipline (`docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s knee-vs-hip split).

## 7. Confidence tiers

- **Baseline ratio (1.412x, 386.77 vs 273.93 %BW)**: **in-vivo-anchored** (OrthoLoad, re-verified LIVE
  this session: 162 trials, 18 distinct subjects, median 273.9309 %BW — matches the established cert
  to 4 decimal places).
- **LP/QP feasibility, envelope, and optimality machinery** (Sec.1-3): **method-only** — exact linear/
  quadratic programming and moment-balance identities, a mathematical result, not an empirical claim;
  either it holds to the stated numerical tolerance or it doesn't (it does — all gates PASS).
- **Gluteus-medius over-recruitment magnitude** (Sec.5): **in-vivo-anchored** for the EMG comparison
  itself (this subject's own real surface EMG, reused verbatim) — inherits the already-disclosed
  EMG-to-force normalization uncertainty from `docs/MECHANISM_EMG_DRIVEN.md`/
  `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md` Sec.9 (not re-derived here).
  The specific *reallocation* used to keep the moment budget valid (Sec.5's LP/QP rebalance) is
  **method-only** (a mathematically valid redistribution, not itself independently measured to be
  what real muscles do) — disclosed, not conflated with the EMG-anchored part.
  The gluteus-maximus-infeasibility explanation is **method-only** (exact geometry) reinforced by the
  independent literature/anomaly flag already established (`docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`).
- **Residual 9-19% over-prediction and the "why does SO's actual solve choose the expensive pattern"
  question**: **diagnosed-gap** — this session localizes and partially closes one mechanism but does
  not fully resolve the remainder; per the task's own framing, closing it further would need either
  an EMG-constrained SO re-solve (a genuinely different redundancy-resolution method) or a
  fundamentally different contact model (deformable/COMAK), neither attempted this session.

## 8. Honest gaps

1. **Single trial, single subject, single instant** (subject2/`walking1`, t=0.55s) — same scope
   caveat as every cert in this family; no claim of generality across subjects/trials/gait phases.
2. **`hip_rotation_r`'s moment budget is reserve-actuator-dominated (122%)**, i.e. near-noise on that
   one axis — a 2-DOF-only robustness check shows this doesn't change the headline LP-min by more than
   4.6 %BW, but a full re-run of every test under a 2-DOF-only budget was not performed (scope limit).
3. **The smooth-QP cross-check is a hip-DECOUPLED approximation** of SO's true whole-body problem — it
   ignores the knee's and other joints' simultaneous moment constraints, which also bind several
   biarticular muscles tested here (`recfem_r`, `sart_r`, the long hamstrings). Its 0.863 correlation
   with SO's actual solution is reassuring but not proof the decoupling is harmless.
4. **The "tight" (activation-extrapolated) Fmax bound is a first-order, linear-in-activation
   approximation** (`active_now/activation + passive`), clamped to never exclude SO's own choice — a
   disclosed approximation, not a re-derivation of the true force-length-velocity-scaled instantaneous
   capacity (which would require re-querying the model's active/velocity multipliers directly).
5. **`recfem_r`/`psoas_r` are geometrically flagged as "expensive" by every LP variant** (consistently
   driven toward zero alongside `glmed1_r`) but have **no real-EMG channel** in this subject's dataset
   — unlike `glmed1_r`, their over-recruitment is a plausible, recurring geometric signal, not an
   independently-verified finding this session. (`recfem_r`'s analogous over-recruitment, via a
   *different* independent channel — vasti's real EMG implying an under-used knee-extensor role
   `recfem_r` may be substituting for — is already flagged in `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`
   Sec.7/9; not re-derived here, cross-referenced only.)
6. **The literature-motivated glmax refill's infeasibility used a restrictive scenario** (glmax free,
   every OTHER non-glmed1_r muscle pinned exactly at its SO baseline) — a more generous scenario
   (glmax + a FEW other abductors free) was not swept; the LP/QP tests in Sec.3/5 already show such
   broader reallocations ARE feasible, so this is a narrow, deliberately strict single-muscle-group
   test, not evidence that no small-group fix exists.
7. **The "smooth" objective is one specific, reasonable choice** (`sum((x/Fmax)^2)`, matching SO's own
   convention) among a broader space of plausible alternative cost functions — not swept further.

## 9. Files

- `scripts/msk/hip_structural_mechanism.py` — the full, re-runnable pipeline (LP + QP cross-checks,
  all validity gates, all three named sub-tests + combined). Real run, exit 0, all machine gates PASS.
  Reuses `validate_joint_force.py`, `static_opt_knee.py`, `validate_hip_force.py`,
  `contact_muscle_decomp.py`, `overprediction_decomp.py` **UNCHANGED** (md5-recorded at run time:
  `e2bb15958cfe8d72f44a3e0ac3d95950` / `d00f962c1a8873606874fffd65e828c8` /
  `a6514b9f081e0b4767fb4b0f94e60200` / `454393a1101b7134a9949a5c7b6d2601` /
  `c43b70d76de2b059216645130f082bc9` — the first three match `docs/MECHANISM_CROSS_SUBJECT.md`'s /
  `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`'s own recorded hashes exactly).
- `data/msk_smoketest/subject2_walking1/hip_structural_mechanism/hip_structural_mechanism_results.json`
  — every number in this document: moment budgets, per-muscle moment arms (all 3 hip DOF), Fmax,
  passive/active split, LP/QP solutions and solution structures (which muscles at ceiling/floor/
  interior) for every test, envelope checks, feasibility flags, condition-number diagnostic.
- Reused, not recomputed: `data/msk_smoketest/subject2_walking1/contact_muscle_decomp/
  contact_muscle_decomp_results.json` (per-muscle tension/activation/k_coef/moment_arm table),
  `data/msk_smoketest/subject2_walking1/emg_driven/emg_hybrid_force.sto` (real-EMG-hybrid glmed1_r
  tension), `data/msk_smoketest/subject2_walking1/static_optimization/so/*.sto` (SO force/activation).
- Prior docs this builds on directly: `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`,
  `docs/MECHANISM_OVERPREDICTION_DECOMP.md`, `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`,
  `docs/MECHANISM_METABOLIC_CALORIMETRY.md`.

Isolation respected throughout: bodytwin only; LabValidation model/session data and OrthoLoad data
read in place on the read-only external drive; `scipy` was added to this repo's own isolated
`.venv-msk` (via `uv pip install scipy`, additive, non-destructive) since it was absent from that
venv's site-packages despite being present under a different Python version in the user environment.
No git commit, no git add, no git push performed. All new files left untracked for the coordinator.
