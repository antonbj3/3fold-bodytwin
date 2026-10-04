# MECHANISM WOBBLING-MASS -> KNEE/HIP JOINT CONTACT FORCE — does soft tissue explain part of the 1.5x gap? (2026-07-21)

Executes the operator's ask: couple the EXISTING segment-level wobbling-mass model
(`docs/MECHANISM_WOBBLING_MASS.md`, `scripts/msk/wobbling_mass.py`) into an estimate of its effect on
this twin's own PEAK knee-r/hip-r **JOINT CONTACT FORCE** — the muscle-force-inclusive,
OrthoLoad-anchored quantity `static_opt_knee.py`/`validate_hip_force.py`/`cmc_second_solve.py`
already certify at 391.11%BW (knee-r, SO) / 387.04%BW (hip-r, SO), 1.41-1.52x the OrthoLoad in-vivo
anchor — to test whether rigid-body-assumption error from unmodeled soft-tissue dynamics can explain
part of that over-prediction. Script: `scripts/msk/wobbling_joint_force_coupling.py` (pure
JSON-to-JSON propagation/arithmetic on already-committed cert outputs — no new OpenSim solve, no
re-optimization). Evidence: `data/msk_smoketest/wobbling_joint_force/wobbling_joint_force_results.json`.

**Headline: CLEAN, FORCED NEGATIVE for walking. Propagating the wobbling-mass segment-force
correction into the knee/hip joint through the standard recursive Newton-Euler free-body-cut
relation (the SAME relation this repo's own `R_old`/`pure_reaction` terms already use) gives a
peak-contact-force change of at most +0.021% (knee) / +0.039% (hip) for walking1 — three orders of
magnitude under the pre-registered 2-3% "meaningful contributor" floor, even using the MOST
GENEROUS combination of tier/summary-statistic this repo's own 100-point generic parameter sweep
allows. Wobbling mass is RULED OUT as a meaningful contributor to walking1's ~1.4-1.6x
over-prediction of OrthoLoad. For DJ1 (real drop-landing): the underlying SEGMENT-level mechanism
is confirmed present and 2.6-4.4x LARGER than in walking1 (matches the literature's small-for-gait/
large-for-impact task-dependence), but propagating it to an actual knee/hip JOINT CONTACT FORCE
number is NOT computable — no Static-Optimization or CMC muscle-force solve exists for DJ1 anywhere
in this repo, so there is no contact-force denominator. This is reported as a precise, disclosed
gap (a PASS under this project's own honest-negative convention), not papered over.**

## 0. Pre-registered claim, adversary, anchor (stated before computing the table below)

**C (walking): wobbling mass reduces peak knee-r/hip-r contact force by <2-3% -> NOT a meaningful
contributor to the 1.5x gap.** ¬C: reduces it by >=2-3% -> a real, if partial, contributor worth
pursuing further. **Adversary (this claim leans negative, so the adversary to force is "the
estimate is unfairly small"):** force the propagation to its most generous form — (1) the
swept-envelope's OWN most-favorable (highest-attenuation) combo per segment, not just one generic
representative point; (2) allow PRIMARY (100Hz IK) and SUPPLEMENTARY (2000Hz real-measured-driver)
tiers, and both segments' own best case, to combine even when that requires DIFFERENT (mu,fn,zeta)
per segment/tier (a strictly more generous stack than requiring one shared parameter point); (3)
ignore the timing mismatch between the wobbling-model's own analysis window and the contact-force
peak instant (ignoring it can only inflate, never deflate, the apparent overlap). **Anchor:** this
repo's own already-certified, OrthoLoad-anchored peak knee-r/hip-r contact force
(`static_opt_knee_results.json`/`hip_force_validation_results.json`'s `jr_summary.peak_pct_bw`, the
SAME numbers the capstone's 1.4-1.6x headline uses) — a real, external, in-vivo-anchored quantity,
never a tautology gate.

## 1. Geometry (derived, not assumed) — how a segment force delta becomes a joint force delta

For a serial open kinematic chain foot-shank-thigh-pelvis, the classical recursive Newton-Euler
DISTAL-TO-PROXIMAL force relation for segment *i* is

```
R_proximal_i(t) = M_i * a_i(t) + R_distal_i(t) - M_i*g - F_muscles_on_i(t)
```

This repo's OWN `static_opt_knee.py`/`validate_hip_force.py` compute exactly this quantity
(`R_old_vec = acc_sum - chain_mass*G_VEC - F_R`, confirmed by direct source read, not assumed) for
the free-body cut at the knee (distal chain = descendants of `walker_knee_r` = shank_r+foot_r) and
at the hip (descendants of `hip_r` = thigh_r+shank_r+foot_r) — this repo's own `R_old`/`pure
reaction` term **IS** `R_proximal_i`, additive across the chain by construction: `R_old(hip) -
R_old(knee) = ` the thigh's own `M*a-g` contribution, exactly (the same recursive-sum structure).

The wobbling-mass model changes ONLY each segment's own `M_i*a_i(t)` term — replacing the rigid
value `F_rigid,i(t) = M_i*a_i(t)` with `F_wobble,i(t) = M_i*a_i(t) + m_w,i*z_i''(t)` (the exact
identity `wobbling_mass.py` already proves to `1.27e-16`). **Holding muscle forces
`F_muscles_on_i(t)` fixed** (a disclosed, first-order, NOT-re-optimized perturbation — see §6 gap
#1), this gives a forced, non-hand-waved propagation:

```
Delta_F_knee_est(t) = F_wobble_shank(t) - F_rigid_shank(t)
Delta_F_hip_est(t)  = Delta_F_knee_est(t) + [F_wobble_thigh(t) - F_rigid_thigh(t)]
```

evaluated at each series' own peak (exactly how `wobbling_mass.py`'s own `attenuation_pct` is
already defined), then expressed as a % of this twin's already-certified peak knee-r/hip-r JOINT
CONTACT FORCE.

**Machine check before trusting this (not assumed):** the propagation above uses each segment's
`peak_rigid_N_representative` together with the SWEPT `attenuation_pct` min/median/max to bound
`Delta_F` across the full 100-combo (mu,fn,zeta) grid, without needing the full per-combo Newton
arrays re-saved. This is licensed by `F_rigid` being INDEPENDENT of (mu,fn,zeta) by construction
(`F_rigid = M_total*a(t)`) — re-verified live this session on a fresh synthetic input (27 combos,
`mu`x`fn`x`zeta`): **max abs diff = 0.0 N, bit-identical, PASS.**

## 2. Denominators — live-loaded, machine-cross-checked (not from doc prose)

| quantity | value | source |
|---|---:|---|
| body weight | 766.88003 N (78.2 kg) | cross-checked identical across all 4 input JSONs |
| knee-r peak contact force, SO (official JointReaction) | 391.11 %BW = 2999.4 N @ t=0.51s | `static_opt_knee_results.json` `jr_summary` |
| hip-r peak contact force, SO (official JointReaction) | 387.04 %BW = 2968.1 N @ t=0.55s | `hip_force_validation_results.json` `jr_summary` |
| knee-r peak contact force, CMC (2nd, decorrelated solve) | 427.54 %BW = 3278.7 N | `cmc_second_solve_results.json` |
| hip-r peak contact force, CMC | 429.26 %BW = 3291.9 N | `cmc_second_solve_results.json` |
| OrthoLoad in-vivo anchor, knee median (n=72) | 258.22 %BW | live-loaded |
| OrthoLoad in-vivo anchor, hip median (n=162) | 273.93 %BW | live-loaded |
| SO/OrthoLoad ratio | knee 1.515x, hip 1.413x | re-derived live |
| CMC/OrthoLoad ratio | knee 1.656x, hip 1.567x | re-derived live |

## 3. Walking1 propagation — every tier x summary-statistic combo (the forced adversary's menu)

| tier | summary stat | Delta_F_shank (N) | Delta_F_thigh (N) | Delta_F_knee_est (N) | % of knee contact force | Delta_F_hip_est (N) | % of hip contact force |
|---|---|---:|---:|---:|---:|---:|---:|
| primary (100Hz IK) | representative (mu=.7,fn=15Hz,zeta=.5) | +0.330 | −0.563 | +0.330 | **+0.0110%** | −0.233 | **−0.0079%** |
| primary (100Hz IK) | sweep best-favorable | +0.636 | +0.520 | +0.636 | **+0.0212%** | +1.156 | **+0.0390%** |
| supplementary (2000Hz measured) | representative | −0.259 | −0.656 | −0.259 | **−0.0086%** | −0.916 | **−0.0309%** |
| supplementary (2000Hz measured) | sweep best-favorable | +0.305 | +0.771 | +0.305 | **+0.0102%** | +1.076 | **+0.0362%** |

**Adversary's strongest shot (max \|%\| across all 4 rows): knee = 0.0212%, hip = 0.0390%.**
**KEY-QUESTION GATE (<2%): PASS — clean negative. Loose-band gate (<3%): PASS.**

Both segments' own PEAK acceleration in walking1 is genuinely gentle (thigh_r 1.61 m/s², shank_r
2.42 m/s² — "a mid-stance-trough-to-loading-response transition... not a sharp heel-strike,"
per `MECHANISM_WOBBLING_MASS.md` §5), so `M*a` for a 3.8-9.7 kg segment tops out at **9-16 Newtons
total** — three orders of magnitude below the ~3000 N muscle-force-dominated contact force at the
same joint. No swept (mu,fn,zeta) combo, no tier, can turn a ~10N correction on a ~10N base into a
multi-percent change of a ~3000N total.

**Timing caveat (disclosed, not hidden):** the knee-r contact-force peak instant (t=0.51s) falls
BEFORE the wobbling-model's own analysis window even starts ([0.5505s, 0.7005s] — different events:
SO's knee peak is mid-stance load acceptance, the wobbling model's own onset detector fires at the
loading-response transition). The hip-r peak (t=0.55s) sits right at the window's own edge. The
propagation above therefore compares each series' OWN peak (same convention `attenuation_pct`
already uses) rather than a simultaneous-instant force balance — ignoring this mismatch can only
inflate, not deflate, the reported %, so it does not threaten the negative conclusion.

## 4. Structural/minority-share cross-check (task-independent, a SECOND, complementary reason)

Wobbling mass, by construction (Gruber 1998 / Pain & Challis's own model class), only ever acts on
the classical rigid-multibody "pure reaction" term (`R_old`/`pure_reaction` above) — never on the
muscle-crossing-force term, which is a separate, independently-solved quantity (SO/CMC) held fixed
in this propagation.

| joint | R_old / pure-reaction (%BW) | muscle-crossing force (%BW) | total contact force (%BW) | R's share of total |
|---|---:|---:|---:|---:|
| knee-r | 102.17 | 306.64 | 391.11 | **26.1%** |
| hip-r | 88.59 | 303.48 | 387.04 | **22.9%** |

This gives a **geometric ceiling independent of spring/damper tuning**: even in the (unrealistic
for gait — see §5) limiting case of wobbling mass fully zeroing out the ENTIRE R-component, the
total contact force could fall by at most ~26%/~23% — and this repo's own measured walking1
kinematics (§3) only ever move a few percent of even THAT smaller term, landing far inside this
ceiling, not near it. This is offered as a structural, not a numerically-independent, cross-check:
it explains *why* the small measured numbers in §3 are physically sensible, it does not re-derive
them from a different data source.

**A specific adversary move explicitly considered and rejected**: could using Pain & Challis
(2006)'s own reported "up to nearly 50% lower" (applied hypothetically to the R-component instead
of this repo's own small swept-envelope max) give a bigger, still-defensible walking estimate
(~13%/~11% of total, i.e. above the 2-3% floor)? **Rejected as an illegitimate adversary**: that 50%
figure was measured for a 43cm DROP LANDING, a high-impact task; importing a landing-scale
attenuation fraction into a gait analysis would contradict the very task-dependence
(small-for-gait/large-for-impact) this document's own DJ1-vs-walking1 contrast (§5) confirms.
The correct forced adversary for THIS task is this repo's own generic (mu,fn,zeta) sweep applied to
walking1's own real (gentle) kinematics — which is what §3 already reports.

## 5. DJ1 (real drop-landing) — segment mechanism confirmed, joint-force coupling an honest gap

**Filesystem+JSON-structural search (fixed once, see below) for any Static-Optimization / CMC /
JointReaction result indexing "DJ1" as an analyzed activity anywhere in `data/msk_smoketest/`:
NONE FOUND.** `docs/MECHANISM_MUSCLE_DAMAGE.md` §8 independently found and disclosed the identical
gap for a different question ("No SO exists for DJ1 in this repo; none was solved here").
Consequently: **part (b) of this task (peak knee/hip CONTACT FORCE reduction for landing) is NOT
computable in this repo as of this session — there is no muscle-resolved contact-force denominator
for DJ1 to propagate into.** Reported as a precise, disclosed gap (a PASS under this project's own
honest-negative convention), not papered over.

**A bug in the first version of this exact search, found and fixed before trusting the "no cert
exists" conclusion**: a raw-text co-occurrence check (`"DJ1" in file_text and "peak_pct_bw" in
file_text`) FALSE-POSITIVED on `cross_activity_validation_results.json` — that file mentions "DJ1"
only as a filename string inside an unrelated trial-availability list
(`kinematics_gap.labvalidation_all_distinct_trial_stems`), while separately containing unrelated
`peak_pct_bw`-keyed force numbers for squat/sit-to-stand elsewhere in the SAME file. **Forced fix**:
a structural JSON-key search (does "DJ1" appear as an actual dict KEY, i.e. an indexed
activity/trial, anywhere?) — re-run, **0 hits**, confirming the negative is real, not a search-scope
artifact. Exactly the class of bug this project's own memory already flags ("whole-module keyword
grep is a scope mismatch for a claim that needs the two things to be structurally related, not
just textually present").

**What IS computable — the segment-level mechanism itself** (from `wobbling_mass.py`, unchanged,
reused directly):

| trial | tier | segment | attenuation % (representative) | validity |
|---|---|---|---:|---|
| walking1 | primary 100Hz IK | shank_r | +3.55% | validity gate PASS |
| walking1 | supplementary 2000Hz | shank_r/thigh_r | −4.50% | not subject to that gate (real-measured driver) |
| DJ1 | primary 100Hz IK | shank_r | +9.37% | validity gate **FAILS** (14.6-18.7%BW residual vs 8%BW ceiling — untrusted, disclosed already in `MECHANISM_WOBBLING_MASS.md`) |
| DJ1 | supplementary 2000Hz | shank_r/thigh_r | **+19.79%** | the most trustworthy real DJ1 number available (real, 2000Hz, no IK differentiation) |

**Task-dependence contrast (pre-registered: DJ1 should exceed walking1, matching Pain & Challis'
small-for-gait/large-for-impact pattern): |attenuation%| ratio DJ1/walking1, shank_r, representative
point — supplementary tier 4.40x, primary tier 2.64x. Both exceed 1.0x: PASS.** This is the walking-
vs-landing CONTRAST the symmetric-QC instruction asked for, and it lands in the literature-predicted
direction even though the absolute joint-force number for DJ1 cannot be produced.

**Illustrative-only scale reference (explicitly NOT a claim, NOT gated, NOT used in any verdict)**:
DJ1's supplementary-tier absolute correction (`Delta_F_knee_est`=25.2N, `Delta_F_hip_est`=88.9N) is
0.84%/3.00% of WALKING1's OWN peak contact force if (invalidly) used as a size reference. DJ1's REAL
contact-force denominator — which does not exist as a cert in this repo — would almost certainly be
far larger than walking1's, since DJ1's own measured peak GRF (465.4%BW) is 3.76x walking1's own
(123.9%BW); a proportionally scaled denominator would make this illustrative % smaller, not larger.
This number exists only to give a size sense, not to substitute for the missing DJ1-specific
estimate.

## 6. Honest gaps (disclosed, not hidden)

1. **First-order, muscle-forces-held-fixed propagation only.** This estimate captures the direct
   segment-translational-inertia pass-through into the joint's Newton-Euler reaction term; it
   structurally CANNOT capture a second, different pathway — a wobbling-driven change to the net
   joint MOMENT that a short-lever-arm muscle could in principle amplify via a full SO/CMC
   re-optimization on a wobbling-augmented model. Closing that gap needs a full ID+SO/CMC re-solve
   on a modified (wobbling-DOF-added) model — out of scope for this first step, consistent with this
   repo's "reuse, don't re-solve" convention; NOT folded into the negative verdict, stated as an open
   question.
2. **Vertical-axis-only numerator vs 3D-resultant denominator.** The wobbling model's segment force
   is a single (vertical) axis quantity; the joint contact force it is compared against is a full 3D
   vector-norm resultant (`np.linalg.norm(F, axis=1)`, confirmed by source read). Treated here as a
   scalar-vs-scalar order-of-magnitude ratio (standard practice for this kind of bound); given the
   gap between the two numbers is ~2 orders of magnitude, this scope mismatch cannot change the
   verdict.
3. **Right-leg thigh/shank only** — same scope limit `wobbling_mass.py` already discloses; foot_r's
   own (unmodeled) wobbling contribution is not included (expected small: foot mass and soft-tissue
   fraction are both smaller than thigh/shank).
4. **DJ1's PRIMARY (100Hz IK) tier fails its own whole-body Newton-residual validity gate** (already
   disclosed in `MECHANISM_WOBBLING_MASS.md`) — only the SUPPLEMENTARY (2000Hz real-measured-driver)
   tier is trusted for DJ1 in this document, consistent with that document's own verdict logic.
5. **Generic, swept (not literature-pinned) mu/f_n/zeta** — inherited limitation from
   `wobbling_mass.py` (its own §2, gap #2); this document's propagation inherits the same envelope,
   not a new one.
6. **Single trial each, right-leg-primary, subject2 only** — no claim of generality across
   subjects/trials/speeds/drop-heights.
7. **The "R_old share of total" structural check (§4) is a plausibility cross-check, not an
   independent data source** — it uses the SAME SO solve's own internal decomposition, not a
   second measurement.

## 7. Files

- `scripts/msk/wobbling_joint_force_coupling.py` — the full propagation + machine checks (grid-
  invariance re-verification, structural DJ1-cert-existence search) + evidence-JSON writer. Run:
  `.venv-msk/bin/python3 scripts/msk/wobbling_joint_force_coupling.py` (~2s wall time, no OpenSim
  solve — pure JSON arithmetic; imports `opensim` only transitively via `wobbling_mass.py` for the
  one grid-invariance self-check).
- `data/msk_smoketest/wobbling_joint_force/wobbling_joint_force_results.json` — every number above,
  traceable back to this one file.
- Read in place, never modified: `data/msk_smoketest/wobbling_mass/wobbling_mass_results.json`,
  `data/msk_smoketest/subject2_walking1/static_optimization/static_opt_knee_results.json`,
  `data/msk_smoketest/subject2_walking1/hip_force_validation/hip_force_validation_results.json`,
  `data/msk_smoketest/subject2_walking1/cmc_second_solve/cmc_second_solve_results.json`.
- Prior work built on: `docs/MECHANISM_WOBBLING_MASS.md` (the segment-level model, unchanged),
  `docs/MECHANISM_CMC_SECOND_SOLVE.md` (the 1.4-1.6x-over-OrthoLoad headline + its CMC
  cross-check), `docs/MECHANISM_MUSCLE_DAMAGE.md` §8 (independently found the same DJ1-has-no-SO gap
  for a different question).

ISOLATION (bodytwin `COORDINATOR.md` Sec.1): bodytwin only; all cert JSONs above read in place, never
mutated; no git commit/push; both files listed above are new.
