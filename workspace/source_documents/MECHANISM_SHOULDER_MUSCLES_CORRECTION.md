# MECHANISM SHOULDER-MUSCLES CORRECTION — closing out the flagged follow-up (2026-07-21)

Closes out `docs/MECHANISM_SIGN_BUG_REMEDIATION.md` Sec.7's most significant flagged-not-fixed
finding: `scripts/msk/validate_shoulder_force_with_muscles.py` defines its OWN, independently-typed
copy of the crossing-muscle sign bug (`arm_crossing_muscles_and_forces`, originally at lines
244-279) —
it never calls the shared, already-fixed `static_opt_knee.knee_crossing_muscles_and_forces`, so
that source fix never reached it. This document: applies the identical one-line fix, verifies it
via a known-answer toy test (mirroring `scripts/msk/audit_sign_bug.py`'s style) plus a falsifier
control confirming the toy test actually discriminates bug from fix, re-derives the corrected
shoulder-with-muscles contact force end to end, adds a new `opensim.JointReaction` cross-check
(the bug-immune referee the original remediation used at knee/hip/ankle, absent from this file
until now), and banners `docs/MECHANISM_ARM_MUSCLES.md`. Every number below is machine-measured
this session. Isolation respected: bodytwin only, no git commit/push (coordinator commits).

## 1. The bug and the fix

`scripts/msk/validate_shoulder_force_with_muscles.py`, `arm_crossing_muscles_and_forces`
(now ~line 271, shifted by an explanatory comment block):

```python
# before (bug): backwards in BOTH branches -- identical pattern to the already-fixed
# static_opt_knee.knee_crossing_muscles_and_forces
inside_idx, outside_idx = (k, k + 1) if in_k1 else (k + 1, k)

# after (fix): when in_k1 is True, body k+1 -- not k -- is the inside point
inside_idx, outside_idx = (k + 1, k) if in_k1 else (k, k + 1)
```

Same geometric derivation as the shared-function fix: `in_k1 = bodies[k+1] in inside_set`; the
branch only executes on a crossing (`in_k != in_k1`). If `in_k1` is `True`, body `k+1` is the one
INSIDE the free-body cut, so `inside_idx` must be `k+1` — the pre-fix code assigned `k` there,
backwards. This function is a genuinely separate code path (confirmed via
`docs/MECHANISM_SIGN_BUG_REMEDIATION.md` Sec.4's own repo-wide grep): it was written standalone
because this graft's donor-1 muscles carry no `_r`/`_l` name suffix, so the shared function's
`name.endswith(f"_{side}")` filter would have silently matched zero of them.

## 2. Toy-test verdict — known-answer unit test, mirroring `audit_sign_bug.py`'s style

New script: `scripts/msk/audit_shoulder_muscles_sign_bug.py`. Adapted to this function's own
4-argument signature (`model, state, muscle_names, inside_bodies` — no `side` filter, unlike the
knee's 5-arg form). Reuses `audit_sign_bug.py`'s duck-typed Fake\* OpenSim stand-ins (read-only
import, no real OpenSim model needed), plus a SECOND, independently-typed local reimplementation
of the fix (case E) as a redundant cross-check, exactly mirroring that script's own 5-case
structure (A/B: path-order diversity; C: crossing mid-path; D: arbitrary 3-D axis; E: independent
second implementation).

Ground truth (identical mechanics to the knee audit): a cable/muscle under tension pulls each
path endpoint TOWARD the other endpoint.

| toy case (post-fix, real function) | verdict | matches +truth | matches −truth |
|---|---|:---:|:---:|
| A: inside point first in path order | **CORRECT** | True | False |
| B: inside point second in path order | **CORRECT** | True | False |
| C: crossing mid-path (3-point path) | **CORRECT** | True | False |
| D: arbitrary 3-D axis orientation | **CORRECT** | True | False |
| E: independent second implementation (sanity) | **CORRECT** | True | False |

**Confirmed exactly as instructed: `np.allclose(computed, +ground_truth)` is True and
`np.allclose(computed, -ground_truth)` is False, for all 5 cases** (JSON:
`data/msk_smoketest/subject2_arm_abduction/shoulder_muscles_sign_bug_audit/shoulder_muscles_sign_bug_toy_test_results.json`).

**Falsifier control (the toy test is not vacuous — it actually detects the bug when present):**
re-ran case A through a standalone copy of the PRE-FIX ternary (not written into any repo file,
a throwaway check): `computed = [0, 1, 0]` vs `truth_A = [0, -1, 0]` → `matches +truth = False`,
`matches -truth = True` — the exact negative, confirming the test discriminates bug from fix
rather than passing regardless of implementation.

## 3. Corrected numbers — full pipeline re-run + a new JointReaction referee

Re-ran `validate_shoulder_force_with_muscles.py` end to end post-fix (rebuilds the locked probe
model, re-runs Static Optimization via Ipopt — 261 frames, 0 NaN, no saturation, reserve usage
≤1.15% of optimal force — then the fixed crossing-detection + subtraction). Added a NEW STEP 6:
an official `opensim.JointReaction` cross-check (this file had none before; the knee/hip/ankle
family's own bug-immune referee, structurally independent — it is fed the SO force.sto directly
via its documented `forces_file` mechanism and never calls `arm_crossing_muscles_and_forces` at
all). Used the LOCKED probe model (not the unlocked graft) for dynamical consistency with the
forces it is fed; targeted `joint_names=acromial_r`, `apply_on_bodies=child`,
`express_in_frame=child` (matching the existing knee/hip JR convention). Output columns
(measured, not guessed): `acromial_r_on_humerus_r_in_humerus_r_f{x,y,z}`.

| | %BW | ratio vs. 72.2 %BW anchor |
|---|---:|---:|
| Pre-graft (0 muscles cross the joint) | 5.306 | 0.0735 |
| **Published (bug artifact)** | **24.177** (t=0.68s) | **0.335** |
| **CORRECTED (self-computed, fixed)** | **19.065** (t=1.700s) | **0.264** |
| **JointReaction (bug-immune referee)** | **19.362** (t=1.700s) | 0.268 |
| corrected-vs-JR relative diff (peak-to-peak) | **1.54%** | — |

**Both peaks land at the IDENTICAL instant, t=1.700s** — not just similar magnitudes at
different times, a stronger agreement signal than magnitude-alone. `n_crossing_muscles=22`,
identical SET (not just count) before and after the fix (`BIC_brevis, BIC_long,
Coracobrachialis, DeltoideusClavicle_A, DeltoideusScapula_M/P, Infraspinatus_I/S,
LatissimusDorsi_I/M/S, PectoralisMajorClavicle_S, PectoralisMajorThorax_I/M, Subscapularis_I/M/S,
Supraspinatus_A/P, TRIlong, TeresMajor, TeresMinor`) — expected, since the fix only changes the
DIRECTION assigned to a detected crossing, never whether one is detected.

**Machine cross-check that the fix touched only the sign-sensitive term** (mirrors the
knee/ankle remediation's own norm-blindness check): every bug-BLIND quantity is bit-identical
before/after, to the full float64 precision on disk —

| quantity | old (buggy run) | new (fixed run) | bit-identical |
|---|---:|---:|:---:|
| `r_old_peak_pctBW` | 5.305774172375711 | 5.305774172375711 | **True** |
| `R_old_mag_pctBW[0]` | 4.917902226001441 | 4.917902226001441 | **True** |
| `R_old_mag_pctBW[-1]` | 4.917902342468531 | 4.917902342468531 | **True** |
| `muscle_force_contrib_mag_pctBW[0]` | 16.976493415726505 | 16.976493415726505 | **True** |
| `muscle_force_contrib_mag_pctBW[-1]` | 13.640497168368695 | 13.640497168368695 | **True** |

— exactly as the mechanism predicts: `R_old` is computed independently of the crossing-fn, and
`np.linalg.norm(v) == np.linalg.norm(-v)`, so any pure-norm-of-crossing-only quantity is
structurally blind to a uniform sign flip. Only the final `F_new = R_old − Σmuscle-crossing`
subtraction (bug-sensitive: `R − M ≠ R − (−M)` in general) moved:

| | old (buggy) | new (fixed) | Δ |
|---|---:|---:|---:|
| `f_new_peak_pctBW` | 24.177169 | 19.064522 | −5.112647 |
| `f_new_at_old_peak_time_pctBW` | 19.961324 | 12.715879 | −7.245445 |

## 4. Why the correction moved the number DOWN here, but UP at the knee/hip — measured, not narrated

`docs/MECHANISM_SIGN_BUG_REMEDIATION.md` Sec.8 explains the knee/hip corrections (233→391,
235→387 %BW, both INCREASES) via antagonist co-contraction: the muscle-crossing pull `M` is
nearly ANTI-aligned with the kinematic reaction `R_old` along the joint's compressive axis, so
correctly computing `R − M` (fixed) is LARGER in magnitude than the erroneous `R + M` (buggy).
The shoulder does the opposite (24.177 → 19.065, a DECREASE). Rather than assume the same
mechanism runs in reverse, the actual angle between `R_old` and `M` was measured directly
(reusing the cached SO output, no new Ipopt solve) at the two frames of interest:

| frame | \|R_old\| (N) | \|M\| (N) | angle(R_old, M) | \|R−M\| (fixed) | \|R+M\| (buggy) |
|---|---:|---:|---:|---:|---:|
| F_new corrected's own peak (t=1.700s) | 35.480 | 154.751 | **69.6°** | 146.202 | 170.406 |
| R_old's own peak (t=1.110s) | 40.689 | 121.720 | **45.3°** | 97.516 | 153.079 |

Both angles are ACUTE (<90°) — i.e. `R_old` and `M` are geometrically more ALIGNED
(constructive) than anti-aligned at this joint/trial, the opposite regime from the knee/hip's
antagonist-compression pattern. Constructive alignment means erroneously ADDING `M` (buggy)
produces a LARGER resultant than correctly SUBTRACTING it (fixed) — exactly the observed
direction. This is plausible on anatomical grounds (a single-DOF abduction sweep with no
antagonist-heavy co-contraction demand, unlike a weight-bearing stance-leg gait cycle) but is
reported here as a MEASURED geometric fact about this specific trial, not asserted as a general
property of the shoulder.

## 5. Symmetric verdict: does the correction move closer to or further from the 72.2 %BW anchor?

**Measured, not assumed: FURTHER.** Ratio moves from 0.335 (buggy) to 0.264 (corrected) — a
−21.1% relative change in the predicted force, and the absolute gap to the anchor widens from
48.02 to 53.14 %BW-points (+5.11 pct-pts further). This is the OPPOSITE direction from the
knee/hip corrections in `docs/MECHANISM_SIGN_BUG_REMEDIATION.md` (which moved past their
anchors, ending up over-predicting in-vivo by ~1.4-1.5×). No hypothesis is offered here that the
correction "should" move the number toward or away from the anchor — Sec.4 above shows the
direction is set by the local `R_old`-vs-`M` geometry of this specific trial/joint, which the
sign bug does not know about and does not consistently push one way.

**This does not change the qualitative graft finding.** Grafting muscles still moves the
prediction substantially up from the zero-muscle baseline:

| | %BW | × pre-graft | % of anchor-minus-pre-graft gap closed |
|---|---:|---:|---:|
| Pre-graft (0 muscles) | 5.306 | 1.00× | 0% |
| Published (buggy) | 24.177 | 4.56× | 28.2% |
| **Corrected** | **19.065** | **3.59×** | **20.6%** |

Zero muscles crossing the joint made Static Optimization structurally impossible (a vacuous,
0-redundancy tautology, per `docs/MECHANISM_SHOULDER_FORCE.md`'s own model audit); grafting
muscles made it solvable and moved the prediction up by a genuine 3.59× — smaller than the
4.56× previously claimed, but the same qualitative story. Only the specific with-muscle
magnitude and its distance from the in-vivo anchor change.

## 6. Every caller — checked, not assumed

`grep -rn "arm_crossing_muscles_and_forces" scripts/msk/*.py`: defined and called only within
`validate_shoulder_force_with_muscles.py` itself (one definition, one call site, inside
`recompute_full_pipeline`'s per-frame loop). No other script in this repo calls this function —
confirmed via repo-wide grep, not assumed from the function's local-sounding name. This closes
the single remaining open item from `docs/MECHANISM_SIGN_BUG_REMEDIATION.md` Sec.7/10 that named
this file specifically; the OTHER two lower-priority items that section flagged
(`validate_elbow_force.py`'s illustrative `arm26` reference comparison, `docs/MECHANISM_ELBOW_FORCE.md`
line 27; `force_scenes_batch.py` + its stale doc) remain out of this task's explicit scope
(fix `validate_shoulder_force_with_muscles.py` + banner `MECHANISM_ARM_MUSCLES.md`) and are
**not** touched here.

## 7. Honest gaps

1. **The residual 1.54% JointReaction-vs-self-computed gap is not fully decomposed.** It is
   plausibly related to this free-body chain's much smaller mass (3.85 kg vs the knee/hip's
   full-leg chains) giving the Savitzky-Golay-differentiated inertial term a smaller
   signal-to-differentiation-noise ratio — but this is a hedged, plausible candidate, not a
   measured decomposition; the gap is small (an order of magnitude tighter than the ~20-70%
   distortions the sign bug itself caused) and does not threaten the corrected headline number.
2. **All pre-existing fidelity limits of the graft itself are unchanged and out of scope**
   (no scapulothoracic joint, no wrapping surfaces ported, brachioradialis absent, single
   synthetic trajectory, EMG-unvalidated activation pattern — full list in
   `docs/MECHANISM_ARM_MUSCLES.md` §8). This remediation is strictly about the sign convention
   of one downstream geometry function, nothing upstream or downstream of it.
3. **The moment-arm-sign verification (47/70 PASS, §6 of the arm-muscles doc) is a different
   mechanism** (`Muscle.computeMomentArm`'s own sign convention, checked in `add_arm_muscles.py`,
   not `arm_crossing_muscles_and_forces`) and was not re-touched or re-verified here.
4. **Single trial, right arm only, synthetic (not measured) kinematics** — same scope caveat
   every cert in this family carries; not re-verified across subjects/trials/loading cases.

## 8. Files

- `scripts/msk/validate_shoulder_force_with_muscles.py` — the one-line fix (§1) with an inline
  comment pointing to this doc; a new `JR_SETUP_TEMPLATE`/`run_joint_reaction`/
  `extract_jr_force_pct_bw` (STEP 6, §3); `main()` updated to run and report the JR cross-check.
- `scripts/msk/audit_shoulder_muscles_sign_bug.py` — **new**, the toy-test verification (§2).
- `docs/MECHANISM_ARM_MUSCLES.md` — top correction banner + inline SUPERSEDED pointers at the
  headline table and §7 (original prose preserved verbatim, not deleted).
- `data/msk_smoketest/subject2_arm_abduction/shoulder_force_with_muscles/shoulder_force_with_muscles_results.json`
  — regenerated by the post-fix end-to-end re-run; now also carries `jointreaction_run`,
  `jointreaction_summary`, `jointreaction_columns_found`,
  `jointreaction_vs_selfcomputed_relative_diff_pct`, `jointreaction_at_selfcheck_peak_time_pctBW`.
- `data/msk_smoketest/subject2_arm_abduction/shoulder_force_with_muscles/jr_output/` — new,
  the JointReaction AnalyzeTool output (setup XML + `*_JointReaction_ReactionLoads.sto`).
- `data/msk_smoketest/subject2_arm_abduction/shoulder_muscles_sign_bug_audit/shoulder_muscles_sign_bug_toy_test_results.json`
  — new, the toy-test machine output (§2).
