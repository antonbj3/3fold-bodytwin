# MECHANISM SIGN BUG AUDIT — independent QC of the knee+hip crossing-muscle sign bug (2026-07-21)

Adjudicates a contradiction between two sibling certs. The spine-cert agent
(`docs/MECHANISM_SPINE_FORCE.md` Sec.5) found a machine-confirmed sign bug in
`scripts/msk/static_opt_knee.py`'s `knee_crossing_muscles_and_forces` and flagged it as unfixed,
still live, in that file's own knee self-computed number (233.20 %BW) — but a global sign flip
should produce a non-physical NEGATIVE result (as it did for the spine, −214 to −284 %BW before
that agent's local fix), yet the knee's published self-computed number is POSITIVE and lands close
to OrthoLoad (233.20/258.22 = 0.90), and the hip's (which imports and calls the EXACT SAME function
object, unmodified) matches too (235.30/273.93 = 0.86). Both cannot be face-value-consistent.
This audit measures, not assumes, which is true. Every number below is machine-measured this
session (`scripts/msk/audit_sign_bug.py`, exit 0), not recalled. Isolation respected: `.venv-msk`
only; `static_opt_knee.py`/`validate_hip_force.py` imported and called UNMODIFIED, never edited;
model/mocap/SO/JR data read in place, never written to; no git operations.

## Headline verdict

| joint | published self-computed | ratio vs OrthoLoad | **CORRECTED** self-computed | **ratio vs OrthoLoad** | independent JointReaction referee | corrected vs JR |
|---|---:|---:|---:|---:|---:|---:|
| **knee** | 233.20 %BW | 0.903 | **391.10 %BW** | **1.515** | 391.11 %BW (unchanged) | **0.0039% apart** |
| **hip** | 235.30 %BW | 0.859 | **386.77 %BW** | **1.412** | 387.04 %BW (unchanged) | **0.07% apart** |

**KNEE: CONFIRMED — CORRECTED-TO-391.10 %BW.**
**HIP: CONFIRMED — CORRECTED-TO-386.77 %BW.**

Neither "UNRESOLVED" nor "the knee/hip code is correct." The bug is real, and it is LIVE, in both
joints' self-computed pipelines — byte-identical code, not merely "the same pattern"
(`validate_hip_force.py` line 325 calls `sok.knee_crossing_muscles_and_forces` directly). The
published 233.20 / 235.30 %BW numbers, and their close-looking 0.90 / 0.86 match to OrthoLoad, are
**artifacts of the bug**, not validated results. There is no contradiction with the spine's finding
— see Sec.3 for exactly why the same bug looks "fine" at the knee/hip and "impossible" at the
spine. The corrected self-computed contact force converges to within a few hundredths of a percent
of the already-published, bug-immune, official `opensim.JointReaction` cross-check at BOTH joints
— the strongest evidence in this document that the correction (not the original) is right.

## 1. Part 1 — known-answer toy-case unit test on the REAL, unmodified function

`scripts/msk/audit_sign_bug.py` imports `static_opt_knee.knee_crossing_muscles_and_forces`
UNMODIFIED (never edited) and calls it directly against minimal duck-typed fake OpenSim path-point
objects (no OpenSim model needed for this part — model-independent, exactly the spine agent's own
falsification style, independently reproduced rather than trusted). Ground truth: a cable/muscle
under tension pulls each of its path endpoints TOWARD the other endpoint — basic, model-independent
mechanics, matching the function's OWN docstring ("the unit pulling direction **toward the outside
anchor**"). Four geometrically diverse cases, plus the audit's locally-defined FIXED function
(branches of the ternary swapped, mirroring the spine agent's own fix — defined standalone in the
audit script, never written into `static_opt_knee.py`):

| case | computed | ground truth | verdict |
|---|---|---|---|
| A: inside point FIRST in path order | (0, +1, 0) | (0, −1, 0) | **SIGN_FLIPPED** (exact negative) |
| B: inside point SECOND in path order | (0, +1, 0) | (0, −1, 0) | **SIGN_FLIPPED** (exact negative) |
| C: crossing in the MIDDLE of a 3-point path | (0, +1, 0) | (0, −1, 0) | **SIGN_FLIPPED** (exact negative) |
| D: arbitrary 3-D axis orientation | (0.498,−0.830,0.249) | (−0.498,0.830,−0.249) | **SIGN_FLIPPED** (exact negative) |
| E: audit's local FIXED function, case A | (0, −1, 0) | (0, −1, 0) | **CORRECT** |

All four real-function cases are the EXACT negative of the correct answer (`np.allclose(computed,
-ground_truth)==True` in every case, machine-checked, not eyeballed) — regardless of path-point
order (A vs B), regardless of whether the crossing is at the path's edge or interior (C), and
regardless of axis orientation (D). This rules out the sign flip being an artifact of any specific
geometric configuration: it is a pure index-logic bug in the ternary `(k, k+1) if in_k1 else (k+1,
k)`, unconditionally backwards for both branches, confirmed independently of — and consistent with
— the spine agent's own finding. The local FIXED function (branches swapped) recovers the correct
answer (case E), confirming the fix itself before using it in Part 2.

## 2. Part 2 — re-derived knee + hip contact force, reusing cached SO/JR output (no Ipopt re-run)

`compute_chain_self_check` in the audit script is a generic re-implementation of
`static_opt_knee.compute_self_cross_check` / `validate_hip_force.compute_self_cross_check_generic`,
parameterized by WHICH crossing-muscle detector to call, so the identical geometry/Newton's-law
pipeline runs once with the PUBLISHED (buggy) detector and once with the audit's FIXED one, without
editing either source file. Static Optimization + JointReaction output already on disk
(`data/msk_smoketest/subject2_walking1/static_optimization/{so,jr}/`, from the original knee-cert
session) is reused directly — no re-run of Ipopt (expensive, unnecessary: this bug is entirely in
the Python-side geometry step downstream of SO, never in SO itself).

**Fidelity gate (this audit's own re-implementation of the PUBLISHED pipeline must reproduce the
committed numbers before anything else here can be trusted):**

| joint | published %BW | this audit's re-run of the SAME (buggy) code | relative diff | verdict |
|---|---:|---:|---:|---|
| knee | 233.20106467 | 233.20106467 | 0.000000% | **PASS** |
| hip | 235.29567885 | 235.29567885 | 0.000000% | **PASS** |

**Exact algebraic identity gate.** Because the fix changes ONLY the sign of the per-crossing unit
vector (linear superposition, nothing else touched), the corrected vector must satisfy
`F_fixed(t) = 2·R_old(t) − F_buggy(t)` at every frame, every axis — a forced consequence of
`F = R − M` with `M_buggy = −M_true` exactly. Measured directly on the full 158-frame vector series
(not just the peak): max abs deviation from this identity = **4.5×10⁻¹³ N** for both joints (14
orders of magnitude below the ~1000-3000 N signal) — **PASS**, both joints. (An earlier version of
this script shared one mutable OpenSim `state` object across all 4 sequential calls and found a
small, 5×10⁻⁶ N, non-zero `R_old` diff between the knee's own buggy/fixed runs, which should be
IDENTICAL by construction since `R_old` never depends on the crossing-muscle detector at all. Forced
via OODA rather than shrugged off as noise: diagnosed as cross-call warm-start hysteresis in
OpenSim's iterative assembly of this model's polynomial-coupled "rolling knee" `CustomJoint`
[the same joint flagged elsewhere in this repo, `docs/MECHANISM_MSK_ELASTIC_BAND.md` Sec.4, as a
numerically sensitive component] — consistent with the hip chain, whose simpler joint showed
EXACTLY 0.000e+00 diff under the same shared-state version. Fix: give every one of the 4 calls its
own fresh `osim.Model`/`state`; the diff dropped to exact floating-point epsilon for both joints,
confirming the diagnosis.)

**Corrected numbers:**

| | knee | hip |
|---|---:|---:|
| Published self-computed (buggy) | 233.20 %BW @ t=0.51s | 235.30 %BW @ t=0.67s |
| **Corrected self-computed (audit_fixed)** | **391.10 %BW @ t=0.51s** | **386.77 %BW @ t=0.55s** |
| Official JointReaction (independent, untouched by this bug) | 391.11 %BW @ t=0.51s | 387.04 %BW @ t=0.55s |
| OrthoLoad in-vivo anchor | 258.22 %BW | 273.93 %BW |
| Corrected ratio vs OrthoLoad | **1.515** (overshoot +51.5%) | **1.412** (overshoot +41.2%) |

The corrected knee number matches JointReaction's peak time exactly (both 0.51s, unchanged — this
agreement already existed even with the bug). The corrected HIP number's peak time **moves from
0.67s to 0.55s** — snapping onto JointReaction's own peak time exactly, whereas the published buggy
hip number peaked at a different instant (0.67s) than JointReaction (0.55s). Both the magnitude
convergence and this independent peak-time convergence are strong, non-tautological, externally
anchored confirmations that the CORRECTION (not the original) is the physically right answer:
JointReaction is a genuinely separate OpenSim/Simbody code path (`calcMobilizerReactionForces`-style
internals fed by the SO `force.sto` via its documented `forces_file` mechanism) that never calls
`knee_crossing_muscles_and_forces` — it cannot be circularly "confirming itself."

## 3. Why the SAME bug is invisible at the knee/hip but screamed "impossible" at the spine

The premise of "the knee/hip code is correct, only the spine's reuse is wrong" is **refuted** by
Secs.1-2 above — the code is equally wrong in both places. What differs is not the bug, but **what
scalar gets extracted from the (buggy or corrected) vector**:

- **Spine** (`validate_spine_force.py` line 750): `compressive = float(np.dot(F_bone_contact,
  axis))` — a **signed projection** onto one fixed anatomical axis (torso-local +Y). A sign error
  in a subtracted term can and did flip this scalar's SIGN, producing an immediately-flagged
  impossible negative (tensile) number.
- **Knee/hip** (`static_opt_knee.py` line 542, `validate_hip_force.py`'s equivalent):
  `Fbc_mag = np.linalg.norm(cross["F_bone_contact_new_vec"], axis=1)` — an **unsigned vector
  magnitude**. `np.linalg.norm(...)` is `>=0` by construction, REGARDLESS of any internal sign
  error. A norm cannot expose a sign bug as "negative" — it can only expose it as "wrong magnitude,"
  which may coincidentally still look plausible.

This audit measured, rather than assumed, that this is the whole explanation — by computing the
analogous SIGNED ground-vertical (y) component at each method's own peak instant (a rough proxy for
"compressive axis" for an upright stance leg, offered as illustrative context, not a formal gate the
way the spine's carefully-verified torso-local-Y check was):

| | knee y-component | hip y-component |
|---|---:|---:|
| Published (buggy) | **+1101.6 N** (+143.6 %BW) — looks like a plausible push | **+1473.8 N** (+192.2 %BW) — looks plausible |
| Corrected (fixed) | **−2605.1 N** (−339.7 %BW) — an "impossible" pull, same signature as the spine's bug | **−2815.3 N** (−367.1 %BW) — same signature |

Had the knee/hip cert reported a signed axial projection (as the spine does) instead of a vector
norm, it would have hit the EXACT same "non-physical negative" trip-wire that tipped off the spine
agent — and the bug would very likely have been caught then, not now. The buggy and corrected
peak-instant vectors point in nearly (not exactly, since `R_old ≠ 0`) opposite directions: 157.7°
apart at the knee, 162.6° at the hip — the expected signature of `F_buggy = R + M_true` vs
`F_true = R − M_true` when `|M_true|` (306.6 %BW at the knee's peak instant) substantially exceeds
`|R_old|` (102.2 %BW) — consistent with, not merely asserted alongside, the algebraic identity in
Sec.2.

**This is a genuine methodological lesson, stated symmetrically (not to relitigate the original
certs, which are otherwise carefully forced):** the ~1.6-1.7x self-computed/JointReaction
disagreement WAS investigated as a forced adversary in both `docs/MECHANISM_STATIC_OPT.md` Sec.5 and
`docs/MECHANISM_HIP_FORCE.md` Sec.5 — but the adversary tested was "did JointReaction accidentally
reread the wrong joint's columns" (correctly cleared) and "is this a differentiation-scheme
artifact" (accepted as the explanation). The self-computed method's OWN sign convention was never
on the table as a suspect, because its output (a positive-looking, close-to-anchor ratio) did not
trigger the kind of alarm a negative number does. A result that looks clean is not the same as a
result that has been forced against its own strongest adversary — this is exactly the
self-flattering-pass pattern the "kills are auditable, an unforced adversary is a false pass" rule
exists to catch.

## 4. The 1.6-1.7× self-computed/JointReaction gap: re-attributed, not merely "closed"

Both original certs found and disclosed a recurring ~1.6-1.7× disagreement between self-computed and
JointReaction, diagnosed as "most likely... different numerical differentiation schemes," and
treated its recurrence across two decorrelated joints as corroboration of that diagnosis:

| joint | JR / published-buggy ratio | JR / corrected ratio |
|---|---:|---:|
| knee | 1.677 (matches the original doc's own "1.677") | **1.00004** |
| hip | 1.645 (matches the original doc's own "1.645") | **1.00070** |

Fixing an entirely different thing (a Python-side muscle-direction sign convention, nothing about
differentiation) collapses the gap from 65-68% to 0.004-0.07% at BOTH joints. This is strong
evidence the differentiation-scheme explanation, though reasoned carefully, was **not the actual
mechanism** — recurrence across two joints corroborated the shared-CODE-bug hypothesis just as well
as the shared-differentiation-scheme hypothesis (both predict a recurring ratio across joints that
share the same pipeline); the two hypotheses were never actually distinguished until this direct
correction. Any genuinely remaining differentiation-scheme effect is now bounded at ≤0.07% — noise
level, not the dominant term it was believed to be.

## 5. Machine-checked gates (PASS/FAIL, not eyeballed)

| gate | pre-registered threshold | measured | verdict |
|---|---|---:|---|
| Toy case: real function sign-flipped, all 4 diverse geometries | `np.allclose(computed, -truth)` | 4/4 | **PASS** |
| Toy case: local FIXED function recovers ground truth | `np.allclose(computed, truth)` | yes | **PASS** |
| Fidelity: audit's buggy re-run reproduces published knee/hip numbers | rel. diff < 0.1% | 0.000000% both | **PASS** |
| Algebraic identity: F_fixed == 2·R_old − F_buggy (per-frame, 3 axes) | max abs diff < 1e-6 N | 4.5e-13 N both | **PASS** |
| R_old identical regardless of crossing-fn (buggy vs fixed) | max abs diff < 1e-6 N | 0.0 N both (fresh-state version) | **PASS** |
| JointReaction re-extraction reproduces published numbers | exact | 391.11478 / 387.03722, bit-identical | **PASS** |
| Corrected self-computed vs JointReaction agreement | — (reported, not gated) | 0.0039% knee, 0.070% hip | over-determined |
| Corrected-vs-OrthoLoad classification (pre-existing gate, ratio>0.6) | PASS-substantial | 1.515 knee, 1.412 hip | **PASS-substantial**, both (regardless of over/undershoot) |
| **Overall (Part1 ∧ Part2 fidelity ∧ identity)** | all PASS | | **PASS**, script exit 0 |

## 6. Verdict per joint (crisp, per the task's own format)

- **KNEE: CONFIRMED — CORRECTED-TO-391.10 %BW.** The published 233.20 %BW (ratio 0.90 vs OrthoLoad)
  was a bug artifact. The corrected number (391.10 %BW, ratio 1.515) matches the independent,
  bug-immune JointReaction cross-check (391.11 %BW) to 0.0039% — effectively the same answer via two
  structurally different code paths. The "0.90, lands just under the anchor" framing does not
  stand; the twin's true muscle-driven knee estimate overshoots OrthoLoad by ~51%, the same
  direction and magnitude the JointReaction cross-check ALREADY independently reported and the
  original cert already pre-registered as PASS-substantial regardless of direction.
- **HIP: CONFIRMED — CORRECTED-TO-386.77 %BW.** The published 235.30 %BW (ratio 0.86) was likewise
  a bug artifact. The corrected number (386.77 %BW, ratio 1.412) matches JointReaction (387.04 %BW)
  to 0.07%, AND its peak time snaps from 0.67s onto JointReaction's own 0.55s. Same conclusion as
  the knee: the "0.86, close match" framing does not stand; the true estimate overshoots by ~41%.
- Neither verdict is UNRESOLVED: fidelity, identity, and external-JointReaction agreement all
  triangulate on the same corrected numbers from independent angles (toy-case ground truth, algebraic
  necessity, and a separately-implemented OpenSim analysis class).

## 7. Honest gaps (full list)

1. **`static_opt_knee.py` and `validate_hip_force.py` were NOT edited** (explicit task instruction;
   read-only audit) — the bug remains LIVE in the committed code right now. This document is
   diagnostic, not a fix-in-place.
2. **`docs/MECHANISM_STATIC_OPT.md` and `docs/MECHANISM_HIP_FORCE.md` (and their `_results.json`
   files) were likewise not edited** — their headline self-computed numbers (233.20/235.30 %BW) and
   "matches OrthoLoad" framing are now known-stale and should be corrected by a future editor
   (flagged here, not fixed, same scope discipline the spine agent itself used on this exact file).
3. **The ground-vertical y-component check (Sec.3) is a rough illustrative proxy**, not a rigorously
   verified "compressive axis" the way the spine cert verified torso-local-Y ≈ world-Y at its
   specific pose — offered as supporting context for WHY the symptom differs, not as a formal gate.
4. **This audit does not relitigate whether a 41-52% OrthoLoad overshoot is "acceptable"** — the
   ratio>0.6-is-PASS-substantial rule (direction-agnostic) was pre-registered by the original cert
   authors before this audit existed; this audit only recomputes the number and applies that
   pre-existing rule mechanically.
5. **OrthoLoad anchors (258.22 / 273.93 %BW) were reused verbatim**, not re-parsed or re-verified
   here — out of scope (the question was about the twin's own computation, not the anchor corpus).
6. **Single trial, right side only** (subject2 `walking1`) — same scope caveat as every upstream
   cert in this family; not re-verified across subjects/trials/left leg.
7. **The spine cert's own already-applied fix and its 296.59/389.17 %BW Tier-2 numbers were not
   re-derived or re-checked here** — out of scope; this audit is knee+hip only, per the task.
8. **SO is an effort-minimizing solution, not measured EMG** — inherited, unchanged limitation from
   every upstream cert in this family; this audit changes the sign of a downstream geometry term,
   not the underlying Static Optimization solve itself (SO `force.sto`/`activation.sto` reused
   verbatim, byte-identical, from the original knee-cert session; re-confirmed via the fidelity gate,
   not re-run).

## Files

- `scripts/msk/audit_sign_bug.py` — the full audit (self-contained, re-runnable; imports
  `validate_joint_force.py`/`static_opt_knee.py`/`validate_hip_force.py` UNMODIFIED for proven
  parse/BFS/SO-reuse/JR-extraction code and the function under audit itself; defines its own local
  FIXED crossing-detector and a generic multi-variant re-implementation of the self-cross-check
  pipeline, never edits any of the three imported files).
- `data/msk_smoketest/subject2_walking1/sign_bug_audit/sign_bug_audit_results.json` — every number
  in this document, machine-written (toy-case results, per-joint per-variant peak vectors, fidelity/
  identity gate results, JointReaction re-extraction, final verdicts).
- Reused, not written by this script: `data/msk_smoketest/subject2_walking1/static_optimization/
  {so,jr}/` (Static Optimization + JointReaction output, already existing from the original knee-cert
  session — no Ipopt re-run).
