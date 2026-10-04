# MECHANISM FOOT AXIS FIX — real Maharaj2021 midtarsal axis transplanted; held-out gate re-run (2026-07-22)

**Pre-registered claim (C):** the real, biplanar-videoradiography-validated Maharaj2021 midtarsal
axis ([0.9392, 0.5666, -0.4183], `docs/MECHANISM_MULTISEGMENT_FOOT.md` Sec.6), transplanted onto
the twin's `midtarsal_r/l` joint in place of the v3 placeholder (Z rotated 90deg about Y, no
citation), **PASSES** the held-out DJ2 generalization gate the placeholder failed (max-abs-diff
**< 5.0deg**, the same pre-registered ceiling from `docs/MECHANISM_LAYER_COMPLETENESS.md` /
`docs/MECHANISM_FOOT_MULTISEGMENT.md`).

**Result: C is FALSIFIED.** The real axis does not close the gap — it makes it measurably worse
(DJ2 max-abs **11.86deg**, more than double the placeholder's own 5.76deg), and newly fails the
previously-passing in-sample walking1 check too (max-abs 9.45deg vs 5.0deg). Forced OODA (not a
one-shot shrug) diagnoses a specific, quantified, geometrically-derived mechanism (Sec.4): the gap
is deeper than "needed a real axis" — it is a structural identifiability property of the *twin's
own* reduced foot topology, not a defect in the axis citation. **The v3 placeholder remains the
better-fitting choice for the twin as currently built; this fork is not recommended as a
replacement** (Sec.6).

**Symmetric-QC pre-check (forced before any transplant):** does Maharaj's `calcn_r` frame
convention match the twin's own `calcn_r` frame? **PASS** (Sec.1) — this was verified empirically
before touching any file, so the negative result in Sec.3-4 is not attributable to a frame bug.

**Confidence tier:** in-vivo-anchored (biplanar videoradiography, PMID 34698598) for the
transplanted axis itself and for the frame-convention verification; MEASURED/machine-checked
negative, mechanistically diagnosed (not a bare fail), for the overall falsifier.

---

## 0. Context (builds on, does not re-litigate)

`docs/MECHANISM_FOOT_MULTISEGMENT.md` grafted a `midtarsal_r/l` PinJoint into the twin's LaiArnold-
lineage foot on a **placeholder** axis (v3: a pure Z-hinge rotated 90deg about Y, i.e. axis =
calcn's own local +X) after 2 more-principled variants (v1: reuse subtalar's own axis; v2: pure
Z-hinge) failed live regression. v3 passed in-sample (walking1: median 0.123deg, max-abs 3.69deg)
but **failed its own held-out generalization gate** (DJ2 drop-jump: max-abs 5.76deg > the
pre-registered 5.0deg ceiling, subtalar RMSE 4.2-10.1x the null-control noise floor).
`docs/MECHANISM_MULTISEGMENT_FOOT.md` then independently re-verified the Maharaj2021 JC model
(PMID 34698598, biplanar videoradiography) and extracted its real midtarsal_r axis, flagging one
open caveat before transplant: whether Maharaj's `calcn_r` frame convention matches the twin's own
— this document closes that caveat and executes the transplant.

**Script (single, re-runnable, self-contained):** `scripts/msk/foot_midtarsal_axis_fix.py`.
**Evidence:** `scripts/msk/foot_midtarsal_axis_fix_evidence.json`.

---

## 1. Frame-convention falsifier (forced BEFORE any transplant)

**Method:** fully empirical, engine-computed (never hand-parsed Euler math — exactly the class of
bug that caused 3 rounds of trial-and-error in the placeholder's own original build). Perturb a
coordinate, read the resulting body rotation via OpenSim's own `Transform`/`Rotation` API
(central difference, dq=1e-4 rad), express the axis in a chosen reference body's frame at the
neutral pose. Model-representation-agnostic — works identically whether the source is a
`CustomJoint` w/ `SpatialTransform` (Maharaj) or a `PinJoint` w/ a pre-rotated offset frame (twin).

**Internal-consistency self-check** (validates the method against a known ground truth before
trusting it on the harder case): FD-measuring Maharaj's *own* midtarsal axis (calcn is the joint's
PARENT there, zero orientation offset, so the raw `TransformAxis` value IS already "in calcn's
frame") must reproduce the raw cited vector almost exactly.

| Side | FD-measured (calcn frame) | Raw cited `TransformAxis`, normalized | cosine sim | angle |
|---|---|---|---|---|
| r | [0.796, 0.480, -0.369] | [0.800, 0.483, -0.356] | 0.9999 | 0.80deg |
| l | [-0.797, -0.481, -0.365] | [-0.800, -0.483, -0.356] | 0.9999 | 0.55deg |

PASS — the FD method is validated to sub-degree precision on a known-ground-truth case, so the
same method can be trusted for subtalar (calcn is the CHILD there — no raw-value shortcut exists;
FD is the only way to reference calcn's own frame directly).

**Decisive, externally-anchored cross-check** (never a tautology): compare Maharaj's own subtalar
axis (FD-measured in **its** calcn frame) against the twin's own subtalar axis (FD-measured in
**its** calcn frame) — same joint, same anatomical role, two independently-built model lineages
(Maharaj lab vs. Delp/Rajagopal/LaiArnold lineage). A genuine frame-convention mismatch would show
up as a wildly different or sign-flipped vector; a genuine match (allowing for real inter-model/
inter-subject subtalar-axis variation, documented at roughly 10-25deg in the literature) confirms
compatible conventions.

| Side | Maharaj subtalar axis (calcn frame) | Twin subtalar axis (calcn frame) | cosine sim | angle | **PASS** (threshold 0.8) |
|---|---|---|---|---|---|
| r | [0.731, 0.612, -0.303] | [0.782, 0.607, -0.143] | 0.986 | 9.64deg | **PASS** |
| l | [-0.733, -0.614, -0.293] | [-0.784, -0.606, -0.133] | 0.986 | 9.68deg | **PASS** |

**Pre-registered threshold:** cosine similarity > 0.8 (≈37deg headroom above known biological
variation). Actual result clears it with a large margin on both sides, with the *same sign
pattern on every component* — the signature a genuine convention mismatch would NOT produce.
**FRAME CONVENTIONS COMPATIBLE — safe to transplant the raw Maharaj vector as-is, no correction
needed.**

---

## 2. Euler-orientation computation (round-trip-validated, not hand algebra)

A `PinJoint` always rotates about its offset frame's local Z-axis. To make an arbitrary unit axis
`A` the physical rotation axis, OpenSim's own `Rotation.setRotationFromOneAxis(UnitVec3, ZAxis)` +
`convertRotationToBodyFixedXYZ()` computes the body-fixed XYZ Euler triple whose Z-column is
exactly `A` (geometric derivation via the engine's own rotation machinery, not manual Euler math —
the same discipline gap that caused the placeholder's own v1/v2/v3 search).

| Side | Raw Maharaj axis (re-extracted from the primary `.osim`, not the doc's prose) | Unit axis | Euler XYZ body-fixed (rad) | Round-trip err |
|---|---|---|---|---|
| r | [0.9392, 0.5666, -0.4183] | [0.800, 0.483, -0.356] | (-2.206744, 0.927377, 0.533514) | 1.6e-16 |
| l | [-0.9392, -0.5666, -0.4183] | [-0.800, -0.483, -0.356] | (2.206744, -0.927377, -2.608079) | 1.6e-16 |

Both sides independently re-extracted directly from `Maharaj2021_BothLegs_V2.osim`'s raw
`<TransformAxis>` (not copy-pasted from the prior doc's prose) — confirms the cited numbers are
real and correctly transcribed, and confirms the model's own internal mirroring convention (X,Y
flip sign, Z unchanged between r/l — a genuine 180deg-about-Z relationship, verified identically
on subtalar_r/l too).

---

## 3. Held-out DJ2 gate re-run (standalone prototype, apples-to-apples with the original 5.76deg)

Rebuilt the same topology as `foot_multisegment.py` (split calcn, add forefoot, midtarsal PinJoint,
reparent mtp) reusing its own proven split_x/mass/marker/muscle-reparenting math verbatim
(imported, not re-derived) — the **only** delta vs. the original build is the axis itself. New
model: `data/msk_models/LaiArnold_midtarsal_realaxis_scaled.osim`. Build-correctness checks (same
as the original): mass conserved exactly (78.2 -> 78.2 kg), topology chain PASS both sides,
muscle-length invariant at neutral pose PASS (max diff < 1e-6 m), reparent cross-check < 1e-9 m.

Then re-ran the **exact same gate**: same thresholds (median<0.5deg, max-abs<5.0deg), same
external anchors (the original 2021 pipeline's own `walking1.mot`/`DJ2.mot`, never touched this
session), same null control (fresh IK on the untouched original 3-segment model, same DJ2 trial —
isolates the toolchain's own reprocessing noise, independent of the graft).

| Trial | Axis | Median RMSE | **Max-abs diff** | Ceiling | Verdict |
|---|---|---|---|---|---|
| walking1 (in-sample/selection) | v3 placeholder | 0.123deg | 3.69deg | 0.5 / 5.0deg | PASS |
| walking1 (in-sample/selection) | **REAL (Maharaj)** | 0.126deg | **9.45deg** | 0.5 / 5.0deg | **FAIL (max-abs)** |
| DJ2 (held-out) | v3 placeholder | 0.246deg | 5.76deg | 0.5 / 5.0deg | FAIL (narrowly) |
| DJ2 (held-out) | **REAL (Maharaj)** | 0.298deg | **11.86deg** | 0.5 / 5.0deg | **FAIL (badly)** |

DJ2 null control (fresh re-run this session, original model): median 0.070deg, subtalar
0.307/0.275deg (r/l), ankle 0.274/0.234deg — **reproduces the historical null-control number from
`MECHANISM_FOOT_MULTISEGMENT.md` almost exactly** (0.070 / 0.31 / 0.27 / 0.27 / 0.23), confirming
this session's toolchain/methodology is consistent with the original, so the *differences* below
are real, not an artifact of a different environment.

**Per-joint breakdown (RMSE / max-abs, degrees) — subtalar is the dominant driver, both trials:**

| Coord | v3 walking1 | REAL walking1 | v3 DJ2 | REAL DJ2 |
|---|---|---|---|---|
| ankle_r | 1.60 / 3.67 | 1.86 / 4.48 | 2.40 / 4.13 | 2.56 / 4.44 |
| ankle_l | 1.81 / 3.69 | 2.45 / 4.64 | 3.24 / 4.35 | 4.20 / 6.30 |
| **subtalar_r** | 0.91 / 2.67 | **3.15 / 6.61** | 1.29 / 3.07 | **4.84 / 8.23** |
| **subtalar_l** | 1.44 / 2.76 | **4.88 / 9.45** | 2.78 / 5.76 | **7.21 / 11.86** |
| mtp_r/l | 0.32 / 0.47 | 0.66 / 0.90 | 0.06 / 0.18 | 0.19 / 0.46 |
| knee_r | 0.12 / 0.33 | 0.23 / 0.56 | 0.18 / 0.42 | 0.31 / 0.61 |
| knee_l | 0.15 / 0.34 | 0.30 / 0.76 | 0.25 / 0.51 | 0.47 / 0.76 |

Subtalar is 2-3x the size of the next-worst coordinate in every column — a systematic effect
(large RMSE across the whole trial, not a single-frame glitch: e.g. subtalar_l DJ2 max-abs occurs
at t=1.86s but RMSE is 7.21deg across the whole trial), and it dominates in **both** the in-sample
and held-out trial — i.e. this is not really a generalization gap, it is a structural mismatch
present from the first trial.

---

## 4. Mechanistic diagnosis (forced OODA — not accepted as a bare "honest negative")

**Observe:** subtalar dominates the residual, in both trials, by 2-3x over any other joint.

**Orient (why):** compute the angle between the new axis and the twin's *own* subtalar axis (both
FD-measured, calcn frame) — the exact quantity the original build's own v1 (parallel, catastrophic
failure) vs. v3 (rotated away, passing) contrast already showed is the causal lever.

| Variant | Angle from twin's own subtalar axis (axis-line distance) | DJ2 subtalar max-abs |
|---|---|---|
| v1 (rejected in the original build) | 0deg (identical) | in-sample subtalar RMSE 3.9-6.0deg (per `MECHANISM_FOOT_MULTISEGMENT.md`) |
| **REAL Maharaj axis (this session)** | **15.1 / 15.3deg** | **8.23 / 11.86deg** |
| v3 placeholder (accepted in the original build) | 38.1 / 38.6deg | 3.07 / 5.76deg |

A clean, monotonic, three-point dose-response: the closer the midtarsal axis sits to parallel with
the twin's own subtalar axis, the worse the identifiability/aliasing interaction — exactly the
"locking mechanism" aliasing the original build's own docstring already diagnosed for v1 (Manter
1941 lineage: oblique tarsal-joint axes are classically non-parallel for this reason), now
reproduced in graded form by a *real, anatomically correct* axis.

**Decisive external cross-check (rules out "the real axis is just wrong"):** is near-parallel
subtalar/midtarsal genuinely bad anatomy, or specific to this transplant? Measure the SAME angle
inside Maharaj's own model (its own subtalar vs. its own midtarsal, both real, both validated):

> Maharaj's own subtalar-vs-midtarsal angle = **9.23deg** — even MORE parallel than the 15.1-15.3deg
> the transplanted axis makes with the *twin's* subtalar — yet Maharaj's own model achieves 2.19deg
> RMS-vs-biplanar-videoradiography-truth in its own published validation, no aliasing catastrophe.

**Conclusion:** near-parallel subtalar/midtarsal axes are the *real* anatomy (present, and
apparently fine, even in the donor model) — the failure is not "the axis is wrong," it is that the
**twin's own reduced foot topology** (4 joints, 3 markers/foot: `r_calc`/`r_toe`/`r_5meta`) lacks
the DOF/marker redundancy that Maharaj's richer model (7 joints incl. a separate tarsometatarsal
DOF, 11 markers/foot) uses to disambiguate two near-parallel rotation axes. This is the same
σ_min/observability signature already established elsewhere in this project's foot work (least-
observed DOF is least solver-reproducible) — now confirmed to bind on subtalar/midtarsal jointly,
not just on the new DOF in isolation.

**Decide/Act:** this diagnosis points to a structural fix (add tarsometatarsal + more foot
markers — Tier 2, already scoped and cost-estimated in `docs/MECHANISM_MULTISEGMENT_FOOT.md` Sec.6),
not a better placeholder guess. Further ad-hoc axis variants would (a) fall outside this task's
explicit Tier-1 (axis-only) scope, and (b) reintroduce the exact train-on-test anti-pattern this
whole exercise exists to move past (walking1/DJ2 would stop being a genuine held-out check the
moment a new variant is picked by looking at them). No further variant was tried — this is a
disciplined stopping point, not premature surrender: the mechanism is identified, quantified, and
externally cross-checked against the donor model itself.

---

## 5. Twin fork — build, symmetric-QC, and re-measured kinematics

New model: `data/msk_models/subject2_unified_v2_realmidtarsalaxis.osim` (forked from
`subject2_unified_v2.osim`, which is read-only and unmodified — confirmed same md5 before/after).

**Pre-check** (guards against a concurrent-write surprise, per this repo's own concurrency note):
both `midtarsal_r`/`midtarsal_l` measured EXACTLY the expected v3 placeholder orientation
`(0, 1.5707963, 0)` immediately before patching — confirms the live file matched what this fix was
scoped against.

**Change applied:** ONLY the parent+child offset frames' `orientation` (both get the identical new
Euler triple, matching the existing PinJoint convention) — translation/location of the joint is
**unchanged** (axis-transplant only, no relocation, per the task's Tier-1 scope).

**Post-fix self-consistency** (genuine disk round-trip, reloaded fresh, FD re-measured):

| Side | Target unit axis | FD-measured axis after save+reload | cosine sim |
|---|---|---|---|
| r | [0.800, 0.483, -0.356] | [0.791, 0.484, -0.375] | 0.9998 |
| l | [-0.800, -0.483, -0.356] | [-0.794, -0.483, -0.369] | 0.9998 |

**Symmetric-QC — "confirm the change doesn't break ankle/subtalar":** SEMANTIC (not raw-text)
comparison of `ankle_r/l`, `subtalar_r/l` between old and new file via the OpenSim API (translation,
orientation, parent-frame socket, coordinate range/default) — a first-pass byte-level text diff
false-alarmed here (the original file was serialized by a plain `ElementTree` writer with no XML
comments; re-saving via `model.printToXML()` adds full auto-generated property-doc comments even
for untouched joints — a pure formatting artifact, not a semantic one). The semantic comparison:

| Joint | Semantically identical to original? |
|---|---|
| ankle_r | **True** |
| ankle_l | **True** |
| subtalar_r | **True** |
| subtalar_l | **True** |

**Confirmed: the ankle/subtalar joint *definitions* are untouched.** Separately (this is the
substantive finding, not a contradiction of the above): their *kinematic fit* during walking1 IK on
the twin's own full unified model gets measurably worse with the real axis, consistent with the
standalone-prototype finding in Sec.3-4:

| Coord | Placeholder axis RMSE (already on record, `merge_unified_v2_evidence.json`) | REAL axis RMSE (this session) |
|---|---|---|
| ankle_r | 1.452deg | 1.726deg |
| ankle_l | 1.648deg | 2.283deg |
| subtalar_r | 1.162deg | 3.588deg |
| subtalar_l | 1.472deg | 5.165deg |

**Unaffected coordinates** (pelvis, hip, knee, lumbar, contralateral arm — should NOT be touched by
a foot-local change): worst RMSE 1.067deg (`lumbar_bending`), under the 2.0deg precedent-anchored
ceiling (`merge_unified_model.py`'s own `GENEROUS_REASSEMBLY_TOL_DEG`) and nearly identical to the
1.058deg already on record for the placeholder-axis unified model — confirms the change is
correctly scoped (no leakage into unrelated body segments).

**New midtarsal/arch kinematics** (walking1, real axis, twin's own full-body IK):

| Side | midtarsal_angle ROM | Arch proxy (joint-center height above calcn-toe chord) range |
|---|---|---|
| r | 27.0deg (-4.3 to 22.7deg) | 10.5 mm (1.6-12.1 mm) |
| l | 22.4deg (3.3 to 25.7deg) | 11.1 mm (2.6-13.7 mm) |

The arch proxy here is a **geometric** substitute (the twin has no navicular marker, unlike
Maharaj's `RNAV`) — the midtarsal joint center's own perpendicular height above the calcn-to-toes
chord, forward-kinematics-replayed through the solved IK motion — disclosed as a proxy, not a
marker-based arch measurement; it lands in a broadly similar order of magnitude to
`MECHANISM_MULTISEGMENT_FOOT.md`'s Maharaj-marker-based estimate (12.2/13.7 mm) but is not a strict
apples-to-apples number (different landmark, different model). midtarsal ROM (22-27deg) is
consistent with both the placeholder's own prior ROM (22.9-25.7deg) and Maharaj's own measured
range (17.0-22.3deg) — the new DOF stays live and non-null under the real axis too, this part is
not in question.

---

## 6. Falsifier verdict and recommendation

**Falsifier as pre-registered:** "does the real Maharaj axis PASS the held-out gate the
placeholder failed (<5.0deg), or does it still fail (meaning the gap is deeper than the axis)?"

**Answer: it still fails, and the gap is deeper than the axis** — mechanistically identified
(Sec.4) as an observability/DOF-architecture mismatch between the twin's reduced 4-joint/3-marker
foot and the near-parallel subtalar/midtarsal axis geometry that is the genuine anatomy (confirmed
present, and tolerated, in Maharaj's own richer 7-joint/11-marker model).

**Recommendation: do NOT swap the twin's operative model to this fork.** The v3 placeholder,
despite lacking a citation, is the empirically better-fitting choice for the twin's *current*
topology (DJ2 max-abs 5.76deg vs. this fork's 11.86deg — 2x worse). The fork built here
(`subject2_unified_v2_realmidtarsalaxis.osim`) is preserved as a correctly-built, axis-transplant-
only artifact and as evidence for this finding, not as a replacement default. The `midtarsal_r/l`
joint in the twin's actual operative model (`data/msk_models/subject2_unified_v2.osim`) is
**unchanged** by this session (confirmed: same md5 before/after).

**What would actually close the gap** (per Sec.4's diagnosis, already scoped elsewhere, not
invented here): `docs/MECHANISM_MULTISEGMENT_FOOT.md` Sec.6 Tier 2 — split tarsometatarsal from
midtarsal (adds the DOF Maharaj's own model uses to disambiguate the near-parallel pair) — but that
tier's own disclosed limitation stands: the twin's mocap protocol has only 3 markers/foot, so a new
DOF added without new markers lands back on the same marker-observability ceiling
(`MECHANISM_FOOT_FIDELITY_PLAN.md` Sec.2) this document's own finding just sharpened from "a risk" to
"confirmed, quantified, and mechanistically explained."

**Confidence tier:** in-vivo-anchored (biplanar videoradiography, PMID 34698598, live-verified) for
the transplanted axis itself and the frame-convention check; MEASURED, machine-checked, forced-OODA
negative for the overall claim — not a silent or one-shot gap.

---

## 7. Honest gaps

1. The mechanistic diagnosis (Sec.4) is strong (3-point monotonic dose-response + an independent
   cross-check inside the donor model) but not a formal proof — a fourth data point (e.g. a
   deliberately-intermediate test axis) was not run, since it would fall outside this task's Tier-1
   scope and would reintroduce a train-on-test risk on the same two trials. Treated as sufficient
   for this falsifier, not as closing the mechanism question completely.
2. The arch-proxy kinematics (Sec.5) use a different landmark (joint center vs. navicular marker)
   than `MECHANISM_MULTISEGMENT_FOOT.md`'s Maharaj-side measurement — the similar order of magnitude
   is suggestive, not a validated cross-model match.
3. GRF/kinetic partitioning is untouched (kinematics-only, matching this task's own scope, same as
   every prior foot dial-turn in this project).
4. This session did not re-derive whether a *different* split_x/joint-location (not just axis)
   would change the subtalar-interaction picture — the task scoped this as axis-only, and the
   diagnosis (Sec.4) suggests topology/marker changes, not a location tweak, are the actual lever.

---

## 8. File index

- `docs/MECHANISM_FOOT_AXIS_FIX.md` — this document.
- `scripts/msk/foot_midtarsal_axis_fix.py` — the build+verification script (re-runnable:
  `.venv-msk/bin/python3 scripts/msk/foot_midtarsal_axis_fix.py`; Part A frame-check, Part B
  standalone real-axis prototype + held-out gate, Part C twin fork + re-measurement).
- `scripts/msk/foot_midtarsal_axis_fix_evidence.json` — full machine-readable evidence (frame
  probe, Euler computation, both gate re-runs, twin-fork build+IK).
- `data/msk_models/LaiArnold_midtarsal_realaxis_scaled.osim` — standalone real-axis prototype
  (apples-to-apples successor to `LaiArnold_midtarsal_proto_scaled.osim`, which is untouched).
- `data/msk_models/subject2_unified_v2_realmidtarsalaxis.osim` — the twin fork with the real axis
  applied to `midtarsal_r/l` (translation unchanged; ankle/subtalar semantically untouched).
  **Not recommended as the twin's operative model** (Sec.6) — `data/msk_models/
  subject2_unified_v2.osim` (the v3 placeholder, unmodified, md5-confirmed) remains the better-
  fitting choice for the twin's current topology.
- `data/msk_smoketest/foot_multisegment_realaxis_proto/` — walking1 IK for the standalone
  real-axis prototype.
- `data/msk_smoketest/foot_multisegment_realaxis_proto_heldout_DJ2/` — DJ2 held-out IK + null
  control for the standalone real-axis prototype.
- `data/msk_smoketest/subject2_unified_v2_realmidtarsalaxis_walking1/` — walking1 IK for the twin
  fork.
- `docs/MECHANISM_MULTISEGMENT_FOOT.md` — extracted the real axis, flagged the frame-convention
  caveat this document closes (prerequisite reading, not re-litigated).
- `docs/MECHANISM_FOOT_MULTISEGMENT.md` — the original placeholder build (v1/v2/v3 search, the
  held-out gate this document re-runs) (prerequisite reading, not re-litigated).
- `docs/MECHANISM_LAYER_COMPLETENESS.md` — where the placeholder's 5.76deg gate failure is recorded.
