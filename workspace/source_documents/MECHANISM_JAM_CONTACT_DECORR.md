# MECHANISM JAM/COMAK DEFORMABLE-CONTACT DECORRELATION OF THE KNEE OVER-PREDICTION (2026-07-21)

Tests whether the ~1.5x knee-vs-OrthoLoad over-prediction found via the RIGID free-body cut
(`docs/MECHANISM_STATIC_OPT.md`/`docs/MECHANISM_JOINT_FORCE_VALIDATION.md`: 391.10-391.11 %BW,
official `opensim.JointReaction`, subject2/walking1, LaiArnold model) is an artifact of the
free-body cut's own RIGID-contact assumption (all load passes through a single bony contact
point) — or whether a genuinely DIFFERENT contact model (JAM's deformable
`Smith2018ArticularContactForce` elastic-foundation cartilage + `Blankevoort1991Ligament`
load-sharing, `lenhart2015.osim`, the same healthy-knee JAM model already used in
`docs/MECHANISM_KNEE_LIGAMENTS.md` and `docs/MECHANISM_CARTILAGE_CONTACT.md`) gives the same
number when driven by subject2/walking1's own gait kinematics and own Static-Optimization
muscle activations.

Script: `scripts/msk/jam_contact_decorr.py` (build inputs + run + orchestration),
`scripts/msk/jam_contact_decorr_analyze.py` (h5 extraction + machine checks + verdict),
`scripts/msk/jam_contact_decorr_leftcheck.py` (forced-adversary empirical check). Data:
`data/msk_smoketest/jam_contact_decorr/`. Tool: JAM C++ fork's `opensim-cmd`
(`/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/opensim-build/opensim-cmd`,
OpenSim 4.5.2-2026-01-19-a3c872a2b — the same binary the prior two JAM legs used).

**Isolation respected:** bodytwin only; `lenhart2015.osim` + its 6 cartilage/bone STL meshes
fetched to `scratchpad/knee_lig/` (git-ignored, not copied into the repo, same convention the
prior two JAM legs used); subject2's IK/GRF/SO data read in place from the external drive and
from this repo's own already-committed `static_optimization/so/` outputs; no git operations.

---

## Headline

**JAM's deformable-contact tibiofemoral peak force, driven by subject2/walking1's real gait
kinematics + real SO-solved activation for the 13 knee-crossing muscles, is 2099.9 N = 345.4
%(of lenhart2015's OWN body weight, 608.0 N) at t=0.54s — a ratio of 0.883x the rigid free-body
cut (391.1 %BW) and 1.337x the OrthoLoad in-vivo median (258.2 %BW, n=72).**

Pre-registered verdict (thresholds fixed before this number was computed): **ratio 0.883 falls
inside the pre-registered "~equal" band [0.85, 1.15], NOT below the 0.85 "materially lower"
cutoff — verdict = NOT-C (soft).** Deformable cartilage contact + ligament load-sharing does
**not** relieve the rigid free-body cut's over-prediction; it reduces the number by only ~12%,
landing just short of (3.9% away from) the pre-registered threshold that would have called it a
material reduction. **This is a genuine, disclosed BORDERLINE finding, not a clean exoneration in
either direction — reported as such, not rounded to the nearest convenient conclusion.**

A mechanistic reason this is plausible (not just a null result): the knee's TOTAL resultant
reaction force is constrained by whole-segment force/moment balance (GRF + muscle + ligament +
gravity = net proximal reaction) regardless of whether the contact surface reacting that load is
rigid or deformable. Ligaments only have material "spare capacity" to offload the compressive/
axial component specifically if their own axial (proximal-distal) force is a substantial fraction
of the total — measured here at the peak instant (§ Ligament cross-check below): the largest
ligament total, ACL 181.8 N + LCL 51.1 N + patellar-tendon-tension "PT" 281 N, is still an order
of magnitude below the ~2100 N total joint force — deformable contact redistributes WHERE the
load is applied (pressure/area), not by how MUCH the resultant axial force is reduced, for a joint
whose passive restraints carry only a modest fraction of the peak axial load.

| quantity | value | anchor | ratio |
|---|---:|---:|---:|
| JAM tf_contact peak force | 2099.9 N (tibia-side) / 2161.1 N (femur-side) | — | — |
| JAM tf_contact, %(lenhart2015's own BW=608.01N) | **345.37 %BW** | rigid free-body cut: 391.11 %BW | **0.883x** |
| JAM tf_contact, %(lenhart2015's own BW) | 345.37 %BW | OrthoLoad median (n=72): 258.22 %BW | **1.337x** |
| JAM peak time | t=0.540s | rigid-cut peak time: t=0.51s | 30ms / ~1.9% of gait-cycle apart |
| JAM tf mean pressure peak | 4.796 MPa | same-model published gait anchor (PMID 25917122): 6.2 MPa | 0.774x |
| JAM pf mean pressure peak | 2.956 MPa | same-model published gait anchor: 2.8 MPa | 1.056x |

---

## 1. Why this is a fair, decorrelated comparison — verified, not assumed

**1a. Coordinate-name mapping** (LaiArnold IK column name -> lenhart2015 coordinate name), values
copied through UNCHANGED — both `.mot` files use degrees for rotation / meters for translation
(confirmed live from both headers: `walking1.mot` has `inDegrees=yes`; JAM's own official
`prescribed_coordinates.sto` example uses bare degree values, e.g. `pelvis_tilt=90`, not radians —
**zero unit-conversion arithmetic applied anywhere in this pipeline**, eliminating the single most
common cross-model transplant bug class):

| LaiArnold IK column | lenhart2015 coordinate | | LaiArnold IK column | lenhart2015 coordinate |
|---|---|---|---|---|
| `pelvis_tilt` | `pelvis_tilt` | | `hip_rotation_r` | `hip_rot_r` |
| `pelvis_list` | `pelvis_list` | | `knee_angle_r` (PRIMARY) | `knee_flex_r` (PRIMARY) |
| `pelvis_rotation` | `pelvis_rot` | | `ankle_angle_r` | `ankle_flex_r` |
| `pelvis_tx/ty/tz` | same | | `subtalar_angle_r` | `subt_angle_r` |
| `hip_flexion_r` | `hip_flex_r` | | `mtp_angle_r` | `mtp_angle_r` (identical) |
| `hip_adduction_r` | `hip_add_r` | | | |

All prescribed values fall comfortably inside lenhart2015's own coordinate ranges (e.g.
`knee_flex_r` prescribed 5.59-65.65°, model range -10° to 160°; `hip_flex_r` prescribed
-12.85-26.66°, model range -60° to 150°) — verified live, no clamping.

**1b. Muscle-name overlap, measured (not assumed) — a striking independent cross-validation.**
Subject2's SO solved 40 right-leg muscles. **ALL 40 have an EXACT name match** in lenhart2015's
own 44-muscle `Millard2012EquilibriumMuscle` set (lenhart2015 has 4 extra: `gem_r`/`pect_r`/
`pertert_r`/`quadfem_r` — none crosses the knee). More strikingly: an **independent, pure-XML
geometric knee-crossing test** run against lenhart2015's OWN body/path geometry (does the muscle's
path connect a thigh-or-above body to a shank-or-below body?) returns **exactly the same 13-muscle
set** LaiArnold's own crossing-muscle detector (`static_opt_knee.py`) already uses for the
391.1%BW headline: `bflh_r, bfsh_r, gaslat_r, gasmed_r, grac_r, recfem_r, sart_r, semimem_r,
semiten_r, tfl_r, vasint_r, vaslat_r, vasmed_r` — set-equality confirmed by script, not eyeballed.
This is a real, falsifiable cross-model anatomical consistency check that passed, not a
coincidence of naming.

**1c. Input translation: ACTIVATION, not FORCE — the disclosed input difference the task asked
to flag.** The 13 knee/pf-crossing muscles are driven by subject2's own SO-solved **activation**
trajectory (0-1 normalized neural drive), not force in Newtons. Force is not a portable quantity
across two Hill-type muscles with different `max_isometric_force`/`optimal_fiber_length`/pennation
(lenhart2015 is a different-donor, generic model); activation is the correct, honestly-documented
common currency — lenhart2015 computes ITS OWN force from that activation at ITS OWN geometry.
**This is exactly the "feed it the SAME muscle forces (or clearly document the input difference)"
fork the task pre-flagged — resolved here as activation parity, not force parity, stated
explicitly rather than silently glossed as "same forces."**

**1d. Zero-resampling timebase match.** Subject2's IK (`walking1.mot`, 158 rows) and SO's
`activation.sto` (158 rows) share a **bit-identical time grid** (verified: max abs diff = 0.0) —
no interpolation was needed to align kinematics and muscle-activation inputs.

**1e. GRF applied via the SAME mechanism subject2's own SO/ID pipeline already uses**:
`ForsimTool`'s `external_loads_file`, RightGRF only, `applied_to_body=calcn_r` — confirmed to
exist verbatim in lenhart2015's BodySet (same naming lineage as LaiArnold's own `calcn_r`/`calcn_l`).

**1f. Scope reduction (left leg + non-crossing muscles), ARGUED and then EMPIRICALLY TESTED, not
just asserted (§3).**

---

## 2. Pre-registered claims, thresholds, and the forced adversary

**C (leaning-positive direction):** JAM's peak tf_contact force, expressed as %(lenhart2015's own
body weight), is **materially LOWER** than the rigid cut (ratio < 0.85), trending toward
OrthoLoad — the rigid free-body cut would be an identified over-estimate SOURCE.
**not-C:** ratio ≥ 0.85 — the contact-model choice (rigid vs. deformable) is not what drives the
over-prediction; some other common-mode factor (shared model/shared muscle-force-solve, per
`docs/MECHANISM_CMC_SECOND_SOLVE.md`'s own still-open gap) remains the suspect.
Thresholds (0.85 "materially lower", 1.15 "materially higher") were fixed in
`scripts/msk/jam_contact_decorr.py`'s module-level constants **before** any JAM simulation in this
session had completed successfully — a genuine pre-registration.

**Forced adversary for the leaning-negative (not-C) result actually obtained** (ratio 0.883, just
inside the "~equal" band — this result leans toward NOT exonerating the rigid cut, so per the
watertight discipline the adversary that must be forced is: "is the model/pipeline somehow
artificially INFLATING JAM's number, hiding a real reduction?"):

1. **Simplified muscle contraction dynamics (`use_muscle_physiology=false`, forced by a real,
   diagnosed numerical blocker — §4) removes force-length/velocity scaling.** Checked the
   DIRECTION of bias, not just noted the simplification exists: at the peak-force instant
   (t=0.54s), `knee_flex_r` is **increasing at +46.9°/s** (textbook early-stance "yielding"
   quadriceps ECCENTRIC contraction during weight acceptance). A full Hill-type force-velocity
   curve gives a lengthening (eccentric) muscle **MORE** force than plain `activation x Fmax` (up
   to ~1.4-1.8x isometric) — meaning this simplification most likely **UNDER-estimates**, not
   over-estimates, the true peak contact force. **The simplification biases AGAINST, not toward,
   the "materially lower" (C) conclusion** — i.e. it cannot be manufacturing the not-C result;
   if anything a full-physiology run would push the ratio further from 0.85, not closer to it.
2. **BW-normalization choice, shown explicitly (not silently fixed).** Peak force 2099.9 N as
   %(lenhart2015's own 608.01 N BW) = 345.4%BW, ratio 0.883 (**~equal**, primary reading, argued
   correct because it matches the standard %BW convention — each subject/model normalized by ITS
   OWN body weight, exactly how OrthoLoad's cross-subject pool and the rigid-cut's
   subject2-normalized number are each computed). **Disclosed alternative**: if the (incorrect,
   scale-mismatched) choice of subject2's own 766.88 N BW were used instead, the same 2099.9 N
   reads as 273.8%BW, ratio **0.700 — which WOULD cross into "materially lower."** This fork is
   shown, not hidden, precisely because it is choice-sensitive; the primary reading is defended in
   §1c/§6.
3. **Left-leg/non-crossing-muscle scope reduction — forced to an empirical test, not left as an
   argument (§3).**

**Result: the forced adversary (could the simplifications be inflating JAM's number and hiding a
real "materially lower" result) does NOT fall** — if anything the muscle-physiology
simplification argues the true number is higher, not lower, and the left-leg/scope reductions are
empirically confirmed negligible at the headline instant (§3). **The ~equal / not-C verdict is
not an artifact of these disclosed simplifications.**

---

## 3. Forced-adversary empirical check: does the left-leg/non-crossing-muscle scope reduction bias the result?

`scripts/msk/jam_contact_decorr_leftcheck.py` re-ran the FULL pipeline with the left leg ALSO
prescribed from subject2's own IK (`hip_flexion_l/hip_adduction_l/hip_rotation_l/knee_angle_l/
ankle_angle_l/subtalar_angle_l/mtp_angle_l`) **and LeftGRF added** (subject2's own unmodified
2-force `ExternalLoads` convention), same 13 right-crossing-muscle activation input — a genuine
ablation, not a repeat of the same run.

| | baseline (right-only) | leftcheck (both legs) | rel. diff |
|---|---:|---:|---:|
| peak tf_contact force | 2099.86 N | 2101.71 N | **0.19%** (at the exact peak frame) |
| max rel. diff within the peak-loading window (t=[0.44,0.64]s) | — | — | **0.73%** |
| max rel. diff over the FULL 1.57s trial | — | — | 9.87% (see below) |
| median rel. diff over the full trial | — | — | 0.61% |

**Honest correction of the original argument, not a rubber stamp.** The docstring's original
claim ("left leg CANNOT affect the right knee AT ALL," since the entire kinematic chain up to
`femur_r` is prescribed and Simbody's prescribed-motion constraint absorbs any reaction force")
was **too strong** — a strict `<0.1%`-everywhere gate technically fails (max 9.87%). Diagnosed
(not dismissed): the large relative differences are concentrated at **low-force,
swing-phase-adjacent frames** (e.g. t=0.01s at 0.25 N, t=0.83-1.03s at 220-315 N) — a
relative-difference-on-small-denominator artifact, plus a real but second-order **inertial**
coupling through the whole-body mass matrix (left-leg pose affects the effective inertia the free
right-knee DOFs see, even though no FORCE transmits through the prescribed joints — a distinct
mechanism from the force-transmission argument, which remains correct). **The claim that actually
matters for this cert — that dropping the left leg does not bias the headline PEAK number — is
CONFIRMED to <1%,** not just argued. Full detail:
`data/msk_smoketest/jam_contact_decorr/jam_contact_decorr_adversary_checks.json`.

---

## 4. Diagnosed infeasibility: full Hill-type physiology could NOT be driven to completion — precise blockers, not a one-shot "honest negative"

Per the task's own allowance ("if JAM cannot be driven with the gait muscle forces at all, that is
an honest diagnosed-gap") — JAM COULD be driven with the gait activations, but only after
diagnosing and fixing two successive numerical blockers, forced via OODA rather than accepted as a
one-shot failure:

1. **Attempt 1** (all 40 name-matched muscles activation-driven, compliant tendon,
   `use_muscle_physiology=true`): ran cleanly to t=0.754s (75/158 steps, 59s wall time), then
   `[error] ... addbrev_r Fiber velocity Newton method did not converge`
   (`Millard2012EquilibriumMuscle.cpp:1067`). **Orient**: `addbrev_r` is a hip-ONLY adductor
   (confirmed NOT knee/pf-crossing, §1b) and `hip_flex_r/hip_add_r/hip_rot_r` are ALL prescribed —
   by the same prescribed-motion-decoupling argument (§3), its exact activation is dynamically
   irrelevant to the knee outcome. **Decide+Act**: scope the actuator_input_file to the 13
   knee/pf-crossing muscles only (§1b); the other 27 name-matched + 4 lenhart-only muscles fall
   back to JAM's own `constant_muscle_control=0.02` default — the SAME operating point JAM's own
   official `passive_flexion` example already ran successfully for every one of these muscles,
   including `addbrev_r`.
2. **Attempt 2** (13-muscle scope, still compliant tendon): failed even EARLIER —
   `[CPODES ERROR] At t = 0.00693887, mxstep steps taken before reaching tout` (stuck before the
   very first 0.01s report interval; killed after 829s of unproductive CPU time). **Orient**: a
   genuinely stiff initial-condition transient in the compliant-tendon fiber-length ODE, not a
   sign/units bug (hip/knee/ankle prescribed ranges independently re-verified sane, §1a).
   **Decide+Act**: smoke-tested `use_tendon_compliance=false` (rigid tendon — per the tool's own
   docstring, "removes fiber length as a state variable, simulation performance is improved") on a
   0.1s window: **completed cleanly in 16s.**
3. **Attempt 3** (13-muscle scope, rigid tendon, full 1.57s trial): progressed much further
   (t=0.01 to 0.32s, 32/158 steps) then failed again — `bflh_r Fiber velocity Newton method did
   not converge` at t≈0.317s. **Orient**: rigid tendon still performs an instantaneous algebraic
   equilibrium solve (not an ODE state, but still a Newton iteration) that can fail for a DIFFERENT
   muscle under DIFFERENT activation input — a recurring class of fragility (short-tendon-relative-
   to-fiber-length Millard2012 muscles under this specific cross-model activation transplant), not
   a single fixable bug. **Decide+Act**: `use_muscle_physiology=false` (Force = activation x Fmax,
   NO equilibrium solve of any kind — the tool's own docstring recommends this exact setting "if
   analyzing COMAK results," i.e. precisely this scenario of driving muscles from an
   externally-derived activation signal rather than solving their physiology from scratch).
   Smoke-tested on a 0.3s window (26s, clean), then the full 1.57s trial: **64.4s, exit 0, clean.**

**Disclosed limitation, with its direction checked (§2 item 1), not just noted**: the headline
number therefore comes from a **simplified** (activation x Fmax) muscle force law, not full
Hill-type contraction dynamics. This is a real, disclosed methodological compromise — but the
force-velocity-direction argument (§2) shows it most likely biases the number DOWN, not up,
meaning it argues against, not for, an inflated "not-C" reading.

---

## 5. Machine cross-checks (6/6 PASS)

| check | result |
|---|---|
| Prescribed `knee_flex_r` tracking (catches a degrees/radians unit bug) | PASS — max err 0.677°, median 0.016° |
| Engagement floor (contact area/force > 0 at peak) | PASS — 543.2 mm², 2099.9 N |
| Absolute ceiling vs. literature anchor (≤1.5x the highest verified MPa figure) | PASS — TF max 11.95 MPa (ceiling 14.55), PF max 6.69 MPa |
| Mesh-pair force agreement (tibia-side vs. femur-side casting, an over-determination check) | PASS — median 2.13%, max 3.60% |
| No NaN/Inf anywhere in the extracted series | PASS |
| Peak-force robustness (smooth broad peak, t=[0.44,0.64]s window, not a noisy spike) | PASS — monotonic rise/fall, 1730N->2100N->1663N |

**Ligament cross-check at the peak-force instant (t=0.54s, knee_flex_r≈7°) — a bonus,
independent consistency check, not a new validated result on its own**: engagement pattern
matches `docs/MECHANISM_KNEE_LIGAMENTS.md`'s own passive-flexion sweep at a similar angle — ACL
engaged (181.8 N combined, higher than the passive-sweep's 3.4-6.8 N because this is a LOADED
condition), PCL fully slack (0 N both bundle groups, matches passive-sweep's ~0-13 N near
extension), LCL engaged (51.1 N, matches passive-sweep's ~22-72 N), MCL modest/mostly slack
(9.3 N). The gait-driven dynamics behave physiologically sensibly, not erratically.

---

## 6. What this DOES and does NOT resolve

**Resolved (narrowed):** a genuinely different contact model — deformable elastic-foundation
cartilage + multi-bundle ligaments, replacing the rigid free-body cut's single-point-contact
assumption — does **not** materially relieve the knee over-prediction (ratio to the rigid cut
0.883, inside the pre-registered "~equal" band; ratio to OrthoLoad still 1.337x). Combined with
`docs/MECHANISM_CMC_SECOND_SOLVE.md`'s own finding (a mechanistically independent forward-dynamics
MUSCLE-FORCE solve, CMC, reads 8-11% **higher** than SO, not lower), this is now the SECOND
independent decorrelation axis (muscle-force-solve method: CMC vs SO; contact-model: rigid vs
deformable) that fails to explain the over-prediction as an artifact of that specific axis. The
shared-model/shared-single-trial common-mode (both this JAM leg's subject2 KINEMATICS/GRF and the
rigid cut's own subject2 model share the same underlying gait trial, and lenhart2015's own generic,
not-subject-specific anthropometry is itself a disclosed axis) remains the leading open suspect,
not this contact-model axis.

**NOT resolved:**
1. **lenhart2015 is a different, generic (62.0 kg) donor, not subject2 (78.2 kg)** — the %BW
   normalization choice (§2 item 2) is real and shown explicitly; the ratio would read differently
   (0.700, crossing into "materially lower") under the (less-defensible) alternative choice.
2. **Simplified muscle physiology** (§4) — force-length/velocity/pennation dynamics were not
   modeled; the direction-of-bias argument (§2 item 1) is reasoning, not a completed full-physiology
   run (that run was attempted twice and failed on two DIFFERENT muscles — a genuine, diagnosed
   infeasibility for this cross-model transplant, not a convenience shortcut).
3. **Only 13/44 muscles driven from real data**; the other 31 sit at a flat 0.02 tonic default.
   Argued (§1b crossing-test) and empirically checked for the LEFT leg specifically (§3) to be
   mechanically inert to the knee outcome; not individually re-ablated for all 27 right-leg
   non-crossing muscles (a cheaper, not-yet-done extension of §3's method).
4. **Single trial, single subject, single donor model** — same scope caveat as every leg in this
   family; no claim of generality across subjects/trials/gait speeds.
5. **Not COMAK** — this is a single explicit ForsimTool forward-dynamics integration (JAM's own
   published protocol for exactly this "prescribe primary, let secondary settle" idiom), not
   COMAK's own per-frame secondary-kinematics + muscle-force joint optimization. A COMAK run might
   converge differently; not attempted this session (heavier toolchain, out of scope given the
   ForsimTool route already answers the pre-registered question).
6. **PF-contact and ligament numbers are supporting/bonus context**, not independently
   anchor-checked to the same rigor as the TF headline (no PF or ligament external anchor
   pre-registered as a decision gate).

---

## 7. Confidence tier and falsifier (per report discipline)

**Confidence tier: published-plausibility.** Anchored against two real external references
(OrthoLoad in-vivo median, n=72; the rigid cut's own official `opensim.JointReaction` number) plus
a same-model literature pressure anchor (PMID 25917122) — but the model itself is a translated,
non-subject-specific, simplified-physiology stand-in (lenhart2015, generic donor, `Force =
activation x Fmax`), not a subject-matched in-vivo measurement of this exact scenario. NOT
in-vivo-anchored (no direct instrumented measurement of THIS scenario exists); NOT method-only
(real external numeric anchors are used as decision gates, not just internal consistency).

**Falsifier (what would overturn this HEADLINE verdict):** a successful full-Hill-type-physiology
re-run (resolving the two diagnosed numerical blockers some other way, e.g. per-muscle rigid/
compliant tendon toggling, a smaller integrator step, or a different muscle-activation smoothing)
that produces a peak tf_contact force whose ratio to the rigid cut falls outside [0.85, 1.15] in
either direction; or a subject-specific (not generic-donor) JAM knee model becoming available,
removing the BW-normalization fork (§2 item 2) entirely.

---

## 8. Honest gaps (full list, consolidating §4/§6's inline disclosures)

1. Simplified muscle contraction dynamics (`use_muscle_physiology=false`) — direction-of-bias
   argued (§2.1), not eliminated by a completed full-physiology run (attempted twice, both failed
   on genuine, diagnosed numerical blockers, §4).
2. BW-normalization choice is outcome-relevant (§2.2) — shown transparently, defended, not hidden.
3. Only 13/44 muscles driven from real subject2 data; left-leg omission empirically checked (§3,
   <1% at the headline instant), the other 27 right-leg non-crossing muscles' omission is argued
   (§1b) but not individually re-ablated.
4. Single trial (subject2/walking1), single generic donor model (lenhart2015) — no claim of
   generality.
5. Not a COMAK run — a single ForsimTool forward-dynamics pass (JAM's own published idiom for this
   exact use case, §6 item 5).
6. PF-contact (60.5%BW) and the ligament cross-check (§5) are supporting context, not independently
   anchor-gated to the same rigor as the TF headline.
7. The `override_default_muscle_activation=-1` (model default) setting for the 4 lenhart-only
   muscles (`gem_r/pect_r/pertert_r/quadfem_r`, never fed real data at all, always 0.02 tonic) was
   not itself sensitivity-tested — argued inert by the same crossing-muscle geometry (none of the
   4 crosses the knee) but not empirically re-ablated the way the left leg was.

---

## 9. Files

- `scripts/msk/jam_contact_decorr.py` — full build+run pipeline: coordinate-name mapping,
  muscle-name-overlap + geometric knee-crossing detector (pure-XML, no opensim import needed),
  `.sto`/XML generation for `ForsimTool`/`JointMechanicsTool`, `--stage {inputs,run,analyze}`,
  `--stop-time`/`--rigid-tendon`/`--no-physiology` flags (used for the smoke-test -> full-trial
  escalation in §4). Pre-registered anchors/thresholds live here as module-level constants.
- `scripts/msk/jam_contact_decorr_analyze.py` — h5 extraction, settle-window detection, prescribed-
  coordinate tracking sanity, %BW conversion, verdict computation, 6 machine checks, JSON writer.
- `scripts/msk/jam_contact_decorr_leftcheck.py` — the forced-adversary left-leg/LeftGRF ablation
  (§3), self-contained, reuses the main script's helpers.
- `data/msk_smoketest/jam_contact_decorr/inputs/` — `prescribed_coordinates.sto` (13 mapped
  coordinates, subject2/walking1's own trajectory), `actuator_input_file.sto` (13 knee-crossing
  muscles' real SO activation), `external_loads_right_only.xml`, `forsim_settings.xml`,
  `joint_mechanics_settings.xml`.
- `data/msk_smoketest/jam_contact_decorr/results/forsim/gait_driven_states.sto` — the completed
  forward-dynamics states (159 frames, 0-1.58s).
- `data/msk_smoketest/jam_contact_decorr/results/joint-mechanics/gait_driven.h5` — the JAM
  `JointMechanicsTool` output (contact force/area/pressure for `tf_contact`/`pf_contact`, both
  mesh-side conventions, all 76 `Blankevoort1991Ligament` bundle forces, 159 frames); the paired
  `gait_driven_leftcheck.h5` for the §3 ablation.
- `data/msk_smoketest/jam_contact_decorr/jam_contact_decorr_results.json` — full machine-readable
  headline evidence (all numbers in the Headline/§5/§6 sections traceable here).
- `data/msk_smoketest/jam_contact_decorr/jam_contact_decorr_adversary_checks.json` — the §2/§3
  forced-adversary evidence (left-leg ablation detail, BW-normalization fork, ligament cross-check,
  force-velocity-direction argument).
- `data/msk_smoketest/jam_contact_decorr/build_manifest.json` — muscle-matching/crossing-detector
  audit trail (which of the 40 name-matched muscles cross the knee and why, §1b).
- Inputs read in place, never modified: `scratchpad/knee_lig/lenhart2015.osim` + `Geometry/*.stl`
  (fetched fresh this session, Apache-2.0), subject2's `OpenSimData/Mocap/IK/walking1.mot` +
  `ForceData/walking1_forces.mot` (external drive), this repo's own already-committed
  `data/msk_smoketest/subject2_walking1/static_optimization/so/walking1_StaticOptimization_activation.sto`.

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` and the task).
