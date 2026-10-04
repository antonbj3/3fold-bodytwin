# MECHANISM KNEE LIGAMENTS — Tier-2 nonlinear passive elastic restraints (2026-07-21)

Extends `docs/MECHANISM_MSK_ELASTIC_BAND.md`'s Tier-1 tension-only `PathSpring` work with the
build charter's named Tier-2 model, `Blankevoort1991Ligament`, for the four major knee ligaments
(ACL, PCL, MCL, LCL). Script: `scripts/msk/add_knee_ligaments.py` + `scripts/msk/knee_ligament_config.json`
(the hotswappable dial, mirrors `band_config.json`'s convention). Model built on the same
Step-1-verified scaled `LaiArnoldModified2017_poly_withArms_weldHand` subject2 model used throughout
this MSK build. Run: `source_repository/.venv-msk/bin/python3 scripts/msk/add_knee_ligaments.py`.

**Isolation respected:** read-only against the OpenSim install and the two external drives (source
`.osim` mtime independently reverified unchanged: `2021-12-06 10:49:36`); nothing pushed/committed;
all new output is a NEW model file, never an in-place edit.

---

## 0. Headline honesty flag — the task's own two citation PMIDs were WRONG

Before anything else: the task asked to cite Lenhart 2015 (given as PMID 26126944) and Blankevoort
& Huiskes 1991 (given as PMID 1748601). **Both PMIDs were verified LIVE (NCBI eutils `esummary`,
not recalled) to be wrong**:

| given PMID | what it actually is (verified via NCBI eutils) |
|---|---|
| 26126944 | Caleffi ER et al., "Isolation and prebiotic activity of inulin-type fructan extracted from *Pfaffia glomerata*..." — *Int J Biol Macromol*, 2015. An unrelated food-science/prebiotics paper. |
| 1748601 | "Drug and alcohol abuse can affect your practice." — *J Am Vet Med Assoc*, 1991. A veterinary-practice-management news item. |

Correct PMIDs were found by NCBI ESearch on author+topic and confirmed via ESummary (title/journal/
authors match exactly):

- **Lenhart RL, Kaiser J, Smith CR, Thelen DG.** "Prediction and Validation of Load-Dependent
  Behavior of the Tibiofemoral and Patellofemoral Joints During Movement." *Ann Biomed Eng.*
  2015;43(11):2675-2685. doi:10.1007/s10439-015-1326-3. **PMID 25917122.**
- **Blankevoort L, Huiskes R.** "Ligament-bone interaction in a three-dimensional model of the
  knee." *J Biomech Eng.* 1991;113(3):263-9. doi:10.1115/1.2894883. **PMID 1921352.**

This is reported prominently, not silently corrected, per the watertight discipline: a citation is
a claim like any other and gets the same "verify, don't propagate" treatment as a numeric result.

---

## 1. Ligament force class used

**`opensim.Blankevoort1991Ligament`** — confirmed LIVE in this repo's stock pip `opensim==4.6` wheel
(`.venv-msk`), the same class already flagged present (but unused) in `docs/MECHANISM_MSK_ENV.md`.
**No spring-proxy fallback was needed** — this is the real Tier-2 nonlinear model, not a workaround.
It is not an `opensim-jam`-plugin-only class; the JAM plugin's own C++-only extras
(`COMAKTool`/`ForsimTool`/`JointMechanicsTool`/`Smith2018ArticularContactForce`) were used **only**
as a read-only data source (§2), never as this build's runtime.

Live-introspected class properties (`.venv-msk`, opensim 4.6-2026-06-22-85aaf64):
`linear_stiffness` (N per unit strain, default 1), `transition_strain` (dimensionless, default
0.06 — the strain at which the force-strain curve switches from quadratic "toe" region to linear),
`damping_coefficient` (default 0.003), `slack_length` (m, default 0), plus a `GeometryPath`.
Confirmed `Force`, **not** `Actuator`/`PathActuator` — same passive-element category as the
elastic-band `PathSpring`, automatically included in system dynamics.

**Machine cross-check of the force law** (not trusting the library black-box): the class's own
documented law was independently re-derived and compared against `getTotalForce()`:
- strain ≤ 0: F = 0
- 0 < strain < transition_strain (εₜ): F = k·strain²/(2εₜ) (quadratic toe region)
- strain ≥ εₜ: F = k·(strain − εₜ/2) (linear region, C¹-continuous with the toe region at strain=εₜ)

Tested on 8 bundles at `knee_angle_r=90°`: **max abs error 0.00e+00 N** (exact match to floating-point
precision) — the library's black-box output matches the published closed form exactly.

---

## 2. Data source — real vendored model, not hand-transcribed from the PDF

Rather than re-typing attachment coordinates from a paper figure/table (transcription risk), the
**actual companion model file** the Lenhart 2015 paper's own authors (Smith/Thelen labs) ship was
fetched live and parsed:

`https://raw.githubusercontent.com/opensim-jam-org/jam-resources/main/models/knee_healthy/lenhart2015/lenhart2015.osim`
(Apache-2.0). `models/knee_healthy/README.txt` (fetched, quoted in full in
`knee_ligament_config.json`) explicitly states: *"Lenhart 2015 [2] — Model first introduced and
validated against dynamic MRI"*, with reference [2] being the exact citation above — **confirming
this is the real Lenhart2015 model**, not a lookalike.

**Cross-corroboration found while searching this environment's own disk**: the JAM plugin's C++
build directory (`/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/`) already contains
real COMAK run outputs (`run_DM_ngait_og1/results/comak-inverse-kinematics/ik_constrained_model.osim`)
whose embedded `model_file` path resolves to
`.../jam-resources/models/knee_tka/grand_challenge/DM/DM.osim` — **the same jam-resources repo,
same Grand-Challenge-in-vivo-knee-loads dataset this repo's own KNEE-CELL flagship anchors against**
(`COORDINATOR.md` §5). That sibling model uses the identical `Blankevoort1991Ligament`
multi-bundle-per-ligament architecture but has **zero ACL bundles** (grep-verified across all 395
`.osim` files reachable on that drive) — consistent with ACL resection being standard in the TKA
procedure that "DM" subject underwent (it retains PCL bundles, consistent with a
cruciate-retaining implant). This is why the **healthy** `lenhart2015.osim` sibling (not the TKA one)
was used as the single, consistent source for all 4 ligaments.

**42 real fiber bundles parsed** (multi-bundle non-uniform-strain representation, the standard
approach since Blankevoort & Huiskes 1991 — ported in full, not collapsed to one bundle per
ligament, to avoid inventing a simplification the verified source doesn't use):

| bundle group | ligament | n fibers | linear_stiffness (N) | reference strain range @ LaiArnold ext. |
|---|---|---:|---:|---|
| ACLam (anteromedial) | ACL | 6 | 820 | −0.14 to −0.05 |
| ACLpl (posterolateral) | ACL | 6 | 820 | −0.12 to +0.03 |
| PCLal (anterolateral) | PCL | 5 | 900 | −0.10 to +0.03 |
| PCLpm (posteromedial) | PCL | 5 | 300 | −0.12 to −0.05 |
| MCLd (deep) | MCL | 5 | 500 | −0.04 to +0.04 |
| MCLs (superficial) | MCL | 6 | 500 | +0.029 to +0.047 |
| MCLp (posterior oblique) | MCL | 5 | 500 | −0.02 to +0.05 |
| LCL | LCL | 4 | 600 | +0.06 (all 4) |

Patellofemoral/capsular/ITB structures also present in the source (PT, PFL, lPFL, mPFL, pCAP, ITB)
were **excluded** — out of scope for "the four major knee ligaments" the task asked for.

### Geometric pipeline (attachment points → LaiArnold)

1. Source points are defined on `femur_distal_r`/`tibia_proximal_r` frames (rigid, zero-rotation
   `WeldJoint` offsets from `femur_r`/`tibia_r` in `lenhart2015.osim`). Converted into `femur_r`/
   `tibia_r`-local coordinates via OpenSim's own `Frame.findStationLocationInAnotherFrame` — loaded
   after stripping 2 JAM-only force classes not registered in stock OpenSim
   (`Smith2018ArticularContactForce`, `Smith2018ContactMesh`; the ligaments themselves needed no
   stripping). **Cross-validated**: manual offset arithmetic from the raw `WeldJoint` XML vs. the
   live API transform matched to 6 decimal places.
2. **Segment scale factors** (isotropic, per-segment): ratio of LaiArnold's own hip-to-knee (femur)
   and knee-to-ankle (tibia) joint-center airline distances to the same quantities measured in
   `lenhart2015.osim`, each at its own default/neutral pose:

   | segment | Lenhart2015 (m) | LaiArnold r (m) | LaiArnold l (m) | scale r | scale l |
   |---|---:|---:|---:|---:|---:|
   | femur (hip→knee) | 0.37422 | 0.49112 | 0.49184 | **1.3124** | **1.3143** |
   | tibia (knee→ankle) | 0.40333 | 0.49896 | 0.49120 | **1.2371** | **1.2179** |

   LaiArnold's subject2 has a notably larger frame than Lenhart2015's MRI subject ("a young healthy
   adult female" per the source README) — a real, substantial (~22-31%) rescale, not a rounding
   tweak. Small L/R asymmetry (femur 1.3124 vs 1.3143, tibia 1.2371 vs 1.2179) is genuine subject
   anthropometry (same phenomenon `MECHANISM_MSK_ELASTIC_BAND.md` found in the heel markers).
3. **Slack length is NOT naively rescaled** by a segment factor (a ligament spans two differently-
   scaled segments, so no single segment's scalar is correct). Instead the **dimensionless
   reference strain** measured in the source model at its own reference pose is preserved, and
   `Blankevoort1991Ligament.setSlackLengthFromReferenceStrain()` (OpenSim's own API, not a
   hand-rolled formula) re-derives `slack_length` in LaiArnold's geometry. Cross-checked against an
   independently-precomputed value (plain numpy, done during config generation): **max abs error
   4.23e-06 m across all 84 bundles** — PASS.
4. `linear_stiffness` (N per unit strain, a material property) is left **unscaled** — no
   cross-sectional-area data available to justify rescaling it. Flagged as a known first-step
   simplification.
5. Left-side (`_l`) bundles mirror the right-side source by negating the local Z coordinate —
   verified live as the correct convention for this model (`femur_r` origin
   (−0.05869, 0.84814, +0.09874) vs `femur_l` (−0.05869, 0.84814, −0.09874): exact Z mirror, X/Y
   identical).

Result: **84 `Blankevoort1991Ligament` bundles** (42 fibers × 2 sides) added to the model. Force
count before/after: 93 → **177** (93 = the exact muscle/actuator count independently reported in
`docs/MECHANISM_MSK_ENV.md`, confirming nothing else changed).

---

## 3. Tension-only verification (task requirement — MEASURED, not assumed)

Every sample in the full flexion sweep (§4) — **2982 (bundle × angle) checks** across all 42 right-side
bundles × 71 angles — asserts: if `strain ≤ 0` then `getTotalForce() == 0.0` exactly (not
"approximately"; the class's own quadratic/linear law is identically zero for non-positive strain,
and this was checked bit-for-bit).

**Result: 0 violations out of 2982 checks → PASS.** Repeated on the left-knee spot-check (10 angles
× 42 bundles = 420 more checks): also 0 violations.

A live, isolated gotcha was found and fixed while building this: calling
`setSlackLengthFromReferenceStrain()` changes a `Property`, but the previously-built `State`'s
cached length/strain/force values are **not** invalidated by a bare `model.realizePosition(state)`
afterward — reproducing this bug and then fixing it (rebuild the state via `model.initSystem()`
again) was done explicitly before trusting any number in this document; the same re-init pattern
was needed a second time after `model.finalizeConnections()`/`printToXML()`. Both are called out
in `add_knee_ligaments.py`'s comments so a future editor doesn't reintroduce them.

---

## 4. Passive knee flexion sweep — force-vs-flexion curve per ligament

Right knee, `knee_angle_r` swept across its full model range `[0°, 140°]` (0°=full extension) in
2° steps, all other coordinates held at 0 (neutral standing reference). Values are the **summed
tension of all bundles in that ligament group** (N):

| angle (deg) | ACL (N) | PCL (N) | MCL (N) | LCL (N) |
|---:|---:|---:|---:|---:|
| 0 | 6.83 | 13.50 | 91.42 | 72.00 |
| 10 | 3.40 | 3.27 | 21.91 | 22.19 |
| 20 | 1.84 | 0.91 | 10.54 | 4.98 |
| 30 | 0.68 | 0.39 | 5.29 | 0.86 |
| 40 | 0.01 | 0.40 | 2.52 | 0.00 |
| 50 | 0.00 | 9.19 | 3.89 | 0.00 |
| 60 | 0.00 | 38.78 | 46.56 | 0.00 |
| 70 | 0.00 | 96.29 | 158.42 | 0.00 |
| 80 | 0.00 | 171.03 | 302.20 | 0.00 |
| 90 | 0.00 | 246.50 | 454.46 | 0.00 |
| 100 | 0.00 | 310.14 | 607.82 | 0.00 |
| 110 | 0.00 | 353.41 | 756.17 | 0.00 |
| 120 | 0.00 | 370.73 | 895.85 | 0.00 |
| 130 | 1.10 | 361.34 | 1027.00 | 0.00 |
| 140 | 5.91 | 334.87 | 1155.34 | 0.00 |

### Pre-registered thresholds (set BEFORE reading the sweep result)

- **ACL tighter at extension**: mean tension over `[0°,15°]` vs `[125°,140°]`.
  **4.318 N vs 2.143 N → PASS.**
- **PCL tighter at flexion**: mean tension over the same bands, reversed direction.
  **6.180 N (extension) vs 355.735 N (flexion) → PASS.**

Both engagement patterns match the textbook expectation the task named (ACL near extension, PCL
near deep flexion) and were checked against a threshold fixed before the numbers were read, not
fit after the fact. LCL is essentially fully engaged at extension (72.0 N, all 4 bundles at their
published +6% reference strain) and fully slack from ~40° onward — the well-established, least
controversial pattern in the knee-ligament literature (LCL is a near-pure extension stabilizer).
Left-knee spot check (10 angles) reproduces the same qualitative pattern (e.g. ACL 6.83→0→4.02 N,
PCL 13.50→381.57 N, confirms bilateral symmetry of the build, not just the right side).

---

## 5. Honest quantitative caveat — deep-flexion magnitudes are qualitative, not validated

Symmetric QC on the above table: MCL and PCL grow very large at deep flexion (MCL 1155 N at 140°).
Rather than presenting this as a validated prediction, the actual per-bundle **strain** (not just
force) was measured across the sweep to check whether this is a broad physiological trend or a
kinematic-mismatch artifact:

| angle (deg) | max\|strain\| ACL | max\|strain\| PCL | max\|strain\| MCL | max\|strain\| LCL |
|---:|---:|---:|---:|---:|
| 0 | 0.140 | 0.120 | 0.050 | 0.060 |
| 40 | 0.138 | 0.041 | 0.132 | 0.035 |
| 80 | 0.132 | 0.133 | 0.238 | 0.118 |
| 100 | 0.145 | 0.182 | 0.329 | 0.156 |
| 120 | 0.134 | 0.209 | 0.463 | 0.184 |
| 140 | 0.165 | 0.230 | **0.582** | 0.188 |

**Reading this honestly**: ACL and LCL stay in a modest, plausible strain band (≤0.19) across the
*entire* range. **MCL exceeds commonly-cited collateral-ligament physiological/pre-injury strain
ranges (~15-20%) beyond roughly 90-100° of flexion, reaching 58% strain (non-physiological — real
MCL tissue would fail well before this) at 140°.** PCL is borderline (23% at 140°, high but less
extreme). The single worst MCL bundle carries up to 276 N alone at 140° (top-2-bundles = 45-49% of
the group total — a broadly-elevated group with some concentration, not one rogue outlier).

**Mechanism, not a bug**: this is the expected consequence of the §2 honesty flag — Lenhart2015's
own knee joint is a **free 6-DOF `knee_r`** (independent `knee_flex_r`/`knee_add_r`/`knee_rot_r` +
3 translations, confirmed live in its XML) whose secondary-DOF path *emerges* from ligament+contact
force balance. LaiArnold's `walker_knee_r` is instead a **prescribed 1-DOF** `CustomJoint`
(confirmed live: `nCoords=1`, range `[0°,140°]`) with ab-adduction/rotation/all 3 translations
riding fixed polynomial functions of `knee_angle_r` alone — a kinematic surrogate tuned for
LaiArnold's own (ligament-free) generic skeleton, not for these specific ported ligaments. Grafting
real attachment points onto a *different* prescribed coupling has no guarantee the two stay
consistent at the extremes of range of motion; they visibly diverge beyond ~100-110°. The
force-law and slack-length cross-checks (§1, §3) rule out an implementation bug (both PASS to
<1e-4); this is a geometric/kinematic-coupling limitation, not a code defect.

**Practical read**: treat the ENGAGEMENT PATTERN (which ligament dominates at which flexion angle,
§4's PASS/PASS thresholds) as the validated result. Treat ABSOLUTE force/strain magnitudes
quantitatively only below ≈100° flexion; beyond that, read qualitatively (order-of-magnitude,
direction of change) only.

---

## 6. Valgus/varus (MCL/LCL) structural + finite-difference check

**Honesty flag restated**: `walker_knee_r` has no free ab-adduction coordinate to sweep (§5), so a
direct free-DOF valgus/varus stress test is not mechanically possible in this joint formulation
(the Lenhart2015/DM source models this data comes from — with their free 6-DOF knees — *could* do
this directly; that is a documented advantage of the COMAK/JAM approach over LaiArnold's prescribed
surrogate). Two alternative, geometry-derived checks were run instead:

**(a) Structural medial/lateral sign** (right knee, femur-local frame): MCL femoral-attachment
centroid Z = **−0.0512 m**, LCL femoral-attachment centroid Z = **+0.0457 m** → opposite sides of
the femur's own centerline (right leg: −Z is medial/toward midline, +Z is lateral, confirmed via
the femur_r/femur_l mirror check in §2) → **PASS**, the geometric precondition for MCL/LCL to
resist opposite frontal-plane perturbations is satisfied.

**(b) Finite-difference virtual-perturbation test** (a "shadow" calculation outside the OpenSim
state, since ab-adduction isn't a free coordinate): at each test angle, a Grood-Suntay-style
floating axis is derived **from the model's own measured geometry** (flexion axis via small-angle
Rodrigues extraction of the actual femur↔tibia coupled rotation; medio-lateral axis from the
MCL↔LCL femoral-attachment separation, Gram-Schmidt-orthogonalized against the flexion axis;
perturbation axis = their cross product — zero hardcoded anatomical convention assumed). The tibia
is rigidly rotated ±2° about this axis through the centroid of the 4 groups' tibial attachments,
and each bundle's length is recomputed by direct 3D point geometry (baseline length cross-checked
against the OpenSim path's own `getLength()`: exact match to 1e-6 m before trusting the perturbed
values):

| angle (deg) | dLen/dΔ ACL (m/rad) | dLen/dΔ PCL | dLen/dΔ MCL | dLen/dΔ LCL |
|---:|---:|---:|---:|---:|
| 0 | +0.0018 | +0.0548 | **−0.4070** | **+0.2346** |
| 45 | −0.0538 | +0.0811 | **−0.3793** | **+0.2352** |
| 90 | −0.0293 | +0.0927 | **−0.3504** | **+0.2406** |

**MCL and LCL have opposite-sign slopes at all 3 test angles → PASS.** Under this virtual
perturbation, one side of the joint always lengthens while the other shortens — the structural
signature of "one resists opening on side A, the other resists opening on side B" — i.e. MCL and
LCL are confirmed to be frontal-plane antagonists by direct geometric measurement, even though the
specific real-world valgus/varus SIGN convention was not independently derived here (that would
require anchoring the axis to a named anatomical direction, out of scope for this structural
check).

---

## 7. Stability / knee-contact-force contribution

**Scope honesty**: LaiArnold has no articular-contact force element (unlike the JAM/COMAK source
models, which pair these same ligaments with `Smith2018ArticularContactForce` cartilage contact).
This section reports the **ligament-side net restraint force vector** — the load a real contact
model would need to react against — not a coupled contact-mechanics simulation.

At each representative pose, the summed bundle-tension vectors (only bundles with positive tension
counted) are projected onto the tibia's own instantaneous axes (proximal-distal ≈ compression/
distraction; anterior-posterior ≈ drawer shear; medial-lateral ≈ frontal-plane):

| knee_angle_r | group | ΣT (N) | proximal-distal (N) | AP shear (N) | ML (N) |
|---:|---|---:|---:|---:|---:|
| 0° | ACL | 6.83 | +4.49 | −4.85 | +1.72 |
| 0° | PCL | 13.50 | +10.77 | +5.49 | −5.94 |
| 0° | MCL | 91.41 | +83.45 | −2.61 | −33.53 |
| 0° | LCL | 72.00 | +70.38 | +10.69 | −10.72 |
| 30° | ACL | 0.68 | +0.35 | −0.56 | +0.17 |
| 30° | MCL | 5.29 | +4.75 | −1.24 | −1.97 |
| 60° | PCL | 38.78 | +31.92 | +20.44 | −7.25 |
| 60° | MCL | 46.56 | +30.42 | +30.73 | −17.15 |
| 90° | PCL | 246.50 | +216.58 | +103.37 | −47.79 |
| 90° | MCL | 454.46 | +291.51 | +321.76 | −131.15 |

**Reading this**: at full extension, all four ligaments contribute a net *compressive* (proximal-
distal, i.e. pulling the tibia toward the femur) component — consistent with the classic
"screw-home"/full-extension "locked knee" stability mechanism where all major restraints
co-tighten. Through flexion, ACL/LCL unload (per §4/§6) while PCL and MCL take over as the dominant
AP-shear and compressive restraints respectively — i.e., **this build reproduces the well-known
handoff from extension-dominant (ACL+LCL) to flexion-dominant (PCL+MCL) restraint**, purely as a
consequence of the ported geometry + LaiArnold's own prescribed kinematics, not asserted by
construction. A properly coupled contact model (COMAK/JAM, already available read-only in this
environment at `/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/opensim-build/opensim-cmd`)
would let tibiofemoral contact force respond to this restraint load self-consistently — flagged as
the natural next step, not attempted here (out of scope; would require the JAM C++ toolchain, not
this venv's stock Python `opensim`).

---

## 8. Honest gaps (full list)

1. **PMID correction (§0)** — the task's supplied citations were wrong; corrected and flagged, not
   silently fixed.
2. **Literature-scaled, not subject-specific** — LaiArnold subject2 has no knee MRI/CT. Attachment
   points are Lenhart2015's real healthy-subject geometry, isotropically rescaled per segment
   (§2) — a first-step per `docs/MECHANISM_MSK_BUILD_CHARTER.md`'s own doctrine, not the ceiling.
3. **Deep-flexion magnitudes are not quantitatively validated** (§5) — MCL strain becomes
   non-physiological (>50%) beyond ~120° due to a genuine kinematic-coupling mismatch between
   LaiArnold's prescribed 1-DOF knee and Lenhart2015's free 6-DOF source knee. Quantitative
   trust bound: ≲100° flexion; qualitative (pattern-only) beyond that.
4. **No free ab-adduction DOF in LaiArnold's knee** (§6) — valgus/varus verified via a geometric
   shadow-perturbation construction, not a free-coordinate stress test; the sign convention for
   "which direction is true anatomical valgus" was not independently anchored.
5. **No articular contact element** (§7) — ligament restraint forces are reported standalone;
   actual tibiofemoral contact-force redistribution would need the JAM/COMAK C++ toolchain
   (available read-only in this environment, not invoked here).
6. **linear_stiffness left unscaled** (§2) — a material property; no cross-sectional-area data
   existed to justify a physically-motivated rescale onto LaiArnold's larger frame.
7. **Single source subject, single source model family** — all 4 ligaments come from one healthy
   subject's model (Lenhart2015); no attempt was made to blend/average across multiple published
   knee-ligament parameter sets (e.g. a second, independent Blankevoort & Huiskes 1991-style
   parameter table) — a decorrelated second source is a natural hardening step, not done here.
8. **Static/passive sweep only** — no muscle activation, external load, or dynamic (velocity-
   dependent damping) case was exercised; `damping_coefficient` (0.003, source-model default) is
   present in every bundle but inactive in this zero-velocity sweep.

---

## 9. Files

- `scripts/msk/knee_ligament_config.json` — the hotswappable dial: full provenance block (citations,
  PMID-correction note, scaling method, scale factors) + all 84 bundle records (name, ligament
  group, side, attachment bodies/local coordinates, stiffness, transition strain, damping,
  reference strain, cross-check slack length).
- `scripts/msk/add_knee_ligaments.py` — build + verification script (model build, slack-length and
  force-law cross-checks, tension-only sweep, engagement-threshold checks, valgus/varus
  perturbation test, stability-contribution report, JSON summary). Two live OpenSim API gotchas
  found and fixed are documented inline (stale cached strain/force after
  `setSlackLengthFromReferenceStrain`; stale state after `finalizeConnections`/`printToXML`) —
  both required rebuilding the `State` via `model.initSystem()` again.
- `data/msk_smoketest/knee_ligaments/model_with_knee_ligaments.osim` — the new model (LaiArnold +
  84 ligament bundles; original source file untouched, mtime-verified).
- `data/msk_smoketest/knee_ligaments/knee_ligament_results.json` — machine-readable results (full
  sweep table, cross-checks, thresholds, valgus/varus test, stability contribution).
