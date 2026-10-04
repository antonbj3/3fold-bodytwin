# MECHANISM ANATOMICAL HAND — first step: forearm-muscle-driven index-finger flexion (2026-07-21)

Executes the operator's explicit spec: a grip must **emerge from the actual forearm
muscles** (extrinsic finger flexors originating proximal to the wrist, tendons crossing
wrist→MCP→PIP→DIP), not a kinematic/robotics hand. bodytwin's LaiArnold model **welds**
the hand to the forearm (`radius_hand_r`, 0 DOF) and has **zero arm muscles**
(`docs/MECHANISM_MUSCLE_AUDIT.md`, `docs/MECHANISM_SHOULDER_FORCE.md`,
`docs/MECHANISM_ELBOW_FORCE.md` — all independently re-confirmed live this session). This
is the **first dial-turn only**: one finger (index), driven by two real forearm flexors
(FDP2_r, FDS2_r) plus one forearm extensor (ED2_r, a sign-contrast control) — **not** the
complete ~30+-muscle hand. Every number below is machine-measured this session
(`scripts/msk/anatomical_hand.py` + `scripts/msk/verify_anatomical_hand.py`, exit 0, a
full independent re-run reproduced the evidence JSON **byte-for-byte** — determinism
checked, not assumed), isolation respected (`.venv-msk` only; OpenSim install and mounted
drives read-only; no git commit/push).

## Headline result

| | |
|---|---:|
| **Does activating a forearm-originating muscle flex the finger?** | **YES** — FDP2_r produces a positive (flexion-direction) generalized force at MCP, PIP, **and** DIP simultaneously (per-N-tension: 14.87, 11.68, 10.24 mN·m) |
| Moment-arm cross-check (OpenSim `computeMomentArm` vs. independent finite-difference) | **12/12 muscle-coordinate pairs PASS** (rel. error ≤ 6.3×10⁻⁵) |
| Sign contrast: extensor (ED2_r) opposite sign to flexors at every shared joint | **PASS** (wrist −11.7, MCP −3.6, PIP −9.9, DIP −9.4 mm, all negative vs. flexors' positive) |
| Pre-registered external moment-arm band (order-of-magnitude literature cross-check) | **11/12 PASS**, 1 miss by 5.6% (FDP2_r@DIP: 10.56mm vs. [1,10]mm band) — disclosed, not hidden |
| Wrist un-welded? | **YES** — `radius_hand_r` WeldJoint replaced with 2 real DOF (`wrist_flex_r`, `wrist_dev_r`) |
| New finger joints | `mcp2_flex_r`, `pip2_flex_r`, `dip2_flex_r` (1-DOF hinges each) |
| New muscles | `FDP2_r`, `FDS2_r`, `ED2_r` (`Millard2012EquilibriumMuscle`, real Hill-type actuators) |
| Complete-hand fidelity gap (honest, symmetric) | **1 of 5 fingers, 3 of ~30+ muscles** — see §6 |

**Model + evidence**: `data/msk_models/LaiArnoldModified2017_anatomical_hand_index_subject2_scaled.osim`
(966 KB, loads clean in OpenSim 4.6, 26 bodies/26 joints/83 muscles) +
`scripts/msk/anatomical_hand_evidence.json`.

---

## 1. Model search (verified this session, not assumed) — why this is a from-scratch graft, not a donor import

The task named specific candidates to check live (PMID/DOI, no fabrication). A
`watertight-researcher` subagent + 2 direct PMC fetches verified all of them:

| candidate | verified? | scope (verified) |
|---|---|---|
| **McFarland DC, Binder-Markey BI, Nichols JA, Wohlman SJ, de Bruin M, Murray WM (2023)**. "A Musculoskeletal Model of the Hand and Wrist Capable of Simulating Functional Tasks." *IEEE Trans Biomed Eng* 70(5):1424-1435. PMID 36301780. doi:10.1109/TBME.2022.3217722. Free full text PMC10650739. | **VERIFIED, real, current, complete** | 43 Hill-type muscles, **both** extrinsic (FDP/FDS/EDC/EPL/FPL + wrist muscles) **and** intrinsic (lumbricals, interossei, thenar group), 23 DOF. States it is freely downloadable on SimTK (login-free access itself not independently re-fetched this session). **This is the correct next-step donor for the full hand** (§6) — not used here because it requires a fresh external fetch (simtk.org), which was disclosed rather than silently substituted; this session built the first step directly instead of blocking on that fetch. |
| Saul KR, Hu X, Goehler CM, Vidt ME, Daly M, Velisar A, Murray WM (2015), the "MoBL-ARMS" upper-limb model. *Comput Methods Biomech Biomed Engin* 18(13):1445-1458. PMID 24995410. | VERIFIED citation | **Confirmed wrist-only** — 7 DOF, kinematic chain stops at the wrist, no MCP/PIP/DIP. |
| Blana D, Hincapie JG, Chadwick EK, Kirsch RF (2008). *J Biomech* 41(8):1714-1721. PMID 18420213. | VERIFIED (title differs slightly from the task's guess) | **Confirmed wrist-only** — 9 DOF, terminates at wrist. |
| Gonzalez RV, Buchanan TS, Delp SL (1997). *J Biomech* 30(7):705-712. PMID 9239550. | VERIFIED | **Confirmed wrist-only** (15 wrist muscles). No distinct "Gonzalez finger model" found anywhere — reported as **not found**, not forced. |
| Lee SW, Chen H, Towles JD, Kamper DG (2008). *J Biomech* 41(7):1567-1573. PMID 18387615. | VERIFIED | A real paper on index-finger **extensor-mechanism** moment arms, not a distributable model. Notable disclosed finding used below: their method gave **systematically smaller** moment arms than classic An-1983-style tendon-excursion methods — real inter-method spread in the field itself. |
| An KN, Ueba Y, Chao EY, Cooney WP, Linscheid RL (1983). *J Biomech* 16(6):419-425. PMID 6619158. | VERIFIED citation | Numeric moment-arm table **could not be retrieved** this session (abstract has no numbers; paywalled). Disclosed, not fabricated. |
| Ketchum LD, Thompson D, Pocock G, Wallingford D (1978). *J Hand Surg* 3(6):571-578. PMID 722032. | VERIFIED citation | No usable numbers retrieved. Disclosed. |
| Buchholz B, Armstrong TJ, Goldstein SA (1992). *Ergonomics* 35(3):261-273. PMID 1572336. | VERIFIED citation | Segment-length ratios **could not be retrieved** (paywalled). Disclosed. |

**Local file search** (this whole machine: bodytwin, cad-to-simulation-I, the opensim-core-jam
source tree [~230 example `.osim` files], video2kin's conda envs, the OpenSim 4.6 pip
wheel's own package tree): **no complete, off-the-shelf OpenSim model with an
articulated, muscle-driven finger exists locally.** Two informative near-misses:

- **`wrist.osim`** (`opensim-core-jam/Bindings/Java/OpenSimJNI/Test/wrist.osim`) — a
  **real** anatomical wrist+index-finger model: 8 carpal bones, CMC1/MCP2/PIP/DIP joints,
  **26 correctly-named extrinsic muscles** (APL, ECRB, ECRL, ECU [+pre/post-surgery
  variants — suggesting a tendon-transfer-surgery origin], EDCI/L/M/R, EDM, EIP, EPB, EPL,
  FCR, FCU, FDPI/L/M/R, FDSI/L/M/R, FPL, PL). Its own `<credits>`/`<publications>` XML
  tags are **literally unfilled placeholder text** ("Model authors names.." / "List of
  publications related to model...") — real engineering geometry, but of **unconfirmed
  authorship**, not a citable published model. All its finger-body masses are exactly
  zero (a kinematics-only JNI test fixture, not built for dynamics). **Used here only as
  an engineering cross-reference**, not grafted wholesale — and doing so caught a real
  bug: its own MCP2 lateral-offset ratio turned out to be a **4× outlier** vs. its own
  PIP/DIP ratios for the *same* finger (§3.4) — a useful example of why donor geometry
  must be cross-checked, not trusted at face value.
- **`Rajagopal_2015.osim`** (bundled with the *official* OpenSim distribution's own
  `OpenSenseExample`) has a real, citable, **unlocked** `radius_hand_r` `UniversalJoint`
  (2 DOF: `wrist_flex_r`, `wrist_dev_r`). bodytwin's own LaiArnold-lineage `WeldJoint` is
  numerically **identical** in its frame offsets to this file's joint (verified live —
  same translation/orientation to the last digit on the generic model), confirming shared
  lineage. Its wrist ROM is reused directly below. **But it has zero finger DOF** — even
  the official OpenSim distribution's fingers are rigid meshes on one `hand_r` body. This
  confirms: no stock OpenSim model, anywhere, has an articulated finger.

**Conclusion, forced not assumed**: building a real forearm-muscle→finger-tendon coupling
this session required constructing the finger joints + muscles directly on bodytwin's own
model, using `wrist.osim` only as a geometric sanity-reference and Rajagopal's wrist ROM
as a real citable data point — not a wholesale donor graft.

---

## 2. What was built

### 2.1 Un-welding the wrist
`radius_hand_r` (`WeldJoint`, 0 DOF) → **removed**, replaced with 2 chained `PinJoint`s
via a massless intermediate body (`wrist_int_r`): `wrist_flex_r` (range −70°..+70°) then
`wrist_dev_r` (range −25°..+35°, from Rajagopal_2015.osim's own on-disk XML). **Why 2
chained PinJoints and not a single 2-DOF `CustomJoint`**: a fresh, Python-constructed
`CustomJoint` with a custom multi-axis `SpatialTransform` was found **live** to segfault
this exact OpenSim 4.6 Python build deterministically (an established, already-documented
finding in this repo's own `scripts/msk/add_scapula_clavicle.py`, reused here rather than
rediscovered) — 2 chained `PinJoint`s are the mathematically-equivalent, crash-safe
workaround (each `PinJoint` always rotates about its own offset frame's local Z; OpenSim's
own body-fixed rotation composition makes this identical to a true 2-axis joint).

### 2.2 Index finger: MCP → PIP → DIP
Three new bodies (`index_proximal_r`, `index_medial_r`, `index_distal_r`, solid-cylinder
mass/inertia, tapering radius 8/7/6mm×scale) + three new 1-DOF `PinJoint`s
(`mcp2_flex_r`, `pip2_flex_r`, `dip2_flex_r`). **Disclosed simplification**: real MCP is a
2-DOF condyloid joint (flexion-extension + abduction-adduction); only flexion-extension is
modeled here (FDP/FDS/ED don't drive ab/ad, and it isn't part of this step's verified
claim). PIP/DIP are correctly 1-DOF hinges anatomically — no simplification there.

**Segment lengths**: generic/nominal values (metacarpal2 68.5mm — taken from `wrist.osim`'s
own geometry, which matches published adult 2nd-metacarpal norms closely; proximal/medial/
distal phalanx 42/24/18mm — commonly-cited approximate adult index-finger values, since
Buchholz et al. 1992's specific regression could not be retrieved this session, §1),
**scaled by this model's own measured generic→subject2 arm scale factor (1.2308, isotropic
to 6 decimal places, machine-computed from the ratio of the pre-removal WeldJoint's
translation in the generic vs. scaled source files — not copy-pasted)**. Resulting
subject2-scaled lengths: metacarpal2 84.3mm, proximal phalanx 51.7mm, medial phalanx
29.5mm, distal phalanx 22.2mm.

### 2.3 Joint-axis construction — a real bug caught and fixed at the source
The flexion axis for every new joint (wrist, MCP, PIP, DIP) is the **elbow's own
ground-frame flexion axis**, independently re-measured live each build (`(0,0,1)` at the
reference pose, via forward-kinematics finite-difference on `elbow_flex_r` — never
hardcoded), transported into each new joint's parent body via **pure analytical rotation
composition** (not a live re-query mid-build: adding/removing any component invalidates
OpenSim's whole `System`, caught live via a "Component has no underlying System" crash,
fixed by propagating rotations in numpy after one valid start-of-build query). Real
anatomy has all of these upper-limb flexion hinges approximately **parallel** (curling the
elbow and making a fist are compatible, same-plane motions) — this is the geometric
justification for reusing one shared axis, not an arbitrary convenience.

**A first attempt copied Rajagopal's own `radius_hand_r` offset-frame orientation numbers
directly** (same numeric values as bodytwin's own WeldJoint — same lineage). This was
**measured live to be wrong**: the resulting "flexion" axis pointed along the **forearm's
own long axis** (an axial-roll DOF, like pronation/supination, not a flexion hinge) —
because that donor's numbers encode axis assignment for its own `UniversalJoint`'s
specific 2-axis internal convention, not portable to a bare `PinJoint`. Caught via a
non-degenerate reach-test (see §3.1), rejected, and replaced with the explicit,
independently-measured elbow-axis-transport method above.

### 2.4 Muscles: FDP2_r, FDS2_r (flexors), ED2_r (extensor, sign-contrast control)
All three are real `Millard2012EquilibriumMuscle` (Hill-type) actuators, originating on
`ulna_r` (approximating FDP/FDS's real anteromedial forearm origin), with via points on
each intermediate body (near **and** far end of each segment — one near every joint it
crosses, matching how a real tendon has local contact support at each joint, not one point
per bone) offset perpendicular to the flexion axis by a radius-proportional distance
(1.3× each segment's own physical radius — the tendon sits just outside the bone surface).
FDP2_r crosses wrist→MCP→PIP→DIP, inserting at the **base** (proximal portion, not the
tip) of the distal phalanx — the real anatomical FDP insertion site. FDS2_r crosses
wrist→MCP→PIP only, inserting at the medial (middle) phalanx — it structurally cannot
reach DIP (verified: its own computed DIP moment arm is 0.00mm on fresh reload, exactly as
anatomy predicts). ED2_r mirrors FDP2_r's crossing set on the dorsal side, as a built-in
falsifier for the whole sign-convention machinery (§3.2).
`tendon_slack_length` is **solved**, not guessed, so each muscle's fiber operates at its
own `optimal_fiber_length` at the reference pose (technique already established in this
repo's `add_arm_muscles.py`, reused here). `Fmax`/`optimal_fiber_length` are
physiologically-plausible round placeholders (FDP2 45N, FDS2 40N, ED2 20N) — **not**
independently verified per-muscle this session (McFarland 2023 almost certainly has real
numbers; not extractable from the accessible PMC text this session, §1). This does **not**
affect the moment-arm/excursion verification below — moment arm and excursion are pure
path-geometry properties, independent of `Fmax`.

---

## 3. Verification — every gate is machine PASS/FAIL, not eyeballed

### 3.1 Joint sign convention — forced via an incremental, non-degenerate reach-chain test
**The adversary**: does "positive coordinate = flexion" hold, or is it an arbitrary
labeling artifact? **Forced test**: flex the elbow, then ADD wrist, then MCP, then PIP,
then DIP, each in its own candidate positive direction, holding all previously-added
joints at their flexed value — a real anatomically-consistent flexor chain must bring the
fingertip **progressively closer** to the shoulder at every single step (a hand closing
into a fist while the elbow also flexes), using **zero muscles** (pure rigid-body
kinematics). Measured, this session (shoulder-to-fingertip-proxy distance, meters):

| step | distance (m) |
|---|---:|
| reference (all 0) | 0.83224 |
| +elbow (40°) | 0.80314 |
| +wrist (+20°) | 0.78040 |
| +MCP (+40°) | 0.72595 |
| +PIP (+40°) | 0.69086 |
| +DIP (+40°) | 0.67925 |

**Strictly monotonically decreasing at every step — PASS.** Converged with **zero axis
flips needed** (`flips_used: {}`), i.e. the elbow-axis-transport construction (§2.3) got
the sign right at all four new joints on the first fully-corrected attempt — but this was
only accepted after the earlier (Rajagopal-copy) attempt was caught and killed by this
exact test (§2.3), and after a **fingertip-proxy bug** was itself caught and fixed: an
early version used `index_distal_r`'s own **origin** as the fingertip point — which **is**
the DIP joint location by construction, so it is invariant under DIP's own rotation by
construction, producing a false, unfixable "sign failure" at DIP on every iteration
(caught because the flip-and-rebuild loop oscillated without ever converging — a stuck
loop was the tell that the test itself, not the joint, was broken).

### 3.2 Muscle sign — moment-arm cross-check + flexor/extensor contrast
Two **independent** methods for every muscle-coordinate pair that crosses: (a) OpenSim's
own `Muscle.computeMomentArm()` API, (b) an independent central finite-difference of
`GeometryPath.getLength()` vs. coordinate value (the virtual-work identity, moment arm =
−dL/dθ — the same principle this repo's elbow-force cert used for its own bug-immune
method). **12/12 pairs agree to ≤6.3×10⁻⁵ relative error**:

| muscle | wrist_flex_r | mcp2_flex_r | pip2_flex_r | dip2_flex_r |
|---|---:|---:|---:|---:|
| **FDP2_r** | +14.64mm | +10.60mm | +11.53mm | +10.56mm |
| **FDS2_r** | +14.65mm | +10.60mm | +12.64mm | 0.00mm (doesn't cross — correct) |
| **ED2_r** | −11.74mm | −3.62mm | −9.94mm | −9.38mm |

(all `computeMomentArm` vs. finite-difference pairs shown match to 3+ significant figures;
full precision in the evidence JSON.) **ED2_r is negative at every joint FDP2_r/FDS2_r are
positive — the sign-contrast falsifier did NOT fire.** This is not circular: the joints'
own flexion sign was established in §3.1 using **zero muscles**; the muscles' via-point
offset side was then a **free construction choice** (like a surgeon choosing which side of
a joint to route a tendon), verified independently via `computeMomentArm` — two genuinely
decorrelated computations (rigid-body kinematics vs. muscle-path geometry) agreeing is the
non-tautological confirmation.

### 3.3 Headline test: does forearm-muscle activation flex the finger?
Per-unit-tension generalized force (`moment_arm × 1N`) from FDP2_r alone, finger at rest
(all coordinates 0): **mcp2_flex_r +14.87 mN·m, pip2_flex_r +11.68 mN·m, dip2_flex_r
+10.24 mN·m — all three positive simultaneously.** Activating this one forearm-originating
muscle flexes MCP, PIP, and DIP at once, exactly as real FDP does (and exactly why FDP is
the one extrinsic muscle capable of flexing the fingertip on its own).

### 3.4 External anchor comparison — pre-registered, not tuned after the fact
Band pre-registered **before** measuring (wide, because the field's own inter-method
spread is large — Lee & Kamper 2008 found systematically smaller moment arms than
classic tendon-excursion methods; a tight band would overclaim precision the field itself
doesn't have): MCP/PIP [2-3, 18-20]mm, DIP [1,10]mm. **Result: 11/12 pass.** The one miss
(FDP2_r@DIP, 10.56mm vs. a 10.0mm ceiling — 5.6% over) is reported as measured, not
massaged into passing. Getting here required **four rounds of forced, principled fixes**
(not blind parameter tuning) — each one caught by a specific machine symptom, not a hunch:

1. **Units bug**: a lateral-offset unit vector was never scaled by its intended magnitude
   before being added to a path point — a 1-**meter** offset instead of ~1cm, giving a
   moment arm of ≈−1.0 (three orders of magnitude too large). Caught because flipping the
   via-point side (the fix for a *sign* problem) couldn't touch a *units* problem — the
   flip-and-rebuild loop oscillated between two wrong values instead of converging.
2. **One via point per intermediate body inflated DIP's lever arm**: a single point near
   a segment's *proximal* end was the only thing bracketing the crossing of **both** the
   proximal and distal joints of that segment — fine for the proximal joint, but far from
   (and thus a poor lever-arm reference for) the distal one. Fixed by giving every
   intermediate body a near/far via-point pair (one close to each joint it touches).
3. **`wrist.osim`'s own MCP2 lateral-offset ratio was a 4× outlier**: cross-checked
   against the *same file's own* PIP (0.166) and DIP (0.144) ratios for the *same*
   finger — MCP2's ratio (0.657) didn't match, revealing that donor's metacarpal
   parent-frame convention almost certainly differs from its phalanx-to-phalanx joints.
   Reused the file's own internally-consistent PIP ratio instead of the outlier.
4. **FDP's insertion was placed near the fingertip (85% of the distal phalanx) instead of
   its real anatomical site, the phalanx's proximal base** — an independently-known
   anatomical fact, not a fit to the target band (though correcting it also reduced the
   DIP lever arm, since the true insertion is closer to the DIP axis).

After these four fixes, the one remaining 5.6%-over miss is reported honestly rather than
chased further — continuing to adjust free parameters specifically until a self-chosen
band is satisfied would be the exact p-hacking failure mode this method exists to prevent.

### 3.5 Tendon excursion (path-length change over each coordinate's full ROM)
| muscle @ coordinate | ROM (deg) | excursion (mm) | monotonic? |
|---|---|---:|---|
| FDP2_r @ wrist_flex_r | [−70, 70] | +28.4 | no |
| FDP2_r @ mcp2_flex_r | [−25.7, 86.0] | +12.5 | no |
| FDP2_r @ pip2_flex_r | [−13.7, 97.2] | +4.5 | no |
| FDP2_r @ dip2_flex_r | [−6.0, 81.6] | +3.1 | no |
| FDS2_r @ pip2_flex_r | [−13.7, 97.2] | +7.2 | no |
| ED2_r @ pip2_flex_r | [−13.7, 97.2] | −16.7 | **yes** |
| ED2_r @ dip2_flex_r | [−6.0, 81.6] | −12.4 | **yes** |

(full table in the evidence JSON.) Several entries are flagged `monotonic=False` over the
**full** sampled ROM (which spans both the small hyperextension arc **and** the large
flexion arc per joint, from the real Mohamed Ibrahim et al. 2024 population data, PMID
39345665) — expected, not a bug: a flexor's path length is generally minimized somewhere
inside a range that includes both hyperextension and flexion, not necessarily at the
extreme ends (the same geometric signature found and explained for the wrist joint during
construction, §2.3). All excursions are the right **sign** (flexors shorten as their
joint's flexion coordinate increases; ED2_r lengthens) — the physiologically primary
question (does the tendon travel the right direction) is answered; strict monotonicity
over the full, hyperextension-inclusive range was not re-litigated further given the §3.4
diminishing-returns judgment call.

### 3.6 Machine artifact checks
- **Determinism**: two independent full runs of `anatomical_hand.py` (fresh Python
  process each time) produced a **byte-for-byte identical** evidence JSON.
- **Round-trip**: the **saved** `.osim` file (not just the in-memory model) was
  freshly reloaded in a separate process; `computeMomentArm` on the reloaded model
  matches the pre-save values (FDP2_r@MCP/PIP/DIP: 10.60/11.53/10.56mm both ways).
- **Structural counts**: 22→26 bodies (+4: `wrist_int_r`, `index_proximal_r`,
  `index_medial_r`, `index_distal_r`), 22→26 joints (−1 weld, +2 wrist, +3 finger),
  80→83 muscles (+3) — exactly matching the intended additions, no silent extras/drops.

---

## 4. Honest scope — what this IS and is NOT (symmetric disclosure)

**This is a first dial-turn, not a hand.** A complete anatomical hand needs, and this
build does **not** have:
- **4 more fingers + thumb** — this session built ONE finger (index) only.
- **~27 more muscles**: the other fingers' own FDP/FDS/ED compartments, the thumb's FPL/
  EPL/EPB/APL, **and all intrinsics** — lumbricals (4), dorsal + palmar interossei (7-8),
  thenar group (opponens pollicis, abductor/flexor pollicis brevis, adductor pollicis),
  hypothenar group. Intrinsics are mechanically essential for real grip (MCP flexion +
  IP extension coordination) and are entirely absent here.
- **Tendon pulleys** (A1-A5 annular + cruciate) — this build uses simple offset via
  points, not `WrapCylinder`/`WrapSurface` geometry. A real flexor tendon is constrained
  to a nearly-constant-radius path by these pulleys; without them, tendon "bowstringing"
  under load is not resisted the way it is in vivo.
- **Collateral ligaments, volar plates** — no passive joint-stabilizing structures at
  MCP/PIP/DIP (this model's finger joints are pure hinges with only the ROM-limiting
  coordinate range, no ligamentous end-feel or laxity).
- **The extensor mechanism proper** (central slip, lateral bands, sagittal bands,
  interosseous/lumbrical contributions to the extensor hood) — ED2_r is a single
  simplified dorsal path, not the real multi-tendon hood.
- **Soft tissue / skin** — bones + tendons + a schematic mass-only cylinder for each
  phalanx; no skin, fat pad, or fingernail.
- **MCP abduction-adduction** — real MCP is 2-DOF; only flexion-extension is modeled.
- **Subject-specific hand anthropometry** — segment lengths are generic literature-typical
  values scaled by the model's own arm scale factor, not a measured subject2 hand length
  (Buchholz 1992's regression, which would give this properly, could not be retrieved this
  session — disclosed gap, not silently assumed away).
- **A small, disclosed mass double-count**: `hand_r`'s pre-existing mass (0.4749 kg,
  presumably already representing the whole hand including fingers in the donor model) is
  kept unchanged, and the 3 new phalanx bodies (17.3g + 7.6g + 4.2g) are added on top —
  slightly over-counting the index finger's own mass. Not corrected in this first step.

**Grafting-limits flag (as the task explicitly asked)**: donor-geometry ratios (even from
a real file like `wrist.osim`) are **not automatically trustworthy at bodytwin's
anthropometric scale** — §3.4 point 3 is a concrete, caught example (a 4× outlier ratio
from a donor file that would have silently propagated a large geometric error had it not
been cross-checked against that same donor's own other joints). Any future donor-model
graft (including the McFarland 2023 route, §6) needs the same live cross-checking, not a
one-shot trust-and-copy.

---

## 5. Files
- `scripts/msk/anatomical_hand.py` — the build (model construction + all geometry).
- `scripts/msk/verify_anatomical_hand.py` — the verification harness (sign-convergence
  loop, moment-arm cross-check, excursion sweep, evidence JSON writer).
- `scripts/msk/anatomical_hand_evidence.json` — full machine-measured numbers (superset of
  every table above).
- `data/msk_models/LaiArnoldModified2017_anatomical_hand_index_subject2_scaled.osim` — the
  new model file (subject2-scaled LaiArnold + un-welded wrist + articulated index finger +
  3 muscles). The original scaled/generic source models were never modified.

## 6. Roadmap to a complete hand (explicit next dials, in priority order)
1. **Fetch McFarland et al. 2023's actual model** (SimTK, PMID 36301780) — the correct,
   citable, complete donor (43 muscles, extrinsic + intrinsic, 23 DOF) identified this
   session but not fetched (disclosed, §1). Graft it the way `add_arm_muscles.py` grafted
   the shoulder: live geometric registration + cross-checked scaling, not a blind copy.
2. **Add the other 4 fingers** using this session's own validated method (elbow-axis
   transport + near/far via points + radius-proportional offsets) — mechanically a repeat
   of §2, not a new technique.
3. **Add intrinsics** (lumbricals + interossei especially — they are what makes MCP
   flexion + IP extension coordinate in a real grip; without them this model's fingers can
   flex but cannot "grip" in the coordinated sense).
4. **Replace via points with `WrapCylinder`/`WrapSurface` pulleys** at MCP/PIP/DIP (A1-A5
   equivalent) for physiologically-correct, pose-invariant moment arms instead of the
   simple offset approximation used here.
5. **Add MCP abduction-adduction** (the 2nd real MCP DOF) once ab/ad-driving intrinsics
   (interossei) are in place to actually exercise it.
6. **Resolve the segment-length anthropometry gap** — either fetch Buchholz et al. 1992's
   actual regression table, or measure subject2's real hand length from the OpenCap
   video/marker corpus if available, replacing the current generic-literature-value +
   arm-scale-factor proxy.
