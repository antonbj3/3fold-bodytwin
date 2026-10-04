# MECHANISM FIDELITY ARCHITECTURE AUDIT — does the current approach converge toward 100% anatomical correctness, or plateau? (2026-07-21)

**Purpose.** A self-audit, not a sales doc. The operator's challenge: is the current OpenSim/Hill/monocular
stack actually *aligned* to work toward 100% anatomical correctness, or does it plateau? This document
formalizes and stress-tests the coordinator's assessment against the actual repo evidence — every number
below is quoted from an existing, dated (2026-07-21) `docs/MECHANISM_*.md` cert, not re-derived or asserted.
Where the evidence refines or contradicts the assessment, that is reported the same way the underlying certs
themselves report an honest negative: as a PASS for the audit discipline, not a failure to be softened.

**Scope note on freshness.** Nearly all Phase-2 (mechanical-twin) evidence cited here is from **one calendar
day** (2026-07-21), on **one subject** (OpenCap `subject2`, 78.2 kg), mostly **one trial** (`walking1`),
**right-leg-primary**. This is a real, additional, orthogonal caveat to the four substrate ceilings below — a
*breadth-of-validation* ceiling, not a substrate one — and it means every number in this audit should be read
as "true for n=1, this trial" until re-run across subjects/trials, which has not happened yet.

---

## 0. Executive verdict

**The method is genuinely aligned and is being applied with unusual rigor** (§1). **The four substrate
ceilings are real, confirmed, and currently binding** — none of them is closed by more parameter-tuning
within the existing stack (§2). **The pieces are NOT yet mutually composable**: the fidelity increments
built today live in at least **9 separate, un-merged model files**, and at least two of them expose a
genuine cross-substrate incompatibility when forced together (ligaments break down kinematically beyond
~100° flexion; the only cartilage-contact model in the repo lives on a different skeleton/individual
entirely) (§4). **100% is not reachable by this stack; it is asymptotic within it, with a hard ceiling from
outside it.** Three of the four ceilings (subject-specific data, deformable continuum, medical-grade
acquisition) require a genuinely different SUBSTRATE (new data, new physics, new sensing), not a smarter
algorithm on what exists. The fourth (Hill→cross-bridge) is shown *traversable in principle* on one muscle,
but is currently a disconnected proof-of-concept, not a load-bearing layer (§3). The shortest credible path
is **not** "graft more literature-scaled tissue" (the thing most of today's effort did) — it is (1) merge the
forked models into one consistent instance before the fragmentation compounds, and (2) acquire ONE real
subject-specific imaging dataset, because that is the single migration that starts dissolving the ceiling
that every other layer's honest-gap section blames first (§5).

---

## 1. The coordinator's three claims — confirmed / refuted / refined

### 1a. "ALIGNED: the method (measure/falsify/decorrelated-anchor/symmetric-QC)" — **CONFIRMED, strong evidence**

This is the best-supported of the three claims. Concrete, dated instances of the method actually catching
itself, not merely being stated as a value:

- **A self-flattering pass caught and reversed.** The knee/hip self-computed contact-force numbers
  (233.20 %BW / 235.30 %BW) looked plausible — close to the OrthoLoad anchor — precisely *because* a sign
  bug in `knee_crossing_muscles_and_forces` was invisible inside a `np.linalg.norm()` (always ≥0). It was
  caught only because a **sibling cert (spine) used a signed projection** instead and got an
  immediately-impossible negative number, triggering a cross-cert audit (`docs/MECHANISM_SIGN_BUG_AUDIT.md`).
  Corrected values (391.10/386.77 %BW) converge on the independent, bug-immune `opensim.JointReaction` to
  **<0.1%** at both joints — the anchor that proves the *correction*, not the original, is right
  (`docs/MECHANISM_SIGN_BUG_REMEDIATION.md`). The bug had been independently re-implemented (copy-pasted, not
  shared) in `validate_shoulder_force_with_muscles.py` and was still live until a same-day follow-up found
  and fixed it too (`docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md`).
- **An honest-negative that stayed honest under pressure to look like a win.** Running RRA to shrink the
  pelvis residual showed a 4.6× improvement *in RRA's own accounting* — but re-validated through the actual
  muscle-driven model, the improvement is gone (175.32 N vs 172.96 N, statistically unchanged) and the
  moment residual is worse (+85%); knee/hip contact-force estimates moved **further** from the OrthoLoad
  anchor, not closer. This was root-caused (undersized ankle/subtalar reserve actuators in a no-muscle RRA
  context), not shrugged off as noise (`docs/MECHANISM_RRA.md`).
- **A fixed heuristic falsified across the diverse instance space, not accepted on the first clip.** The
  MediaPipe→OpenSim axis remap (`OpenSim_Y = -mediapipe_y`) passed a sign check but failed a magnitude check
  (0.24 m vs an expected ~1.4 m); the *replacement* guess ("X is vertical") was tested against all 5 pilot
  clips and falsified on 2/5 (the correct diverse-instance-space test, not a single confirming clip) before
  the real per-clip geometric fix was derived and re-verified at full 264-clip scale
  (`docs/MECHANISM_POSE_PIPELINE.md` §4.1, §9).
- **A cross-bridge result reported as mixed, not spun.** History-dependence (RFE/FDE) was tested and found
  *provably absent* in the model's own current structure — the honest, symmetric verdict recorded is "adds
  one form of history-dependence Hill cannot have, not the other," not "cross-bridge wins"
  (`docs/MECHANISM_CROSSBRIDGE_MODEL.md` §6.4).

Verdict: **CONFIRM.** The forced-adversary / external-anchor / symmetric-QC discipline is not aspirational
here — it is visibly running and visibly catching real defects, including ones that made the twin look
*better* than it is.

### 1b. "the layer ONTOLOGY (joint→muscle→tendon→ligament→cross-bridge)" — **REFINE**

This ontology names a real, working chain for exactly one tissue axis (the muscle-tendon unit's mechanistic
depth) and should not be read as the anatomical layer stack — it omits cartilage/contact, fascia,
vasculature, nerve, fat, skin entirely (§3 covers all of these; several are at ~zero fidelity, a
categorically different state than "joint/muscle/tendon/ligament"). Within its own stated scope:

- joint → muscle → tendon → ligament is real and independently verified (IK to sub-1.5°, 80-193 muscles
  depending on fork, Millard tendon curves force-anchored to two literature sources, 84 Blankevoort ligament
  bundles pattern-verified) — **but "→ cross-bridge" is not yet a fifth operating layer**. It is a single
  muscle (soleus), explicitly built **decoupled from the full-body simulation** — "a standalone
  sarcomere/half-sarcomere, not wired into the `.osim` musculoskeletal chain or any joint-force cert in this
  repo" (`docs/MECHANISM_CROSSBRIDGE_MODEL.md`, its own scope line). No cert in this audit's evidence base
  uses cross-bridge output for anything.

Verdict: **REFINE.** The chain is real for 4 of its 5 named links, operating on real (if literature-scaled)
data; the 5th link is a proof-of-concept sitting next to the chain, not yet in it.

### 1c. "a PROVEN multi-scale link (Hill emerges from cross-bridge kinetics — a coarse-graining ladder, not a dead end)" — **CONFIRM-NARROW, with a major symmetric asterisk**

What is genuinely proven (`docs/MECHANISM_CROSSBRIDGE_MODEL.md`):
- Isometric force-length plateau + descending limb **emerges from pure sarcomere filament-overlap geometry**
  (not fit), reproducing Gordon-Huxley-Julian 1966's plateau bounds exactly on the frog calibration case.
- A Hill force-velocity hyperbola **emerges** (never coded in) from Huxley 1957 two-state kinetics on the
  **concentric** branch: R²=0.9987 against the pre-registered ≥0.95 gate, RMS diff 0.0537 vs the real
  Millard2012 curve on the same muscle.
- A genuine advantage over Hill exists: a finite relaxation time constant (5.4 ms) after a length
  perturbation — Hill/Millard is *provably* incapable of this (RFE_Hill = 0.0000000000%, no internal state to
  relax) — a structural, not fitted, distinction.

What symmetrically qualifies the claim, measured not asserted:
- The **eccentric** (lengthening) branch is substantially wrong: 2.57× over-prediction of peak force vs
  Millard's own physiologically-capped 1.4×, checked against a bounded-rate variant and found to get
  *worse*, not better — a genuine, non-fixable-by-this-model-class limitation (motivating why real
  multi-state models like Razumova/Campbell 1999 exist, cited but not built).
- **RFE/FDE — the most commonly-cited reason to prefer cross-bridge over Hill — is *provably absent*** in
  this specific model construction (uniform half-sarcomere, no titin): −0.085%/+0.005%, shown to shrink
  monotonically toward exactly zero with hold duration. On this specific, famous discriminator, the current
  cross-bridge model behaves **like Hill**, not better than it.
- Parameters are canonical/illustrative throughout (Huxley's own f1:g1:g2 ratio, calcium kinetics, bare-zone
  geometry) — not fit to this subject, not fit to soleus specifically, several flagged as "not independently
  re-verified against primary sources this session" (paywalled originals).
- **Zero coupling to the twin.** It runs standalone; nothing in this audit's evidence base is downstream of it.

Verdict: **CONFIRM the physics, REFINE the claim's scope.** The coarse-graining ladder is real,
non-trivial, and "not a dead end" in the sense that the *method* is sound and extensible — but it currently
contributes **zero measured fidelity** to the actual twin, reproduces Hill only on the physiologically
dominant (concentric) side, and — on the one property (RFE) most people would cite cross-bridge models FOR —
currently offers no improvement over Hill at all. Scaling it to "every muscle, coupled, subject-specific" is
a large, mostly-unstarted undertaking, not a remaining formality.

---

## 2. The four substrate ceilings — confirmed, and refined by measured evidence

| # | Ceiling (as stated) | Confirmed? | Refinement from this session's evidence |
|---|---|---|---|
| 1 | generic → subject-specific (needs THIS person's imaging) | **CONFIRMED, pervasive** | Not one bone, muscle, tendon, or ligament parameter in any fork traces to subject2's own imaging. Bone geometry = a donor cadaver template, linearly rescaled per marker-based `Scale`. Muscle strength (Fmax, pennation) is **bit-identical between the generic template and the "scaled" model for all 80 native muscles and every grafted set** (`docs/MECHANISM_MUSCLE_AUDIT.md` §3) — scaling touches geometry only, never strength. Ligaments are Lenhart2015's healthy-donor geometry, isotropically rescaled (`docs/MECHANISM_KNEE_LIGAMENTS.md`). Tendon compliance is the OpenSim class *default*, identical across all introspected muscles (`docs/MECHANISM_TENDON_ELASTIC.md` §7.1, its own "DOMINANT UNCERTAINTY"). |
| 2 | phenomenological Hill → mechanistic cross-bridge in EVERY muscle | **CONFIRMED as a real gap, but LOWER-URGENCY than the others in the current regime** | §1c: Hill already matches cross-bridge to R²=0.9987 on the concentric branch, and matches real physiology on eccentric (1.4× cap) *better* than the plain cross-bridge model does (2.57×). The gap matters most for the regimes Hill structurally cannot do at all — RFE/fatigue/damage — none of which any current cert exercises. Scaling cross-bridge to "every muscle" also re-imports Ceiling #1 (subject-specific rate constants don't exist) and adds a new numerical-coupling problem (integrating ~100+ stiff ODEs against multibody dynamics) not yet attempted at any scale beyond 1 muscle. |
| 3 | rigid-body + independent line-actuators → deformable continuum multi-physics | **CONFIRMED, and one assumed precedent is REFUTED** | Cartilage contact exists, but only in a *separate* model (JAM/COMAK "DM," a different individual, TKA-implant, zero ACL) never merged with the muscle/ligament model used everywhere else in this audit. Fascia, myofascial lateral transmission, muscle bulging, and skin sliding have **no model of any kind** — not low-fidelity, absent. **The coordinator's implied precedent, CS's "HUM-DEFORMABLE" cell, is NOT a human-tissue deformable model** — verified directly from its own claim/probe JSON (`cad-to-simulation-I/data/claims/I-HUM-DEFORMABLE.json`, `reports/probes/hum_deformable.json`): it is a **PET polarizer-film lamination onto an LCD panel for a pico-projector**, a manufacturing cell using Newton's generic cloth/FEM solver (`SolverVBD`), `regime_gate.open=false`, blocked, with zero biological content. What IS real and useful: the underlying engine capability (Newton `add_cloth_grid`/`add_soft_grid`/`add_tetrahedron`, `SolverVBD/XPBD/Style3D/ImplicitMPM`) is proven to run stably on a real (non-biological) deformable-material problem in the same shared stack — lowering the *engineering* barrier to a future biological port, but no such port exists today, and bodytwin's own `MSK-CONNECTIVE` literature cell explicitly records this coupling as **"blocked"** (`docs/MECHANISM_CELL_INVENTORY.md`). |
| 4 | monocular social-media pose → medical-grade acquisition | **CONFIRMED, and partly a hard information-theoretic limit, not an engineering gap** | MediaPipe's `pose_world_landmarks` are **re-centered on the hip in every single frame independently** — this is not a bug to be fixed by a better model, it is how the representation is defined, and it destroys all whole-body-translation signal (jump height, squat depth) for every clip processed, permanently (`docs/MECHANISM_POSE_PIPELINE.md` §4.2). What IS recoverable — relative joint angles — is good: ~4.7° lower-limb RMSE (matches a ~4.5° multi-camera anchor) on the older, currently-broken-in-place video chain, and 228/264 (86.4%) of the real IG athletic corpus passes a full, independently-cross-checked gate suite this session. 31/264 (11.7%) of the "priority" corpus was found to be mislabeled blank video — a corpus-quality ceiling layered on top of the acquisition ceiling. |

---

## 3. Per-layer fidelity table (top-to-bottom anatomical stack)

| Layer | Current fidelity (cert / number) | Substrate ceiling | Concrete migration | Falsifier |
|---|---|---|---|---|
| **Skeleton** | Donor-lineage rigid bodies (LaiArnold/Rajagopal), linearly rescaled per-segment to subject2's markers; mass preserved exactly (78.2 kg); IK tracks real mocap to median RMSE 0.055°, max 1.328° vs an independent reference (`MSK_ENV.md`). | #1 — geometry is a rescaled population template, not subject2's own bones. | CT/MRI/DXA-derived subject-specific segment geometry + inertial properties. | Overlay subject's own DXA/CT outline on the scaled segments — already known to fail by construction (no imaging was ever used). |
| **Joints** | Verified IK to sub-1.5° (`MSK_ENV.md`); multi-segment foot present out of the box (ankle+subtalar+mtp); a midtarsal-joint prototype passes in-sample (median 0.123°) but **fails held-out** generalization (max-abs 5.76° > the pre-registered 5.0° ceiling, subtalar RMSE 4–10× the noise floor, `MECHANISM_FOOT_MULTISEGMENT.md` §4); knee is a **prescribed 1-DOF** polynomial-coupled joint, not free; wrist is a 0-DOF `WeldJoint` (confirmed by direct XML grep, `MSK_ENV.md`); a real scapulothoracic+sternoclavicular chain was grafted but the acromioclavicular loop **cannot be closed** below a measured ~12.6 cm gap after an exhaustive multi-start search — shipped kinematically locked (`MECHANISM_SCAPULA_CLAVICLE.md` §4). | #1 (subject-specific joint centers/axes) + a live OpenSim 4.6 engine bug: `InverseDynamicsTool` misreports `knee_angle_r/l`/`hip_flexion_r/l` generalized force by **~1000–1800×** at non-trivial poses (`MECHANISM_MSK_ELASTIC_BAND.md` §4) — a software-correctness ceiling every force cert in this repo now has to route around. | Subject-specific functional joint-axis calibration or imaging; a free-DOF knee (JAM/COMAK-style) merged into the SAME model as the muscles/ligaments; real wrist DOF; resolved AC registration; an upstream OpenSim fix. | The ID-tool bug's own falsifier is already true and reproducible: at the model's literal zero pose, gravity only, `knee_angle_r_moment`/`hip_flexion_r_moment` read −1252.6/+1337.1 N·m — anatomically impossible for a hanging leg. |
| **Muscles** | 80 native `Millard2012EquilibriumMuscle` (0/80 physiologically-implausible, `MUSCLE_AUDIT.md`); +88 erector-spinae/multifidus (ceiling 228 N·m vs Arokoski in-vivo 147–170 N·m — sane, `ERECTOR_SPINAE.md`); +25 arm/shoulder (47/70=67% moment-arm-sign checks pass, `ARM_MUSCLES.md`); +11 scapulothoracic stabilizers (`SCAPULA_CLAVICLE.md`). **These three grafts live in three separate, never-merged model files.** Fmax/pennation are **never** subject-rescaled anywhere (ratio exactly 1.000000 for all 80 native muscles, confirmed 3 independent ways). | #1 (cadaver/literature Fmax, not subject2's measured strength) **and** #2 (Hill, not cross-bridge) bind simultaneously here. | Dynamometry or MRI-volume-based subject-specific strength; extend the 1-muscle cross-bridge prototype to all muscles, coupled into the dynamics loop, if RFE/fatigue/damage regimes are ever needed. | This subject's own isokinetic torque-velocity curve vs the model's Fmax-derived prediction — never measured. |
| **Tendons (+pulleys)** | Millard tendon curve, energy accounting self-tested to 0.0000% error, anchored to Lichtwark&Wilson2005 (ratio 0.20 of the hopping ceiling, sane direction) and Fukunaga2001 walking strain (12.13/9.30 mm vs ~7 mm anchor) (`TENDON_ELASTIC.md`). | #1 (curve shape is the OpenSim **default**, identical across every muscle checked, not ultrasound-measured) **and** a #3-adjacent gap: zero hysteresis (real tendon 17–35%), so "29.8% elastic" is an upper bound, not measured recoil efficiency; Achilles modeled as 3 independent tendons, not one confluent structure; wrap surfaces (pulleys) were **not ported** for any grafted muscle, causing measured pose-dependent moment-arm sign flips (deltoid/supraspinatus, `ARM_MUSCLES.md` §6). | Ultrasound-measured subject-specific tendon compliance; ported wrap geometry; a hysteretic/viscoelastic tendon model. | Ultrasound-measure this subject's real Achilles strain under a known load vs the model's prediction at matched force — never done. |
| **Ligaments** | Blankevoort1991Ligament, 84 bundles, ported from a real independent donor (Lenhart2015, PMID 25917122), engagement pattern (ACL-tight-extension/PCL-tight-flexion) PASSES pre-registered thresholds, MCL/LCL confirmed frontal-plane antagonists geometrically (`KNEE_LIGAMENTS.md`). | #1 (donor geometry, isotropically rescaled, not subject2's knee MRI) **and a measured cross-substrate incompatibility with #3**: ligaments parametrized on a free-6-DOF donor knee, grafted onto a prescribed-1-DOF target knee — MCL strain reaches a **non-physiological 58%** at 140° flexion (real tissue fails ~15–20%), verified to be a kinematic-coupling mismatch, not a code bug (force-law and slack-length checks both PASS to <1e-4). No articular-contact element exists in this same model at all. | Either a free-DOF knee formulation merged with the muscle+ligament model, or an explicit validity bound (<100° flexion); ultimately subject-specific ligament geometry from MRI. | **Already fired**: the 58%-strain-at-140°-flexion prediction is itself the falsifier — real MCL would have ruptured well before that strain. |
| **Cartilage (contact)** | Exists only in the separate JAM/COMAK "DM" model: predicted tibiofemoral contact 2.503 BW vs Grand-Challenge measured 2.584 BW on trial og1, ~3% error, r=0.802 — og1 is the best of the 3 trials actually run (og4 3.76%/r=0.741, og3 10.35%/r=0.658; honest range 3-10% peak error). **CORRECTION 2026-07-26**: this cell was headlined here (and in `MECHANISM_UNIFIED_MODEL.md`, `MECHANISM_LAYER_COMPLETENESS.md`) as "the single strongest cert in the repo" — a domain-default audit reverses that superlative: OrthoLoad's independent population-median baseline (zero fit to subject DM) beats the model on mean absolute error across the 3 trials (3.06% vs the model's 5.72%) and wins 2 of the 3 held-out trials (og1, og3); the model wins only og4 (`AUDIT-VOIDFLOOR-DOMAINDEFAULT-KNEE-SARCOMERE-2026-07-26`, `AUDIT-KNEE-BODYWEIGHT-WINDOW-CORRECTION-OG1-SELECTION-2026-07-26`). What the model still demonstrably owns against a fair adversary: it beats a cross-subject-scrambled-activation null by 5-7x on RMS (`VOID-FLOOR-KNEE-GENERIC-WAVEFORM-CLEARS-BAR`) — the most-instrumented mechanistic cert in the repo, no longer its most accurate number. | #3 (only for an artificial/implant joint state) **and** #1 (that model is a different, non-subject2 individual — the Grand-Challenge TKA patient). | Couple a native (non-implant) cartilage-contact model into the SAME model instance carrying subject2's muscles+ligaments; needs MRI-segmented cartilage geometry for a healthy joint. | Run COMAK on subject2's own ligamented model against any native-knee in-vivo contact measurement — not attempted; the 3% flagship accuracy has never been shown to generalize past the one TKA patient it was measured on, and a same-subject domain-default baseline now beats it on 2 of 3 trials (above). |
| **Fascia (myofascial transmission)** | **~zero.** Repo-wide search finds no lateral-force-transmission model; "fascia" hits are a muscle-name substring, an unmodeled video-caption keyword, and one RE-SCOPED literature leg explicitly flagged qualitative/bovine-ex-vivo/non-quantitative (`MECHANISM_CELL_INVENTORY.md`, MSK-CONNECTIVE). | #3, total absence. | A continuum/membrane FEM myofascial layer — the generic engine capability exists and runs stably elsewhere in the shared stack (§2 row 3) but has never been pointed at biology. | Any lateral force-transmission or fasciotomy redistribution measurement — moot; no mechanism exists to even attempt a prediction. |
| **Vasculature** | **~zero** in the mechanical model. Present only as population/epidemiology cells (`ORG-CARDIAC-MECHANICS`, `ORG-CARDIAC-PLAQUE`) and one fat-identifiability leg (NIRS perfusion), fully decoupled from the MSK sim. | #3, substrate absent (no fluid/network representation of any kind). | A lumped-parameter or 1D circulatory/perfusion model coupled to muscle metabolic demand — a new substrate, not a parameter upgrade. | Any exercise NIRS/Doppler blood-flow measurement vs the twin — moot, no mechanism predicts it. |
| **Nerves / neural drive** | Proxy only: Static Optimization's activation-minimizing solution. Timing validated vs literature EMG: **0/8 PASS, 3/8 PARTIAL, 5/8 FAIL** on the pre-registered Jaccard-vs-circular-null gate (mean precision 0.809, mean recall 0.631) — but the failures independently **reproduce** the model lineage's own published CMC-vs-EMG mismatch (glut max/medius) and a specific literature-predicted defect (mid-swing hamstring) was checked for and correctly found **absent** (`EMG_TIMING.md`). | #2 (an optimization proxy, not a neuromuscular control model — no reflex arc, no motor-unit recruitment dynamics) **and** #1 (no EMG was ever recorded for this subject; windows are population-averaged Perry-1992 convention). | Record this subject's own EMG during captured trials (not currently in the corpus) to constrain or replace SO; eventually a real recruitment/motor-unit model. | **Already fired**: the pre-registered gate already returned 0/8 PASS. |
| **Fat** | Advanced on the **literature/identifiability track only**: fusing surface geometry + bioimpedance cuts regional-fat error 6.3–14.4% vs best single modality (pre-registered, CI excludes 0, 6/6 region×cycle). **Zero presence in the mechanical model** — rigid segments are single lumped bone+muscle+fat+skin masses; fat cannot redistribute load or change segment inertia dynamically. | #3 — no soft-tissue compartment substrate exists in the rigid-body sim, independent of how good the *identifiability* estimate is. | Couple the identifiability track's fat-fraction estimate into a genuine per-segment soft-tissue or FEM compartment. | Compare the rigid model's fixed per-segment inertia against the identifiability track's own fat-distribution estimate for the same subject — never attempted; the two tracks have never cross-checked each other despite living in the same repo. |
| **Dermis / epidermis** | **~zero mechanically.** A literature-only cell (`ORG-SKIN-BARRIER`) covers dermoscopic lesion classification and population wound-healing — a different question (dermatology) entirely. | #3, total absence. | Skin-sliding/compliance FEM layer (same unclaimed engine precedent as fascia). | Skin-marker-slippage-during-motion measurement — moot; notably, every kinematic marker used everywhere in this audit rides on unmodeled skin, so this is a live, silently-uncorrected error source baked into every IK number cited above, not a purely hypothetical future gap. |
| **Follicles / glands** | **Zero.** At most a listed entry in the source catalog (`WAVE_PLAN.md`); never instantiated as a built cell. | N/A — not attempted; arguably out of the MSK mission's scope (enters via the separate thermoregulation/calorimetry cert-twin, not the musculoskeletal one). | Would only enter via the metabolism/calorimetry coupling (`COORDINATOR.md` §5), not directly relevant to load-bearing fidelity. | N/A. |

---

## 4. Composability verdict — NOT YET COMPOSABLE

Every individual layer above is verified in isolation with real rigor. The pieces do **not** yet compose into
one consistent instance. This is not a hypothetical risk — it is the audited, current state:

**1. At least 9 separate model-file forks exist, none merged:**
`..._scaled.osim` (base, 80 muscles) · `..._erector_spinae_subject2_scaled.osim` (+88 muscles, 168 total,
used only for the spine cert) · `..._arm_muscles_subject2_scaled.osim` (+25 muscles, 105 total, built on the
*base*, not on the erector-spinae fork) · `..._scapula_clavicle_subject2_scaled.osim` (+2 bodies, 116
muscles, built on the arm-muscles fork) · `knee_ligaments/model_with_knee_ligaments.osim` (+84 ligament
bundles on the *base*, no erector spinae, no arm muscles, no scapula) · `LaiArnold_midtarsal_proto_scaled.osim`
(+1 DOF/side, on the *base*) · the hallux-split prototype (dial-turn 1, likewise base-only) · the elastic-band
model pair · and the wholly separate JAM/COMAK "DM" model (different skeleton lineage, TKA implant, zero
ACL, used only for the KNEE-CELL flagship). **No single file has erector spinae + arm/shoulder muscles +
scapula/clavicle + knee ligaments + the midtarsal foot together.** Each fidelity increment was deliberately
built as an isolated fork "so the evidence trail stays legible" — a sound verification practice, but it means
today's gains do not compound; they sit next to each other.

**2. A genuine, measured cross-substrate incompatibility, not a hypothetical one**: the knee-ligament layer
(parametrized on a free-6-DOF donor joint) and the joint layer it was grafted onto (a prescribed-1-DOF
`CustomJoint`) visibly disagree beyond ~100–110° of flexion (MCL strain 58%, non-physiological) — precisely
the kind of "OpenSim Hill hand vs continuum-FEM hand" incompatibility the operator's framing anticipated,
except it shows up **one level down**, between two rigid-body kinematic conventions, before continuum
substrates even enter the picture.

**3. The cartilage/contact substrate and the muscle/ligament substrate are on different individuals
entirely** — the JAM/COMAK "DM" model (Grand-Challenge TKA patient) has never been merged with, or even
run against, subject2's LaiArnold-lineage muscle+ligament model. The repo's most-instrumented cert (og1's
3% error vs Grand-Challenge — no longer the unqualified "single strongest," see the Cartilage row above and
`AUDIT-KNEE-BODYWEIGHT-WINDOW-CORRECTION-OG1-SELECTION-2026-07-26`) and its most detailed ligament work (84
Blankevoort bundles) have never touched the same model file.

**4. The claimed cross-domain deformable-tissue precedent does not exist for biology** (§2 row 3) — CS's
`HUM-DEFORMABLE` is a projector-manufacturing cell, `regime_gate.open=false`, and bodytwin's own graph
already records the intended coupling as `blocked` — self-consistent with this audit's independent finding.

**5. The shared substrate (OpenSim's own standard tool) is not internally trustworthy for this model
family**: the `InverseDynamicsTool` bug forces *every* force cert (knee, hip, ankle, elbow, shoulder, spine,
elastic-band) to build its own bespoke Newton's-law/virtual-work workaround rather than use the tool as
shipped — and the workarounds themselves have then had to be independently re-audited for a **second**,
unrelated bug (the crossing-muscle sign error), which was found live and unfixed in at least one further
file (`validate_shoulder_force_with_muscles.py`, fixed same-day) and remains unfixed in two peripheral docs'
saved numbers (`FORCE_SCENES_BATCH`, `reconcile_diff_scheme.py`'s JSON) as of this writing. Composability
failures are propagating at the *code* level, not only the *model* level.

**Net:** the individual layers are each real, mostly good-to-excellent verification work. The **system** is
currently a set of well-tested parts, not a tested whole.

---

## 5. Is 100% reachable, and what is the shortest path?

**Not reachable by this stack; asymptotic, with a hard external ceiling.** Three of the four substrate
ceilings (#1 subject-specific data, #3 deformable continuum, #4 medical-grade acquisition) require a
genuinely different substrate — new imaging data, a new physics representation, new sensing — not a smarter
use of what already exists. No amount of further literature-scaled muscle/ligament grafting closes them; the
muscle audit already proved the current scaling pipeline **structurally never touches strength**, only
geometry (§3, Muscles row) — an algorithmic ceiling, not a data-quality one, inside the current pipeline. The
fourth ceiling (#2, Hill→cross-bridge) is shown traversable in principle on one muscle, but is the
*lowest-urgency* of the four given the current regime (Hill already matches the concentric branch to
R²=0.9987, and matches real eccentric behavior *better* than the plain cross-bridge model tested so far).

**Ranking migrations by leverage, given everything already built:**

1. **Highest immediate leverage, lowest cost: merge the forked models and pay down the sign-bug/tool-bug
   debt (§4).** This is pure engineering, needs no new data or physics, and is a *prerequisite* for the other
   three migrations to compound rather than dead-end as isolated demos — right now, closing Ceiling #1 for
   the skeleton would still leave the ligaments, muscles, and cartilage each on their own un-merged fork.
2. **Highest ceiling-clearing value: ONE real subject-specific 3D imaging dataset** (even a single
   consenting subject's CT/MRI/DXA). This directly starts dissolving Ceiling #1, and — unlike the other two
   remaining ceilings — gives the ligament/muscle grafts a real registration target instead of another
   population-average donor (the AC-joint 12.6 cm gap and the >100° MCL blow-up are both, at root,
   cross-donor registration problems that subject-specific geometry would bound, though not eliminate, since
   material properties still need separate measurement).
3. **Ceiling #4 (monocular acquisition) is the one most clearly requiring a different substrate outright**:
   MediaPipe's hip-centered-per-frame representation discards translation *by construction* — this needs a
   different perception algorithm or a return to multi-view/calibrated capture, not a better monocular model.
   It is also the ceiling that gates whether this can ever run on the 893-clip social-media corpus the
   project's own mission depends on, so it has the largest data-scale multiplier even though its per-subject
   fidelity payoff is smaller than #1's.
4. **Ceiling #2 (cross-bridge everywhere) is the correct lowest priority right now**: it mainly buys
   correctness in regimes (fatigue, damage, RFE) that no current cert exercises, it re-imports Ceiling #1
   the moment it's taken seriously (subject-specific rate constants), and scaling from 1 decoupled muscle to
   ~100+ coupled ones is a nontrivial numerical-integration project in its own right.

**One line:** the fastest way to make today's real progress compound, rather than plateau as a pile of
individually-verified but mutually-incompatible forks, is to unify the model before extending it further,
and to trade one more literature-scaled muscle graft for one real scan.

---

## Evidence index (files read/verified for this audit)

`docs/MECHANISM_MISSION.md` · `MECHANISM_STATE.md` · `MECHANISM_HARDENED_CONVENTIONS.md` ·
`MECHANISM_CELL_INVENTORY.md` · `MECHANISM_SCALING_DOCTRINE.md` · `MECHANISM_MSK_BUILD_CHARTER.md` ·
`MECHANISM_MSK_BUILD_PLAN.md` · `MECHANISM_MSK_ENV.md` · `MECHANISM_MUSCLE_AUDIT.md` ·
`MECHANISM_TENDON_ELASTIC.md` · `MECHANISM_KNEE_LIGAMENTS.md` · `MECHANISM_CROSSBRIDGE_MODEL.md` ·
`MECHANISM_JOINT_FORCE_VALIDATION.md` · `MECHANISM_STATIC_OPT.md` · `MECHANISM_SIGN_BUG_AUDIT.md` ·
`MECHANISM_SIGN_BUG_REMEDIATION.md` · `MECHANISM_ERECTOR_SPINAE.md` · `MECHANISM_MSK_ELASTIC_BAND.md` ·
`MECHANISM_SPINE_FORCE.md` · `MECHANISM_ARM_MUSCLES.md` · `MECHANISM_SCAPULA_CLAVICLE.md` ·
`MECHANISM_ELBOW_FORCE.md` · `MECHANISM_EMG_TIMING.md` · `MECHANISM_FOOT_MULTISEGMENT.md` ·
`MECHANISM_POSE_PIPELINE.md` · `MECHANISM_RRA.md` — all under `source_documents/`.
External (read-only, per COORDINATOR.md §1 isolation): `~/projects/cad-to-simulation-I/data/claims/I-HUM-DEFORMABLE.json`,
`~/projects/cad-to-simulation-I/reports/probes/hum_deformable.json`.

No writes were made outside this document. No git commit or push performed (isolation-respected, per task).
