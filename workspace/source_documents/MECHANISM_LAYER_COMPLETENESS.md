# MECHANISM TWIN LAYER-COMPLETENESS MAP (2026-07-21, snapshot cutoff ~19:36)

**Purpose.** `docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md` (written ~17:44) mapped the anatomical stack.
`docs/MECHANISM_TRUST_LEDGER.md` (~18:25, 66 docs/72 rows) then machine-cross-checked every headline claim
into 5 tiers. This repo is **live and concurrently written** (two mechanism instances, per `COORDINATOR.md` §0) —
in the ~70 minutes since the ledger, **19 more scoreable docs** landed, taking the doc count from 66 to **85**.
Three separate freshness re-checks were run while writing this map, each catching a doc that appeared mid-write
(`MECHANISM_CARDIAC.md`, then `MECHANISM_HIP_ANKLE_LIGAMENTS.md`/`MECHANISM_FULL_HAND.md`/`MECHANISM_EMG_DRIVEN.md`,
then `MECHANISM_SPINE_LIGAMENTS.md` — the last one specifically because it falsified an earlier draft of this
map's own "spine ligaments not built" line). **Hard cutoff: 85 docs, re-confirmed stable immediately before
finalizing.** This doc folds those 19 into the ledger's own per-layer stack, using the ledger's own tier
definitions, not a re-badging of prose.

**Method.** The original 72 ledger rows are reused as-is (not re-verified end-to-end here; 2 of its most
load-bearing numbers were spot-checked directly against raw JSON below and held — knee `self_computed` stale
233.20%BW + new `_STALE_PRE_SIGN_FIX_NOTICE` key, and bone_stress's Taylor-1996 contradiction). The 19 new
rows were each independently pulled from on-disk JSON (11 read directly, 8 via two QC'd research passes,
each spot-checked against its most surprising claim before being trusted — same discipline the ledger itself
used for its own delegated rows). Isolation respected: bodytwin only, read-only, no git commit, no push.

---

## HEADLINE FRAME — read this before the table

**Coverage got broader. Validation depth did not get proportionally deeper, and the dominant failure mode
is unchanged: one shared upstream solve, reused downstream and mistaken for independent confirmation.**

- Of the 85 docs, the plurality tier is still **method-only-no-external-anchor (32/91 rows, 35%)** — internally
  rigorous, never compared to anything outside this twin's own pipeline. A layer *existing* is not a layer
  *validated*.
- **A concrete, newly-observed common-mode chain**: `metabolic_cost.py` (itself downstream of the single
  subject2/`walking1` Static-Optimization activation solve) is now the **trunk** for three "independent-sounding"
  new layers — thermoregulation, respiratory gas-exchange, and cardiac output all read `metabolic_cost_results.json`
  directly, none re-solves anything. `docs/MECHANISM_RESPIRATORY.md` says this about itself, and it is the single
  best sentence in this entire audit: *"the twin-internal number itself is one shared instrument chain across
  all three layers... reported as '3 external anchors agree,' never as '3 independent measurements.'"*
- **A concrete, newly-observed error-propagation instance** (not hypothetical): `MECHANISM_BONE_STRESS.md`'s own
  headline number is **in tension with** (its words) the Taylor 1996 literature anchor — bending-dominated where
  Taylor found compression-dominated. `MECHANISM_BONE_REMODELING.md` then classifies bone state ("formation, not
  maintenance") purely from that same number, and its own §4 admits the classification "would flip to resorption"
  if the upstream contradiction were resolved — a downstream layer inheriting, not independently testing, an
  upstream defect.
- **One genuine bright spot against the common-mode pattern**: `MECHANISM_EMG_DRIVEN.md` is the first cert in this
  family to drive real muscles from **measured surface EMG** instead of the SO solve, for 7/80 muscles. It finds
  real co-contraction SO structurally cannot see (gated CCI 0.39 vs 0.19), and a joint-specific, non-uniform force
  effect (+48% at the knee, -12% at the hip) — but only 8.75% of muscles are covered, and against the OrthoLoad
  in-vivo anchor the EMG-driven estimate is *further* off at the knee (only marginally closer at the hip) — a
  real, decorrelated check, and a mixed, not flattering, result.
- **Fragmentation is mixed, not resolved.** `MECHANISM_UNIFIED_MODEL.md` successfully merged 6 forks into one
  232-muscle/84-ligament/27-body instance — genuine engineering progress — and, encouragingly, the new spine
  ligaments were built **directly on that merged model** rather than beside it (the first real sign the
  unified model is accumulating, not just a one-time merge-and-abandon artifact). But the hand build (index
  finger, then +middle finger/+2 intrinsics) and the new hip/ankle ligaments both still live in **further,
  separate, unmerged** fork(s) built the same evening — confirmed directly: `MECHANISM_UNIFIED_MODEL.md`
  contains zero mentions of hand/wrist/finger. Net: the absolute number of live, mutually-incompatible model
  files is still larger now than at audit time, even though one branch (spine) is now compounding correctly.
- **"Unified model" is not a hand claim.** The task-framing "hand (anatomical index-finger / unified model)" reads
  as one bullet but these are two **unrelated** builds — flagged explicitly so it isn't silently conflated.

---

## 1. Per-layer table

Tier legend (verbatim from the ledger — reused, not redefined): **in-vivo-anchored** (real instrumented
human telemetry) > **cadaveric-or-published-plausibility** (real external/decorrelated literature figure,
not a direct in-vivo force measurement of this twin's own quantity) > **method-only-no-external-anchor**
(rigorous internally, no external comparator — a real PMID for donor *geometry* does not count) >
**diagnosed-gap** (the doc itself reports a structural limitation/honest negative — a PASS for discipline,
not a claim to trust yet).

### Skeleton

| Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|
| Donor template (LaiArnold/Rajagopal) linearly rescaled to subject2's own markers; mass preserved exactly (78.2–78.976 kg); IK tracks the subject's own captured mocap to median 0.055°/max 1.328°. **Ceiling #1 re-confirmed structural, not incidental**: `MECHANISM_SUBJECT_SPECIFIC.md` traced OpenSim's C++ `Model::scale()` directly — geometry + gross mass are rescaled; muscle/ligament **strength** is referenced nowhere in the scaling code path, on ANY tissue class (0/80 muscles differ, 0/480 cross-method mismatches). | method-only-no-external-anchor (IK self-consistency against the same subject's own mocap, not an external population) | N/A (geometry pipeline, pre-dates any SO solve) | Real subject imaging (CT/MRI/DXA) — zero exists. Cheapest real step: Handsfield-et-al-style anthropometry muscle-volume regression (no new acquisition needed). |

### Joints

| Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|
| Sub-1.5° IK; multi-segment foot (ankle+subtalar+mtp) out of the box; midtarsal prototype **fails held-out generalization** (max 5.76° > pre-registered 5.0° ceiling); knee is a prescribed 1-DOF joint; wrist is 0-DOF in the base model (2-DOF only in the separate hand-graft fork); scapulothoracic+SC chain grafted but the AC loop cannot close (~12.6 cm gap after exhaustive search); **new finding**: the unified model's first full marker-driven IK trial through the scapula fork surfaced a pre-existing **~40–55° rotational registration offset** at the shoulder (right arm only, confirmed not a merge artifact); the `InverseDynamicsTool` bug (misreports knee/hip generalized force by ~1000–1800× at the literal zero pose) is now independently re-persisted (`neutral_pose_idbug_results.json`: knee 1826×, hip 1018×). | diagnosed-gap (multiple structural/software-correctness findings, several self-fired falsifiers) | N/A (kinematics-only; independent of SO) | The AC/shoulder registration needs re-derivation at the source (`add_scapula_clavicle.py`) before right-arm IK from the unified model is usable for any force analysis. |

### Muscles

| Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|
| 80 native + grafts (erector spinae +88, arm/shoulder +25, scapula stabilizers +11, trunk flexors +28 → 232 total in the merged unified model). 0/80 native flagged outside physiological bounds. Force-generation validated only through the joint-force family: knee/hip over-predict OrthoLoad in-vivo median by 1.412–1.515× (see Joint-force family below). Fmax/pennation confirmed **bit-identical** generic-vs-scaled (never subject-rescaled) two independent ways. Two Fmax "subject-specific" proxies prototyped and tested against the knee/hip in-vivo anchor: neither is a validated improvement — the "mass-informed" proxy's apparent knee improvement (1.515×→1.348×) is coincidental (the proxy itself carries no real per-segment signal for this subject); the geometrically-defensible proxy makes the knee estimate worse (1.515×→1.668×). | diagnosed-gap (subject-specific strength remains unsolved; prototype explicitly non-validated) for the strength question; in-vivo-anchored for the underlying force magnitude (shared with Joint-force family) | YES — reruns the full SO+dual-method contact-force pipeline | Ultrasound fascicle-length/pennation measurement per muscle (cheap, no MRI) — the ranked next step named in `MECHANISM_SUBJECT_SPECIFIC.md` §5. |

### Tendons

| Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|
| Millard tendon curve; Achilles peak stored elastic energy 7.76 J = 20% of Lichtwark & Wilson 2005's measured one-legged-hopping ceiling (38 J, cross-activity not same-task); Fukunaga 2001 walking-strain anchor (12.13/9.30mm vs ~7mm). Curve shape is the OpenSim class **default**, identical across every muscle checked (never ultrasound-measured); zero hysteresis; wrap surfaces not ported for any grafted muscle (moment-arm sign flips found in shoulder). | cadaveric-or-published-plausibility (ledger, unchanged) | YES (walking1 SO) | Ultrasound-measured subject-specific tendon compliance under a known load. |

### Ligaments (knee / hip / ankle / spine)

| Joint | Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|---|
| **Knee** | 84 `Blankevoort1991Ligament` bundles, real digitized donor (Lenhart2015, PMID 25917122); engagement pattern (ACL-tight-extension/PCL-tight-flexion) PASSES; MCL/LCL confirmed antagonists. **Structural mismatch, reproduces on the merged unified model**: MCL strain reaches non-physiological 58.2% at 140° flexion — the free-6-DOF donor vs prescribed-1-DOF target-knee incompatibility. Absolute stiffness never independently scaled. | method-only-no-external-anchor (ledger, unchanged) | NO (static geometric sweep) | Scale stiffness by cross-sectional area, compare vs a cadaveric/Grand-Challenge dataset. |
| **Hip** (NEW) | 4 bundles/side (iliofemoral ×2 arms, pubofemoral, ischiofemoral), same force class. Engagement pattern vs Martin 2008 (real cadaveric/arthroscopic study, PMID 18237703): ILFL tighter in extension (300N vs 17N, PASS), ISFL controls internal rotation both flexion levels (PASS), PFL controls external rotation in extension (PASS on the fairer normalized metric, FAILS on raw-Newton gap — both reported). 3377 tension-only checks, 0 violations. **Explicitly flagged weaker than the knee**: no digitized cadaveric geometry exists for the hip anywhere found this session — attachment points are geometric constructions (real joint centers + literature-named directions), not digitized coordinates; no absolute stiffness table recovered (ordinal-only literature). | cadaveric-or-published-plausibility (qualitative engagement pattern only) | NO (static sweep) | A digitized hip-ligament donor model (none located; searched live, genuinely absent from the JAM/opensim-core-jam ecosystem). |
| **Ankle** (NEW) | 7 bundles/side (ATFL/CFL/PTFL, deltoid ×4 bands). ATFL tighter in plantarflexion (PASS), PTFL tighter in dorsiflexion (PASS), CFL carries the inversion-restraint verification since ATFL/PTFL structurally cannot respond to isolated subtalar rotation (a genuine, pre-registered kinematic fact, not a bug). Two real bugs caught+fixed during the build (inversion sign backwards on first attempt; 3 deltoid bands 90–127mm, ~3–5× too long from a wrong reference point) — both diagnosed via the model's own checks, not asserted. | cadaveric-or-published-plausibility | NO (static sweep) | Same as hip: no digitized ankle-ligament donor exists; navicular body is missing entirely (2 deltoid bands proxied onto calcaneus, the weakest single approximation in this build). |
| **Spine** (NEW — appeared mid-write, see cutoff note) | 5 lumped bundles (ALL/PLL/LF/ISL/SSL), built directly on **`subject2_unified.osim`** (the merged model — a genuine positive: this is the first new fork built ON the merge rather than beside it). Attachment points derived from the donor's real bundled lumbar1/sacrum STL mesh vertices (real geometry, not literature-recalled dimensions) — geometrically stronger than the hip/ankle build in this one respect. 0/455 tension-only violations; all 5 PASS their engagement-direction threshold (ALL tighter in extension, PLL/LF/ISL/SSL tighter in flexion) — but PLL only after a disclosed, 2-round, live-forced OODA correction (a first geometric placement and a first reference-strain both failed before the fix). Same gap as hip/ankle: no absolute stiffness table recoverable (Pintar 1992 paywalled, 8 independent fetch attempts exhausted); generic representative values used, ordered only by relative robustness. | method-only-no-external-anchor — same style of self-consistency check as the knee-ligament build (engagement DIRECTION verified, absolute magnitude not); unlike hip/ankle, no specific external study's own quantitative/directional finding is tested against here (the cited PMIDs are provenance/corroboration-of-disagreement, not a functional anchor) | NO (static geometric sweep) | Same as hip/ankle: no digitized lumbar-ligament donor model exists anywhere searched; a real per-vertebra attachment source would replace the current single-lumped-bundle-per-ligament simplification. |

### Cartilage (knee)

| Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|
| Two separate models, still unmerged with anything else and with each other: **(1)** JAM/COMAK "DM" implant model — predicted vs Grand-Challenge in-vivo contact force, ~3% error on its best-of-3 trial (og1; honest range 3-10% across the 3 trials run) — previously headlined here as the single strongest cert in the repo, a superlative a 2026-07-26 domain-default audit reverses: OrthoLoad's independent population-median baseline beats the model on mean absolute error across the 3 trials (3.06% vs the model's 5.72%) and wins 2 of 3 held-out trials (og1, og3), the model winning only og4 (`AUDIT-VOIDFLOOR-DOMAINDEFAULT-KNEE-SARCOMERE-2026-07-26`, `AUDIT-KNEE-BODYWEIGHT-WINDOW-CORRECTION-OG1-SELECTION-2026-07-26`) — still on a TKA-implant, zero-ACL, non-subject2 individual. **(2)** NEW elastic-foundation contact (`MECHANISM_CARTILAGE_CONTACT.md`, ledger-scored): TF 2.11–4.78 MPa / PF 1.43–6.13 MPa, 16/16 machine cross-checks PASS, anchored against the SAME donor's (Lenhart2015) own real gait/running MPa figures (one-sided physical-bound test: both ratios <1, correctly). This model uses the **same lenhart2015.osim source** as the knee-ligament donor (same lineage, not merged into one file), and is **still not subject2** — passive-flexion only, not loaded gait. | cadaveric-or-published-plausibility (ledger, unchanged) | NO | Couple a native (non-implant) cartilage-contact model into subject2's own muscle+ligament instance; needs MRI-segmented healthy-joint geometry. |

### Neural

| Sub-layer | Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|---|
| **Motor-unit** | Fuglevand-style recruitment pool (soleus, N=120): all pre-registered gates PASS (size-principle order holds above a theory-derived resolution floor, 0/140,000+ pairs inverted; force-curve convexity 8–11σ above a 200-shuffle null). Anti-onion-skin firing direction matches Oya 2009 (Spearman −0.965). Standalone, no `opensim` import. | method-only-no-external-anchor (ledger, unchanged — Oya anchors only the qualitative direction of a design choice) | NO | Wire into the twin's actual activation output (currently fully decoupled). |
| **EMD** (electromechanical delay, NEW) | EMD_classical (activation+coupling) 5.85ms mean, EMD_full_chain 30.67ms — anchored to Nordez et al. 2009 (real human EMG-evoked EMD, PMID 19359617). Forced OODA: Gate B failed on 2/3 muscles at face value, re-diagnosed and fixed (Gate B′, 3/3 PASS). | cadaveric-or-published-plausibility | NO (isolated single-muscle ablation) | Subject-specific tendon-compliance measurement (the largest lever on the delay). |
| **Proprioception / spindle** (NEW) | Ia afferent firing rate (mean 96.2, peak 130.2 imp/s) via Prochazka & Gorassini 1998's published regression (r²=0.91) applied to the twin's simulated fiber length/velocity — more than a bare rate signal, but **not mechanistic**: no intrafusal bag1/bag2/chain fibers, gamma-drive static/fixed, **open-loop by the doc's own statement** ("this script only senses, it does not close any loop"). Real anchor is cross-species AND cross-muscle-group (cat hamstring → human soleus). 5/7 gates PASS. | cadaveric-or-published-plausibility | YES (reads SO's `activation.sto`) | Implement Mileusnic et al. 2006's intrafusal model (already calibrated on this muscle group per the doc's own next step). |
| **EMG-driven muscle force** (NEW — the one partial break from the SO common-mode) | Real surface EMG (7/80 muscles, subject2/walking1) drives actual contraction dynamics, compared against SO 3 ways (method-controlled). Real co-contraction found that SO structurally cannot see: gated CCI 0.391 (EMG) vs 0.190 (SO), diff exceeds the pre-registered 0.10 threshold; 32/158 frames EMG-active/SO-silent vs 6/158 reverse. Joint-force effect: **+48.4% at the knee (material), −12.2% at the hip (modest, opposite direction)**. Vs OrthoLoad: EMG-hybrid is FURTHER from the anchor at the knee, marginally CLOSER at the hip than SO — a genuine, non-uniform, not-uniformly-flattering result. Two real bugs found+fixed first (a numerical-integrator pathology; a 5-frame EMG lead-in padding artifact in the raw data, confirmed across all 3 trials). | in-vivo-anchored (via the same OrthoLoad comparison) | PARTIAL — 73/80 muscles (91.2%) still fall back to SO-activation-as-excitation even in this "hybrid" construction | Extending real EMG coverage past 7/80 muscles; a true MVC calibration trial (none exists — normalization is the doc's own disclosed dominant uncertainty). |

### Sensory (nociception)

| Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|
| Three channels: ligament-strain (176 events, 96 trusted/80 flagged-artifact, threshold anchored to Provenzano 2002 subfailure-strain data, PMID 11744679 — real); bone-stress (54.968/190 MPa = 0.289, correctly silent, no false trigger vs Torzilli 1999/Reilly 1975); joint-force (reuses knee 391.115%BW vs OrthoLoad 258.221%BW). **The single flashiest number (MCL 58.2% strain at 140°) is explicitly self-flagged by the doc as a known kinematic artifact, not physiology** — an honest disclosure, not a hidden flaw. | cadaveric-or-published-plausibility for the non-artifact ligament crossings + correctly-silent channels; diagnosed-gap for the 140° number specifically | PARTIAL (bone+joint-force channels trace to the SO solve; ligament channel is an independent live kinematic sweep; cartilage channel uses a separate passive-flexion trial) | A fitted firing-rate/recruitment curve (currently a linear overshoot ratio only). |

### Energetics (metabolic / thermal / muscle-perfusion / respiratory / cardiac)

| Sub-layer | Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|---|
| **Metabolic** | Umberger2010 7.683 / Bhargava2004 6.198 J/kg/m; twin reads 1.57–1.95× higher than Koelewijn's published level-walking cost of transport. | cadaveric-or-published-plausibility (ledger, unchanged) | YES (SO activation, trunk of the whole energetics chain) | Independent confirmation the gap is real physiology, not a model artifact. |
| **Metabolic, cross-activity** (NEW) | Squat/STS metabolic rate vs Nakagata et al. 2019/2020 (real, PMID 30934628/32379233): squat 1.545–1.603× above walking, STS 0.834–0.852× below; per-rep ratios 0.25–0.76×. **Doc's own governing verdict overrides the anchor comparison: "NOT independently trustworthy"** — squat/STS SO solves are separately flagged implausible upstream. | diagnosed-gap (self-downgraded) | PARTIAL (own SO solves per activity; walking baseline reuses the shared solve) | Re-solve squat/STS SO fixing muscle-pinning/reserve-saturation before trusting the anchor comparison. |
| **Thermal** (NEW) | Combined-corrected metabolic-heat config (400.4W) is the only one inside Malchaire 2006's own measured validity domain (100–450W, n=672 real lab dataset, PMID 16922181); sweat-rate/core-temp derived from it. | cadaveric-or-published-plausibility | YES (via `metabolic_cost.py`) | Real environmental+perfusion heat-dissipation term (currently a zero-dissipation adiabatic bound). Doc's own words: "physiologically unsustainable if taken literally" without it. |
| **Muscle-perfusion** (NEW) | (a) Magnitude: all 80 muscles stay below the 250 ml/min/100g ceiling (max 131.7, median 27.5), anchored to Andersen & Saltin 1985 / Joyner & Casey 2015 (real). (b) Fine-grained activation-vs-flow phase-locking: **overall_pass = False** — a 3-attempt forced OODA still fails the stringent per-muscle timing gate (only 41.2% of muscles pass a shuffle-null margin). Explicitly no cardiac-output constraint on its own (each muscle modeled independently) — this gap is what the new CARDIAC doc partially addresses. | (a) cadaveric-or-published-plausibility; (b) diagnosed-gap | YES (reads the same `activation.sto` directly) | Subject-matched NIRS/Doppler flow measurement (the doc's own named, non-existent falsifier). |
| **Respiratory** (NEW) | Combined-corrected VO2 = 14.698 mL/kg/min PASSES Ainsworth 2011's published range [12,15]; both primary metabolic configs FAIL (1.5–1.8× over). | cadaveric-or-published-plausibility | YES, explicitly self-disclosed (quoted in Headline Frame above) | Dead-space/alveolar decomposition + breath-by-breath VO2 kinetics. |
| **Cardiac** (NEW — appeared during this audit, confirmed genuine by mtime, not on any prior scan) | Fick-chain (Q=VO2/a-vO2diff, HR=Q/SV), explicitly 0-D steady-state, no pulsatile/valve/Frank-Starling mechanics. **REST vs Higginbotham et al. 1986 (real right-heart catheterization + radionuclide angiography, n=24, PMID 3948345): Q 5.65 vs 5.70 L/min (0.8% diff), HR 72.6 vs 73 bpm (0.6% diff)** — tight, genuine in-vivo agreement. WALKING: Q robustly clears the task anchor across the full sweep; HR clears only at the upper edge of the pre-registered a-vO2diff range (117.7 vs [90,110] at the midpoint) — disclosed as open modeling uncertainty, not hidden. Void-floor forced test: pinning a-vO2diff at its resting value demands HR ≥219–401 bpm at mere walking effort — physiologically impossible at any config, confirming the widening mechanism is load-bearing, not decorative. | in-vivo-anchored (rest tight; walking HR partial, disclosed) | YES — reads `metabolic_cost_results.json` directly, no re-solve | A subject-specific CPET or echocardiography measurement (stroke volume is still a generic population interpolation). |

### Tissue (bone-stress / bone-remodeling / wobbling-mass)

| Sub-layer | Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|---|
| **Bone-stress** (NEW) | Peak stress −54.968 MPa (axial −3.13, bending 51.84). Passes only the **task-given, non-literature, explicitly-non-gating** 40–60 MPa band. The two REAL literature anchors fare worse: Duda 1998's converted band (~34–51 MPa) is near-missed; **Taylor 1996's compression-dominated finding is directly contradicted** — this model is bending-dominated (doc's own words: "IN TENSION with Taylor 1996's finding, not agreement" — independently confirmed by reading the doc directly, line 154). Own §6.1 quantifies a distal/proximal bending-moment discrepancy (148 N·m) **larger than the reported bending moment itself** (93 N·m). | cadaveric-or-published-plausibility, but weak/mixed (one anchor contradicted, one near-missed) | YES, maximal (reuses SO `force.sto` + JR `ReactionLoads.sto` directly) | A rotational-inertia correction (I·α + ω×Iω) — the doc's own named next step, and the most likely source of the bending-moment discrepancy. |
| **Bone-remodeling** (NEW) | Medial/lateral/axial-only microstrain (3233/2865/184 µε) classified against real Frost 1987 thresholds (1500–3000 µε, PMID 3688455) as "formation, not maintenance." **Doc's own §4, titled "the critical caveat": this verdict is likely biased HIGH, not settled** — would flip to resorption if bone-stress's own disclosed bending/Taylor-1996 defect were resolved. The doc states outright it "cannot adjudicate... inherits, rather than resolves" the upstream defect. | diagnosed-gap (self-undermined by its own §4) | YES, 100% (pure classifier on bone-stress's output, recomputes nothing upstream) | Identical to bone-stress's own fix — this layer's fidelity is fully gated by that one upstream correction. |
| **Wobbling-mass** (NEW) | Identity/void-floor self-tests pass exactly. Anchored to Pain & Challis 2006 (real drop-landing FE data, PMID 16271595, "~50% lower" — same direction as the surviving number). **Only the supplementary DJ1 path (+19.79%) clears both its own validity gate and the external anchor** — the PRIMARY 100Hz-IK path fails its own validity gate (2–3.5× over ceiling); the walking1 primary result is sub-floor/sign-unstable. Doc's own verdict string: "MECHANISM-VERIFIED-SUPPLEMENTARY-DATA-ONLY-PARAMETER-SENSITIVE" — not a clean win. | cadaveric-or-published-plausibility (weakest member of this tier — only the supplementary path actually clears) | NO — driven by IK/GRF kinematics, not SO activations (a genuine decorrelation from the dominant common-mode chain) | Real per-segment ultrasound-measured soft-tissue mass/natural-frequency/damping (currently generic/swept). |

### Hand

| Sub-layer | Current fidelity | Tier | Common-mode? | One thing to advance it |
|---|---|---|---|---|
| **Anatomical hand** (index finger, then +middle finger/+2 intrinsics, NEW extension) | Index: 2 real forearm flexors + 1 extensor, wrist un-welded (0→2 DOF); 12/12 moment-arm cross-checks PASS; 11/12 vs a pre-registered literature band. Extension (`MECHANISM_FULL_HAND.md`): +middle finger (mirrors index method exactly, 3 more muscles) + bodytwin's **first 2 real intrinsic muscles** (1st lumbrical, 1st dorsal interosseous) — machine-verifies the textbook lumbrical action (MCP flexion **+** PIP extension from the SAME muscle, +8.891mm/−2.257mm). External anchor Nakajima et al. 2022 (real per-muscle table, PMC9345905): lumbrical 7.4% off, interosseous 44% off, both inside a wide pre-registered band. 17/17 moment-arm cross-checks PASS. The correct full-hand donor (McFarland 2023, 43 muscles) was **genuinely fetch-blocked** this session (JS-gated SimTK search, forced through a real multi-channel OODA search, not a one-shot give-up) — scoped from the paper's own accessible text instead. **Progress: 2/5 fingers have extrinsics, 2/~30+ intrinsics (both serving the index finger only), thumb entirely absent, no MCP ab/ad, no pulleys/ligaments/extensor mechanism.** | cadaveric-or-published-plausibility | NO (static geometric moment-arm verification) | Re-attempt the McFarland 2023 fetch if a JS-capable fetch tool becomes available — would supersede much of the current scaffolding at once. |
| ("Unified model" — NOT a hand claim, flagged to prevent conflation) | `MECHANISM_UNIFIED_MODEL.md` merges 6 OTHER forks (erector spinae, trunk flexors, arm muscles, scapula/clavicle, knee ligaments, foot multisegment) into one 232-muscle/84-ligament/27-body/329-force-element instance, all counts exact. Contains **zero** mentions of hand/wrist/finger (confirmed by direct grep) — the hand build is not merged into it and was never in scope for it. | method-only-no-external-anchor (ledger, unchanged) | N/A | Merging the hand fork into the unified model is unstarted work, not a near-term byproduct of either effort. |

---

## 2. Updated tier counts

| Tier | Ledger (66 docs/72 rows, ~18:25) | + newly scored this pass (19 rows, 18:25→~19:36) | **Combined (85 docs / 91 rows)** |
|---|---:|---:|---:|
| in-vivo-anchored | 8 (11%) | 2 (CARDIAC-rest, EMG_DRIVEN) | **10 (11%)** |
| cadaveric-or-published-plausibility | 11 (15%) | 10 | **21 (23%)** |
| method-only-no-external-anchor | 31 (43%) | 1 (spine ligaments) | **32 (35%)** |
| diagnosed-gap | 12 (17%) | 6 | **18 (20%)** |
| N/A-process | 10 (14%) | 0 | **10 (11%)** |
| **Total** | **72** | **19** | **91** |

Reading this honestly: method-only is still the single **largest** tier in absolute terms (32 rows) even
though its *share* fell from 43% to 35%, because this wave of new docs happened to find real literature
anchors more often than not (thermal/respiratory/cardiac/EMD/proprioception/hip-ankle-ligament/hand-extension
all landed on genuine, functionally-comparable PMIDs). That is a real, measured broadening of external
anchoring — **not** the same thing as validation getting deeper: most of those same new anchors are
population-literature comparisons (cadaveric-or-published-plausibility, one step below in-vivo), several
explicitly self-downgrade (diagnosed-gap rose from 17% to 20%: vasculature timing, metabolic cross-activity,
nociception's flashy number, bone-remodeling, RRA-task-gains, subject-specific-strength), and the two new
in-vivo-anchored rows are themselves partial (EMG-driven covers 7/80 muscles; cardiac's HR only partially
clears its own pre-registered walking band). Notably, the spine-ligament build — despite citing 6 real PMIDs —
still landed method-only-no-external-anchor, the same tier as the knee-ligament precedent it was built to
match: citing real literature for donor provenance or corroboration is not the same as testing against a
specific external study's own quantitative finding (which is exactly what separated the hip/ankle-ligament
build's cadaveric-or-published-plausibility tier, via Martin 2008/Siegler 1988's own named results).

## 3. The honest framing, stated plainly

**Broad, not deep. One shared solve, not many independent confirmations.** The twin now has SOME
representation touching nearly every named layer — skeleton, joints, muscles, tendons, knee+hip+ankle
ligaments, knee cartilage, motor-unit/EMD/proprioception, nociception, metabolic/thermal/perfusion/
respiratory/cardiac energetics, bone-stress/remodeling/wobbling-mass, and a two-finger hand. That breadth is
real and newly measured here. But:

1. **43%→35% method-only is still the modal tier** — most of what exists has never been compared to anything
   outside this twin's own pipeline.
2. **The energetics cluster is a tree, not four branches**: metabolic_cost.py is the trunk; thermal,
   respiratory, and cardiac are downstream consumers of its one number, not three separate measurements
   (the respiratory doc says this about itself).
3. **Errors propagate downstream, concretely, not hypothetically**: bone-stress's disclosed contradiction with
   Taylor 1996 flows straight into bone-remodeling's self-graded "not settled" verdict.
4. **Composability did not resolve**: one successful 6-fork merge coexists with a growing number of further,
   still-separate forks (hip/ankle ligaments, the 2-stage hand build) built the same evening.
5. **The one real crack in the common-mode wall (EMG-driven, 7/80 muscles) produced a mixed, not a flattering,
   result** — exactly the kind of finding this audit exists to surface, not soften.

A map that read "the twin is nearly complete" would be false. The accurate read is: **wide coverage, shallow
and heavily common-mode validation, with the project's own discipline (forced adversaries, disclosed
contradictions, honest self-downgrades) doing real, machine-checked work catching this — which is itself the
strongest asset here, stronger than any single layer's fidelity number.**

## Evidence index

Read directly this session: `docs/MECHANISM_TRUST_LEDGER.md`, `docs/MECHANISM_TRUST_LEDGER_REMEDIATION.md`,
`docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md`, `docs/MECHANISM_UNIFIED_MODEL.md`,
`docs/MECHANISM_SUBJECT_SPECIFIC.md`, `docs/MECHANISM_CARTILAGE_CONTACT.md`, `docs/MECHANISM_MOTOR_UNIT.md`,
`docs/MECHANISM_JOINT_FORCE_SCORECARD.md`, `docs/MECHANISM_ANATOMICAL_HAND.md`, `docs/MECHANISM_CARDIAC.md`,
`docs/MECHANISM_HIP_ANKLE_LIGAMENTS.md`, `docs/MECHANISM_FULL_HAND.md`, `docs/MECHANISM_EMG_DRIVEN.md`,
`docs/MECHANISM_SPINE_LIGAMENTS.md`, plus direct spot-checks against
`data/msk_smoketest/subject2_walking1/static_optimization/static_opt_knee_results.json` and
`data/msk_smoketest/subject2_walking1/bone_stress/bone_stress_results.json`.
Cross-checked via two independent, QC'd research passes (each verified against raw on-disk JSON, spot-checked
against their most decisive claims before being trusted here): `docs/MECHANISM_THERMOREGULATION.md`,
`docs/MECHANISM_VASCULATURE.md`, `docs/MECHANISM_RESPIRATORY.md`, `docs/MECHANISM_METABOLIC_CROSS_ACTIVITY.md`,
`docs/MECHANISM_RRA_TASK_GAINS.md`, `docs/MECHANISM_EMD.md`, `docs/MECHANISM_PROPRIOCEPTION.md`,
`docs/MECHANISM_NOCICEPTION.md`, `docs/MECHANISM_WOBBLING_MASS.md`, `docs/MECHANISM_BONE_STRESS.md`,
`docs/MECHANISM_BONE_REMODELING.md`.

No git commit, no git push performed (isolation-respected, per `COORDINATOR.md` §1 and the task).
