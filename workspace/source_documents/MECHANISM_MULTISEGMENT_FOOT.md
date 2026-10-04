# MECHANISM MULTISEGMENT FOOT — Maharaj2021 JC model: characterization + validation + coupling scope (2026-07-21)

Characterizes the just-fetched `data/external/multisegment_foot/` corpus (Maharaj2021_BothLegs_V2.osim +
walk.trc/walk_ik.mot + SampleData_SCALED.osim, biplanar-videoradiography-validated, via operator SimTK
login), independently re-runs its IK, extracts multisegment foot kinematics, anchors to the model's own
live-verified in-vivo validation, and scopes (not builds) a coupling path into the subject2 mechanism.
Builds on, and **corrects a stale premise in**, the existing foot-fidelity sub-project
(`docs/MECHANISM_FOOT_FIDELITY_PLAN.md`, `docs/MECHANISM_FOOT_MULTISEGMENT.md`) — read first for context, not
re-litigated here.

**Correction up front:** the twin's foot is *not* currently a single rigid segment — a prior dial-turn
(`docs/MECHANISM_FOOT_MULTISEGMENT.md`) already merged a `midtarsal_r/l` joint into
`data/msk_models/subject2_unified_v2.osim` (confirmed live this session: `ankle→subtalar→midtarsal→mtp`, 4
joints/side). That joint runs on a **placeholder axis** (no cited numeric angle was available at build time)
and **failed its own held-out generalization gate** (DJ2 max-abs 5.76° > 5.0° ceiling, per
`docs/MECHANISM_LAYER_COMPLETENESS.md`). The new Maharaj corpus matters because it supplies exactly the input
that gap flagged as missing — a **real, cited, oblique midtarsal axis** (extracted below) — not because the
twin currently has zero foot detail.

---

## 1. Model characterization — machine-enumerated via OpenSim 4.6 Python API

**V1-vs-V2 bug found in transit**: `Maharaj2021_BothLegs.osim` (V1) has `femur_r/l`/`tibia_r/l` mass+inertia
**zeroed** (total mass 49.31 kg vs V2's 74.44 kg; the 25.13 kg gap = exactly 2×(femur 8.984 kg + tibia 3.581
kg)). V2 is the corrected file — used throughout, per the task's own instruction to load V2.

**Foot topology, per side** (talus→calcn→midfoot→forefoot→digits, all `CustomJoint`):

| Joint | DOF | Coordinates (range) |
|---|---|---|
| ankle | 3 | ankle_flexion [−60,60]°, ankle_in_ev [−5,5]°, ankle_rot [−5,5]° |
| subtalar | 1 | subtalar_angle [−30,30]° |
| midtarsal | 1 | midtarsal_angle [−30,30]° |
| tarsometatarsal | 1 | tarsometatarsal_angle [−30,30]° |
| digits (mtp) | 1 | mtp_angle [−60,30]° |

**Total: 5 bodies, 7 DOF/side** — matches the paper's own abstract ("6-segments, and 7 degrees of freedom")
**exactly** (6 segments = tibia + 5 foot bones), confirming this file is genuinely the JC model the paper
describes, not just a same-named download.

**Markers actually driving this session's IK** (`SampleData_SCALED.osim` + `walk.trc`, 42 markers total —
**not** the 56-marker generic template, see correction below): calcn {RDC,RLC,RMC}, midfoot
{RNAV,RMB1,RMB2,RMB5}, forefoot {RMH1,RMH2,RMH5}, digits {RD1} = **11 markers/foot**.

**Self-caught correction**: an initial pass counted 14 markers/foot by introspecting the *generic*
`Maharaj2021_BothLegs_V2.osim` (56 markers total, including `RPC/RD2/RD5/LPC/LD2/LD5/RS1-4/LS1-4` that are
**not** present in the actual scaled+walk-trial dataset). Re-derived from the file actually used for IK — 11,
not 14 — before it went into any downstream number.

### Twin vs. Maharaj — the concrete resolution gain

| | Twin (`subject2_unified_v2.osim`) | Maharaj2021 JC model |
|---|---|---|
| Bodies/side | 4 (talus, calcn, forefoot, toes) | 5 (talus, calcn, **midfoot**, forefoot, digits) |
| DOF/side | 4 (ankle×1, subtalar, midtarsal, mtp) | 7 (ankle×**3**, subtalar, midtarsal, **tarsometatarsal**, mtp) |
| Markers/foot | 3 (r_calc, r_toe, r_5meta — all on calcn; toes carries 0) | **11** (3.7×) |
| Midtarsal axis | placeholder (v3: Z rotated 90° about Y), no citation, **failed held-out gate** | real, oblique, extracted below |
| mtp/digits marker support | 0 (architecturally null, per FOOT_FIDELITY_PLAN.md §2) | 1 (RD1) — non-null but the weakest of the 4 new segments |

The 3.7× marker gain directly targets what `FOOT_FIDELITY_PLAN.md` §2 diagnosed as *the* limiting factor
("adding more foot segments...without adding markers...does not add resolvable detail" ) — this corpus is the
first candidate in this project's search that actually clears that bar for most (not all — see §4) of its new DOF.

---

## 2. Literature anchor — live-verified this session, not recalled

**PMID 34698598 / DOI 10.1080/10255842.2021.1968844** — Maharaj JN, Rainbow MJ, Cresswell AG, Kessler S,
Konow N, Gehring D, Lichtwark GA. "Modelling the complexity of the foot and ankle during human locomotion:
the development and validation of a multi-segment foot model using biplanar videoradiography." *Comput
Methods Biomech Biomed Engin*, epub **2021-10-26**, print 2022;25(5):554-565 (the epub-2021 date is why the
task calls this "Maharaj 2021"). Fetched live via PubMed this session (title/authors/journal/dates
cross-checked, not from training-data recall).

**Verbatim abstract quote (the model's own validation number):**
> "The model has 6-segments, and 7 degrees of freedom; motion between foot segments were constrained with a
> single oblique axis...[Joint Constrained (JC) model]... Compared to the 6-DoF model, our JC model
> demonstrated significantly smaller RMS differences [**JC: 2.19° (1.43–2.73)**; 6-DoF: 3.25° (1.37–5.89)]
> across walking and running" — i.e. JC-model-vs-biplanar-videoradiography-truth RMS error.

**Full text: paywalled** — independently reconfirmed this session (HTTP 403, `tandfonline.com/doi/full/...`),
matching (not just copying) the prior session's own finding. **No per-joint ROM/arch-mm numbers are
accessible live** — only the RMS-accuracy figures above. Also checked: PMID 28783556 (Maharaj/Cresswell/
Lichtwark, *Gait Posture* 2017, the earlier model-description paper) — confirms "a custom multi-segment foot
model in OpenSim" existed, but its own abstract states "no model could be created for any of the kinematic
measures analysed" — not a ROM source either.

**The task brief's own ballpark** ("arch deformation ~few mm / ~5-10° midfoot motion") could **not** be
independently confirmed in any accessible text this session (WebSearch quota exhausted, full text paywalled)
— treated below as an **unverified prior**, compared against (not assumed to validate) the measured numbers.

---

## 3. IK re-solve — forced OODA on a real discrepancy, not a one-shot accept

Ran `InverseKinematicsTool` independently on `walk.trc` using the model's own scaled instance
(`SampleData_SCALED.osim`, 90 kg) and its own marker task/weight set (`Maharaj2021_BothLegs_IK_Tasks.xml`).

**Observed**: default accuracy (1e-5) hits an IPOPT max-iteration failure near the trial's end (t>2.765s of
2.925s) — a **gradual** RLC marker-residual climb over ~10 consecutive frames (not a single-frame glitch),
diagnosed as genuine late-trial marker-quality decline.

**First pass (relaxed accuracy=1e-4, full trial) vs. the provided `walk_ik.mot`**: large cross-check RMSE
(mtp 7.6°, midtarsal 4.2°, subtalar 3.2–3.9°) — suspicious enough to force a controlled test rather than
accept it.

**Forced diagnostic (OODA, not a shrug)**: re-ran at matched tight tolerance (1e-5, the tool's true default)
over the identical converged window (t=0–2.70s) and compared to the *same* window of the provided file.
Result: RMSE roughly **halves** (e.g. mtp 7.8°→5.2°, subtalar 3.07°→1.47°) and signal-roughness ratios
(mine/provided) collapse from **0.04–0.41** to **0.99–1.48** — most of the initial mismatch was my own
tolerance choice, not a real discrepancy. **All numbers below use accuracy=1e-5**, truncated to the
confirmed-converged window (554/586 frames, 94.5% of the trial).

**Residual cross-check vs. provided `walk_ik.mot`, matched tight tolerance (RMSE, degrees):**

| Coordinate | RMSE (r/l) | vs. paper's own RMS-vs-truth ceiling (1.43–2.73°) |
|---|---|---|
| ankle_flexion | 0.79 / 0.66 | well under |
| ankle_in_ev | 0.35 / 0.33 | well under |
| ankle_rot | 0.41 / 0.51 | well under |
| subtalar_angle | 1.47 / 1.43 | **at the floor of the band** |
| midtarsal_angle | 2.08 / 2.11 | inside the band |
| tarsometatarsal_angle | 1.89 / 1.69 | inside the band |
| **mtp_angle** | **5.13 / 2.54** | **right side exceeds the worst-case (2.73°); left inside** |

**Mechanistic (not coincidental) explanation**: mtp/digits carries only **1** informing marker (RD1) vs 3–4
for every other new segment — the least-observed new DOF is the least solver-reproducible. This is the same
σ_min/observability signature already established for the twin's own architecturally-null mtp (§1's marker
comparison) — just far less severe here because a real marker exists at all.

---

## 4. Kinematics results (machine-computed, 554/586 frames, both legs)

**Whole-window ROM (°):**

| Coordinate | Right | Left |
|---|---|---|
| ankle_flexion | 28.9 (−7.4 to 21.4) | 29.0 (−4.8 to 24.1) |
| ankle_in_ev | 3.8 | 5.7 |
| ankle_rot | 5.4 | 6.7 |
| subtalar_angle | 23.3 | 30.1 |
| **midtarsal_angle** | **17.0** (−2.7 to 14.2) | **22.3** (−3.8 to 18.4) |
| tarsometatarsal_angle | 13.1 | 15.4 |
| mtp_angle (raw, offset not anatomical zero) | 33.5 | 30.6 |

**Medial longitudinal arch proxy** (perpendicular height of navicular marker above the
calcaneus–1st-metatarsal-head chord, sagittal plane, via forward-kinematics marker replay of the resolved
IK motion — a standard navicular-height/arch-drop construction): whole-trial range **12.2 mm (right)**,
**13.7 mm (left)**; min ~2.8–4.7 mm, max ~17.0 mm.

**Gait events**: Zeni-style relative-marker (foot-marker-x minus pelvis-x local extrema) heel-strike/toe-off
detection, self-validated (not tuned to any external number) via stance-fraction convergence to **63–66%**
(canonical normal-walking ≈60–62%) and tight stride-time repeatability. **A boundary-extrema bug was caught
and fixed in this same step**: the naive extrema-finder treated trial-boundary frames as automatic local
extrema, producing a spurious 29%-stance pseudo-stride at t=0 on the left side (right side was coincidentally
unaffected). Fixed version requires a genuine interior turning point; all numbers below use the fixed version.

**Windlass MTP dorsiflexion excursion, heel-strike→push-off peak** (offset-invariant — the raw coordinate's
zero is not confirmed to be an anatomical neutral, so the *excursion* is reported, not the absolute angle):

| Side | n strides | Excursion (°) |
|---|---|---|
| Right | 1 | 11.9 |
| Left | 2 | 3.2, 2.1 |
| **Pooled** | **3** | **2.1 – 11.9** |

Sample-size caveat (stated, not buried): only 1 (right) / 2 (left) strides survive the tight-tolerance
truncation + bug fix. The right-vs-left difference is **not** asserted as a genuine physiological asymmetry
— too few strides to separate from noise; more strides (untruncated trial or additional trials) would be
needed. **Stance-phase arch compression is far more bilaterally consistent** in the same strides: right 10.1
mm, left 10.5–11.0 mm.

**Sign-convention independence check**: a marker-geometry-only (model-coordinate-free) sagittal toe-vs-
forefoot angle (RD1−RMH1 vector vs. RMH1−RDC vector) correlates **r=0.78** with the model's own `mtp_angle_r`
— confirms the coordinate tracks real toe-relative-to-forefoot rotation, not a decoupled/meaningless DOF, and
independently fixes the sign used to define "dorsiflexion" above.

---

## 5. Falsifier verdict

**Pre-registered threshold**: solver-reproducibility RMSE (a lower-bound proxy for true identifiability) at
or below the paper's own published worst-case RMS-vs-biplanar-videoradiography-truth (2.73°) counts as
"usable at the model's own designed precision" for that DOF.

**Forced adversary** (leaning positive — the tempting skip was "it loaded, produced smooth periodic curves,
ship it"): forced via (a) two tolerance settings compared under a controlled same-window test rather than
accepting the first run, (b) cross-check against a pre-existing file I did not produce, (c) a
convention-independent second instrument for mtp specifically, (d) an internally-consistent, untuned
gait-event detector, (e) a caught-and-fixed bug in that same detector rather than reporting its first output.

**Result: 6 of 7 per-side DOF PASS** (ankle×3, subtalar, midtarsal, tarsometatarsal — RMSE 0.33–2.11°, at or
under the paper's own accuracy ceiling). **mtp is a diagnosed, bounded exception** (5.13° right / 2.54° left
— mechanistically explained by having only 1 informing marker, not an unexplained anomaly).

**What was NOT tested** (disclosed, not glossed over): the paper's specific per-joint ROM numbers are
inaccessible (paywalled) — I cannot directly confirm my measured 17.0–22.3° midtarsal ROM or 2.1–11.9°
windlass excursion against the paper's own reported ranges. What I confirmed instead is that this model
**functions correctly, is genuinely the described JC model** (7-DOF match), and its **kinematic
reproducibility sits at-or-inside its own published in-vivo accuracy band** for the well-marker-observed
joints. The task brief's own ROM ballpark (~5–10° midfoot) is an unverified prior my measured midtarsal ROM
(17–22°) sits 2–3× above — flagged as a discrepancy against an *unverified* number, not a failure against a
confirmed citation.

**CONFIDENCE: in-vivo-anchored (biplanar videoradiography, live-verified PMID/DOI) for the model's own
published accuracy ceiling and the 6/7-DOF reproducibility-vs-that-ceiling check; OPEN for direct per-joint
ROM literature matching (inaccessible); diagnosed-bounded-exception for mtp, not a silent gap.**

---

## 6. Coupling scope — a path, not a build

Three tiers, cheapest/highest-value first. None executed this session (scoping only, per the task).

**Tier 1 — axis-transplant only (cheapest, directly retires the #1 disclosed gap in `FOOT_MULTISEGMENT.md`).**
Extracted the real midtarsal joint axis from the Maharaj model: raw vector **[0.9392, 0.5666, −0.4183]**
(normalized **[0.800, 0.483, −0.356]**), expressed directly in `calcn_r`'s own body-local frame (both
`PhysicalOffsetFrame`s have zero orientation offset — no extra rotation to compose). This is a **real,
oblique, non-coordinate-aligned axis** from a peer-reviewed, biplanar-videoradiography-validated model —
independently corroborating `FOOT_MULTISEGMENT.md`'s own recalled-but-unverified claim that the true
transverse-tarsal axis is classically oblique, and concretely different from that doc's own placeholder v3
choice (pure Z rotated 90° about Y). This is exactly the artifact `FOOT_MULTISEGMENT.md` §5 flagged as
highest-value-missing ("a real, cited midtarsal axis angle would very plausibly close most of the §4 gap").
**One caveat before transplanting as-is**: whether Maharaj's `calcn_r` local-frame convention matches the
twin's own `calcn_r` local frame is not yet checked (different labs, different digitization pipelines, same
body name ≠ same local axes) — a single, cheap, geometric verification step, not a re-derivation. **Next
step**: substitute this axis into the twin's existing `midtarsal_r/l` joint and re-run
`scripts/msk/foot_multisegment.py`'s own already-built held-out DJ2 regression gate (previously FAILED at
5.76° vs the 5.0° ceiling) — a decisive, cheap, re-runnable test of whether a real axis closes that gap.

**Tier 2 — full skeletal graft** (same method as the existing 2 dial-turns, one step further): split
`forefoot` into `midfoot`+`forefoot`, upgrade `ankle` from 1-DOF to 3-DOF, add `tarsometatarsal`. Mechanically
proven feasible (the existing dial-turn scripts already did comparable topology surgery twice). **Does not
by itself solve the marker problem**: the twin's own mocap protocol only has 3 markers/foot
(`r_calc`/`r_toe`/`r_5meta`) vs. the 11 this model uses — new DOF beyond what's already graftable would
remain marker-starved on the twin's *own* capture data, landing back on the exact ceiling
`FOOT_FIDELITY_PLAN.md` §2 already diagnosed, unless a marker crosswalk or richer capture is added. Scoped
effort: comparable to the 2 existing dial-turns combined (mass-splitting, muscle re-routing, topology
surgery, held-out validation each).

**Tier 3 — full lower-limb re-scale/swap-in** (the `FOOT_FIDELITY_PLAN.md` §3 "Option A" already assessed
and deprioritized — that reasoning is unchanged, just re-pointed at a now-downloadable file instead of a
gated one): scale the *entire* Maharaj2021_BothLegs model to subject2's anthropometry via its own
`Scale_MeasurementSet`/`Scale_Tasks` XMLs, either as a standalone foot-specific analysis tool or as a full
hip-down replacement for the twin's LaiArnold lineage. Still hits the same marker-crosswalk problem (subject2's
own mocap/video doesn't carry `RDC/RNAV/RMB*/RMH*/RD1` either), plus the previously-documented GRF-partitioning
problem for any kinetic (not just kinematic) use, plus the Kim & Kipp 2019 diminishing-returns finding already
in `FOOT_FIDELITY_PLAN.md` (3-seg vs 5-seg mostly don't differ once you're past 2-seg) — none of these change
by having the file in hand. Not recommended as the next step; Tier 1 is.

---

## 7. Honest gaps

1. Per-joint ROM literature numbers from the Maharaj paper itself are inaccessible (paywalled) — Tier-1's
   axis-transplant recommendation rests on the abstract's RMS-accuracy figures + the model's own internal
   consistency, not a direct ROM-to-ROM literature match.
2. mtp/digits kinematics carry a diagnosed, larger (2.5–5.1°) solver-path uncertainty than the other 4 new
   DOF — usable as a topology/feasibility signal, not yet as a precise angle.
3. Windlass-excursion bilateral asymmetry (11.9° right vs. 2.1–3.2° left) is explicitly **not** claimed as
   real physiology — n=1/2 strides is too small; would need the untruncated trial or more trials.
4. The local-axis-frame-convention match between Maharaj's and the twin's `calcn_r` (needed before Tier 1's
   axis number can be transplanted as literal numbers, not just as directional corroboration) is unverified.
5. `JointKinematics.tif` ships with the data — not used for any numeric claim here (figures are forensic-only
   per this project's method discipline); all numbers above are machine-computed from `.mot`/`.trc` directly.
6. GRF/kinetic partitioning across the (now 5-segment) foot is untouched — this entire session is
   kinematics-only, matching the task's own scope.

---

## 8. File index

- `docs/MECHANISM_MULTISEGMENT_FOOT.md` — this document.
- `data/msk_smoketest/maharaj2021_foot_probe/MAHARAJ_FOOT_EVIDENCE.json` — consolidated, final evidence
  (provenance, model characterization, IK method + tolerance-diagnosis trail, kinematics results, coupling
  findings incl. the extracted axis vector).
- `data/msk_smoketest/maharaj2021_foot_probe/foot_kinematics_evidence.json` — raw/intermediate evidence
  (full audit trail incl. the superseded pre-bugfix gait-event numbers, kept for transparency).
- `data/msk_smoketest/maharaj2021_foot_probe/run_ik.py`, `run_ik_tight_truncated.py`, `run_ik_final.py` —
  re-runnable IK scripts (the three tolerance/window variants used to force the diagnosis in §3).
- `data/msk_smoketest/maharaj2021_foot_probe/analyze_foot_kinematics.py` — re-runnable kinematics analysis
  (ROM, arch proxy, gait events, push-off windlass excursion).
- `data/msk_smoketest/maharaj2021_foot_probe/walk_ik_FINAL.mot` — the trusted IK solve (554/586 frames,
  tight tolerance) all numbers in this doc derive from.
- `data/msk_smoketest/maharaj2021_foot_probe/marker_positions_ground.npy` — forward-kinematics marker
  positions in the ground frame (feeds the arch proxy + gait events).
- `data/external/multisegment_foot/Maharaj2021_BothLegs_V2.osim`, `SampleData_SCALED.osim`, `walk.trc`,
  `walk_ik.mot`, `Maharaj2021_BothLegs_IK_Tasks.xml` — the source corpus (read-only, unmodified).
- `docs/MECHANISM_FOOT_FIDELITY_PLAN.md`, `docs/MECHANISM_FOOT_MULTISEGMENT.md` — prerequisite reading, the
  existing dial-turns this session builds on (not re-litigated).
