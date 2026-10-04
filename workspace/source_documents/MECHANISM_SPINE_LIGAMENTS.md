# MECHANISM SPINE LIGAMENTS — lumbar tension-only restraints on the unified model (2026-07-21)

Extends the ligament tier (`docs/MECHANISM_KNEE_LIGAMENTS.md`, hip/ankle in progress) to the **lumbar
spine**: the 5 major lumbar ligaments the task named — anterior longitudinal (ALL), posterior
longitudinal (PLL), ligamentum flavum (LF), interspinous (ISL), supraspinous (SSL) — added as
`opensim.Blankevoort1991Ligament` tension-only passive restraints (same Tier-2 nonlinear force class
already validated for the 4 knee ligaments, MSK build charter). Script:
`scripts/msk/add_spine_ligaments.py` + `scripts/msk/spine_ligament_config.json` (the hotswappable
dial, regenerated fresh each run). Built on **`data/msk_models/subject2_unified.osim`** (the merged
model: 27 bodies/27 joints/232 muscles/84 knee-ligament bundles/13 reserves/3 constraints/49
markers — re-verified live at the top of the script, not assumed from the doc). Run:
`source_repository/.venv-msk/bin/python3 scripts/msk/add_spine_ligaments.py`.

**Isolation respected:** `subject2_unified.osim` is read-only input (never modified in place); output
is a NEW model file (`data/msk_smoketest/spine_ligaments/model_with_spine_ligaments.osim`); nothing
pushed/committed (bodytwin isolation, task instruction).

---

## 0. Headline

**5 `Blankevoort1991Ligament` bundles added (ALL, PLL, LF, ISL, SSL — one lumped bundle each, midline/
unpaired, unlike the bilateral knee ligaments). All 5 PASS tension-only verification (0 violations /
455 checks) and all 5 PASS their pre-registered flexion/extension engagement-direction threshold —
but PLL only passes after a real, two-round, live-forced OODA correction (§5), not on the first
attempt.** Force count 329→334 (232 muscles unchanged + 84 pre-existing knee ligaments unchanged + 5
new). Structural counts, slack-length derivation, and the force law were independently cross-checked
against the class's own documented closed form (max error ≤1.39e-17 N, floating-point noise).
**Honest headline caveat**: exact stiffness/reference-strain numbers could NOT be extracted from
Pintar et al. 1992's own table this session (paywalled, no accessible full text — 8 independent
live-fetch attempts documented in §6); generic, clearly-flagged representative values are used
instead. The qualitative engagement-direction claim — the task's actual falsifiable ask — is shown
to be governed by geometry (attachment point placement), not by the stiffness magnitude, so this gap
does not undermine the verified result.

---

## 1. Ligament force class used

**`opensim.Blankevoort1991Ligament`** — the same class already validated in this repo for the 4 knee
ligaments (`docs/MECHANISM_KNEE_LIGAMENTS.md`), reused here without modification. Force law
(quadratic toe region below `transition_strain`, linear above, zero for non-positive strain) was
independently re-derived and checked against the live `getTotalForce()` output for all 5 new bundles
at `lumbar_extension=30°` — **max abs error 1.39e-17 N** (floating-point exact) — the same
machine-cross-check discipline as the knee build, re-applied to this specific instance/parameter
regime rather than assumed to carry over from the class-level knee verification.

---

## 2. Data source and geometry pipeline

### 2a. Attachment points — real mesh geometry, not literature-recalled dimensions

The lumbar spine in `subject2_unified.osim` (same as every fork that touches it) is collapsed to
**one joint**, `back` (pelvis → torso), 3 DOF (`lumbar_extension`/`bending`/`rotation`) — confirmed
live: `pelvis_offset` translation = `(-0.10502, 0.08500, 0)` m, `torso_offset` translation =
`(0, 0, 0)` m, `lumbar_extension` range = ±90.0000001836389° (a permissive coordinate limit, not a
physiological ROM claim — see §6). This is the model's **own fidelity limit** (real spine has 5
lumbar vertebrae + sacrum) — the same single-hinge reduction `docs/MECHANISM_ERECTOR_SPINAE.md` and
`docs/MECHANISM_TRUNK_FLEXORS.md` already apply to this exact joint for their muscles; ligaments
inherit it here, not a new simplification.

The **same donor** already verified+used for erector spinae/trunk flexors
(`Model_Pose2Sim_muscles_flex.osim`, Beaucage-Gauvreau et al. 2019 lineage, PMID 30714401) has real
per-vertebra bodies `lumbar1`..`lumbar5` + `sacrum`, each shipping a **real bundled STL surface
mesh** (`Geometry/lumbar1.stl`, `Geometry/sacrum.stl`, ASCII format, 1272 and 1140 vertices
respectively). This build:

1. Reads the **real mesh vertices** for `lumbar1` (top of the lumbar column, borders T12 → torso-side
   reference) and `sacrum` (bottom, S1 borders L5 → pelvis-side reference).
2. Filters each to its own **superior (top-Y) 20% slice** (the endplate-level cross-section relevant
   to inter-vertebral ligament attachment) — measured live: `lumbar1` top slice X range
   `(-0.0361, +0.0084)` m (150 verts), `sacrum` top slice X range `(-0.1435, -0.0875)` m (188 verts).
3. Places each of the 5 ligaments at a generic, anatomically-ordered fraction along that **real
   measured** AP extent: ALL at the anterior extreme (vertebral-body front wall, fraction 1.00), PLL
   just posterior of it (0.55, see §5 for why not the first guess of 0.65), ligamentum flavum at
   lamina depth (0.45), interspinous along the spinous-process shaft (0.20), supraspinous at the
   posterior extreme (spinous-process tip, 0.00). **This ordering is uncontroversial gross anatomy
   (same epistemic status as "MCL is medial, LCL is lateral" in `add_knee_ligaments.py`); only the
   specific proportions are a generic first-step approximation** (HONESTY FLAG, not derived from
   segmenting the mesh into named sub-regions).
4. Converts each point into LaiArnold subject2's own scaled pelvis/torso local frames via the
   **exact same forward-kinematics + per-axis-scale pipeline** `add_erector_spinae.py` already built
   and validated for muscle points (imported and reused, not reimplemented): donor-body-local →
   donor pelvis-frame (`point_in_frame`) → pelvis-side: × measured `pelvis_scale`
   `(1.042901, 1.042901, 1.277999)`; torso-side: donor-pelvis-frame → generic-target-torso-frame →
   × measured `torso_scale` `(0.987945, 0.987945, 0.987945)` (same scale factors
   `MECHANISM_ERECTOR_SPINAE.md` measured, re-measured fresh here, matching to displayed precision).
5. **Live sanity check** (not assumed): donor pelvis-frame Y of `sacrum`=0.049835 m <
   `split_y`(back-joint threshold)=0.081500 m < `lumbar1`=0.216929 m → **straddles=True** — the two
   reference vertebrae really do straddle the model's single lumbar hinge, the geometric precondition
   for a non-degenerate (nonzero) moment arm (derived below, §2b).

**Donor cross-check**: `Model_Pose2Sim_muscles_flex.osim` was directly checked for any pre-existing
ligament content (its Force elements were enumerated by class) — **zero** `Blankevoort1991Ligament`,
`PathSpring`, or ligament-keyword-named objects exist in it (it ships `DeGrooteFregly2016Muscle`,
`SmoothSphereHalfSpaceForce` contact, and `CoordinateCouplerConstraint`s only). Both external drives
were also searched for any `.osim` file containing ligament-name keywords
(`flavum|supraspinous|interspinous|anterior longitudinal|posterior longitudinal`) across all 395
`.osim` files reachable — **zero hits**. Unlike the knee build (a real vendored `Blankevoort1991Ligament`
model existed to port directly), **no ready-made lumbar-ligament model exists anywhere in this
environment** — attachment points here are a geometric construction from real vertebra meshes, not a
ported real ligament model.

### 2b. Geometric derivation — why this placement scheme can give a correctly-signed engagement

The `back` joint's `pelvis_offset`/`torso_offset` frame origins **coincide exactly** (a pure-rotation
joint — confirmed live: no translation term couples to the 3 rotational coordinates). A point offset
purely by `(dx, 0, 0)` at the **same** small radius on both sides of such a joint gives
`length²(θ) = 2·dx²·(1−cos θ)` — **symmetric in θ → −θ**, i.e. it cannot discriminate flexion from
extension at all (a degenerate placement, derived and rejected before use, not found by trial and
error). The fix used here: give the pelvis-side point a real **inferior** offset (`−h_lower`, at the
sacrum's own real height) and the torso-side point a real **superior** offset (`+h_upper`, at
lumbar1's own real height), in addition to the anterior/posterior (`dx`) offset. With
`P = (dx, −h_lower)` and `T(θ) = R_z(θ)·(dx, h_upper)`, a first-order expansion at `θ=0` gives
`d(Length)/dθ|₀ ∝ dx·(h_lower + h_upper)`. Since `h_lower, h_upper > 0` by construction (§2a's live
straddle check), **the sign of the length response to `lumbar_extension` is exactly the sign of
`dx`**: anterior (`dx>0`) lengthens in extension, posterior (`dx<0`) lengthens in flexion. This is
why the mesh-derived vertical (Y) separation is not cosmetic — it is the geometric feature that makes
engagement direction well-defined at all.

---

## 3. Tension-only verification (task requirement — MEASURED)

Every sample in the full flexion/extension sweep (§4) — **455 (ligament × angle) checks** across all
5 ligaments × 91 angles (`lumbar_extension` from −90° to +90° in 2° steps, `lumbar_bending` =
`lumbar_rotation` = 0) — asserts: if `strain ≤ 0` then `getTotalForce() == 0.0` exactly.

**Result: 0 violations out of 455 checks → PASS.**

Slack lengths (via `setSlackLengthFromReferenceStrain` at the neutral pose) were independently
cross-checked against a hand-computed value (`slack_length = L_ref / (1 + reference_strain)`, using
the path's real geometric length `L_ref` measured live at the neutral pose): **max abs error 0.00e+00
m across all 5 bundles → PASS**. The same stale-state gotchas already documented in
`add_knee_ligaments.py` (cached strain/force not invalidated by a bare `realizePosition` after
`setSlackLengthFromReferenceStrain`, or after `finalizeConnections`/`printToXML`) were guarded against
here too (both require rebuilding the `State` via `model.initSystem()` again).

---

## 4. Passive lumbar flexion/extension sweep

`lumbar_extension` swept across its full model range `[-90°, +90°]` (positive = **extension**,
negative = **flexion** — sign convention reused, not re-derived, from
`docs/MECHANISM_ERECTOR_SPINAE.md` §6's already-published "bent-over, 30° flexed
(`lumbar_extension = -30°`)", independently corroborated there by 88/88 erector-spinae muscles
showing positive `computeMomentArm` = true extensors), 2° steps, `lumbar_bending`/`lumbar_rotation`
held at 0. Values are each ligament's own total tension (N):

| angle (deg) | ALL (N) | PLL (N) | LF (N) | ISL (N) | SSL (N) |
|---:|---:|---:|---:|---:|---:|
| −90 | 0.000 | 0.000 | 0.000 | 0.000 | 3.806 |
| −70 | 0.000 | 0.000 | 0.003 | 0.286 | 6.422 |
| −50 | 0.000 | 0.000 | 1.252 | 0.869 | 6.513 |
| −30 | 0.000 | 0.018 | 2.862 | 0.694 | 4.075 |
| −10 | 0.000 | 0.525 | 2.998 | 0.060 | 0.417 |
| +10 | 4.811 | 0.178 | 1.526 | 0.000 | 0.000 |
| +30 | 12.388 | 0.000 | 0.055 | 0.000 | 0.000 |
| +50 | 13.014 | 0.000 | 0.000 | 0.000 | 0.000 |
| +70 | 6.141 | 0.000 | 0.000 | 0.000 | 0.000 |
| +90 | 0.022 | 0.000 | 0.000 | 0.000 | 0.000 |

### Pre-registered thresholds (set BEFORE reading the sweep result)

Extension band = mean over `[+20°, +35°]`; flexion band = mean over `[−35°, −20°]` — generic, round
bands chosen well inside the model's permissive ±90° coordinate limit, **not** claimed to equal a
specific measured lumbar ROM value (see §6 honesty note).

| ligament | expected tighter in | mean force, extension band (N) | mean force, flexion band (N) | verdict |
|---|---|---:|---:|---|
| ALL | extension | **11.718** | 0.000 | **PASS** |
| PLL | flexion | 0.000 | **0.108** | **PASS** |
| LF | flexion | 0.177 | **2.955** | **PASS** |
| ISL | flexion | 0.000 | **0.591** | **PASS** |
| SSL | flexion | 0.000 | **3.450** | **PASS** |

**All 5 ligaments PASS** — ligamentum flavum, supraspinous, interspinous, and PLL all tighten in
flexion; ALL tightens in extension — exactly the pattern the task named.

---

## 5. OODA-forced fix — PLL failed twice before it passed (not a one-shot)

Per the "honest-negative is not a free pass" discipline, this section documents the real, live
correction process rather than presenting only the final passing numbers.

**Round 1 — Observe**: the first attempt placed PLL at AP fraction 0.65 (anatomically motivated: a
lumbar vertebral body occupies roughly the anterior 35-40% of the total sagittal footprint, so its
posterior wall sits near fraction 0.60-0.65) with `reference_strain_at_neutral = 0.00`. The full
sweep showed PLL strain **peaking at neutral and decaying in both directions** (a symmetric-ish
hump), and the pre-registered band check read **0.000 N in both bands → FAIL**. This was not
accepted as an honest negative on the first result.

**Round 1 — Orient**: a direct computational scan of AP fraction over `[0.00, 1.00]` (holding
`reference_strain=0` for pure shape diagnosis) found the real zero-moment-arm crossover for this
specific geometry falls **between fraction 0.60** (flexion-favoring: `strain(−30°)=−0.0204 >
strain(+30°)=−0.0226`) **and 0.65** (flips to extension-favoring: `strain(−30°)=−0.0290 <
strain(+30°)=−0.0142`) — 0.65 landed just past the crossover, on the wrong side.

**Round 1 — Decide/Act**: moved PLL to fraction 0.55, comfortably inside the flexion-favoring
region, while preserving the anatomical ordering `ALL(1.00) > PLL(0.55) > LF(0.45) > ISL(0.20) >
SSL(0.00)`.

**Round 2 — Observe**: re-running with fraction 0.55 **still gave 0.000 N / 0.000 N** in the two
bands.

**Round 2 — Orient (root cause #2)**: with `reference_strain_at_neutral` still 0.00, PLL's length is
maximized AT neutral and decreases moving to **either** side (confirmed: `strain(−30°)=−0.012`,
`strain(+30°)=−0.031`, both negative/slack). The fraction fix corrected the *direction* of the
asymmetry (flexion decays slower than extension) but strain never actually crosses zero (goes taut)
within the pre-registered band at `reference_strain=0`.

**Round 2 — Decide/Act**: directly scanned candidate `reference_strain_at_neutral` values
(0.000→0.030) at the fixed fraction-0.55 geometry, measuring the actual pre-registered-band force
means for each (not a guess):

| R (reference_strain_at_neutral) | ext band mean (N) | flex band mean (N) |
|---:|---:|---:|
| 0.000 | 0.0000 | 0.0000 |
| 0.010 | 0.0000 | 0.0218 |
| 0.012 | 0.0000 | 0.0440 |
| **0.015** | **0.0000** | **0.0981** |
| 0.020 | 0.0025 | 0.2573 |
| 0.030 | 0.0864 | 0.8784 |

`R=0.015` was chosen: a clean, comfortable margin (extension band still cleanly slack at 0.000 N,
flexion band clearly engaged), smaller than ALL's `+0.02` (consistent with PLL being generally
described as a thinner/lower-tension structure than ALL at rest). Re-verified after **both** fixes:
all 5 ligaments pass (§4 table).

This is reported as a load-bearing part of the result, not hidden: a false-positive first "it looks
plausible" pass on the wrong geometry (fraction 0.65 with no sweep check) would have shipped a PLL
that doesn't actually discriminate flexion from extension at all.

---

## 6. Honest gaps (full list)

1. **Exact stiffness/reference-strain numbers are NOT from Pintar et al. 1992's own table** — this
   was searched for, hard, and forced through multiple angles before accepting the gap: Pintar FA,
   Yoganandan N, Myers T, Elhagediab A, Sances A Jr. "Biomechanical properties of human lumbar spine
   ligaments." *J Biomech.* 1992;25(11):1351-6. **PMID 1400536**, DOI
   `10.1016/0021-9290(92)90290-h` — verified LIVE via Europe PMC (title/authors/journal/year/volume/
   issue/pages match exactly; this is the task's own named source). The paper is paywalled with **no
   indexed PubMed abstract** (pre-dates routine PubMed abstracting), Europe PMC's `fullTextXML`
   endpoint for this PMID returns **HTTP 404**, Semantic Scholar confirms
   `openAccessPdf.status: CLOSED`, and the DOI redirect chain (`doi.org` → `linkinghub.elsevier.com`)
   surfaces no visible content. **Three** different open-access secondary papers were live-fetched
   looking for a reproduced numeric table: Wiczenbach et al. 2023 (PeerJ, PMID 37583909, PMCID
   PMC10424670 — gives ligament **thickness only**, 0.5-5.0mm across 7 ligaments × 3 levels, citing
   Chazal et al. 1985 + Pintar 1992, but no stiffness/modulus table); Liu et al. 2019 (Front Bioeng
   Biotechnol, PMID 31921829, PMCID PMC6928040 — defers its own ligament property table to an
   earlier, unretrieved paper); Naserkhaki et al. 2018 (J Biomech, PMID 28549604 — full abstract
   retrieved, confirming **8 independently published ligament property datasets disagree
   substantially**, "None of the datasets yielded results in agreement with all reported
   measurements" — itself external, published corroboration that a single generic representative
   value is not an unusual simplification for a first-step model — but its own full text is
   paywalled). **Generic, clearly-flagged representative values are used instead**
   (`linear_stiffness`: ALL 350N > PLL 250N > LF 200N > SSL 100N > ISL 60N; `transition_strain`:
   0.06 default, 0.15 for LF), ordered by uncontroversial relative ligament robustness, never
   presented as extracted-from-Pintar numbers.
2. **The qualitative engagement-direction claim is robust to gap #1** — §2b's derivation shows the
   *sign* of the length response to `lumbar_extension` is governed by attachment-point placement
   (geometry), not by `linear_stiffness` magnitude; a stiffer or softer spring at the same geometric
   placement still engages/disengages at the same angle, only the force magnitude scales. The
   verified result (§4) is therefore the falsifiable claim the task actually asked for; the stiffness
   gap affects absolute force magnitude realism, not the tighten-in-flexion/extension pattern.
3. **Attachment-point proportions are a generic anatomical ordering, not a mesh segmentation** — the
   absolute scale (real measured AP extent of the lumbar1/sacrum top-slice) is real; the specific
   fractions (1.00/0.55/0.45/0.20/0.00) placing each ligament along that extent are a first-step
   approximation (uncontroversial gross anatomy, not derived from labeling individual mesh vertices
   as "lamina" vs "spinous process").
4. **Single lumped bundle per ligament, matching the model's own single-hinge fidelity** — real
   ALL/PLL/LF/ISL/SSL each span up to 6 segments (5 lumbar vertebrae + sacrum); this build represents
   each as ONE pelvis↔torso band, the same reduction `MECHANISM_ERECTOR_SPINAE.md`/
   `MECHANISM_TRUNK_FLEXORS.md` already apply to this joint for muscles. Multi-bundling (as
   `add_knee_ligaments.py` did for the knee, which had a real 42-fiber-bundle donor source) was
   considered and rejected: no real per-level ligament attachment source exists anywhere in this
   environment (§2a), so multi-bundling here would add spurious precision, not real fidelity.
5. **Torso Z-axis scale factor is assumed, not measured** — same pre-existing, disclosed limitation
   `add_erector_spinae.py` found (torso has zero native muscle points to derive a clean per-axis
   scale from); inherited here via reuse of `measure_torso_scale()`, not re-litigated.
6. **Deep-range magnitudes (beyond a generic ±35° band) should be read qualitatively, not
   quantitatively** — ALL's force is non-monotonic across the full ±90° sweep (peaks ~13N at +50°,
   falls to 0.02N at +90°); this mirrors the same single-hinge/prescribed-kinematics artifact
   `MECHANISM_KNEE_LIGAMENTS.md` §5 found for MCL at deep flexion. The model's ±90° coordinate limit
   is a permissive software range, not a physiological lumbar ROM claim (real lumbar flexion-extension
   arc is well inside that range) — no specific ROM citation is claimed here since it was not
   re-verified live this session; the pre-registered ±20-35° bands were chosen as a generic,
   comfortably-sub-physiological test range specifically to avoid relying on an unverified ROM number.
7. **capsular (CL) and intertransverse (ITL) ligaments are out of scope** — present in the standard
   7-ligament set confirmed via Liu et al. 2019, but not named in the task.
8. **Static/passive sweep only** — no muscle activation, external load, or velocity-dependent damping
   case was exercised (`damping_coefficient=0.003`, present but inactive at zero sweep velocity), same
   scope limitation as the knee build.
9. **Single donor, single donor subject** — all 5 ligaments' geometric scale derives from one donor's
   vertebra meshes; no second independent geometric source was cross-checked (the literature-property
   gap in #1 is the dimension that WAS cross-checked against multiple sources; the geometry dimension
   was not).

---

## 7. Citations (all live-verified this session, not recalled)

| citation | PMID | DOI | role |
|---|---|---|---|
| Pintar FA, Yoganandan N, Myers T, Elhagediab A, Sances A Jr. Biomechanical properties of human lumbar spine ligaments. *J Biomech.* 1992;25(11):1351-6. | 1400536 | 10.1016/0021-9290(92)90290-h | task's named source; identifies the standard ligament set (numeric table not retrievable, §6.1) |
| Christophy M, Faruk Senan NA, Lotz JC, O'Reilly OM. A musculoskeletal model for the lumbar spine. *Biomech Model Mechanobiol.* 2012;11(1-2):19-34. | 21318374 | 10.1007/s10237-011-0290-6 | task's "validated lumbar FE model" alternative; structural lineage of this repo's donor |
| Nachemson AL, Evans JH. Some mechanical properties of the third human lumbar interlaminar ligament (ligamentum flavum). *J Biomech.* 1968. | 16329292 | 10.1016/0021-9290(68)90006-7 | justifies LF's elevated `transition_strain`+positive pre-strain modeling choice |
| Naserkhaki S, Arjmand N, Shirazi-Adl A, Farahmand F, El-Rich M. Effects of eight different ligament property datasets on biomechanics of a lumbar L4-L5 finite element model. *J Biomech.* 2018;70:33-42. | 28549604 | 10.1016/j.jbiomech.2017.05.003 | corroborates that the field's own literature disagrees on ligament property values |
| Wiczenbach T, Pachocki L, Daszkiewicz K, Łuczkiewicz P, Witkowski W. Development and validation of lumbar spine finite element model. *PeerJ.* 2023. | 37583909 (PMCID PMC10424670) | 10.7717/peerj.15805 | open-access; confirmed ligament thickness range, citing Chazal1985+Pintar1992 |
| Liu T, Khalaf K, Adeeb S, El-Rich M. Numerical Investigation of Intra-abdominal Pressure Effects on Spinal Loads and Load-Sharing in Forward Flexion. *Front Bioeng Biotechnol.* 2019. | 31921829 (PMCID PMC6928040) | 10.3389/fbioe.2019.00428 | open-access; confirms standard 7-ligament naming set |
| Beaucage-Gauvreau E, Robert-Lachaine X, Brandon SCE, Larivière C. Validation of an OpenSim full-body model with detailed lumbar spine... *Comput Methods Biomech Biomed Engin.* 2019;22(5):451-64. | 30714401 | 10.1080/10255842.2018.1564819 | donor model provenance (already verified in `MECHANISM_ERECTOR_SPINAE.md`, reused for mesh geometry only) |

---

## 8. Files

- `scripts/msk/spine_ligament_config.json` — the hotswappable dial: full provenance block
  (citations, both honesty flags, scaling method, measured scale factors, donor mesh profiles) + all
  5 bundle records (attachment points in scaled pelvis/torso local coordinates, stiffness, transition
  strain, damping, reference strain). Regenerated fresh on every run (idempotent, matches
  `add_erector_spinae.py`'s own rebuild-from-scratch convention).
- `scripts/msk/add_spine_ligaments.py` — build + verification script: pre-flight structural
  re-verification of the unified model, real-mesh geometry extraction, model build, slack-length and
  force-law cross-checks, tension-only sweep, engagement-threshold checks, JSON summary. The two live
  OpenSim state-staleness gotchas already documented in `add_knee_ligaments.py` are guarded against
  identically (inline comments point to the exact lines).
- `data/msk_smoketest/spine_ligaments/model_with_spine_ligaments.osim` — the new model (unified
  model + 5 ligament bundles; source `subject2_unified.osim` untouched, verified fresh from disk:
  27 bodies, 27 joints, 232 muscles unchanged, 89 `Blankevoort1991Ligament` = 84 pre-existing knee
  bundles unchanged + 5 new, 49 markers, 3 constraints, forward-dynamics 5ms integration runs without
  exception).
- `data/msk_smoketest/spine_ligaments/spine_ligament_results.json` — machine-readable results (full
  sweep table, both cross-checks, engagement thresholds, full config).
