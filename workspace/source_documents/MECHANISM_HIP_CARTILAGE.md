# MECHANISM HIP CARTILAGE CONTACT — extending the knee COMAK/JAM leg to the hip (2026-07-21)

Extends `docs/MECHANISM_CARTILAGE_CONTACT.md` (tibiofemoral/patellofemoral cartilage contact
pressure/area via OpenSim-JAM's `Smith2018ArticularContactForce`) to the **hip** (acetabular/
femoral-head). Per `COORDINATOR.md` §5's knee closing-cell precedent, this is a candidate 5th
decorrelated leg for a future hip closing-cell — but this build's headline result is a **forced,
machine-verified HONEST NEGATIVE**: no ready-to-run hip cartilage-contact model exists, locally or
in the public opensim-jam ecosystem, as of this session. Section 1 is that search trail. Section 2
is what would be needed to close the gap. Section 3 is a **synthetic idealized toolchain-readiness
smoke test** built to de-risk that future work (NOT a hip result — labeled as such throughout).
Section 4 is the literature anchor. Section 5 is honest gaps, including the two caveats the task
explicitly asked to flag (different-individual composability; passive-vs-loaded).

**Isolation respected:** read-only against the OpenSim/JAM install (`git status` there before/after
this session shows zero tracked-file changes — only pre-existing stray `opensim.log` test-run
artifacts, untouched by this session); nothing pushed/committed; only new files added to bodytwin
(`scripts/msk/hip_cartilage_contact.py`, this doc, and the small — 33KB — synthetic-geometry output
under the already-gitignored `data/msk_smoketest/`).

---

## 1. The search — forced, not a one-shot miss (all machine/HTTP-verified this session)

Before declaring a negative, every plausible source was checked directly (raw HTTP status codes /
exhaustive repo-tree greps, not narration — full detail also embedded as `SEARCH_LOG` in
`scripts/msk/hip_cartilage_contact.py`):

| # | channel | method | result |
|---|---|---|---|
| 1 | Local disk (JAM install + broad filesystem search) | `find` for hip/acetab/femoral_head | **0 hits** |
| 2 | `opensim-jam-org` GitHub org | `GET /orgs/opensim-jam-org/repos` | 3 repos total: `opensim-core` (fork), `jam-plugin`, `jam-resources` — no hip-named repo |
| 3 | `opensim-jam-org/jam-resources` full tree | `GET .../git/trees/main?recursive=1` | **1001 paths, not truncated, 0 hip/acetab hits.** `models/` has ONLY `knee_healthy`, `knee_tka` |
| 4 | `clnsmith/opensim-jam` (lead author's actively-maintained personal continuation — same Colin R. Smith as the `Smith2018*` namesake) | same tree API | **4760 paths, not truncated, 0 hip/acetab hits.** `models/` has ONLY `lenhart2015` (knee) + `cams` (verified = **CAMS-Knee** TKA implant dataset, NOT hip/FAI-cam morphology despite the name — checked its file list: `femur_implant`/`tibia_implant`/`tibia_inlay` STLs) |
| 5 | GitHub repo-search API, 4 query variants | `hip+cartilage+STL+acetabular`, `acetabular+cartilage+mesh`, `Harris+hip+cartilage+contact`, `SCIInstitute+hip` | **0 hits**, all 4 |
| 6 | SimTK `simtk.org/projects/opensim-jam` | live HTTP GET | 200 OK — same knee-only codebase as #2-4 |
| 7 | SimTK `simtk.org/home/hip_muscles` ("Hip Musculoskeletal Model") | live HTTP GET | 200 OK — **real project**, but MUSCLE actuator geometry only (`Hip_wNotes.osim`, `Hip_LaiArnold2017_hiptorso_OS40.osim`, `Hip_23_72_Baseline_3.osim`) — 0 mentions of cartilage/acetabular/contact-mesh on the page |
| 8 | The one literature hit that DOES do hip JAM/COMAK cartilage contact (Gaffney et al. 2022, §4) — its own cited data link | live HTTP GET on `mrl.sci.utah.edu/software/hip-data/` (href confirmed byte-for-byte from the paper's own PMC HTML, not a transcription guess) | **HTTP 404**, dead this session — and even if live, was only CT images, not derived contact meshes/a redistributable `.osim` |
| 9 | Web search (the `WebSearch` tool's own budget was exhausted mid-session — confirmed via its own returned error string, not assumed; fell back to direct `curl` against `duckduckgo.com`) | keyword search | Surfaced the same 4 sources above (#2, #4, #8, plus a 2024 *J Biomech* OpenSim-JAM "personalization framework" paper, knee-lineage, not hip) — no new hip-model hit |
| 10 | OpenSim/JAM C++ source itself (`Smith2018ArticularContactForce`/`Smith2018ContactMesh` `.h`/`.cpp`) | `grep -i "knee\|tibia\|femur\|patella"` | 0 hits in actual code logic — the only matches are doxygen **example** comments ("For example, in a knee joint…"). **The engine is joint-agnostic by construction** (generic `OpenSim_DECLARE_SOCKET` to any `Smith2018ContactMesh`/`PhysicalFrame`) |

**Conclusion (C, pre-registered threshold: at least one redistributable hip acetabular+femoral-head
`Smith2018ContactMesh`-ready `.osim` or mesh set found in 10 independent channels spanning the
official org, the actively-maintained fork, general code search, and the one relevant paper's own
citations): NOT MET, 0/10.** This is a forced negative (ten decorrelated channels, not one grep),
matching the discipline that a one-shot miss is not license to declare failure.

---

## 2. What's needed to close this gap — the hip analog of "the knee needed lenhart2015"

The knee leg became possible by fetching ONE asset family: `lenhart2015.osim` + 6 cartilage/bone
STLs (Apache-2.0, `opensim-jam-org/jam-resources`). The hip needs the **same shape of asset**,
which does not yet exist publicly in redistributable form:

1. **Acetabular + femoral-head cartilage surface meshes** (triangulated, `.stl`/`.obj`/`.vtp`),
   ideally with a bone backing mesh per side for `use_variable_thickness=true` (the knee's mode) —
   generic (population-representative) or at minimum one open subject.
2. **Material properties**: elastic modulus + Poisson's ratio for hip cartilage. NOT the same as
   knee — Gaffney et al. 2022 (§4) found **E=15 MPa, ν=0.475** gives the closest elastic-foundation-
   vs-FE match for the hip, vs. the knee's own `lenhart2015` default of E=5 MPa, ν=0.45 — a real,
   cited, joint-specific difference to carry forward, not reuse blindly.
3. **A `Smith2018ArticularContactForce` wiring it into a pelvis+femur `.osim`** — mechanically
   trivial once #1-2 exist (Section 3 proves the wiring/tooling works end-to-end on synthetic
   geometry with zero source changes).
4. **Kinematics to drive it**: either a passive hip-flexion sweep (this repo's own JAM install can
   generate one directly, as demonstrated in Section 3) or gait kinematics (already have subject2
   `walking1` IK in this repo, per `docs/MECHANISM_HIP_FORCE.md`).
5. **A path to get #1**: (a) contact the Gaffney/Harris/Weiss group (Univ. Utah SCI Institute /
   WashU) directly — their paper's own data link is dead, and per the FE-idealization caveat below,
   even their published CT images would need re-segmentation, not a direct STL; or (b) build cartilage
   as an idealized surface offset from an existing generic bone model (e.g. this repo's own
   `LaiArnoldModified2017`/`gait2392`-family pelvis+femur bones already used elsewhere in this repo)
   — but Anderson et al. 2010 (PMID 20176359, §4) already show idealized (spherical/conchoid)
   geometry **underestimates peak and average contact pressure by ~50% and ~25%** respectively vs.
   subject-specific geometry — a known, documented, and now doubly-relevant bias (Section 3's own
   smoke test uses exactly this kind of idealized sphere, and inherits the same caveat).

---

## 3. What WAS built instead — a synthetic idealized toolchain-readiness smoke test

**Read this section knowing it is NOT a hip cartilage-contact result.** Its only purpose: prove the
already-installed JAM engine (same `opensim-cmd` binary the knee leg used,
`/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/opensim-build/opensim-cmd`, OpenSim
4.5.2-2026-01-19-a3c872a2b) runs a **complete ForsimTool → JointMechanicsTool contact-pressure
pipeline on novel, non-knee, hip-scaled geometry with ZERO source modifications** — i.e., that the
Section 1/2 gap is purely a **data asset** gap, the same shape of gap the knee build closed, not a
tooling gap. Script: `scripts/msk/hip_cartilage_contact.py` (self-contained: generates geometry,
writes the `.osim` + JAM settings XML, runs both tools, analyzes the result, all in one file, staged
via `--stage {mesh,model,inputs,run,analyze,all}`).

### 3.1 Geometry (derived, not decorative)

Two triangulated spherical caps (half-angle 75°, ~2832 triangles each), sharing a center at a 1-DOF
`hip_flex_synth` `PinJoint`'s rotation axis: a **femoral head** (radius 25mm) and an **acetabulum**
(radius 26mm) — a uniform 1mm nominal radial gap at perfect pole alignment (`hip_flex=0`). Constant-
thickness `Smith2018ContactMesh` (2.0mm femoral / 2.5mm acetabular — no `mesh_back_file` needed),
`elastic_modulus=15 MPa`, `poissons_ratio=0.475` — **Gaffney et al. 2022's own best-fit hip values**
(§4), not the knee's. `hip_flex_synth` is directly **prescribed** 0→90° over 0.90s (91 steps,
1°/0.01s) — there are zero free dynamical coordinates in the whole model (pelvis welded to ground),
so this is a pure kinematic sweep, not a settling-dynamics simulation.

### 3.2 A real bug, found and fixed via OODA (not swept under "honest negative")

**First attempt failed**: the anatomically "obvious" convention (femoral-head mesh normals pointing
OUTWARD, away from center, toward the socket; acetabulum normals pointing INWARD, toward center,
toward the ball) gave **zero contact area at every one of the 92 frames**, including `hip_flex=0`
where a 1mm gap vs. 4.5mm total cartilage thickness should trivially engage. Orient: read
`Smith2018ArticularContactForce.cpp::computeMeshProximity()` — the ray-cast direction used is the
**negative** of the casting mesh's own triangle normal, not the normal itself. Rather than layer a
second guess about what that implies, ran a cheap, decisive 4-way empirical test at a single static
pose (`scratchpad/hip_orient_test/probe.py`, ~1-2s per combo): only **casting (femoral head) =
INWARD-declared + target (acetabulum) = OUTWARD-declared** engaged both mesh-sides consistently
(2905mm²/3142mm², 26.8 MPa); the "obvious" combination gave 0/0 at every combination that included
it. Fixed at the source (`stage_mesh()`'s `invert_normals` arguments swapped), documented in the
script so a future reader doesn't re-derive it wrong a second time.

### 3.3 Results + machine checks (5/5 PASS — two genuine FAILs surfaced first, diagnosed via closed-form geometry, not loosened blindly)

| quantity | value |
|---|---:|
| contact area | 1014.6 – 2905.1 mm² (shrinks monotonically as `hip_flex` sweeps 0→90°) |
| mean pressure | 26.68 – 27.02 MPa |
| max (single-triangle) pressure | up to 28.35 MPa |
| contact force | up to 49.0 kN |

| check | result |
|---|---|
| engagement floor (area>1mm², meanP>0.1MPa) | PASS |
| pathology ceiling (not a units/blow-up bug, <1000 MPa) | PASS (28.35 MPa max) |
| **geometric falsifier**: contact area must shrink monotonically as `hip_flex` sweeps 0→90° (derived independently from spherical-cap-intersection geometry — two 75°-half-angle caps' pole-separation angle equals the flexion angle itself, so overlap must shrink from 0° and only fully vanish past 150°) | **PASS — Pearson r(time, area) = −1.0000**, i.e. essentially perfect monotonic decrease |
| force ≈ mean_pressure×area identity, within the curvature-cancellation ceiling | PASS (35.64% mean rel. err, inside the closed-form bound `(1+cos75°)/2` → 58.9% ceiling for a uniform-pressure 75°-cap; NOT the knee's flat-patch 10%, see §3.4) |
| casting/target mesh-pair area agreement, within the radius-ratio bound | PASS (7.54% median rel. diff, inside `(26/25)²−1+3pp` = 11.16% ceiling) |

Full detail: `data/msk_smoketest/hip_cartilage_contact/hip_cartilage_contact_results.json`.

### 3.4 The two diagnosed FAILs — geometric derivation, not threshold-shopping

Both checks above FAILED under the **knee script's own thresholds** on the first pass; both were
diagnosed via closed-form geometry from THIS synthetic shape (a deliberately curved 75°-half-angle
cap and a deliberate 25/26mm radius mismatch — neither applies to the knee's much shallower, nearly
radius-matched articular patches), not by loosening a number to make data pass (same discipline the
knee script itself models):

- **Curvature cancellation**: per-triangle pressure vectors on a curved cap point in different
  local-radial directions, so their vector sum (`|F|`) is less than the scalar product
  `mean_pressure × area` (which ignores direction). Closed form for uniform pressure over a
  half-angle-θ cap: `|F|/(meanP·A) = (1+cos θ)/2` → for θ=75°, a ceiling of 58.9% relative error.
  Measured: 35.64%, comfortably inside — and at `t=0` specifically (where the two concentric
  spheres' gap is *exactly* uniform by construction, not just approximately), the measured ratio
  (0.62924) matches the closed-form uniform-pressure prediction (0.62941) to **0.03%** — a precise,
  falsifiable, unplanned confirmation that both the mesh generator and the tool's own contact
  computation are behaving exactly as the geometry predicts at that frame.
- **Radius-ratio mesh-pair mismatch**: the two meshes share the same angular grid but different
  radii, so the same solid-angle patch maps to different physical areas: `(26/25)²−1 = 8.16%`
  predicted vs. **7.54% measured** — within 0.6 percentage points.

### 3.5 What this smoke test does NOT prove

It says nothing about real hip cartilage contact magnitudes, locations, or gait behavior — no
muscles, no ground-reaction force, no subject-specific geometry, fixed (not load-dependent)
interference. The 26.7-27.0 MPa mean pressure is **same order of magnitude** as the literature (§4:
5.1-15.65 MPa range) but is not compared as a validation — a categorically different loading regime
(idealized fixed interference vs. walking/stair-climbing) makes a tighter comparison meaningless, the
same "one-sided bound, not equality" discipline the knee's passive-flexion leg already used.

---

## 4. Literature anchor (all PMID/DOI live-verified via NCBI eutils + PMC full-text HTML this session)

- **Anderson AE, Ellis BJ, Maas SA, Peters CL, Weiss JA.** "Validation of finite element predictions
  of cartilage contact pressure in the human hip joint." *J Biomech Eng.* 2008;130(5):051008.
  doi:[10.1115/1.2953472](https://doi.org/10.1115/1.2953472). **PMID 19045515.** Cadaveric
  pressure-**film** experiment (not just simulation) + subject-specific FE, walking/stairs/stair-
  climbing: experimental peak **10.0 MPa** (film detection limit), average **4.4-5.0 MPa**, area
  **321.9-425.1 mm²**; FE-predicted peak **10.8-12.7 MPa**, average **5.1-6.2 MPa**, area
  **304.2-366.1 mm²**.
- **Anderson AE, Ellis BJ, Maas SA, Weiss JA.** "Effects of idealized joint geometry on finite
  element predictions of cartilage contact stresses in the hip." *J Biomech.* 2010;43(7):1351-7.
  doi:[10.1016/j.jbiomech.2010.01.010](https://doi.org/10.1016/j.jbiomech.2010.01.010). **PMID
  20176359.** Directly relevant caveat for Section 3's own idealized-sphere approach: spherical/
  conchoid bone geometry with smoothed cartilage **underestimated peak/average contact pressure by
  ~50%/~25%** vs. the subject-specific model.
- **Harris MD, Anderson AE, Henak CR, Ellis BJ, Peters CL, Weiss JA.** "Finite element prediction of
  cartilage contact stresses in normal human hips." *J Orthop Res.* 2012;30(7):1133-9.
  doi:[10.1002/jor.22040](https://doi.org/10.1002/jor.22040). **PMID 22213112.** 10 normal
  volunteers, subject-specific FE: peak stress **7.52±2.11 MPa** (walking heel-strike, 233% BW) to
  **8.66±3.01 MPa** (descending-stairs heel-strike, 261% BW).
- **Gaffney BMM, Williams ST, Todd JN, Weiss JA, Harris MD.** "A Musculoskeletal Model for
  Estimating Hip Contact Pressure During Walking." *Ann Biomed Eng.* 2022;50(12):1954-1963.
  doi:[10.1007/s10439-022-03016-w](https://doi.org/10.1007/s10439-022-03016-w). **PMID 35864367**,
  PMCID PMC9797423. The closest methodological analog found: **COMAK kinematics + an elastic-
  foundation (Smith-family ray-cast+OBB) contact model, IN OpenSim**, on subject-specific
  acetabular+femoral cartilage meshes, 5 subjects, walking, explicitly citing
  `simtk.org/projects/opensim-jam` as the tool. Table 3 (E=15 MPa, ν=0.475 — our exact material
  choice above): EF peak pressure **5.23-15.65 MPa** (subject range), contact area **23.3-39.9%** of
  acetabular surface, cross-validated against subject-matched FEA (**11.98-15.58 MPa** peak) to a
  mean **39.96±24.64%** relative difference — a useful outside calibration that this field's own
  EF-vs-FEA agreement is already looser than the knee build's own <5% mesh-pair agreement.

---

## 5. Honest gaps (including the two caveats flagged in the task)

1. **No real hip cartilage-contact result exists in this repo** — Section 1's negative is the
   headline finding, not a side note. Section 3 is a toolchain de-risking exercise, not a
   substitute.
2. **★ Different-individual composability caveat (explicitly flagged per task).** Any future hip
   cert assembled from today's available pieces would necessarily mix sources from DIFFERENT
   individuals/labs: hip cartilage geometry (if obtained at all, e.g. from Utah SCI/Weiss-lab
   subjects) + this repo's own gait kinematics (subject2, a different, unrelated individual) + a
   generic muscle geometry model (SimTK `hip_muscles`, yet another source lineage) + material
   constants fit to Gaffney et al.'s 5 subjects. Each individual pairing is a real, documented
   modeling convention (the knee build already composes `lenhart2015`'s cartilage with this repo's
   own gait data across the SAME composability seam) — but stacking 3-4 independent-individual
   sources compounds anatomical-mismatch uncertainty (bone size/cartilage-thickness/joint-
   congruence scale differently per person) in a way a single-donor model does not. Any future hip
   cert must state this compounding explicitly, not silently inherit the knee leg's single-swap
   framing.
3. **★ Passive-vs-loaded caveat (explicitly flagged per task).** Section 3's smoke test is neither
   "passive" (no ligaments, no tonic muscle activation, no gravity-driven settling — unlike the
   knee's passive-flexion leg) nor "loaded" (no ground-reaction force, no muscle-driven joint
   compression) — it is a **third, distinct regime**: pure prescribed-kinematic geometric
   interference with a fixed radius mismatch. A real future hip leg should be explicit about which
   of these three regimes (passive/kinematic-interference/loaded-gait) any reported number comes
   from, exactly as `docs/MECHANISM_CARTILAGE_CONTACT.md` §0/§6 already had to do for the knee.
4. **Idealized geometry is a known-biased shortcut, not a free equivalent** (§2 item 5, §4's
   Anderson 2010 citation) — even if a "quick" hip model were built from idealized spherical
   surfaces rather than real segmented cartilage, the literature's own FE community has already
   shown this underestimates peak/average pressure by ~50%/~25%. Section 3's numbers should not be
   read as "what a simplified hip model would show" for this reason — they are a pure toolchain
   check, deliberately not offered as a cheap substitute for real geometry.
5. **Synthetic mesh material/thickness values are placeholders**, not fit to any specific
   individual or even a population average of the SHAPE parameters (only the elastic
   modulus/Poisson's ratio are literature-sourced, from Gaffney et al.'s fitted values — the radii,
   half-angle, and thickness are chosen for a clean, falsifiable geometric test, not anatomical
   fidelity).
6. **Single-frame vs. swept-frame diagnostic asymmetry**: the orientation bug (§3.2) was found and
   fixed using a single static pose; the full 91-frame sweep was only re-run after the fix, so there
   is no direct record of what the WRONG orientation's area/pressure would have looked like across
   the full sweep (not needed for the fix, but noted for completeness).
7. **No regional (medial/lateral or anterior/posterior) contact-patch split** — same class of gap
   already flagged for the knee's tibiofemoral compartment in `docs/MECHANISM_CARTILAGE_CONTACT.md`
   §7.3; `Smith2018ArticularContactForce`'s `*_regional_*` outputs were not exercised here either.

---

## 6. Files

- `scripts/msk/hip_cartilage_contact.py` — self-contained: `SEARCH_LOG` (Section 1's machine-
  checkable record), the literature `ANCHOR` dict (Section 4), synthetic mesh generation +
  self-check (Section 3.1-3.2), `.osim`/JAM-settings XML generation, pipeline runner, and the
  geometrically-derived analysis/checks (Section 3.3-3.4). Staged via `--stage {mesh,model,inputs,
  run,analyze,all}` (default `all`).
- `data/msk_smoketest/hip_cartilage_contact/` (gitignored, like all `data/msk_smoketest/`) —
  `model/synthetic_hip.osim` + `model/Geometry/*.stl` (the two synthetic caps, ~230KB each,
  wholly-generated-by-this-script, not third-party), `inputs/` (JAM settings XML + prescribed
  coordinates), `results/` (ForsimTool states + JointMechanicsTool `.h5`, 92 frames), and
  `hip_cartilage_contact_results.json` (full machine-readable evidence: search log, synth params,
  mesh self-checks, run status, analysis + all 5 checks).
- `/tmp/.../scratchpad/hip_orient_test/probe.py` — the throwaway 4-way orientation diagnostic
  (Section 3.2); not part of the deliverable, kept only for this session's own reproducibility.
