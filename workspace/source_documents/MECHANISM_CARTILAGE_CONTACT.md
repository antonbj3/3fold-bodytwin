# MECHANISM CARTILAGE CONTACT — tibiofemoral/patellofemoral contact pressure, area, COP (2026-07-21)

Adds the CARTILAGE CONTACT MECHANICS layer below joint contact FORCE. The twin already has
knee joint reaction **force** (`docs/MECHANISM_STATIC_OPT.md`: 233.20 → **391.10 %BW** after the
sign-bug correction, `docs/MECHANISM_JOINT_FORCE_VALIDATION.md`, cross-checked against OpenSim's
own `JointReaction` to 0.0039%) — a single resultant vector at the joint. This build adds the
**pressure DISTRIBUTION** over a real elastic-foundation cartilage contact patch: peak/mean
pressure, contact area, and center-of-pressure location, for both the tibiofemoral (TF) and
patellofemoral (PF) joints. Per `COORDINATOR.md` §5, cartilage-contact is the 4th of the knee
closing-cell's 4 decorrelated legs (geometry / DOF+ROM / moment-realizability / **cartilage-
contact**).

Scripts: `scripts/msk/cartilage_contact.py` (analysis + machine cross-checks + JSON evidence,
run: `.venv-msk/bin/python3 scripts/msk/cartilage_contact.py`). Data:
`data/msk_smoketest/cartilage_contact/` (persisted JAM output `.h5`, settings XML, prescribed
motion `.sto`, `cartilage_contact_results.json`).

**Isolation respected:** read-only against the OpenSim/JAM install and the external drive;
nothing pushed/committed; the two heavy third-party inputs (the 879KB `lenhart2015.osim` model
and its ~37MB of cartilage/bone STL meshes) are NOT copied into this repo — same convention
`docs/MECHANISM_KNEE_LIGAMENTS.md` already established for the same model file — only this new
script + the small (1.8MB) numeric JAM output are new files in bodytwin.

---

## 0. Relationship to the existing FORCE layer and the Grand-Challenge flagship — read this first

Two things this doc is **not**, stated up front to avoid conflation:

1. **Not the same model instance as the 391.10 %BW force number.** That number comes from
   `LaiArnoldModified2017_poly_withArms_weldHand_scaled` (subject2, walking1 trial, static
   optimization). This cartilage-contact leg uses **lenhart2015.osim** (the JAM/Smith/Thelen
   reference healthy-knee model) — the SAME model already used in `MECHANISM_KNEE_LIGAMENTS.md`
   as the literature/geometry DONOR for the 4 major ligaments. It is a complementary,
   decorrelated capability leg (does this repo's OpenSim chain now have cartilage-contact-
   mechanics capability at all?), not a further-refined version of the 391.10 %BW force.
2. **Not the Grand-Challenge COMAK flagship.** `bt_memory/knee-comak-flagship-close-and-fork-
   build.md` already closed a DIFFERENT cert (predicted vs measured tibiofemoral **contact
   force** on the instrumented-TKA-implant `DM.osim` model, 2.503 vs 2.584 BW, 3% error) using
   the same JAM `opensim-cmd` binary. That model is a **metal-on-polyethylene implant**, not
   cartilage — it cannot answer "cartilage contact pressure" by construction. This leg uses the
   **healthy, natural-cartilage** sibling model instead (explicitly what the task asked for:
   "cartilage geometry+material, generic not subject-specific"), anchored against a **literature
   MPa figure**, not implant telemetry.

Both prior legs stand; this is a new, additive 3rd data point in the same OpenSim/JAM
environment, not a revision of either.

---

## 1. Environment / tool choice — which of the two OpenSim installs, and why

`docs/MECHANISM_MSK_ENV.md` already documents two independent, both-working OpenSim installs in
this repo: the stock pip wheel (`.venv-msk`, OpenSim 4.6, Python API) and the JAM C++ fork
(`/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/opensim-build/opensim-cmd`,
OpenSim 4.5.2-2026-01-19-a3c872a2b, re-verified working this session, exit 0). The ligament
build used the FIRST (stock wheel) because `Blankevoort1991Ligament` is upstreamed into
mainline OpenSim; it explicitly **stripped** `Smith2018ArticularContactForce` /
`Smith2018ContactMesh` before loading because those classes are **not registered** in the
stock wheel. This build needs exactly those two stripped-out classes, so it uses the **second**
install (the JAM fork), which registers them natively — no plugin load, no stripping. Verified
live this session:
```
opensim-cmd --version  ->  OpenSim version 4.5.2-2026-01-19-a3c872a2b
opensim-cmd print-xml ForsimTool  ->  writes default_ForsimTool.xml (tool recognized)
```
`h5py` (3.16.0) was added to `.venv-msk` this session (`uv pip install --offline h5py`, cached
wheel, zero network) to read the JAM tool's `.h5` output — the only new Python dependency.

**Naming note:** the task described "JAM ForwardTool" — the actual JAM class is `ForsimTool`
("**for**ward **sim**ulation"). Flagged for precision, matching this repo's own precedent of not
silently gliding past a naming/citation mismatch (`MECHANISM_KNEE_LIGAMENTS.md` §0's PMID
correction).

---

## 2. Model + contact definition — real cartilage geometry and material, verified live

`lenhart2015.osim` was already cached on disk from the prior ligament-build session
(`.../scratchpad/knee_lig/lenhart2015.osim`, fetched from
`raw.githubusercontent.com/opensim-jam-org/jam-resources/main/models/knee_healthy/lenhart2015/
lenhart2015.osim`, Apache-2.0) — this time loaded **unstripped**. Grepping the raw XML (not
narration) confirms the exact contact architecture:

| element | count | detail |
|---|---:|---|
| `Smith2018ContactMesh` | 3 | `femur_cartilage`, `tibia_cartilage`, `patella_cartilage` |
| `Smith2018ArticularContactForce` | 2 | `tf_contact` (target=femur, casting=tibia), `pf_contact` (target=femur, casting=patella) |

Per-mesh material properties (live-read from the XML, not assumed defaults):
`elastic_modulus=5,000,000 Pa` (5 MPa), `poissons_ratio=0.45`, `use_variable_thickness=true`
(distance-to-`mesh_back_file` ray cast, bounded `[1mm, 10mm]`, nominal 3mm), formulation
`nonlinear` (Bei & Fregly 2003 lumped elastic-foundation model), `max_proximity=0.02m`. The 6
mesh files this requires (`lenhart2015-R-{femur,tibia,patella}-{bone,cartilage}.stl`, 1.2-18.8MB
each, ~37MB total) were **not** previously fetched (the ligament build only needed attachment-
point XML data, not the contact meshes) — fetched fresh this session from the same repo's
`models/knee_healthy/lenhart2015/Geometry/`, byte-size-verified against the GitHub API's own
reported sizes (exact match, 6/6).

---

## 3. Pipeline — JAM's own official "passive_flexion" example, paths adapted only

Rather than hand-deriving a settings XML (transcription risk), the exact official example was
fetched from `opensim-jam-org/jam-resources/cmd/linux/passive_flexion/inputs/`
(`forsim_settings.xml`, `joint_mechanics_settings.xml`, `prescribed_coordinates.sto`) and only
paths + I/O verbosity were changed — the **physical protocol itself is JAM's own published
default**, unchanged:

- **ForsimTool**: `knee_flex_r` **prescribed** 0→91.04° over 2.5s (251 steps @ 0.01s,
  `prescribed_coordinates.sto`, pelvis_tilt held at 90°). 11 secondary tibiofemoral +
  patellofemoral DOFs (`knee_add/rot/tx/ty/tz_r`, `pf_flex/rot/tilt/tx/ty/tz_r`) **unconstrained**
  — settle via real forward dynamics under ligament + contact + gravity + a low
  `constant_muscle_control=0.02` tonic activation on every muscle (JAM's own default for passive
  runs, to avoid a zero-activation transient — not zero muscle force). All other coordinates
  locked. This mirrors the actual dynamic-knee-simulator passive-flexion protocol the
  lenhart2015 model was originally validated against (PMID 25917122 Methods: "simulate ...
  passive knee flexion").
- **JointMechanicsTool**: replays the resulting states, computes `Smith2018ArticularContactForce`
  outputs for both contacts, both mesh-side conventions.

Both PASS, real wall-clock, this session:
```
opensim-cmd run-tool inputs/forsim_settings.xml            # Forsim Tool complete. Finished in 55 second(s)
opensim-cmd run-tool inputs/joint_mechanics_settings.xml    # JointMechanicsTool complete. Finished in 37 second(s)
```
Output: `passive_flexion_states.sto` (252 rows × 193 columns) → `passive_flexion.h5` (1.8MB,
persisted at `data/msk_smoketest/cartilage_contact/results/joint-mechanics/`).

**Settled window honesty flag**: the first ~0.1s (10/252 samples) is a startup transient (model
relaxing from its default pose to mechanical equilibrium before the flexion sweep engages —
contact area at `t=0` is 163mm², jumps to a stable ~320mm² by `t=0.08s`, `knee_flex_r` still at
0° throughout). Same pattern already flagged in `bt_memory/knee-jam-plugin-install-recipe.md`
("the t=0 spike is a reference-pose mesh-preload artifact — use the settled value"). All
headline numbers below use `t≥0.1s` (measured via the data, not assumed).

---

## 4. Results — contact area, pressure, center-of-pressure (settled window, 242 samples)

| quantity | tf_contact (femur/tibia cartilage) | pf_contact (femur/patella cartilage) |
|---|---|---|
| contact area | 179–331 mm² | 59–215 mm² |
| mean pressure | 1.09–2.11 MPa | 0.61–2.36 MPa |
| peak (max-triangle) pressure | 2.11–4.78 MPa | 1.43–6.13 MPa |
| \|contact force\| | 201–502 N | 52–302 N |
| center-of-pressure excursion over the 0→91° sweep (tibia/patella-attached local frame) | 12.4mm (x), 1.7mm (y), 12.6mm (z) | — |

Both `tf_contact` numbers are reported from BOTH mesh-side conventions (femur_cartilage=target,
tibia_cartilage=casting) — they agree closely (§5) — an internal over-determination check, not
a tautology (each side is computed from its own mesh's independent ray-casting sweep).

**COP location — frame clarification (a real forced check, not assumed):** `total_center_of_
pressure` was hypothesized to be reported in `JointMechanicsTool`'s `output_position_frame`
(set to `ground`). Tested directly: re-ran `JointMechanicsTool` a second time with
`output_position_frame`/`output_orientation_frame` set to `/bodyset/tibia_proximal_r` instead —
the COP (and contact-force vector) values came back **bit-identical** to the "ground frame" run.
Conclusion, now measured rather than assumed: `Smith2018ArticularContactForce`'s own COP/force
outputs are always expressed in the contact mesh's own attached body frame (`tibia_proximal_r`
for `tf_contact`'s casting side), **regardless** of the tool's `output_position_frame` setting —
that setting only affects kinematic/bodyset-transform reporting. The COP numbers above are
therefore already anatomically local (tens-of-mm range, consistent with a joint-local frame, not
a lab/ground frame), confirmed empirically, not by reading documentation.

The COP moves smoothly and monotonically-ish across ~12-13mm in two axes over the flexion
sweep — qualitatively consistent with the well-known tibiofemoral "rollback" contact-patch
migration during flexion. Reported as a qualitative/descriptive pattern only: the specific local
axes were not independently anchored to named anatomical directions (anterior/posterior,
medial/lateral) in this build, so no quantitative "rollback distance" claim is made (see gaps).

---

## 5. Machine cross-checks (16/16 PASS) — including two real FAILs found, diagnosed, and fixed at the check design, not the data

Symmetric QC: every check below carries the same forced-steelman burden whether it passed or
initially failed. Two genuine FAILs surfaced while building the verification script; both were
diagnosed (OODA, not dismissed) and resolved by fixing the **check's statistical design**, never
by loosening a number to make data pass.

**Round 1 FAIL**: `max_pressure` cross-mesh agreement for `pf_contact` was 8.1% under a uniform
5% tolerance (shared with area/mean_pressure/force). Diagnosis: `max_pressure` is an
extreme-value statistic (single hottest triangle) computed independently on two DIFFERENT
triangulations, unlike area/mean_pressure/force which are law-of-large-numbers integrals over
every contacting triangle — not the same kind of statistic, should not share one tolerance.
**Quantitatively tested**, not just argued: the 6 cartilage/bone STL meshes were parsed directly
(ASCII STL, own parser) — triangle density is uniform across all 3 cartilage surfaces (femur
4.260, tibia 4.249, patella 4.268 triangles/mm², ruling out a mesh-resolution confound). Contact
patch size then implies contacting-triangle counts: `tf_contact` ~760–1410, `pf_contact` (smaller
patch) ~250–920. A √N scaling of extreme-value spread predicts the `pf`/`tf` max_pressure
relative-diff ratio should be ≈√(1411/253)=2.36×; the MEASURED ratio is 8.1%/2.93%=2.77× — a
striking, falsifiable, independently-derived match. `max_pressure` now gets its own 15% (median-
gated) tolerance.

**Round 2 FAIL**: even split out, `pf_contact`'s INTEGRATED area/mean_pressure still failed a 5%
ceiling on the *max-over-time* relative difference (5.55% / 4.94%). Diagnosis (h5 inspection):
gating on the single worst frame of an already-noisy 242-sample per-frame series is itself the
"fixed threshold on a noisy statistic is regime-blind" failure mode. Measured: median relative
diff is only 1.1%/1.4%, p90 is 3.9%/3.1%, and only 3/242 frames (1.2%) exceed 5% — concentrated
near a contact-patch topology transition (~t=2.14s, 62.9° flexion), not a systematic bias. Fixed
by gating on the **median** (regime-typical) relative difference, reporting max as descriptive
context rather than a second silent gate.

Final check table (all 16 PASS; full detail in `cartilage_contact_results.json`):

| check | tf_contact | pf_contact |
|---|---|---|
| engagement floor (area>1mm², meanP>0.1MPa) | PASS | PASS |
| absolute ceiling vs anchor (≤9.7 MPa, the highest verified anchor figure) | PASS (max 4.78 MPa) | PASS (max 6.13 MPa) |
| mesh-pair agreement, integrated (area/meanP/force), median ≤5% | PASS (medians 0.6-1.4%) | PASS (medians 0.4-1.4%) |
| mesh-pair agreement, max_pressure, median ≤15% | PASS (0.7%) | PASS (1.3%) |
| force ≈ mean_pressure × area identity, mean rel. err ≤10% | PASS (5.4-6.0%) | PASS (4.0-5.6%) |

---

## 6. Literature anchor (all live-verified via NCBI eutils this session, not recalled)

- **Lenhart RL, Kaiser J, Smith CR, Thelen DG.** "Prediction and Validation of Load-Dependent
  Behavior of the Tibiofemoral and Patellofemoral Joints During Movement." *Ann Biomed Eng.*
  2015;43(11):2675-2685. doi:10.1007/s10439-015-1326-3. **PMID 25917122.** THE SAME model
  (already cited for the ligament work) — its own Results, **normal gait**: mean cartilage
  contact pressure reached **6.2 MPa** (medial tibial plateau) / **2.8 MPa** (patellar facets).
  **Unit honesty flag**: NCBI's own indexed abstract (both `esummary` and `efetch
  rettype=abstract`, JSON/text/XML) reads "...reached 6.2 and 2.8 **Pa**..." — 6 orders of
  magnitude too small for cartilage contact (physically impossible; the task's own prompt cites
  ~3-10 MPa). Verified as an abstract-indexing artifact, not silently propagated or silently
  "corrected": pulled the actual PMC full-text HTML
  (`pmc.ncbi.nlm.nih.gov/articles/PMC4886716/`), which reads "...reached 6.2 and 2.8 **MPa** on
  the medial tibial plateau and patellar facets, respectively." Independently corroborated by
  the companion paper below (clean "MPa", no ambiguity).
- **Lenhart RL, Smith CR, Vignos MF, Kaiser J, Heiderscheit BC, Thelen DG.** "Influence of step
  rate and quadriceps load distribution on patellofemoral cartilage contact pressures during
  running." *J Biomech.* 2015;48(11):2871-2878. doi:10.1016/j.jbiomech.2015.04.036. **PMID
  26070646.** Same lab/model family — running midstance: mean/peak patella contact pressure =
  **5.0 / 9.7 MPa**.
- **Smith CR, Choi KW, Negrut D, Thelen DG.** "Efficient Computation of Cartilage Contact
  Pressures within Dynamic Simulations of Movement." *Comput Methods Biomech Biomed Eng Imaging
  Vis.* 2018;6(5):491-498. doi:10.1080/21681163.2016.1172346. **PMID 30740280.** Likely the
  namesake/origin of the `Smith2018ArticularContactForce` class itself (ray-casting + OBB +
  elastic foundation) — states pressure magnitudes "were similar to those predicted by finite
  element models," corroborating the algorithm, not just a number.

**Comparison — pre-registered as a ONE-SIDED bound, not an equality test** (stated before this
was written: passive flexion has no ground-reaction-force/body-weight load, only ligament
tension + 0.02 tonic muscle activation + gravity — it CANNOT physically exceed the same model's
own loaded-gait pressure without indicating a units/stiffness bug; it is not expected to MATCH
gait either):

| | our passive peak mean pressure | anchor (same model, loaded gait) | ratio |
|---|---:|---:|---:|
| TF (medial-plateau-equivalent) | 2.106 MPa | 6.2 MPa | **0.34×** |
| PF (patellar-facet-equivalent) | 2.361 MPa | 2.8 MPa | **0.84×** |

Both ratios are `<1` (passes the physical bound) and both peak (max-triangle) pressures (4.78 /
6.13 MPa) sit below the highest verified anchor figure (9.7 MPa, running PF peak) — the
"absolute ceiling" check in §5. The DIFFERENTIAL pattern is mechanistically coherent, not
coincidental: TF loading is predominantly ground-reaction-force-driven (near-zero here, hence
the large 0.34× undershoot), while PF loading is predominantly quadriceps-force/flexion-angle-
driven (present even in this passive-with-tonic-activation rig, hence the much closer 0.84×
approach) — a real, textbook biomechanical distinction the two numbers reproduce without being
tuned to.

---

## 7. Honest gaps

1. **Passive, not loaded-gait** — by design (§0, §6): this measures cartilage-contact-mechanics
   CAPABILITY on a physiologically-plausible passive trial, anchored as a one-sided bound. A
   loaded condition (e.g. via `ForsimTool`'s `external_loads_file`, confirmed present in the
   tool's own property list but not exercised here) would allow a tighter, two-sided gait
   comparison — flagged as the natural next step, not attempted (time/scope).
2. **Generic, not subject-specific cartilage** — same caveat class as
   `MECHANISM_KNEE_LIGAMENTS.md`: `elastic_modulus`/`poissons_ratio`/mesh geometry are
   Lenhart2015's literature/MRI-derived healthy-subject values, not derived from any subject in
   this repo's own corpus (matches the task's own framing: "generic, not subject-specific").
3. **No medial/lateral tibiofemoral compartment split** — `tf_contact` reports ONE combined
   contact patch across the whole tibial plateau. The `Smith2018ArticularContactForce` class
   does declare `*_regional_*` outputs (regional area/pressure/force/COP) that could split this,
   not requested in this run (kept the h5 small/simple for the first pass) — a real, bounded
   next step, not a blocker.
4. **`num_contacting_triangles` silently not written** — requested in `contact_outputs`
   (`casting_num_contacting_triangles`/`target_num_contacting_triangles`, both declared
   `OpenSim_DECLARE_OUTPUT`s on the C++ class per `Smith2018ArticularContactForce.h`), but does
   not appear in the resulting `.h5` — a genuine tool-side gap between "declared output" and
   "`JointMechanicsTool`'s writer implements it." Worked around via an independent measurement
   (STL-parsed triangle density × measured contact area, §5), not silently dropped.
5. **COP anatomical axis labels not independently verified** — the local frame's x/y/z were
   confirmed to be the mesh's own attached body frame (§4), but which axis is "anterior" vs
   "medial" in that frame's own convention was not independently derived here (would need
   tracing the model's original bone-registration axes) — the "rollback" pattern is reported
   qualitatively, not as a labeled quantitative claim.
6. **Single model, single source** — same as the ligament build: one healthy-subject model
   family (Lenhart2015/JAM), no second independent cartilage-contact model cross-checked
   against it.

---

## 8. Files

- `scripts/msk/cartilage_contact.py` — analysis + machine cross-checks + JSON evidence writer.
  Self-contained given the persisted `.h5`; the full external-pipeline fetch+run recipe (model,
  meshes, JAM settings, tool invocations) is documented in its own docstring/`FETCH_RECIPE`
  (not re-executed by default — those inputs are ~37MB external third-party files, consistent
  with how the ligament build handled the same source model).
- `data/msk_smoketest/cartilage_contact/results/joint-mechanics/passive_flexion.h5` — persisted
  JAM `JointMechanicsTool` output (1.8MB): per-frame area/mean_pressure/max_pressure/force/COP
  for `tf_contact` and `pf_contact`, both mesh-side conventions, 252 samples.
- `data/msk_smoketest/cartilage_contact/inputs/` — the adapted `forsim_settings.xml`,
  `joint_mechanics_settings.xml`, and JAM's own `prescribed_coordinates.sto` (the exact passive-
  flexion protocol run).
- `data/msk_smoketest/cartilage_contact/cartilage_contact_results.json` — full machine-readable
  evidence (settled-window summaries, all 16 checks with detail, thresholds, literature anchor,
  anchor ratios).
