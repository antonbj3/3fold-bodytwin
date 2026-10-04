# MECHANISM FOOT MULTISEGMENT — dial-turn 2: midtarsal (Chopart) joint (2026-07-21)

Operator's explicit priority: "all joints incl. feet" (all joints incl. feet). Executes the
task's 3-part ask on top of the already-verified core (`docs/MECHANISM_FOOT_FIDELITY_PLAN.md`,
`scripts/msk/foot_detail_proto.py`, dial-turn 1 = the hallux split): (1) re-survey + extend the
multi-segment-foot model-candidate search to local disk, (2) assess graftability of a midtarsal
joint specifically, (3) build it. Script: `scripts/msk/foot_multisegment.py` (1000+ lines,
re-runnable, self-contained — one command produces the model, every check, and the held-out
validation). New model: `data/msk_models/LaiArnold_midtarsal_proto_scaled.osim`. Original source
model verified bit-identical before/after (mtime `1638784176`, md5 `0eb06c2a2b235b579b76855c8c8bff88`
— unchanged across every experiment this session).

---

## 1. Model candidate survey — re-verified + extended to local disk

`docs/MECHANISM_FOOT_FIDELITY_PLAN.md` §1 already surveyed the literature live (KU Leuven
Malaquias/Postolka lineage, U Queensland Maharaj "JC" model, Bruening/Saraswat lineage, Oxford
Foot Model) and concluded: **no ready-to-download, confirmed-free, OpenSim-native, more-than-3-
segment foot-ankle model exists.** This session:

- **Independently re-verified 2 of that survey's citations live** (not recalled): PMID 27381808
  (Malaquias TM, Silveira C, Aerts W, De Groote F, Dereymaeker G, Vander Sloten J, Jonkers I.
  "Extended foot-ankle musculoskeletal models for application in movement analysis." *Comput
  Methods Biomech Biomed Engin* 2017;20(2):153-159) — fetched via NCBI eutils `esummary`
  (title/authors/journal/date match exactly) and its PubMed abstract directly: "two detailed 3D
  multibody foot-ankle models generated based on CT scans using a semi-automatic tool... **five
  rigid segments (talus, calcaneus, midfoot, forefoot and toes), connected by five joints (ankle,
  subtalar, midtarsal, tarsometatarsal and metatarsophalangeal)**, one with 15DOF and the other
  with 8DOF." No axis-orientation or segment-mass numbers are in the abstract; full text is
  paywalled (Taylor & Francis, confirmed via direct WebFetch this session — matches the prior
  session's own finding, not contradicted). PMID 38723393 (Postolka et al. 2024 follow-up)
  abstract also re-fetched: covers talocrural+subtalar contact mechanics only (4D-CT validated,
  12 healthy + 4 calcaneal-fracture patients), no midtarsal-specific numbers.
- **WebSearch is quota-exhausted this session** — hit the identical "2000 of 2000" wall the prior
  session already documented. Worked around with WebFetch (separate quota) for the two PMID
  checks above, same method the prior session used.
- **NEW this session — exhaustive local-disk search** (the prior survey was web-only):
  - `grep -rliE "malaquias|bruening|saraswat|midtarsal|chopart|multi-?segment.*foot|footankle"` across
    `~/projects/cad-to-simulation-I` (*.py/*.md/*.json) — **zero hits**. `find` for any `*.osim` in
    that tree — **zero hits**.
  - `find` for every `*.osim` across all mounted drives (`/media/anton/8838D60F38D5FBDE`,
    `/media/anton/183E48713E484A48`, `/mnt/shared_data`) — **~70 files found and inventoried**.
    Every one is either the LaiArnold/Rajagopal lineage (3-segment foot, the model this script
    extends — confirmed present in 10 subjects under `LabValidation_withVideos/`), the
    gait2392/gait2354/gait10dof18musc family, or generic OpenSim example/test models (`arm26`,
    `pendulum`, `ThreeMasses`, `WalkerModel`, ...). **None has a midtarsal/tarsometatarsal joint
    or more than 3 segments distal to the tibia.**
  - Checked the OpenSim JAM C++ build tree (`opensim_jam_build/`, both `opensim-core-jam` source
    and every `run_*` results directory) and the pip-installed `opensim` wheel's own bundled
    tests/examples — same result.

**Bottom line (extends, does not contradict, the prior survey): there is no multi-segment foot
model sitting unused anywhere in this project's reachable storage.** The KU Leuven / U Queensland
models remain gated behind either a from-scratch CT-scan rebuild or a human-identity contact
request to the respective labs (Jonkers, Lichtwark/Maharaj) — unchanged from the prior finding,
now confirmed not to be a local-search miss.

---

## 2. Graftability assessment

**DOF count.** This dial-turn adds exactly 1 DOF/side (`midtarsal_angle_r`/`_l`), taking the foot
from 3 joints (ankle, subtalar, mtp) to 4 (ankle, subtalar, **midtarsal**, mtp) — one joint short
of Malaquias' own reduced 8-DOF variant (which further splits tarsometatarsal from midtarsal and
likely gives ankle+subtalar more than 1 DOF each), far short of the 15-DOF variant. Deliberate:
one increment, per the first-step doctrine and the Fidelity Plan's own staged recommendation.

**Marker/experimental requirements — the harder graft than dial-turn 1, measured directly.**
`calcn_r`/`calcn_l` carry the model's *only* 3 foot markers: `r_calc` (hindfoot), `r_toe` +
`r_5meta` (forefoot) — confirmed via the `MarkerSet` API. Splitting calcn (unlike dial-turn 1's
zero-marker `toes` split) **forces** 2 of 3 markers to move to the new body. A live, easily-missed
naming gap: **the model's marker names are not symmetric between sides** —
`r_calc`/`r_toe`/`r_5meta` (lowercase) on the right, `L_calc`/`L_toe`/`L_5meta` (uppercase `L`) on
the left, verified by grepping every `<Marker name=...>` in the raw XML. A script assuming
lowercase `l_*` would have silently found nothing on the left side; handled explicitly (`MARKER_NAMES`
dict in the script), not assumed symmetric.

**Muscle re-routing — measured, not the 2 muscles dial-turn 1 touched, but 8.** Introspecting
every `Muscle`'s `GeometryPath` for points parented on `calcn_r` found **11 muscles** with any
point there. **8 are extrinsic tendons with multiple points spanning calcn's full
anterior-posterior extent** — `edl`/`ehl`/`fdl`/`fhl` (extensor/flexor digitorum+hallucis longus),
`perbrev`/`perlong` (peroneus brevis+longus), `tibant`/`tibpost` (tibialis anterior+posterior) —
every one of these needed its calcn-side points examined and (for the anterior ones) reassigned to
the new forefoot body. Only **3** (`gaslat`/`gasmed`/`soleus` — triceps surae via the Achilles
tendon) have a single, purely posterior point (local x=0.00514 m, the calcaneal tuberosity) and
needed no reassignment. This is the concrete, measured form of the task's "plantarflexors/toe
flexors would need to span the new joint" concern.

**Existing `mtp` joint required re-parenting** (not just adding a parallel joint, unlike dial-turn
1's `mtp1`): `mtp`'s parent-offset frame sits at local x=0.209 m on calcn — anterior of the split
point — so the joint itself had to move from calcn (hindfoot) to the new forefoot body. Verified
via a pre-tested-in-isolation JointSet remove+re-add mechanic (confirmed working before use in
the real build) — structurally the single biggest change in this script.

**Anthropometry (LaiArnold subject2 scaling already applied).** DOF and marker requirements above
are measured directly on `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` (subject2) —
the same scaled instance the rest of this project's MSK pipeline (IK/ID/muscle work) already uses.

---

## 3. The build — midtarsal joint, methodology + result

### 3.1 Split threshold (geometric, not a heuristic — and a self-correction worth keeping)

Pooled the local-x of every muscle with ≥2 points on calcn (excluding single-point muscles, whose
lone point carries no information about *where* a natural break sits — verified: including them
lets the single Achilles point dominate a naive gap search and pick an anatomically meaningless
gap). **First attempt used "largest consecutive gap"** — picked x=0.170 m. **Forced OODA check
before accepting it**: this data has *at least 3* comparably-large local gaps (0.024 m, 0.019 m,
0.038 m) — not cleanly bimodal — so "largest gap" is itself an arbitrary tie-break dominated by
just 2 adjacent samples. **Replaced with minimum within-cluster variance** (Fisher 1958 optimal
1-D 2-partition, exhaustive over all candidate splits) — a standard statistic using every sample,
not just two. Lands at **x=0.1308 m** (both sides, symmetric with the excluded-muscle-set
identical by construction). **Independently corroborated** (over-determination, not circularity):
this sits only **0.0141 m** from calcn's own real, pre-existing center of mass (0.1168 m) — a
completely independent geometric quantity (the segment's actual mass distribution has nothing to
do with tendon via-points) landing close by.

### 3.2 Mass split (constraint-derived, not a bare guess)

Unlike dial-turn 1's disclosed *assumed* 30/70 hallux split, this fraction is **solved** so that
two point-masses (at the centroid of each sub-body's own AP span) exactly reproduce calcn's real
total mass **and** real center-of-mass location — two genuine invariants, not a free parameter.
Result: **hindfoot 45.6%/49.2%, forefoot 54.4%/50.8%** (r/l) — `m_hindfoot`=0.591/0.639 kg,
`m_forefoot`=0.706/0.659 kg (right/left asymmetry reflects real measured differences between the
subject's two calcn segments, not a bug). Total model mass conserved exactly: 78.2 kg → 78.2 kg
(machine-asserted, not eyeballed). Still a simplification (assumes uniform density within each
piece's own span, ignores the parallel-axis shift for inertia) — disclosed, not hidden.

### 3.3 Axis orientation — forced through 3 live variants (the real OODA work of this session)

No verified numeric citation for the true midtarsal axis exists this session (full text
paywalled). Three variants tested, in order:

| Variant | Choice | Result on walking1 |
|---|---|---|
| v1 | Reuse subtalar's own oblique axis exactly (same-family placeholder) | **REJECTED** — passed topology/moment-arm checks but subtalar_angle RMSE vs external baseline jumped to 3.9–6.0° (from ~0.2–0.27° in the unmodified model). Diagnosed: parallel axes at different pivot points create a near-redundant direction in the marker-residual Jacobian (both DOFs can trade off to explain the same r_toe/r_5meta residuals) — classical foot biomechanics describes oblique tarsal-joint axes as deliberately *non-parallel* for exactly this reason (the transverse-tarsal "locking mechanism", Manter 1941 lineage — not independently re-verified this session). |
| v2 | Pure Z-axis hinge (deliberately non-parallel to subtalar) | **REJECTED** — fixed subtalar (RMSE down to 0.71–0.84°) but broke ankle instead (RMSE up to 2.3–3.3°, from ~0.14–0.17°) — this model's own ankle axis is itself close to pure-Z, so the aliasing simply moved to the other neighbor. Proves the issue is structural (inserting a new oblique DOF between two already-oblique joints sharing a sparse marker set), not a one-off bug fixable by any single simple guess. |
| **v3** | Z-axis rotated 90° about Y — a stand-in for the classical **"longitudinal axis"** component of the transverse tarsal joint (distinct from the "oblique axis" component tried in v1) | **ACCEPTED** — passed topology, moment-arm, muscle-length-invariant, AND the walking1 regression gate (median RMSE 0.123°, max-abs 3.69°, both under the 0.5°/5.0° pre-registered ceiling). Subtalar RMSE 0.91–1.44°, ankle RMSE 1.6–1.8° — both elevated vs. baseline but within ceiling. |

**Still a placeholder** (no verified numeric angle — the anatomical motivation for "try a
longitudinal-axis-like direction" is a recalled, NOT independently re-verified this session,
classical two-axis description of the Chopart joint). Because v3 was *selected* using the same
walking1 trial the regression gate checks, this is exactly a **train-on-test risk** — addressed
in §4.

### 3.4 Build result (machine-checked, reload-from-disk)

| Check | Result |
|---|---|
| Model round-trips through XML | PASS |
| Bodies/joints/coordinates | 22→24 bodies, 22→24 joints (net: +2 `midtarsal_r/l`; `mtp_r/l` removed+re-added = net 0), 35→37 coordinates |
| Total mass conserved | 78.19999999999997 kg → 78.19999999999997 kg (exact) |
| Topology: exact chain talus→[subtalar]→calcn(hindfoot)→[midtarsal]→forefoot→[mtp]→toes | PASS, both sides |
| New DOF mechanically live | forefoot moves 1.9–9.1 mm for 15° midtarsal perturbation (depending on axis variant tested) |
| Downstream coupling (toes moves WITH forefoot) | PASS — correctly a *serial chain* effect, not a bug (unlike dial-turn 1's parallel-branch independence requirement) |
| Upstream unaffected (calcn/talus don't move when midtarsal is perturbed) | PASS, after a forced recalibration: first attempt used a 1e-6 m tolerance and failed (calcn "moved" 4.9e-6 m) — an angle-sweep diagnostic (0.001°→60°, same method dial-turn 1's own docstring already used) showed this is a **flat ~5e-6 to 2.6e-5 m double-precision reassembly floor, not linear-in-angle coupling** — same class of floor dial-turn 1 already documented (~1e-5 m there). Retuned to the same 1e-4 m precedent tolerance, with the sub-linearity itself now machine-verified (ratio 60°-displacement/0.001°-displacement ≈ 0.5, i.e. *less* motion at the larger angle — pure noise). |
| Muscle-length invariant at neutral pose (all 11 calcn-crossing muscles, grafted vs. original) | max abs diff = **1.11e-16 m** (machine floating-point zero) — proves the split is geometry-EXACT at the reference pose, not approximate |
| Moment arms about `midtarsal_angle` at neutral | crossing muscles: 6/8 nonzero and finite (edl −0.0070, ehl +0.0121, fdl ~0, fhl +0.0053, perlong +0.0057, tibant +0.0338 m); Achilles group (negative control, should NOT cross the new joint): **all exactly 0.0** — correct |
| Real IK, walking1.trc | 158 frames, 0 NaN, monotonic time, regression median 0.123°/max-abs 3.69° (ceiling 0.5°/5.0°) — PASS |
| New DOF ROM (walking1) | `midtarsal_angle_r/l` range 22.9–25.7°, non-null (unlike dial-turn 1's mtp1 — this DOF keeps 2 real markers on forefoot) |

---

## 4. Held-out validation (never train-on-test) — the honest generalization finding

The v3 axis was *selected* using walking1 — trusting it required checking a trial that had zero
influence on that selection. **DJ2** (drop-jump 2, subject2) qualifies: an independently
pre-computed reference (`DJ2.mot`, produced by the original 2021 pipeline on the untouched
3-segment model) already exists, and per the Fidelity Plan's own cross-trial table, drop-jump is
a genuinely different, higher-ROM movement regime (e.g. ankle range 70–80° vs. walking's 27–33°).

**Null control, run first**: fresh IK on the *unmodified original* model, same DJ2 trial, vs. the
pre-existing `DJ2.mot` — this isolates the toolchain's own reprocessing noise floor, independent
of anything this script changed. Result: median RMSE 0.070°, subtalar RMSE 0.31°/0.27° (r/l),
ankle RMSE 0.27°/0.23° — **essentially identical in magnitude to walking1's own established noise
floor** (0.271°/0.203° subtalar, per the ENV doc) — ruling out "DJ2 is just inherently noisier."

**Grafted model on DJ2**: regression median 0.246° (under the 0.5° ceiling — PASS), but max-abs
5.76° (**narrowly exceeds** the 5.0° ceiling — FAIL), driven by `subtalar_angle_l`. Subtalar RMSE
1.29°/2.78° (r/l) — **4.2×/10.1× the null-control floor**, ankle RMSE 2.40°/3.24° (also elevated).
`midtarsal_angle` itself stays large and non-null (17–25° range, consistent with walking1).

**Conclusion (measured, not assumed): the v3 placeholder axis reduces but does not fully
eliminate a real identifiability interaction with the neighboring ankle/subtalar joints, and the
residual grows on a more dynamic movement regime.** This is a genuine, quantified, disclosed
limitation of *this specific placeholder choice* — not a bug (the skeletal build itself is
geometry-exact, per §3.4's muscle-length-invariant check) and not fixable by more guessing within
this session's scope (3 variants were already forced). The script's own exit code is
**conservative**: it reflects this DJ2 gap (exit 1), while a separate, explicitly-labeled
`OVERALL_CORE_mechanical_topological_walking1` line reports PASS — so a reader (or an automated
caller) sees both "the skeletal/geometric graft is sound" and "the specific axis angle chosen has
a measured cross-regime generalization gap" as two distinct, un-conflated facts.

---

## 5. Honest gaps / what's needed to go further

1. **A real, cited midtarsal axis angle** would very plausibly close most of the §4 gap — this
   session forced 3 principled variants but could not access the one thing that would settle it:
   the Malaquias/Postolka full text (paywalled) or a from-scratch CT-derived axis. Next step: a
   human-identity request to the KU Leuven (Jonkers/Malaquias/Postolka) or U Queensland
   (Lichtwark/Maharaj) groups (same pattern as the standing OrthoLoad request), or full-text
   access via a library/institutional route.
2. **True per-bone (~26-segment) fidelity is explicitly OUT of scope** and needs subject-specific
   CT/MRI imaging — no amount of further placeholder-tuning of the *existing* mocap-derived model
   substitutes for real bone geometry. This dial-turn (and dial-turn 1) are first steps on that
   ladder, not the ceiling.
3. **The core marker-observability limit from the Fidelity Plan (§2) still applies and now shows
   a second face**: dial-turn 1's zero-marker hallux DOF was *architecturally* guaranteed null.
   This dial-turn's midtarsal DOF keeps 2 real markers and is *not* null (20–29° ROM in both
   trials) — but a marker-informed new DOF inserted between two already marker-informed
   neighbors is exposed to a *different* risk dial-turn 1 never faced: axis-orientation-dependent
   aliasing. Going further (tarsometatarsal split, hallux+midtarsal composed together) will need
   either real axis data (item 1) or additional markers/keypoints that specifically discriminate
   the new segment from its neighbors (e.g. an actual midfoot marker) — the mocap protocol used
   here has none, and MediaPipe's video pipeline has even fewer foot landmarks (per the Fidelity
   Plan's own finding).
4. **Composing dial-turn 1 (hallux) and dial-turn 2 (midtarsal) into one model** is a natural next
   increment, explicitly deferred: it would require re-parenting dial-turn 1's `mtp1` joint from
   calcn onto the new forefoot body too (anatomically, the 1st metatarsal is forefoot, not
   hindfoot) — not done here to keep this dial-turn's evidence trail isolated and legible.
5. **GRF/kinetic partitioning across the new segment** (needed for inverse dynamics once a
   dynamic trial is used) is a separate, actively-researched sub-problem (Bruening et al.,
   already flagged in the Fidelity Plan) — not addressed by this kinematic-topology increment.
6. **Inertia is not parallel-axis corrected** (moments scaled by mass fraction off the original
   tensor, same simplification dial-turn 1 used) — fine for a load/solve/moment-arm feasibility
   test, not for a dynamics-accuracy claim.

---

## 6. File index

- `scripts/msk/foot_multisegment.py` — the build+verification script (re-runnable:
  `.venv-msk/bin/python3 scripts/msk/foot_multisegment.py`; exit code reflects ALL checks
  including the DJ2 held-out gap).
- `data/msk_models/LaiArnold_midtarsal_proto_scaled.osim` — the new 4-segment-per-foot model
  (original 3-segment source untouched, verified via mtime+md5 before/after every run this
  session).
- `data/msk_smoketest/foot_multisegment_proto/` — walking1 IK setup/output +
  `foot_multisegment_summary.json` (full machine-readable evidence, all phases).
- `data/msk_smoketest/foot_multisegment_proto_heldout_DJ2/` — DJ2 held-out IK setup/output +
  the null-control (original-model) rerun.
- `docs/MECHANISM_FOOT_FIDELITY_PLAN.md` — dial-turn 1's survey + prototype (prerequisite reading,
  not re-litigated here).
- `docs/MECHANISM_FOOT_MULTISEGMENT.md` — this document.
