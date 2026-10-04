# MECHANISM FOOT FIDELITY PLAN — survey + recommendation + prototype (2026-07-21)

Answers `docs/MECHANISM_MSK_BUILD_PLAN.md` §6 honest-gap #4 ("detailed foot beyond subtalar+mtp — no
specific model verified this session, WebSearch quota exhausted"). Per the FIRST-STEP DOCTRINE, the
current LaiArnold foot (talocrural + subtalar + mtp, 3 segments) is a first step, not the ceiling. This
doc: (1) surveys freely-available detailed multi-segment foot-ankle models with LIVE verification, not
recall, (2) measures — on this project's own real data — what actually limits foot-model fidelity today,
(3) recommends an upgrade path, (4) prototypes the cheapest safe increment on a COPY.

**Method note:** WebSearch hit the same wall the build-plan agent hit (shared-session quota, 0 available
at the start of this task — confirmed by the tool's own "2000 of 2000 used" response, not assumed).
Worked around it with tools on a *different* quota: `WebFetch` (DuckDuckGo HTML result pages, SimTK/GitHub
pages directly) and raw `curl` to NCBI's E-utilities (`eutils.ncbi.nlm.nih.gov`, unauthenticated, returns
real PubMed records — no summarizer in the loop for these, i.e. no hallucination surface). Every citation
below with a PMID/DOI was fetched this way and cross-checked; nothing here is from training-data recall.

---

## 1. Survey of candidate models — verified live, not recalled

| Candidate | Verified real? | Segments / DOF | OpenSim-native? | Freely downloadable? |
|---|---|---|---|---|
| **KU Leuven "extended foot-ankle model"** — Malaquias TM, Silveira C, Aerts W, De Groote F, Dereymaeker G, Vander Sloten J, Jonkers I. "Extended foot-ankle musculoskeletal models for application in movement analysis." *Comput Methods Biomech Biomed Engin* 2017;20(2):153-159. doi:10.1080/10255842.2016.1206533, PMID 27381808. **Follow-up:** Postolka B, Killen BA, Boey H, Malaquias TM, et al. (Jonkers senior author) "Hindfoot kinematics and kinetics — a combined in vivo and in silico analysis approach." *Gait Posture* 2024;112:8-15. doi:10.1016/j.gaitpost.2024.04.023, PMID 38723393. | **YES** — real, peer-reviewed, active 2016→2024 lineage, same lab | **5 rigid segments** (talus, calcaneus, midfoot, forefoot, toes) + **5 joints** (ankle, subtalar, **midtarsal**, **tarsometatarsal**, mtp). Two DOF variants built: 15-DOF (full) and 8-DOF (reduced). 2024 paper adds an elastic-foundation cartilage-contact model at talocrural+subtalar, validated against **4D CT in 12 healthy + 4 calcaneal-fracture patients** (real numbers: talocrural ROM 15.9°±3.9° dorsi/plantarflexion, subtalar ROM 5.9°±3.9° inv/eversion). | Built from CT-scan bone geometry via a "semi-automatic tool" — the modeling *approach* is OpenSim-compatible-in-spirit (rigid segments + pin/custom joints), not confirmed to be literally an `.osim` file. | **NOT CONFIRMED.** No SimTK/GitHub/Zenodo link in either abstract. Needs subject-specific CT to personalize — a much higher data bar than anything else in this project. Contact would be a human-identity ask to the KU Leuven Human Movement Biomechanics group (Jonkers lab), same pattern as the existing OrthoLoad request. |
| **U Queensland "Joint Constrained (JC)" model** — Maharaj JN, Rainbow MJ, Cresswell AG, Kessler S, Konow N, Gehring D, Lichtwark GA. "Modelling the complexity of the foot and ankle during human locomotion... using biplanar videoradiography." *Comput Methods Biomech Biomed Engin* 2022;25(5):554-565. doi:10.1080/10255842.2021.1968844, PMID 34698598. Earlier: Maharaj, Cresswell, Lichtwark, *Gait Posture* 2017;58:159-165, PMID 28783556 ("a custom multi-segment foot model **in OpenSim**"). | **YES** — real, validated against biplanar videoradiography (dual-fluoroscopy, gold-standard dynamic bone motion) | **6 segments, 7 DOF**, single-oblique-axis joint constraints (a middle ground between naive 6-DOF-per-joint and a rigid lumped foot). RMS 2.19° (JC) vs 3.25° (naive 6-DOF model) vs the fluoroscopy truth, across walking+running. | **Confirmed OpenSim** (2017 paper states "in Opensim" explicitly) for this lineage. | **NOT CONFIRMED.** QUT ePrints mirror (`eprints.qut.edu.au/239553/`) hosts only the citation — WebFetch-checked directly: "no downloadable model file... no data-availability statement... copyright owner permission required." |
| **Utah/Shriners "Bruening/Saraswat" lineage** (the task's named candidate) | **YES, but two DIFFERENT artifact classes, both real** | **Kinematic side** (Saraswat P, MacWilliams BA, Davis RB, *Gait Posture* 2012 PMID 22192872 + 2013 PMID 22858244, "technical coordinate systems", pediatric planovalgus feet): a **marker-cluster protocol**, same class as the Oxford/Rizzoli foot models — computes segment angles directly from skin markers, is NOT a loadable dynamic multibody model. **Musculoskeletal side**: Saraswat P, Andersen MS, MacWilliams BA. "A musculoskeletal foot model for clinical gait analysis." *J Biomech* 2010;43(9):1645-52. doi:10.1016/j.jbiomech.2010.03.005, PMID 20385385 — 3-segment, muscle+ligament-driven, EMG-validated. **Kinetic-methodology side**: Bruening DA, Cooney KM, Buczek FL, *J Biomech* 2010 PMID 20825944 + *Gait Posture* 2012 PMID 22421190; Bruening DA, Takahashi KZ, *Gait Posture* 2018 PMID 29544155 — these are about **partitioning one force-plate's GRF across multiple foot segments**, a real, still-nontrivial sub-problem for KINETIC (not just kinematic) multi-segment foot analysis, not an alternative model. | Musculoskeletal-side paper co-authored by M. S. Andersen, a core **the reference model Modeling System** developer — **very likely built in the reference model, NOT OpenSim** (inferred from authorship, not stated in the abstract — flagged as inference, not fact). | Not confirmed for any of the three sub-lineages. |
| **Oxford Foot Model (OFM) → OpenSim port** | Marker protocol itself is real/widely used; an OpenSim port is **evidenced by a working precedent**: Zandbergen MA, Schallig W, Stebbins JA, Harlaar J, van der Krogt MM. *Gait Posture* 2020;77:14-19. doi:10.1016/j.gaitpost.2020.01.010, PMID 31951914 — explicitly built "musculoskeletal lower body models... in OpenSim with either a mono- or multi-SFM based on... the Newington gait model and **Oxford Foot Model**." A 2024 "OpenOFM" open-source-software paper also exists (tandfonline DOI 10.1080/10255842.2024.2448558 — resolves/exists, blocked by a 403 bot-wall on direct fetch, title only confirmed, abstract NOT independently verified this session). | Marker math is open-source (OpenOFM); at least one group's OpenSim *model* built on it is published | Not confirmed as a ready `.osim` download; the marker-protocol math (OpenOFM) is the most likely genuinely-open piece. |
| **Rajagopal "detailed foot" variant** | **Checked and REFUTED**, not just "not found" | Directly inspected this session (both `opensim.Model` API introspection and raw XML grep on `Rajagopal_2015.osim` and `LaiArnoldModified2017...osim`, corroborating the build-plan-agent's independent prior finding): **identical joint/segment topology** to LaiArnold — `ankle, subtalar, mtp`, 3 segments. | — | Rajagopal is not a source of extra foot detail. Drop this candidate. |
| **"footankle_model" on simtk.org** | **UNVERIFIED — treat as FALSE until re-checked.** | A first-pass WebFetch-summarized DuckDuckGo result claimed a SimTK project "footankle_model" with "tibia, talus, calcaneus, midfoot, forefoot, toes... 7 DOF" (suspiciously identical numbers to the real Maharaj JC model above). Forced verification: (a) `simtk.org` is **currently unreachable** — confirmed independently via both direct `curl` (DNS resolves to `171.64.65.11`, TCP connect times out) AND `WebFetch` (`ECONNREFUSED` on the same IP) — i.e. the site itself is down, not a sandbox artifact (general egress to github.com/google.com/ncbi.nlm.nih.gov all worked fine in the same session). (b) A direct phrase search `"footankle_model" simtk` returned **zero results anywhere**. Most likely explanation: the small summarizer model pattern-completed a plausible SimTK URL by mixing the real Maharaj 6-segment/7-DOF fact with an invented project slug. | — | **Do not build on this URL/name.** Re-verify from scratch when simtk.org is reachable again; do not carry the name "footankle_model" forward as if confirmed. |

**Bottom line of the survey:** there is **no ready-to-download, confirmed-free, OpenSim-native, more-than-3-segment foot-ankle `.osim` file** verifiable in this session. Every real multi-segment lineage found is a peer-reviewed *custom research model*, access-by-request or platform-mismatched (the reference model), not a SimTK-hosted turnkey package the way LaiArnold/Rajagopal/gait2392 are. This is the honest, load-bearing conclusion — it closes build-plan gap #4 with much more specificity than "not found," but the specificity is "found and gated," not "found and ready."

---

## 2. The sharper finding: what actually limits foot fidelity right now is not the skeleton — it's the markers

Verified directly against this project's own real data (not literature), via `opensim.Model`/`MarkerSet` API introspection on the exact model+data already in use (`LaiArnoldModified2017_poly_withArms_weldHand_{generic,scaled}.osim`, `LabValidation_withVideos/`):

- The **entire foot distal to the ankle** carries exactly **3 markers** in the real mocap protocol used to scale/drive this model: `r_calc`, `r_toe`, `r_5meta` — **all three parented to the SAME rigid body** (`calcn_r`). The `toes_r`/`toes_l` bodies carry **zero markers**.
- Measured `mtp_angle` range-of-motion across **6 independent real trials** (`walking1`/`walking2`/`walking3`/`DJ2` drop-jump, 4 different subjects, pulled directly from `LabValidation_withVideos/subject{2,4,6,7,8,10}/OpenSimData/Mocap/IK/*.mot`):

  | trial | ankle_angle_r range | subtalar_angle_r range | mtp_angle_r range |
  |---|---|---|---|
  | subject2/walking1 | 32.4° | 23.6° | **2.3°** |
  | subject6/DJ2 (drop jump) | 70.5° | 19.1° | **2.0°** |
  | subject6/walking2 | 27.4° | 24.1° | **1.8°** |
  | subject10/walking1 | 33.1° | 24.8° | **2.0°** |
  | subject4/DJ2 (drop jump) | 80.5° | 31.1° | **2.5°** |
  | subject8/DJ2 (drop jump) | 68.0° | 17.3° | **2.7°** |

  `mtp_angle` stays at 1.8–2.7° range in **every** trial — an order of magnitude below the genuinely marker-observed ankle (27–80°) and subtalar (17–31°) ranges, **including the drop-jump trials where real toe-break/push-off dorsiflexion should be large.** This holds across a diverse instance space (2 movement types × 4 subjects), so it isn't a one-trial fluke.
- **Geometric reason, not a heuristic:** resolving an *additional* rigid segment via optical mocap needs ≥3 non-collinear markers on **that segment specifically** (to fix its 6 DOF). The current protocol has 3 markers total for the *entire* foot-distal-to-ankle and zero on `toes`. `mtp_angle`'s near-null range is the direct, predictable consequence — it is not "solved by IK", it is sitting near its unconstrained default because nothing pulls on it. **Adding more foot segments (midtarsal, per-ray toes) without adding markers/keypoints on each new segment does not add resolvable detail — it adds more unknowns to the same underdetermined system.** This is the exact same rank/observability problem the more-detailed literature models (§1) solve by adding CT scans or dense multi-segment marker sets — inputs this project's current corpus doesn't have.
- **The markerless video pipeline is, if anything, worse on this axis**: MediaPipe Pose's 33-landmark topology gives exactly `heel`, `ankle`, `foot_index` (3 points, no independent midfoot/hallux landmark) per foot — the same or a smaller information budget than the already-insufficient mocap protocol above.

This finding directly gates the recommendation in §3: any added foot DOF is honest-by-default only as a *kinematic-topology* increment until markers/keypoints exist to inform it.

---

## 3. Trade-off: swap in a detailed literature model vs. incrementally augment LaiArnold

**Option A — swap in one of the §1 models.**
- Pro: anatomically correct midtarsal/tarsometatarsal separation (KU Leuven) or a validated oblique-axis joint-constraint scheme (Maharaj JC), each checked against a gold-standard imaging modality in its own paper.
- Con, concrete: (1) none confirmed freely downloadable — every route is a human-identity request with unknown turnaround, or a from-scratch rebuild (KU Leuven needs a real CT scan — a much higher data bar than this project has cleared for any other joint); (2) the marker-observability ceiling in §2 applies with *more* force — a 5–6-segment model needs markers/keypoints on midfoot AND forefoot AND toes independently, which neither the existing mocap protocol nor the video pipeline supplies; (3) **the literature itself flags diminishing returns exactly here** — Kim H, Kipp K, *J Orthop Res* 2019;37(10):2231-2240 (PMID 31206865) directly compared 2-segment vs 3-segment vs 5-segment foot models during hopping: the 2-segment model significantly overestimated ROM/tissue strain vs both 3- and 5-segment, **but 3-segment and 5-segment did not differ from each other in most variables** — i.e. returns flatten out right around where LaiArnold (3 segments) already sits; (4) kinetic (not just kinematic) use of any multi-segment foot model hits the GRF-partitioning problem documented by Bruening DA et al. (*J Biomech* 2010 PMID 20825944, *Gait Posture* 2012 PMID 22421190, *Gait Posture* 2018 PMID 29544155) — splitting one force-plate's measured GRF across ≥2 foot segments is itself a nontrivial, actively-researched sub-problem, an entirely new methodological dependency this project doesn't currently have.
- Contrast, for calibration: Zandbergen et al. 2020 (PMID 31951914, §1) found mono- vs multi-segment foot models differ by 7.5°±1.2° ankle dorsiflexion and ~2–2.6% triceps surae MTC length in gait — **and the difference is largest in *deformed* feet (CP patients with planovalgus/cavovarus), smallest in healthy gait.** This project's corpus (athletic/gymnastic movement, presumably healthy feet) sits in the regime where mono-vs-multi matters *least*, per this same literature.

**Option B — add ONE incremental DOF to a COPY of LaiArnold** (what this doc prototypes, §4).
- Pro: zero new toolchain/model acquisition, reuses the already-proven-and-scaled-to-10-subjects pipeline (IK/ID/muscle paths/COMAK bridge/elastic-band work already built); matches the FIRST-STEP DOCTRINE's "hotswappable dial" requirement — a single reversible increment, measured, before deciding whether to go further; sits exactly in the "3-seg → +1" zone the Kim & Kipp finding says still plausibly matters, without jumping to the "5-6-seg, CT-personalized" zone where the same literature says returns flatten AND where this project has no acquisition path.
- Con: still bounded by the §2 observability ceiling — an added DOF is honest only as a topology/feasibility increment until real markers/keypoints exist to inform it; anatomically a simplification vs. a true CT-personalized axis.

**Recommendation: Option B, staged.** Add incremental DOFs to a LaiArnold copy one at a time, each with its own load/assemble/IK-solve smoke test (cheap, per §4), rather than committing to a full literature-model swap that is presently ungated-access and would still be starved of the markers needed to make its extra segments meaningful. Revisit Option A specifically (contacting the KU Leuven Jonkers lab or the U Queensland Lichtwark/Maharaj lab — human-identity asks, same pattern as the standing OrthoLoad request) if/when either (a) a CT scan or dense foot marker/keypoint capture becomes part of the corpus, or (b) the operator wants the axis-level anatomical fidelity those groups have already validated against imaging, independent of what this project's own IK can currently resolve.

---

## 4. Prototype — built and verified this session (COPY only, original never touched)

**Script:** `scripts/msk/foot_detail_proto.py`. **New model:** `data/msk_models/LaiArnold_hallux_proto_scaled.osim` (935,636 bytes). **Outputs:** `data/msk_smoketest/foot_detail_proto/`.

**What it does:** splits the existing single-rigid `toes_r`/`toes_l` body into `toes_r/l` (kept: lesser toes 2–5, same name, same existing muscle attachments) + a **new** `hallux_r`/`hallux_l` body, connected via a **new** `mtp1_r`/`mtp1_l` PinJoint (1st metatarsophalangeal) parented off `calcn_r`/`calcn_l` in parallel with the existing `mtp_r`/`mtp_l` joint.

**Why this split, not a midtarsal joint, for the first increment (verified, not arbitrary):**
- `toes_r`/`toes_l` carry **zero markers** (§2) — splitting them risks no marker re-parenting. `calcn` carries the model's *only* foot markers — a midtarsal split would force an immediate decision about re-parenting real markers, a materially bigger first step.
- The model's **own existing muscle anatomy already distinguishes hallux from lesser toes**: introspecting every `Muscle`'s `GeometryPath` shows exactly 4 muscles cross `toes_r` — `ehl_r`/`fhl_r` (extensor/flexor **hallucis** longus) and `edl_r`/`fdl_r` (extensor/flexor **digitorum** longus, toes 2–5). The skeleton is currently coarser than its own muscle set; this increment brings it one step closer to what the muscles already imply, rather than an arbitrary cut.

**Explicit placeholders (state before trusting anything downstream of "does it load"):**
- Mass split 30% hallux / 70% lesser-toes is an assumed round number, not cadaver/CT-derived. **Total foot mass is exactly conserved** (asserted in code, machine-verified: 78.19999999999997 kg → 78.2 kg, conserved to <1e-9 kg).
- Inertia scaled by the same mass fraction off the original `toes` tensor (same assumed mass-distribution shape, not independently derived).
- The new joint reuses the *existing* `mtp` joint's location + axis convention (same point on `calcn`, same rotation sense) — real anatomy places the 1st MTP center medial/anterior to the 2–5 break; not modeled here.
- Coordinate range copied from the existing `mtp` joint (−45°..+30°) — literature 1st-MTP push-off ROM is known to run larger in places; not independently verified/cited this session, flagged as a refinement not fabricated as "verified."
- **No muscle re-routing done** — `ehl_r`/`fhl_r` still attach to the old `toes_r` body. Re-routing these two hallux-specific muscles onto the new `hallux_r` body is the natural next step, explicitly deferred, not silently skipped.

**Checks run, all machine-checked (exit code independently re-verified: `0`, not just the script's self-report; output files independently confirmed on disk via `ls`; original source model's mtime/md5 independently confirmed unchanged: `1638784176` / `0eb06c2a2b235b579b76855c8c8bff88`):**

| check | result |
|---|---|
| Model round-trips through XML (fresh reload from disk, not just the in-memory object) | PASS |
| New DOF is mechanically live (perturbing `mtp1_angle_r` by 20° moves `hallux_r`'s mass center by 7.69 mm) | PASS |
| New DOF is topologically independent of the old `mtp` joint (exact tree check: different child-base-body, same parent-base-body `calcn_r`) — the deterministic check; a first attempt at a floating-point motion-threshold check was miscalibrated (1 nanometer — tighter than double-precision assembly noise) and was caught and fixed via an angle-sweep diagnostic (0.001°→60°) showing the residual is a flat ~1e-5 m numerical floor, not linear-in-angle coupling | PASS |
| Real `InverseKinematicsTool` run, same `walking1.trc` trial as the existing smoke test | ran to completion, 158 frames, 0 NaN, monotonic time |
| Regression: do the pre-existing 35 coordinates still solve the same way with the new, marker-disconnected DOF present, vs. the pre-existing smoke-test output on the *unmodified* model? | median RMSE 0.057°, max abs diff 1.344° (pre-registered thresholds: <0.5° / <5.0° — comfortably inside, matches the original smoke test's own 0.055°/1.328° almost exactly) |
| Pre-registered prediction: new `mtp1_angle` stays near-null (unobserved by any marker, same mechanism as §2) | range 2.175°, max abs 7.20° — inside the pre-registered <5° range threshold, and near-identical in magnitude to the *existing* `mtp_angle_r`'s own 2.3° range in this same trial — an internally-consistent cross-check of the §2 mechanism, not a coincidence |

**VERDICT: PASS** (all sub-checks true; full JSON at `data/msk_smoketest/foot_detail_proto/foot_detail_proto_summary.json`).

**Honest scope of this PASS:** it demonstrates the *skeletal topology* increment loads, assembles, and IK-solves without disturbing the existing 35-coordinate solution — it does **not** demonstrate anatomical validity of the placeholder mass/axis numbers, and does **not** yet demonstrate the new DOF carries real signal (§2 predicts, correctly, that it currently cannot — no marker informs it).

---

## 5. Honest gaps / next steps

1. No confirmed-free, ready-to-download, OpenSim-native detailed (>3-segment) foot model exists per this session's live search — re-verify `simtk.org` when it's reachable again (down at time of check) before fully closing this question; do not trust the specific name "footankle_model" without independent re-confirmation.
2. Contacting the KU Leuven (Jonkers/Malaquias/Postolka) or U Queensland (Lichtwark/Maharaj) groups for their actual model files is a human-identity action (mirrors the standing OrthoLoad request pattern) — not something this session can complete.
3. The prototype's mass/inertia split (30/70) and joint-axis reuse are explicit placeholders, not validated anthropometry — fine for a load/solve feasibility test, not for any dynamics/ID accuracy claim.
4. Muscle re-routing (`ehl_r`/`fhl_r` onto the new `hallux_r` body) is identified but not done — the next cheap increment if this path is continued.
5. The core limiter is markers/keypoints, not the skeleton (§2) — before adding more foot segments, the higher-leverage move is acquiring either (a) a few extra physical markers on the foot in a future capture session (midfoot + hallux, ≥3 non-collinear points each), or (b) a foot-specific pose-estimation model with more than MediaPipe's 3 foot landmarks — without one of these, additional skeletal DOFs (this prototype's `mtp1`, or any future midtarsal joint) will predictably sit near-null exactly like the existing `mtp_angle` does now.
6. Bilateral (`_r`/`_l`) symmetry was built and verified for both sides; only the right-side numbers are quoted in tables above for brevity — the left side matches (`mtp1_angle_l` range 2.175°, identical to `_r`, confirmed in the JSON summary).

## 6. File index

- `docs/MECHANISM_FOOT_FIDELITY_PLAN.md` — this document.
- `scripts/msk/foot_detail_proto.py` — the prototype build+check script (re-runnable: `.venv-msk/bin/python3 scripts/msk/foot_detail_proto.py`).
- `data/msk_models/LaiArnold_hallux_proto_scaled.osim` — the new 4-segment-per-foot model copy (original 3-segment source untouched, verified via mtime+md5).
- `data/msk_smoketest/foot_detail_proto/` — IK setup XML, fresh `.mot` output, marker-error `.sto`, and `foot_detail_proto_summary.json` (full machine-readable results).
