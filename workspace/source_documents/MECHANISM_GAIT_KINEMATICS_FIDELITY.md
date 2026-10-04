# MECHANISM GAIT KINEMATICS FIDELITY — is the INPUT to every force number normal human gait? (2026-07-21)

Every joint-force cert in this repo (`docs/MECHANISM_JOINT_FORCE_SCORECARD.md` and its 6 children:
knee, hip, ankle, spine, shoulder, elbow) takes subject2/walking1's IK joint angles as a **given
input** and asks whether the resulting FORCES match in-vivo/cadaveric anchors. None of them
separately asked whether the input itself — the sagittal hip/knee/ankle flexion waveforms the
whole force chain is driven by — looks like normal human gait. This doc is that decorrelated check:
extract the twin's ROM/peak-angle/timing for hip flexion, knee flexion, and ankle dorsi/
plantarflexion over the gait cycle, and test them against normative healthy-gait literature and a
second, independent, real gait dataset.

**Headline: the kinematic-fidelity floor HOLDS. 16/17 pre-registered magnitude+timing gates PASS
directly; the 1 non-passing gate (hip peak-flexion timing, measured at 87%GC vs a pre-registered
88–100%GC window) is machine-diagnosed as an argmax-of-a-flat-plateau artifact, not a kinematic
anomaly** — the curve sits within 0.11° of its own value at the band edge, on a 14-%GC-wide
near-peak plateau that already overlaps the passing window. Cross-validated against a **second,
fully independent real dataset** (OpenSim's own canonical gait2392 example trial, `subject01_walk1`,
fetched live from `github.com/opensim-org/opensim-models` this session): the twin's hip/knee/ankle
ROM lands at 0.84–1.24× that trial's ROM — squarely inside normal inter-subject gait variability, not
a method-specific artifact of the normative bands chosen. **No kinematic anomaly was found that
could explain the joint-force family's 1.4–1.5× knee/hip over-prediction** — the over-prediction
documented in `docs/MECHANISM_JOINT_FORCE_SCORECARD.md` is not a garbage-in-garbage-out symptom of
bad input angles; the angles feeding that pipeline are themselves normal, textbook-consistent
walking kinematics.

**Confidence tier: published-plausibility** (normative gait databases + one real independent
dataset — a population/band-membership check, not a per-subject-matched statistical validation;
same epistemic honesty convention this repo's cadaveric/literature-tier docs already use).

## 0. Pre-registration (stated before any number below was computed)

- **Falsifier:** for each of {hip, knee, ankle} × {peak magnitude, ROM, peak timing} (right leg,
  primary — matches this repo's established single-trial/right-side-primary convention), does the
  twin's measured value fall inside a pre-registered normative band? Any joint with ≥1 magnitude
  value outside its band, or a timing value clearly outside its band, is flagged **anomalous** and
  named, not averaged away.
- **Bands set BEFORE extraction** (`NORMATIVE_BANDS` in the script, degrees / %GC): hip peak
  flexion 18–38° near 0/100%GC, peak extension −25 to −2° at 35–65%GC, ROM 32–58°; knee
  stance-phase (loading-response) peak 8–28° at 5–30%GC, swing-phase peak 45–78° at 60–85%GC, min
  −8 to 12°, ROM 45–80°; ankle peak dorsiflexion 2–18° at 28–55%GC, peak plantarflexion −30 to −8°
  at 50–72%GC, ROM 18–42°. These reproduce the task's own pre-registered envelope (knee ~60° swing
  / ~15–20° stance, ankle ~10° dorsi to ~20° plantar, hip ~30° flex to ~10° ext) widened to a
  literature-consensus population band, not a single-point target.
- **Decorrelation:** the %GC phase clock is derived from the real FORCE-PLATE GRF signal
  (`contact_waveform_analysis.model_pct_gc_mapping`, already established and re-verified by every
  force-family doc: `T_stride=1.327s`, `R_heel_strike_t=1.2425s`) — a different measurement
  modality than the marker-based IK joint angles being phase-mapped onto it. If the two modalities
  disagreed about where in the cycle events happen, the already-passing R/L toe-off plausibility
  gate (inherited, re-run, still PASS) would fail.
- **Bonus over-determination (not the primary falsifier):** does the ankle's peak-plantarflexion
  ANGLE instant land near the plantarflexor peak-FORCE instant from the independently-built SO
  muscle-force pipeline (t=0.59s/50.8%GC, `docs/MECHANISM_PUSHOFF_PLANTARFLEXOR.md`)?

## 1. Headline numbers (right leg, primary; all from `gait_kinematics_fidelity_results.json`)

| joint | quantity | twin (measured) | normative band (pre-registered) | gate |
|---|---|---:|---:|---|
| Hip | peak flexion | 25.90° @ 87%GC | 18–38° @ [0–12 ∪ 88–100]%GC | value **PASS**; timing **borderline** (see §3.1) |
| Hip | peak extension | −12.80° @ 55%GC | −25 to −2° @ 35–65%GC | **PASS** |
| Hip | ROM | 38.71° | 32–58° | **PASS** |
| Knee | stance-phase peak | 20.72° @ 18%GC | 8–28° @ 5–30%GC | **PASS** |
| Knee | swing-phase peak | 65.61° @ 74%GC | 45–78° @ 60–85%GC | **PASS** |
| Knee | midstance min | 5.60° @ 42%GC | −8 to 12° | **PASS** |
| Knee | ROM | 60.01° | 45–80° | **PASS** |
| Ankle | peak dorsiflexion | 15.84° @ 50%GC | 2–18° @ 28–55%GC | **PASS** |
| Ankle | peak plantarflexion | −16.35° @ 68%GC | −30 to −8° @ 50–72%GC | **PASS** |
| Ankle | ROM | 32.19° | 18–42° | **PASS** |

**Gates: 16/17 PASS as literally pre-registered.** The 1 non-passing gate (hip peak-flexion
timing, 87%GC vs the 88–100%GC lobe of the wraparound band) is forced through an OODA diagnostic in
§3.1 and shown to be a metric-conditioning artifact, not a magnitude, shape, or genuine timing
anomaly.

## 2. Method

Reused, unedited, zero new OpenSim run: `scripts/msk/gait_kinematics_fidelity.py` imports
`parse_mot`, `model_pct_gc_mapping`, `verify_phase_alignment` directly from
`scripts/msk/contact_waveform_analysis.py` (same convention as `scripts/msk/ankle_waveform_analysis.py`)
— the SAME already-verified GRF-derived %GC clock every force-family doc uses, so every %GC value
here is directly comparable to theirs. New: reads the IK **coordinates** `.mot`
(`walking1_smoketest.mot`, confirmed `inDegrees=yes`, 158 rows, t∈[0,1.57]s) instead of a
JointReaction `.sto` — extracts `hip_flexion_r`, `knee_angle_r`, `ankle_angle_r` (+ `_l` for a
secondary symmetry check), resamples each onto the identical 101-point %GC grid via linear
interpolation (`period=100`).

**Sign convention confirmed from data, not assumed:** gait2392 convention in this model —
`hip_flexion_r`/`knee_angle_r` positive = flexion; `ankle_angle_r` positive = dorsiflexion,
negative = plantarflexion — consistent with the printed values' own sign pattern (knee rises from
~5° to ~65° through swing; ankle goes positive in stance, negative at push-off).

**Phase-alignment re-verified before trusting anything else:** `R_toe_off_pct_gc=66.88` (band
55–68, PASS), `L_toe_off_pct_gc=17.03` (band 8–20, PASS) — both inherited gates re-run and still
PASS on this run.

**Environment:** `/usr/bin/python3` (numpy 2.2.6 + openpyxl 3.1.5 via the `contact_waveform_analysis`
import), never `.venv-msk`, never imports OpenSim — same convention as the rest of the waveform
family.

## 3. Results in detail

### 3.1 Hip — the one non-passing gate, forced through OODA, not waved away

Raw finding: peak hip flexion measured at 87%GC, 1 percentage point short of the pre-registered
88–100%GC passing lobe. **Observed** (instrumented the right observable): dumped the curve at 1%GC
resolution from 70–100%GC and the 0–15%GC wraparound (see script + JSON `curves.hip_flexion_r`).
**Oriented** (why): the curve rises steeply from 8.91° (70%GC) to 22.45° (80%GC), then goes nearly
FLAT — 23.97° (82%GC) → 25.90° (87%GC, the measured argmax) → 25.79° (88%GC) → 25.70/25.88/25.40°
(89/90/91%GC) → 24.29° (100%GC/0%GC) — a textbook flat-topped maximum. A hard %GC boundary on the
ARGMAX of a nearly-flat curve is ill-conditioned in its *location* (noise-sensitive) even though the
*value* is rock-solid. **Decided/Acted** (machine-measured, not eyeballed): computed the %GC width
where the curve stays ≥95% of its own peak value — **14 %GC wide (7%GC and 84–96%GC)** — and this
plateau **already overlaps the passing timing band** (`plateau_overlaps_passing_timing_band=true`,
since 88–96%GC ⊂ the plateau). The curve's value at the band's own edge (88%GC, 25.795°) differs
from the measured argmax value (87%GC, 25.904°) by **0.109° absolute (0.4% relative)** — negligible
against this pipeline's own Savitzky-Golay/IK-residual noise floor. **Verdict: this is a
metric-conditioning artifact of gating a flat maximum's location against a 1-%GC-wide boundary, not
a kinematic anomaly** — the hip reaches and holds near-peak flexion through the terminal-swing →
initial-contact transition exactly as the textbook gait cycle describes; a 1-%GC miss on which
frame is technically "the max" of an 0.4%-flat plateau is not a finding.

Peak extension (−12.80° @ 55%GC) sits centrally in its band and lands, as expected, in late stance
before this trial's own toe-off (66.88%GC) — textbook timing, PASS outright.

### 3.2 Knee — clean pass, correct two-peak structure

The stance-phase (loading-response) peak (20.72° @ 18%GC) and the swing-phase peak (65.61° @
74%GC, confirmed to land AFTER the trial's own measured toe-off at 66.88%GC —
`swing_peak_lands_in_swing_not_stance=true`) reproduce the textbook double-rise knee shape exactly:
small early-stance flexion for shock absorption, extension through midstance (min 5.60° @ 42%GC,
near but not quite full extension — inside band), then the large swing-phase flexion for foot
clearance. All 7 knee gates PASS with no near-misses (closest margin: swing-phase peak 65.61° is
12.6° inside its 78° upper bound).

### 3.3 Ankle — clean pass, AND a mechanistically-verified over-determination

Peak dorsiflexion (15.84° @ 50%GC) and peak plantarflexion (−16.35° @ 68%GC) both PASS outright,
magnitude and timing. The pre-registered **bonus** over-determination check (angle-peak vs. the
independently-derived plantarflexor peak-FORCE instant, 50.8%GC) misses its 15-%GC tolerance
(offset 17.2%GC) — but forcing this adversary (OODA again, not accepted at face value) shows why:
peak plantarflexOR **FORCE** and peak plantarflexion **ANGLE** are two mechanistically DIFFERENT
events, not the same event measured twice. Force peaks while the foot is still loaded in terminal
stance (driving push-off); the joint continues rotating into plantarflexion under tendon-recoil and
momentum after force has already started to fall, reaching its own extreme at/after toe-off (the
"forefoot rocker" completing). The mechanistically-correct comparison is angle-peak vs. the twin's
own independently-measured (GRF-derived) toe-off instant: **68.0%GC vs. 66.88%GC — offset 1.12%GC,
PASS.** Both numbers are reported, not just the one that passes: the naive 17.2%GC "miss" is
explained, not hidden.

### 3.4 Left-leg bilateral symmetry (secondary context, not a primary gate)

ROM asymmetry (phase-independent, not subject to the mapping's own contralateral-symmetry
construction): hip **1.38%**, knee **1.35%**, ankle **2.23%** — all small, consistent with healthy
symmetric gait. **Disclosed limitation:** the left-leg TIMING gates in the raw JSON use
right-leg-referenced windows without the ~50%GC antiphase relabeling a proper left-leg timing
analysis would need, so those specific numbers are not reported as findings here (would require a
dedicated left-leg phase convention, out of scope for this right-side-primary doc, consistent with
`docs/MECHANISM_JOINT_FORCE_SCORECARD.md`'s own disclosed scope: "left-side numbers, where reported,
are coarse order-of-magnitude symmetry checks only").

### 3.5 Second external anchor: OpenSim's own gait2392 canonical example trial

To go beyond a band-membership check, `scripts/msk/gait2392_reference_cross_check.py` fetched
(live, this session) the actual `subject01_walk1_ik.mot` — the canonical example gait trial
distributed with OpenSim itself, from `github.com/opensim-org/opensim-models`
(`Pipelines/Gait2392_Simbody/OutputReference/subject01_walk1_ik.mot`) — and compared ROM/peak
magnitudes directly (phase-independent quantities only; full %GC re-derivation for this second
trial's own multi-cycle GRF file was out of scope, disclosed, not silently skipped).

**Forced sign-convention adversary (measured, not assumed):** the reference file's raw
`knee_angle_r` ranges **[−70.18°, +1.05°]** — large NEGATIVE values, physically impossible as a
human knee flexion reading if taken at face value (a 70° reading in that direction would be
hyperextension no knee performs). `hip_flexion_r` ([−22.37°,+21.93°]) and `ankle_angle_r`
([−8.90°,+17.12°]) in the SAME file both land in small-magnitude, correctly-signed ranges matching
the twin's own convention — so this is an isolated, knee-specific flipped-sign convention (a known
historical artifact across gait2392/gait2354 distribution vintages), not a general axis mismatch.
Negating the reference knee curve produces a physiologically sane shape (min −1.05°≈0° full
extension, max +70.18° swing flexion) where the raw curve did not — confirmed via an explicit
in-script assertion (`knee_physio_sane` check on the corrected range), not silently applied.

| joint | twin ROM | gait2392-reference ROM (subject01_walk1) | ratio (twin/ref) |
|---|---:|---:|---:|
| Hip | 38.71° | 44.30° | **0.874×** |
| Knee | 60.01° | 71.23° (sign-corrected) | **0.842×** |
| Ankle | 32.19° | 26.03° | **1.237×** |

All three ratios sit within 0.84–1.24× of a fully independent real subject's gait trial — squarely
inside normal inter-subject variability (for contrast, the joint-FORCE family's own knee/hip
over-prediction ratio is 1.4–1.5×, a materially larger and differently-signed deviation). The ankle
ratio (1.24×) is the largest of the three, driven mainly by the twin's larger peak-plantarflexion
magnitude (−16.35° vs. the reference's −8.90°) — both values individually sit inside normal
published ranges (ankle push-off plantarflexion is one of the more speed/subject-variable sagittal
gait quantities in the literature), reported plainly as a real, disclosed difference between two
individual subjects, not a failure.

### 3.6 IK tracking-quality signal (a decorrelated instrumentation check on the input itself)

`walking1_ik_marker_errors.sto` (already-committed IK solver output, read not modified):
`marker_error_RMS` mean **0.0130 m**, max over the trial **0.0142 m**; `marker_error_max`
(worst single marker per frame) mean 0.0340 m, max 0.0358 m. Using this repo's own already-accepted
precedent for the same column (`docs/MECHANISM_SECOND_TRIAL_FORCE.md` §Orient: "`marker_error_RMS`
stays flat at 0.016–0.020 m... max over the whole trial 0.0279 m... ruled out [as a quality
concern]") — subject2/walking1's RMS marker error (0.0130/0.0142 m) is comfortably **better** than
that already-accepted trial's 0.016–0.029 m range. The IK solve underlying every angle in this doc
is itself good-quality, not a degraded/noisy tracking result silently feeding the rest of the
analysis.

## 4. Anchors — literature verified LIVE this session (PMIDs fetched, not recalled)

- **Fukuchi RK, Fukuchi CA, Duarte M (2018).** "A public dataset of overground and treadmill
  walking kinematics and kinetics in healthy individuals." *PeerJ* 6:e4640. DOI
  [10.7717/peerj.4640](https://doi.org/10.7717/peerj.4640), PMID
  [29707431](https://pubmed.ncbi.nlm.nih.gov/29707431/), PMCID PMC5922232. Confirmed real and
  population-matched (n=24 young adults, 27.6±4.4y, overground self-selected speed, 26-marker
  Leardini et al. 2007 set). **Disclosed limitation:** the paper itself reports ensemble-average
  waveforms as FIGURES (Figs 3, 8), not tabulated numeric ROM/peak values in text — per this
  protocol's own rule against eyeballing figures, no number was extracted from its plots; it is
  cited here as population-validity corroboration (a real, matching-methodology public dataset
  exists and is exactly the shape of gait this doc checks against), not as the numeric source of
  the bands in §0.
- **Kadaba MP, Ramakrishnan HK, Wootten ME (1990).** "Measurement of lower extremity kinematics
  during level walking." *J Orthop Res* 8(3):383–392. DOI
  [10.1002/jor.1100080310](https://doi.org/10.1002/jor.1100080310), PMID
  [2324857](https://pubmed.ncbi.nlm.nih.gov/2324857/). Confirmed real (n=40 normal young adults, 3
  sessions each) — the classic normative gait-kinematics reference underlying the standard
  hip/knee/ankle band description used throughout clinical gait analysis. Abstract confirmed live;
  full numeric tables are in the paywalled full text (disclosed, not fabricated).
- **Bovi G, Rabuffetti M, Mazzoleni P, Ferrarin M (2011).** "A multiple-task gait analysis
  approach: kinematic, kinetic and EMG reference data for healthy young and adult subjects." *Gait
  Posture* 33(1):6–13. DOI [10.1016/j.gaitpost.2010.08.009](https://doi.org/10.1016/j.gaitpost.2010.08.009),
  PMID [21123071](https://pubmed.ncbi.nlm.nih.gov/21123071/). Confirmed real (n=40: 20 aged 6–17,
  20 aged 22–72; comprehensive kinematic/kinetic/EMG reference dataset, self-selected speed
  condition included).
- **Winter DA.** *The Biomechanics and Motor Control of Human Gait: Normal, Elderly and
  Pathological*, 2nd ed., University of Waterloo Press, 1991. **Honestly disclosed: this is a
  book monograph, not a journal article — confirmed absent from PubMed under any searched title
  variant (no PMID exists for it), not fabricated.** It remains the standard textbook source for
  the normative sagittal gait-angle description used in §0's bands, corroborated (not
  independently re-derived) by the three PMID-verified papers above and by the OpenSim
  gait2392-reference curve-vs-curve comparison in §3.5.
- **OpenSim gait2392 canonical example** (`subject01_walk1`): not a paper, a real distributed
  dataset — fetched live this session, used as a direct quantitative second anchor (§3.5), not a
  citation.

## 5. Honest gaps

1. **n=1 subject/trial (subject2/walking1), right leg primary** — same scope caveat every doc in
   this family carries; left-leg check is symmetry-context only (§3.4), not independently validated
   against normative bands.
2. **Normative bands are literature-consensus, hand-set widths, not a formally pooled
   meta-analytic CI** — disclosed at pre-registration (§0); reported alongside the raw values so a
   reader can see how close/far each measurement sits from its own boundary (§1 table, §3.1's
   explicit margin numbers), not just a binary pass/fail.
3. **Fukuchi 2018's own numeric tables were not extractable via text** (graphical ensemble averages
   only) — used as population-validity corroboration, not a numeric source (§4), consistent with
   the "never eyeball a figure" rule.
4. **The gait2392-reference cross-check (§3.5) compares ROM/peak MAGNITUDE only** — its own
   %GC timing re-derivation (a separate heel-strike segmentation of that trial's own GRF file) was
   out of scope for this addendum; timing validation rests on §0's literature bands + §3.3's
   internal toe-off over-determination check, not on this second dataset.
5. **The knee sign-convention correction in §3.5 is a disclosed, measured, asserted-in-script
   correction** (not an assumption) — but it means the gait2392-reference knee comparison depends on
   that correction being right; the physiological-sanity check (`knee_physio_sane=true`, corrected
   range [−1.05°,+70.18°]) is the evidence it is.
6. **Savitzky-Golay / IK-solver conventions are inherited unchanged** from the rest of this repo's
   waveform family — a design choice, not the literal raw marker signal, same caveat as every
   sibling doc.
7. **The 87%GC hip-timing near-miss (§3.1) is resolved via a machine-measured flatness diagnostic
   this session, not a pre-existing gate** — a legitimate OODA-forced resolution (measuring WHY,
   not just reporting pass/fail), but it is a new diagnostic invented to explain this specific
   near-miss, disclosed as such rather than presented as if it had been pre-registered alongside
   the original bands.

## 6. Verification (machine-checked, not narrated)

- Phase-alignment: both inherited out-of-sample gates re-run and PASS
  (`phase_alignment_verification.all_pass=true`: R toe-off 66.88%GC ∈ [55,68]; L toe-off
  17.03%GC ∈ [8,20]).
- IK `.mot` confirmed `inDegrees=yes`, 158/158 rows, required columns present (asserted in-script,
  would raise otherwise).
- 16/17 pre-registered right-leg gates PASS by direct threshold comparison (`right_leg.*.gates`,
  JSON); the 1 non-passing gate's flatness diagnostic is machine-computed
  (`peak_flatness_diagnostic`), not eyeballed from a printed table.
- Knee sign-flip in the gait2392 reference detected by an in-script assertion
  (`knee_sign_flip_detected`, requires raw max<10° AND raw min<−50°) before any correction was
  applied, and the correction's physiological sanity independently checked
  (`knee_physio_sane`, requires corrected min∈[−10,15]° AND corrected max∈[45,90]°) — both `true`.
- Bonus over-determination: both the naive (vs. force-peak) and mechanistically-corrected (vs.
  toe-off) comparisons computed and reported (`bonus_overdetermination_check`), not just the
  passing one.
- All numbers in this doc are pulled directly from the two JSON files listed below, not
  transcribed from console prose.

## 7. Files

- This doc: `docs/MECHANISM_GAIT_KINEMATICS_FIDELITY.md` (new, this session).
- `scripts/msk/gait_kinematics_fidelity.py` (new) — primary extraction + normative-band gating.
  Run: `/usr/bin/python3 scripts/msk/gait_kinematics_fidelity.py` (needs numpy+openpyxl
  transitively via its import of `contact_waveform_analysis`; never `.venv-msk`; never imports
  OpenSim). Imports `parse_mot`, `model_pct_gc_mapping`, `verify_phase_alignment` from
  `scripts/msk/contact_waveform_analysis.py` (read-only, not edited) — exit 0.
- `scripts/msk/gait2392_reference_cross_check.py` (new) — second-anchor ROM/peak cross-check
  against the OpenSim gait2392 canonical example trial. Run: `/usr/bin/python3
  scripts/msk/gait2392_reference_cross_check.py` — exit 0.
- `data/msk_smoketest/subject2_walking1/gait_kinematics_fidelity/gait_kinematics_fidelity_results.json`
  (new) — every number in §1–§3.4, §3.6: the %GC mapping + alignment verification, the IK `.mot`
  metadata, the full 101-point %GC curves for all 6 (R+L × hip/knee/ankle) coordinates, per-joint
  feature extraction + gates + the hip flatness diagnostic, the left-leg symmetry context, and the
  bonus over-determination check (both variants).
- `data/msk_smoketest/subject2_walking1/gait_kinematics_fidelity/gait2392_reference_cross_check_results.json`
  (new) — every number in §3.5: pre/post sign-correction ranges, the sign-flip-detection and
  physiological-sanity assertions, and the twin-vs-reference ROM comparison table.
- `data/external/opensim_gait2392_reference/subject01_walk1_ik.mot` +
  `subject01_walk1_grf.mot` (new, external) — fetched live this session from
  `github.com/opensim-org/opensim-models` (`Pipelines/Gait2392_Simbody/`, Apache-2.0-licensed
  OpenSim example data), persisted for reproducibility.
- Inputs read in place, never modified: `data/msk_smoketest/subject2_walking1/walking1_smoketest.mot`
  (IK joint angles), `data/msk_smoketest/subject2_walking1/walking1_ik_marker_errors.sto` (IK
  tracking quality), `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/ForceData/walking1_forces.mot`
  (GRF, for the %GC mapping, via the imported `contact_waveform_analysis` function).
- **Not touched:** `scripts/msk/contact_waveform_analysis.py` (imported read-only),
  `scripts/msk/ankle_waveform_analysis.py`, and any concurrent, independent instance's in-progress
  work.

No git commit, no git push performed (isolation respected). All new files are untracked, for the
coordinator to commit.
