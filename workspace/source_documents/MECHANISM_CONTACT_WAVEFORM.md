# MECHANISM CONTACT-FORCE WAVEFORM — full stance-phase shape, not just the peak (2026-07-21)

Every knee/hip force cert in this repo so far (`docs/MECHANISM_JOINT_FORCE_SCORECARD.md`,
`docs/MECHANISM_CMC_SECOND_SOLVE.md`, `docs/MECHANISM_CROSS_SUBJECT.md`) compares exactly ONE
number — the PEAK contact force over the trial — against the OrthoLoad in-vivo median/peak.
In-vivo knee/hip contact force during gait is a characteristic DOUBLE-HUMP waveform
(early-stance weight-acceptance peak + late-stance push-off peak); collapsing the whole
time-course to a single max cannot tell "matches shape but over-shoots one hump" from
"over-predicts uniformly" from "wrong hump ratio" — very different trust stories. This doc
extracts the FULL subject2/`walking1` SO **and** CMC knee-r/hip-r contact-force TIME SERIES,
phase-aligns it to %gait-cycle (heel-strike to heel-strike), and compares the whole waveform —
not just the peak — to an OrthoLoad in-vivo ensemble waveform built from real per-subject/
per-trial OrthoLoad data (not a digitized figure — SCENE-EYES: numbers, not pixels).

**Headline: the model DOES reproduce the real double-hump shape (Pearson r=0.95 knee, 0.76-0.78
hip, both far above noise). The uniform-vs-concentrated question resolves MIXED, not clean, but
with a consistent, 4-for-4 directional signal that survives every one of the four
solver×joint combinations tested: the push-off hump is ALWAYS over-predicted more than the
weight-acceptance hump.** For SO (this repo's primary, most-used solve method), the effect is
sharp and clean: the knee's weight-acceptance hump (245.4 %BW) is statistically
indistinguishable from real OrthoLoad in-vivo data (224.6±21.6 %BW, Z=0.97, ratio 1.09×) while
its push-off hump (390.3 %BW vs 266.9±46.2 %BW, Z=2.67, ratio 1.46×) is a clear, significant
over-prediction — crossing the pre-registered 25%-relative-spread threshold into
**CONCENTRATED**. CMC and both hip series land just inside the **UNIFORM** side of that same
threshold (16-24% spread, vs the 25% cutoff), but their own Z-scores still show the push-off
hump 1.3-3.8 SDs further from the population mean than the weight-acceptance hump in every
single case — so even the "uniform" verdicts are not phase-blind; they are "elevated everywhere,
elevated more at push-off." **A genuinely new finding surfaced by looking at the full
waveform for the first time: every prior peak-only cert in this repo (391.10/427.54 %BW knee,
etc.) was, without anyone knowing it, always measuring the PUSH-OFF hump (44.8-47.8 %GC) — never
the weight-acceptance hump — because that is where this model's global maximum happens to sit.**

## Headline numbers

| Series | hump1 (weight-accept.) model/ortho ±SD | ratio₁ | Z₁ | hump2 (push-off) model/ortho ±SD | ratio₂ | Z₂ | rel. spread | verdict (±25% band) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| SO knee | 245.4 / 224.6±21.6 %BW | **1.09×** | 0.97 | 390.3 / 266.9±46.2 %BW | **1.46×** | 2.67 | 28.9% | **CONCENTRATED** |
| CMC knee | 295.7 / 224.6±21.6 %BW | 1.32× | 3.30 | 425.1 / 266.9±46.2 %BW | 1.59× | 3.42 | 19.0% | UNIFORM |
| SO hip | 326.3 / 249.9±44.4 %BW | 1.31× | 1.72 | 386.9 / 232.2±34.0 %BW | 1.67× | 4.56 | 24.3% | UNIFORM (borderline) |
| CMC hip | 393.2 / 249.9±44.4 %BW | 1.57× | 3.23 | 429.2 / 232.2±34.0 %BW | 1.85× | 5.80 | 16.1% | UNIFORM |

Waveform-level shape match (whole 0-100%GC curve, not just the two humps): SO knee r=0.951
RMSE=48.0%BW; CMC knee r=0.952 RMSE=73.6%BW; SO hip r=0.756 RMSE=91.8%BW; CMC hip r=0.782
RMSE=124.5%BW (hip RMSE is inflated by the hip curve's own much higher absolute magnitude, not
a proportionally worse fit — see §3.1). Mean ratio₁ across all 4 series = 1.32×; mean ratio₂ =
1.64×; **ratio₂ > ratio₁ in 4/4 series** (n=4 is a small basis for a formal significance claim on
its own, stated honestly, but it converges with the Z-score pattern, the SO-knee threshold-cross,
and a documented gastroc/soleus-driven-push-off mechanism from the literature — an
over-determined, not a lone, signal).

**Confidence tier: in-vivo-anchored for knee** (OrthoLoad's own official, `averDTW`-aligned,
%BW-normalized per-subject curves, 8 subjects) — **in-vivo-anchored-but-self-built for hip** (the
official OrthoLoad "Standardized Loads" hip archive is a corrupted/truncated download already
present in this repo's data folder, §6; the hip ensemble used here is built directly from raw
per-trial OrthoLoad AKF force telemetry — real patient data — using a segmentation method
explicitly calibrated against the official knee curve first, r=0.928, before being trusted on
hip, but it is this session's own derivation, not an OrthoLoad-published curve).

## 1. Pre-registration (stated before the final consolidated hump numbers were computed)

**Phase-alignment validity gate (symmetric-QC, forced BEFORE trusting any hump comparison):**
the model's %GC mapping and the raw-AKF segmentation method (needed for hip, where no official
%GC-normalized curve exists) must each pass an independent, out-of-sample plausibility check —
see §2.1/§2.3 — else the whole comparison is blocked. Both passed.

**Uniform-vs-concentrated decision rule:** let r1, r2 be the model/OrthoLoad ratio AT hump 1 and
hump 2 respectively. `|r1-r2| / mean(r1,r2) < 0.25` → **UNIFORM** (consistent with a scaling /
over-strong-model / contact-model-wide explanation); `>= 0.25` → **CONCENTRATED** (localizes the
error in TIME, implicating the muscles/mechanisms dominant in that specific phase). The 25% band
matches this repo's own established convention (`docs/MECHANISM_CMC_SECOND_SOLVE.md`'s
pre-registered 25% SO-vs-CMC agreement band) rather than being invented for this doc.

**Calibration gate (method-validation, run to decide whether the self-built hip ensemble could be
trusted at all):** the raw-AKF trough+steepest-rise segmentation must reproduce the OFFICIAL
OrthoLoad knee ensemble at Pearson r > 0.8. Measured 0.928 — **PASS**. This is a legitimate
method-validation check run during development, not post-hoc threshold-shopping on the actual
claim (the actual hump-ratio/verdict numbers above were computed only after this gate passed).

## 2. Method

### 2.1 Model heel-strike-to-heel-strike %GC mapping

subject2/`walking1`'s own GRF file (`ForceData/walking1_forces.mot`, 2000 Hz, real force-plate
data) contains only **one** full right-foot heel-strike inside its 1.579s window — the trial
starts already mid-right-stance (`R_ground_force_vy(t=0)=641.0 N`, `L_ground_force_vy(t=0)=0 N`).
Directly measured events (20N threshold-crossing, right and left vertical GRF): R toe-off
t=0.803s, R heel-strike t=1.2425s, L heel-strike t=0.579s, L toe-off t=1.4685s. Because only ONE
right heel-strike falls in-window, the cycle-closing boundary is recovered via the standard
contralateral-symmetry convention (opposite-foot heel-strikes are 50% GC apart in periodic gait):
`T_stride = 2×(1.2425−0.579) = 1.327s`. **Verification (the "verify the alignment" step this
task's own mandate requires, so a phase offset cannot silently manufacture a false shape
mismatch):** two events NOT used to build the mapping are checked against their textbook-expected
bands — ipsilateral (R) toe-off lands at **66.9%GC** (band 55-68%, PASS) and contralateral (L)
toe-off lands at **17.0%GC** (band 8-20%, end of initial double support, PASS). Both independent
checks pass; `phase_alignment_verification.all_pass = true` in the evidence JSON.

### 2.2 Model knee-r/hip-r force time series

Re-extracted directly from the ALREADY-COMMITTED, already-certified `opensim.JointReaction` `.sto`
outputs — SO: `static_optimization/jr/walking1_JointReaction_ReactionLoads.sto`; CMC:
`cmc_second_solve/jr_cmc/walking1_JointReaction_ReactionLoads.sto` — **no new OpenSim run**.
Self-consistency gate (machine-checked before trusting anything else): this script's own
re-derived overall-peak values match `docs/MECHANISM_CMC_SECOND_SOLVE.md`'s already-published
numbers to **<0.01% relative difference** on all 4 series (391.11/387.04/427.54/429.26 %BW) —
confirms the new extraction code is bug-free against already-trusted numbers before anything new
is built on top of it.

### 2.3 OrthoLoad ensembles

**Knee (official):** the per-subject `data/external/orthoload/knee/standard_loads/K*_Walking.xlsx`
"Dat_BW" sheets (8 subjects: K1L,K2L,K3R,K5R,K6L,K7L,K8L,K9L) — OrthoLoad's own `averDTW`
(Dynamic-Time-Warping)-aligned, %BW-normalized, per-subject representative-cycle curve.
**A forced adversary caught here:** the sibling "Dat_AVER75"/"Dat_HIGH100" sheets in the SAME
workbook look like an obvious multi-subject-average alternative, but their own header states
*"The individual peak value of Fres was adapted to the average peak of Fres from all subjects...
All 6 load components were then adapted by the same Multiplication Factor"* — these are
implant-fatigue-TEST-RIG standard loads (peak forcibly harmonized across subjects for mechanical
test rigs), not a physiological ensemble. Caught by reading the sheet's own disclosure text before
using it, not assumed innocent. "Dat_BW" carries no such rescaling note and its 8 subjects'
peaks span a real 202-344 %BW range (not harmonized) — used instead.

**Hip (self-built, no official alternative — see §6):** the equivalent
`StandardLoads-Hip_CompleteData.zip` is a truncated download already sitting in this repo's data
folder (valid `PK\x03\x04` local-file-header magic but no end-of-central-directory record anywhere
in 117MB — confirmed via direct byte search, not a tool-version quirk). Built instead from raw
per-trial OrthoLoad AKF force telemetry: the SAME "primary" trial selection
`validate_hip_force.py` already uses (gen1 "walking free" OR gen2 "level walking", excluding
crutch-assisted) — 162 trials, 155 yielding valid cycles, 757 cycles pooled. Segmentation:
`find_swing_troughs` (greedy global-minimum, value-thresholded at 40% of trial peak) +
steepest-post-trough-rise heel-strike proxy — see the FORCED FIX below.

**Forced adversary #1 (segmentation bug, caught via OODA, not shipped blind):** a first version
used per-point local-minima + temporal non-max-suppression. Observed: the first pooled ensemble
showed only a single broad hump, no push-off hump — contradicting every individual raw trace
inspected by eye. Oriented: dumped one file's own minima list
(`k1l_110108_1_86p.akf`) and found the detector had latched onto the MID-STANCE dip (~190 %BW, a
real but shallow local minimum between the two humps) as a cycle boundary just as often as the
true SWING trough (~20-27 %BW), since the two are >0.4s apart (passes a purely temporal
separation test but not a value-based one). Decided: swing troughs are always far deeper than
mid-stance dips for a limb bearing body weight (measured ~7-9% of trial peak vs ~65-70%) — fixed
to greedily take the GLOBAL minimum first (value priority), excluding a window around each pick,
stopping once the next candidate exceeds 40% of the trial peak. Acted: re-ran, got clean
~1.1-1.2s cycles (mean duration 1.09-1.10s, SD 0.125s both joints — physiologically sane) with the
expected double-hump shape visible in individual traces. Converged, no further loop needed.

**Forced adversary #2 (hump-finder bug, same discipline):** a first "top-2 local maxima anywhere
in [0,62%GC]" hump-finder mis-reported CMC-knee's hump1 at 44%GC (the SAME region as hump2) because
SO/CMC's per-frame solve has real early-stance jitter (several shallow local minima at
12/17/20/23/27%GC, not one clean dip), so "top-2 by value anywhere" sometimes picked two nearby
wiggles of the SAME late broad hump. Fixed to a fixed 30%GC split (measured to sit cleanly in the
24-30%GC dip band shared by all 6 curves in play, model and OrthoLoad alike) + plain argmax on each
side — verified it clips no true peak on any of the 6 curves.

**Environment note (disclosed, not hidden):** this analysis needs numpy + openpyxl and touches NO
OpenSim API (it only re-parses already-computed `.sto`/`.akf`/`.xlsx` files already on disk) — run
with the system `/usr/bin/python3` (numpy 2.2.6, openpyxl 3.1.5 confirmed present), NOT
`.venv-msk` (which lacks both, by this repo's own established scipy/openpyxl-avoidance
convention — `static_opt_knee.py`'s own Savitzky-Golay implementation carries the same "no scipy
in .venv-msk" note). OpenSim is never imported.

## 3. Results in detail

### 3.1 Does the model reproduce the double-hump shape?

Yes, in all 4 series and both OrthoLoad ensembles: a clear early-stance rise to hump 1
(12-24%GC), a dip (24-30%GC, shallower for the model than for OrthoLoad — the model's per-frame
solve is noisier here, see the smoothing-sensitivity check below), a second rise to a TALLER
hump 2 (41-48%GC), then a fall toward a swing-phase trough (~75-85%GC) before rising back into the
next cycle. This qualitative shape is universal (dip-then-taller-second-hump) across SO, CMC, and
both OrthoLoad-derived ensembles — hump2 > hump1 for knee in every series (model AND OrthoLoad);
for hip, hump1 > hump2 in the OrthoLoad ensemble (249.9 vs 232.2 %BW) matching the textbook
"HIP98"-style pattern where the first (weight-acceptance) peak dominates, while BOTH model series
still read hump2 slightly higher than hump1 in absolute terms even though the shape's relative
hump-height ORDER differs from real hip data — itself informative (the model doesn't just scale
real hip loading, it changes which hump dominates). Waveform correlation: knee r=0.95 (both
solvers) is a strong, decisive shape match; hip r=0.76-0.78 is weaker, most plausibly because the
hip OrthoLoad ensemble is this session's own raw-AKF-pooled build (naive linear-phase averaging
blurs sharp features relative to OrthoLoad's own DTW-based averaging, exactly the attenuation the
calibration step already quantified on knee — RMSE 36.2 %BW there, comparable to that ensemble's
own 1-SD band) rather than a sign the model's hip shape is qualitatively wrong.

### 3.2 Is the over-prediction uniform or concentrated?

The pre-registered ratio-spread rule gives a MIXED answer (1 CONCENTRATED, 3 UNIFORM-but-close),
and the most honest single statement is: **not cleanly one or the other — but never phase-blind.**
Two convergent, independently-computed views support this:

1. **Ratio view** (model/OrthoLoad AT each hump): ratio₂ > ratio₁ in **4/4** series (SO-knee
   1.46>1.09; CMC-knee 1.59>1.32; SO-hip 1.67>1.31; CMC-hip 1.85>1.57). Only SO-knee's spread
   (28.9%) crosses the pre-registered 25% line into CONCENTRATED; CMC-knee (19.0%), SO-hip (24.3%,
   essentially AT the line), and CMC-hip (16.1%) stay just inside UNIFORM.
2. **Population-spread (Z-score) view** (model value vs that OrthoLoad ensemble's own mean±SD AT
   that hump's %GC — a variance-normalized complement the ratio alone cannot provide): Z₂ > Z₁ in
   **4/4** series too, more starkly (SO-knee 2.67 vs 0.97 — hump1 is NOT a statistically
   distinguishable outlier at all; CMC-knee 3.42 vs 3.30 — both significant, nearly equal; SO-hip
   4.56 vs 1.72 — a large gap even though the RATIO rule called this one "uniform"; CMC-hip 5.80 vs
   3.23).

Read together: **SO (this repo's primary, default solve method, and the one underlying nearly
every prior joint-force cert) shows the cleanest, most actionable localization** — its
weight-acceptance hump is not distinguishable from real OrthoLoad in-vivo variability at all
(Z=0.97, well inside a typical ±2SD band), while its push-off hump is a real, significant
over-prediction (Z=2.67). CMC (the mechanistically-decorrelated forward-dynamics second solve)
and the hip joint both show a more evenly-elevated pattern across both humps — consistent with
`docs/MECHANISM_CMC_SECOND_SOLVE.md`'s own already-disclosed finding that CMC reads 8-11% higher
than SO EVERYWHERE, not at one instant — but even there, the push-off side is consistently further
from the real population than the weight-acceptance side. **A secondary, non-hump finding in the
same direction:** during SWING (70-90%GC, where real joint loading should be near its minimum),
CMC-knee reads 79.3 %BW vs OrthoLoad's 26.4 %BW (3.01×) while SO-knee reads only 34.3 %BW (1.30×)
— a real, absolute (not small-denominator-inflated) gap suggesting CMC's forward-dynamics tracking
controller carries meaningful residual joint loading even when the limb should be nearly
unloaded, consistent with CMC's own previously-disclosed incomplete pelvis-residual gate
(`docs/MECHANISM_RRA_TASK_GAINS.md`) and boosted-reserve-actuator dependency
(`docs/MECHANISM_ANKLE_RESERVE_FIX.md`).

**Robustness (forced, not assumed): smoothing sensitivity.** A 5-point moving average on each
model curve shifts hump values by only 3-9% and hump LOCATIONS by at most 9%GC (e.g. CMC-knee
hump1: raw 295.7%BW@24%GC → smoothed 282.7%BW@16%GC) — the underlying per-frame solver jitter in
early stance makes the EXACT hump1 location somewhat noise-sensitive, but the hump VALUES (and
therefore the ratios/verdicts above) are stable to this perturbation; no series flips its
verdict under smoothing.

### 3.3 A new finding: every prior "peak" cert was measuring the push-off hump

`model_curves_meta` (evidence JSON) shows the already-established overall-peak times (SO: t=0.51s
knee / 0.55s hip; CMC: t=0.51s knee / 0.54s hip) map to **44.8-47.8%GC** in every case — squarely
inside the push-off-hump window, nowhere near the weight-acceptance hump (12-24%GC). This was not
previously known because no prior cert looked at %GC at all. It means every "1.4-2.0×
over-prediction" headline in `docs/MECHANISM_CROSS_SUBJECT.md`/`docs/MECHANISM_CMC_SECOND_SOLVE.md`
was, unknowingly, always a push-off-hump number — this doc is the first to show that (at least for
SO-knee) the SAME model's weight-acceptance hump tells a substantially better-trust story.

## 4. Honest gaps

1. **Hip has no official OrthoLoad ensemble to compare against** (the standardized-loads archive
   is a corrupted download, §2.3/§6) — the hip anchor here is this session's own raw-AKF
   segmentation, method-validated on knee (r=0.928) but not independently checked on hip itself.
   If the true hip curve is sharper (less blurred) than this session's build — as the knee
   calibration's own ~36 %BW RMSE / peak-attenuation pattern suggests is likely — the hip ratios
   above are, if anything, mild UNDER-estimates of the true over-prediction, not over-estimates.
2. **n=1 subject/trial for the model side** (subject2/`walking1` only) — same scope caveat
   inherited from every cert in this family; no claim of generality across subjects/speeds.
3. **The contralateral-symmetry stride-time estimate (T=1.327s) is a measured-but-indirect
   recovery, not a direct heel-strike-to-heel-strike measurement** — the trial's own single
   directly-observed right heel-strike is real, but the cycle-CLOSING boundary rests on the
   standard 50%-GC contralateral-symmetry convention, cross-checked (§2.1) against two
   independent events landing in-band, not independently re-derived a third way.
4. **The 25%-relative-spread uniform/concentrated threshold is a binary cut on what the data show
   is actually a continuum** — SO-hip sits at 24.3%, one percentage point from flipping verdict;
   treated honestly above (§3.2) via the complementary Z-score view rather than hiding behind the
   binary label.
5. **Hump1's exact %GC location is noise-sensitive for the model curves** (§3.2 smoothing check)
   because SO/CMC are independent per-frame solves with no temporal-smoothness constraint —
   disclosed, not smoothed away in the primary numbers (smoothing is reported only as a robustness
   check, not substituted into the headline table).
6. **No per-muscle/per-phase decomposition was done here** — WHICH muscles drive the push-off-hump
   excess (the task's own hypothesis names the gastrocnemius/soleus plantarflexors as a plausible,
   literature-consistent candidate, given push-off's known ankle-plantarflexor-dominated mechanics)
   is not tested in this doc; a sibling in-progress analysis in this same repo
   (`data/msk_smoketest/subject2_walking1/contact_muscle_decomp/`, `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`)
   looks adjacent but was not read in detail or relied upon here (concurrent, independent work by
   another instance per this repo's own multi-instance convention) — flagged as a natural follow-up
   pointer, not a citation of its conclusions.
7. **RMSE is reported in absolute %BW, not normalized** — hip's larger RMSE (91.8-124.5 vs 48.0-73.6
   for knee) partly reflects hip's own higher absolute force magnitude, not a proportionally worse
   fit (Pearson r, which is scale-invariant, is the primary shape-match statistic used for the
   verdict logic; RMSE is reported for completeness, not as the ranking criterion).
8. **The two corrupted-archive/broken-download findings (hip standardized-loads zip, §2.3/§6) were
   not fixed or re-fetched this session** — flagged for a future session with fetch access.

## 5. Verification (machine-checked, not narrated)

- Self-consistency: this script's own re-derived overall-peak values match
  `docs/MECHANISM_CMC_SECOND_SOLVE.md`'s published numbers to <0.01% relative difference on all 4
  series (asserted in-script; a mismatch would raise, not silently pass).
- Phase-alignment: both out-of-sample cross-checks (R toe-off, L toe-off landing in their
  textbook %GC bands) pass (`phase_alignment_verification.all_pass = true`).
- Calibration: raw-AKF segmentation method vs official OrthoLoad knee ensemble, Pearson r=0.928,
  RMSE=36.2 %BW (gate: r>0.8, PASS) — computed BEFORE the hip ensemble (which has no official
  cross-check available) was trusted.
- Hip archive corruption: confirmed via direct byte search for the ZIP end-of-central-directory
  signature (absent in 117MB of otherwise-valid-looking data) — a forensic fact, not a tool-version
  assumption.
- All numbers in this document are pulled directly from
  `data/msk_smoketest/subject2_walking1/contact_waveform/contact_waveform_results.json`, not
  transcribed from console prose.

## 6. Files

- `scripts/msk/contact_waveform_analysis.py` (new) — the full, re-runnable pipeline. Run with
  `/usr/bin/python3 scripts/msk/contact_waveform_analysis.py` (needs numpy+openpyxl; deliberately
  NOT `.venv-msk`, see §2.3 environment note; never imports OpenSim). Reimplements (standalone,
  not imported, since `validate_joint_force.py`/`validate_hip_force.py` both `import opensim` at
  load time) the same `.mot`/`.sto` parser convention and the same pre-registered AKF trial-selection
  filters those modules already use, for internal consistency — verified identical selection counts
  (72 knee "Level Walking" trials, matching `validate_joint_force.py`'s own established n=72; 162
  hip "primary" trials, matching `validate_hip_force.py`'s own established n=162).
- `data/msk_smoketest/subject2_walking1/contact_waveform/contact_waveform_results.json` (new) — every
  number in this document: the model %GC mapping + alignment verification, both OrthoLoad
  ensembles (full mean/SD curves across the 0-100%GC grid, 101 points), the knee raw-AKF
  calibration, and the 4 full comparison records (RMSE, Pearson r, hump1/hump2 values+locations+
  ratios+Z-scores, verdict, swing-region ratio, smoothing-sensitivity values).
- Inputs read in place, never modified: `data/msk_smoketest/subject2_walking1/static_optimization/jr/`,
  `data/msk_smoketest/subject2_walking1/cmc_second_solve/jr_cmc/` (model force time series),
  `data/external/orthoload/knee/standard_loads/K*_Walking.xlsx` (official knee ensemble),
  `data/external/orthoload/{knee,hip_gen1,hip_gen2}/database_api/akf/*.akf` (raw AKF corpus),
  `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/ForceData/walking1_forces.mot`
  (GRF, for heel-strike detection).
- **Not touched**: `scripts/msk/validate_joint_force.py`, `scripts/msk/static_opt_knee.py`,
  `scripts/msk/validate_hip_force.py`, `scripts/msk/cmc_second_solve.py`, and the sibling
  `contact_muscle_decomp` work-in-progress (concurrent, independent instance).

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task).
