# MECHANISM Fmax-FIX CROSS-SUBJECT — does the selective push-off correction generalize to subject3/subject4? (2026-07-21)

Runs the **full stance-phase waveform pipeline** (`scripts/msk/contact_waveform_analysis.py`, reused
unedited) on the **Fmax-corrected** subject3 and subject4 models, mirroring exactly what
`docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md` already did for subject2 (knee push-off hump
390.3→324.0 %BW, Z 2.67→1.23, closing 53.8% of the OrthoLoad gap, WHILE weight-acceptance stayed
flat at Z 0.97→0.74).

**Headline: the mechanism GENERALIZES IN DIRECTION AND SELECTIVITY-PATTERN to both new subjects, but
WEAKENS SUBSTANTIALLY IN MAGNITUDE relative to subject2 — a real, disclosed, partial replication, not
a clean repeat of subject2's result.** Push-off (hump2) drops toward OrthoLoad for **both** subject3
knee (512.1→459.8 %BW, ratio 1.919x→1.723x, **10.2% raw drop, closing 21.3%** of the gap, Z
5.30→4.17) **and** subject4 knee (433.1→393.5 %BW, ratio 1.623x→1.475x, **9.1% raw drop, closing
23.8%** of the gap, Z 3.59→2.74) — robust to 5-point smoothing (10.05%/9.10%) and self-consistent
with each subject's own independently-computed peak-level test (subject3: −10.2% exactly reproduces
`docs/MECHANISM_FMAX_PCSA_VALIDATION.md` Sec.4; subject4: −9.3%, new this session). **Both remain
large, statistically-significant outliers even after correction** (Z=4.17 and Z=2.74), unlike
subject2 which crossed under Z=2 into "no longer a stark outlier" territory — the fix is **~2-2.5x
less sufficient** for subject3/4 than for subject2, even though the muscle-Fmax cut applied to each
subject's own 13 detected knee-crossing muscles is comparable in size (mean ratio 0.815 for
subject3, 0.814 for subject4, an ~18.5% cut for both). Weight-acceptance (hump1) **barely moves for either subject**
(subject3: +0.2% raw/+1.6% smoothed; subject4: +1.2% raw/+1.4% smoothed) — both **smaller than
subject2's own already-marginal 2.0% hump1 drop**, which that doc itself flagged as sitting inside
the pipeline's established noise band. **Answering the task's explicit question — does the
correction help subject4's already-over weight-acceptance, or only push-off — the honest answer is
overwhelmingly "only push-off"**: subject4's weight-acceptance ratio/Z technically move in the
helpful direction (1.266x→1.251x, Z 2.77→2.61) and this small movement IS robust to smoothing (not
pure numerical noise), but it is an order of magnitude smaller than the push-off effect (1.2-1.6%
vs 9.1-9.4%) and leaves weight-acceptance a clear, large, still-significant over-prediction (Z=2.61,
barely changed from 2.77) — **the fix does not meaningfully close subject4's weight-acceptance
excess.** The **hip specificity control replicates cleanly at both new subjects** (all 4 hip humps
move <2.3%, near the noise floor) — confirming the effect stays knee-crossing-specific, not a
generic re-solve artifact, exactly as subject2's own control predicted.

**A forced, load-bearing correction to this task's own premise**: the task stated subject4's
Fmax-corrected model already existed. **Verified false, not assumed true** (Sec.0): only
subject2's and subject3's corrected models existed on disk; `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s
own Sec.4/Sec.7 explicitly disclose subject4 was never causally re-solved. This session built it
(`scripts/msk/fmax_pcsa_correction_test_subject4.py`, new) before the waveform comparison could run
at all.

## Headline table (machine-read from the evidence JSON, not transcribed from console prose)

| subject | joint/hump | baseline (ratio, Z) | corrected (ratio, Z) | raw drop | gap-closed | smoothing-robust? |
|---|---|---:|---:|---:|---:|---|
| subject3 | knee hump1 (WA, alt reading*) | 225.1 %BW (1.097x, Z=0.74) | 224.6 %BW (1.095x, Z=0.72) | +0.2% | 2.3% | yes (+1.6% smoothed) |
| subject3 | **knee hump2 (push-off)** | 512.1 %BW (1.919x, Z=5.30) | 459.8 %BW (1.723x, Z=4.17) | **+10.2%** | **21.3%** | yes (+10.05% smoothed) |
| subject3 | hip hump1 | 426.2 %BW (1.706x, Z=3.97) | 424.7 %BW (1.700x, Z=3.94) | +0.4% | — (control) | — |
| subject3 | hip hump2 | 455.2 %BW (1.960x, Z=6.57) | 445.7 %BW (1.920x, Z=6.29) | +2.1% | — (control) | — |
| subject4 | knee hump1 (WA, over at baseline) | 284.3 %BW (1.266x, Z=2.77) | 280.9 %BW (1.251x, Z=2.61) | +1.2% | 5.7% | yes (+1.35% smoothed) |
| subject4 | **knee hump2 (push-off)** | 433.1 %BW (1.623x, Z=3.59) | 393.5 %BW (1.475x, Z=2.74) | **+9.1%** | **23.8%** | yes (+9.10% smoothed) |
| subject4 | hip hump1 | 449.1 %BW (1.797x, Z=4.49) | 453.5 %BW (1.815x, Z=4.58) | −1.0% | — (control) | — |
| subject4 | hip hump2 | 418.4 %BW (1.802x, Z=5.48) | 408.7 %BW (1.760x, Z=5.20) | +2.3% | — (control) | — |
| *reference* subject2 | *knee hump1 (WA)* | *245.4 %BW (1.093x, Z=0.97)* | *240.6 %BW (1.071x, Z=0.74)* | *+2.0%* | *23.1%* | *(established)* |
| *reference* subject2 | *knee hump2 (push-off)* | *390.3 %BW (1.462x, Z=2.67)* | *324.0 %BW (1.214x, Z=1.23)* | *+17.0%* | *53.8%* | *(established)* |

\* subject3's official fixed-30%GC-split hump1 is a boundary-clipped ramp value (no real dip either
side, on BOTH the baseline and corrected curve, re-checked fresh in this run — Sec.3.3); the
genuine-local-peak alternate reading (22%GC) is used for both sides consistently, exactly as
`docs/MECHANISM_WAVEFORM_CROSS_SUBJECT.md` Sec.3.3 already established for the baseline side alone.

## 1. Pre-registration (stated before any subject3/4 corrected-model hump number was computed)

**C_generalizes** (restates the task's own literal wording, not a looser subject2-style bar — that
bar presupposes weight-acceptance already matched at baseline, which is false for subject4): (i)
push-off ratio-vs-OrthoLoad moves DOWN for BOTH subject3-knee AND subject4-knee, AND (ii) for
subject4 specifically, the correction also measurably helps weight-acceptance
(ratio1_corrected < ratio1_baseline).

**NOT-C (subject2-specific)**: hump2 fails to move toward 1.0 at subject3-knee OR subject4-knee, OR
subject4's weight-acceptance excess does NOT shrink.

**Falsifier**: any ONE of (ratio2_corrected ≥ ratio2_baseline at subject3-knee), (same at
subject4-knee), (ratio1_corrected ≥ ratio1_baseline at subject4-knee) would falsify C_generalizes.
**None fired** — see Sec.4 for the precise, magnitude-qualified reading of what "generalizes" means
here.

**Secondary, non-decisive check** (kept as a supplement, not a substitute, per this repo's own
"fixed-threshold-on-a-noisy-statistic-is-regime-blind" convention): the Z-score direction/magnitude,
and the STRICTER subject2-calibrated H_selective bar (Z1_corrected<1.0), reported for full
transparency even though it cannot be a fair bar for subject4 (Sec.4.3).

**Hip specificity control** (same structural prediction as subject2/subject3's peak-level test): near
-zero net change at both hip humps in both subjects — corroborates the knee effect is knee-crossing
-specific.

**Convergence/self-consistency symmetric-QC clause** (task's own explicit ask): each corrected
model's SO/JR run must reproduce its OWN already-published overall peak to <0.5% before any
hump-level number is trusted, and must show no degenerate solve.

## 2. Section 0 — a forced correction to the task's own premise (verified live, not assumed)

The task stated subject3 **and subject4's** Fmax-corrected models already existed. `find` across
`data/msk_models/` and `data/msk_smoketest/` (run before writing any new code) showed **only**
`subject2_scaled_handsfield_fmax_corrected.osim` and `subject3_scaled_handsfield_fmax_corrected.osim`
— **zero** subject4 artifacts. Cross-checked against the documentary record:
`docs/MECHANISM_FMAX_PCSA_VALIDATION.md` Sec.4's own decisive-causal-test table lists only subject2
and subject3 rows; its Sec.7 honest-gaps item 5 explicitly states "Only the
specific-tension-corrected knee/hip re-solve was tested for subject2 and subject3 (**not
subject4**...)". Subject4's Sec.3 table row (mean ratio@60=1.454) is a **different, non-causal**
computation (the raw per-muscle PCSA ratio, not a re-solved model).

**Fix**: `scripts/msk/fmax_pcsa_correction_test_subject4.py` (new, this session) built it, reusing
`fmax_pcsa_correction_test.build_corrected_model` **imported unedited** (not re-derived, not
copy-pasted) — the identical PCSA-Fmax-rewrite function subject2/subject3's own corrected models were
built with, differing only in which subject's model/height/mass are passed in. Output:
`data/msk_models/subject4_scaled_handsfield_fmax_corrected.osim` (80 muscles, Fmax_new/Fmax_old ratio
min=0.521 max=1.235 mean=0.816) + a fresh SO+JR re-solve
(`data/msk_smoketest/subject4_walking1/fmax_pcsa_correction_test/{so,jr}/`). Peak-level result (new):
knee_r ratio 1.681x→1.525x (**−9.3%**, just under this repo's own historically-used ≥10% "real
contributor" bar — see Sec.4.4 for why this near-miss should not be over-read); hip_r ratio
1.529x/1.816x (self/JR) → 1.553x/1.834x (**+1.6%/+1.0%, ROSE slightly** — consistent with the
hip-crossing-muscles-are-a-genuine-mix pattern already established for subject2/subject3).
Convergence gates non-degenerate (`muscles_over_fmax_ratio_1p5={}`,
`joint_reserve_leaning_hard={}`, `n_frames=131` exactly matching the baseline run's own frame count).

## 3. Method (reuse discipline, machine-checked)

`scripts/msk/contact_waveform_fmax_corrected_cross_subject.py` (new) imports
`contact_waveform_analysis` (cwa) **and** `contact_waveform_cross_subject` (cwcs) as modules, **zero
edits to either** — md5-verified before/after this run (`cwa`: `9d54780f5146b8fea3821b3f92bac18d`,
matching the value already on record in `docs/MECHANISM_WAVEFORM_CROSS_SUBJECT.md`; `cwcs`:
`b80ed220e0acce5513852ec372c651d8` — both unchanged, asserted in-script, would raise otherwise).

**The def-time-default-binding trap (task's own explicit warning), avoided by construction, not by
discipline alone**: this script never calls `cwa.model_curve_on_grid` (the function whose hardcoded
internal `joint_force_mag_pct_bw(reaction_sto, colprefix)` call binds `bw_n=BW_N` once at
def-time — monkey-patching `cwa.BW_N` afterward is inert, already diagnosed via a standalone toy
closure in `docs/MECHANISM_WAVEFORM_CROSS_SUBJECT.md` Sec.2 Adversary #1). It exclusively uses
`cwcs.build_model_curve(reaction_sto, colprefix, to_pct_gc, bw_n=...)` (imported unedited from the
already-fixed script), threading each subject's own body weight explicitly: **subject3=622.7222750000001N,
subject4=613.89629N** (never subject2's 766.88N).

**Phase mapping reused, not recomputed from scratch, and cross-checked against the published
baseline**: since the Fmax correction changes neither the GRF file nor the model's IK/kinematics, the
%GC mapping is, by construction, identical between baseline and corrected runs. This script
re-derives it fresh (via the same try-strict-then-robust-fallback strategy `contact_waveform_cross_
subject.py` already established: subject3 falls back to `cwcs.robust_model_pct_gc_mapping`, subject4
uses `cwa.model_pct_gc_mapping` directly) and asserts the result matches the already-published
`contact_waveform_cross_subject_results.json` mapping_diag to **exactly 0.0 absolute difference**
(both subjects) — confirming no drift, not assumed.

**Self-consistency gate** (forced before trusting anything else, task's own explicit requirement):
each subject's fresh re-extraction of the corrected model's overall knee_r/hip_r peak reproduces the
already-published `fmax_pcsa_correction_test*_results.json` `jointreaction_peak_pct_bw` values to
**0.0000% relative difference** (all 4: subject3 knee/hip, subject4 knee/hip) — via a genuinely
decorrelated code path (raw `.sto` vector-magnitude parse + %GC grid interpolation) vs the
correction-test scripts' own biomechanical force-decomposition + `opensim.JointReaction` route, not a
tautological re-read of the same number.

**Convergence gates re-surfaced and asserted, not silently trusted**: both corrected runs show
`muscles_over_fmax_ratio_1p5={}` and `joint_reserve_leaning_hard={}` (no degenerate solve); raw
`convergence_pass=False` for both, but this is the **already-diagnosed subject2-specific hardcoded
`t≥1.50s` literal** (`docs/MECHANISM_CROSS_SUBJECT.md` §4) that fails BY CONSTRUCTION for every
subject but subject2 (whose trial is the longest of the ten) — not a real non-convergence; both
subjects' `n_frames` (132, 131) exactly match their own already-established baseline frame counts.

## 4. Results in detail

### 4.1 Push-off generalizes in direction for both subjects — but 2-2.5x weaker than subject2

Both subject3 (10.2%) and subject4 (9.1%) show a real, smoothing-robust push-off drop of the same
order of magnitude as each other, but roughly **half** subject2's 17.0% raw drop. In gap-closed
terms (fraction of the ratio's excess-over-1.0 removed): subject3 closes 21.3%, subject4 closes
23.8%, vs subject2's 53.8% — **subject2's fix was 2.3-2.5x more sufficient**. Both subject3/4 remain
strong statistical outliers post-correction (Z=4.17, Z=2.74) while subject2 crossed under the Z=2
"no longer a stark outlier" line (Z=1.23). This is a real, quantified, generalizing-but-weaker
finding, not a clean replication.

**Consistency check (over-determination, not circular)**: subject3's waveform-level push-off drop
(10.2%) reproduces its own already-published peak-level causal test (`docs/
MECHANISM_FMAX_PCSA_VALIDATION.md` Sec.4: −10.2%) to the same precision, computed by a **different**
script/code path (waveform %GC-hump extraction vs the peak-level self-computed+JR cross-check) —
genuine agreement, not a repeated read of the same stored number. Subject4's waveform-level drop
(9.1%) is consistent with (within 0.2 points of) its own newly-computed peak-level drop (9.3%,
Sec.2) via the same cross-check structure.

### 4.2 Weight-acceptance: technically helped for subject4, but far too small to matter practically

Subject4's hump1 ratio/Z do move in the helpful direction (1.266x→1.251x, Z 2.77→2.61) and this small
movement **survives a 5-point smoothing perturbation** (1.2% raw → 1.35% smoothed — same sign, same
order of magnitude, not an artifact of frame-to-frame jitter). But it is **smaller than subject2's
own hump1 drop** (2.0%), which `docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md` Sec.5 item 3 itself
already flagged as "sitting well within this pipeline's own established smoothing-sensitivity band
... on its own would not be distinguishable from noise." By that **same, already-established**
standard, applied symmetrically here: subject4's weight-acceptance "help" is real in direction but
**practically negligible** — Z=2.61 is barely different from Z=2.77, still a clear, large,
significant over-prediction. **Direct answer to the task's question**: the correction helps almost
exclusively push-off; weight-acceptance is not meaningfully rescued.

Subject3's own hump1 (using the alt genuine-peak reading) shows an even smaller raw movement (+0.2%,
+1.6% smoothed) — consistent with subject2's pattern of an already-good weight-acceptance match
staying flat, not degrading (no "over-correction" signature: neither subject's hump1 ratio fell
below 1.0 or flipped Z negative).

### 4.3 The stricter, subject2-calibrated H_selective bar does not transfer cleanly to subject4 — by construction, not failure

Applying `docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md`'s own strict bar (Z1_corrected < 1.0 AND
ratio1 does not worsen AND hump2 moves toward 1) for full transparency: **subject3-knee clears it**
(Z1=0.72<1.0, ratio1 improved, hump2 moved toward 1 → H_SELECTIVE_CORROBORATED, matching subject2).
**Subject4-knee does not** (Z1=2.61, nowhere near <1.0) — but this bar implicitly assumes
weight-acceptance was ALREADY a statistical match at baseline (true for subject2/subject3, Z=0.97/
0.74-1.22, false for subject4, Z=2.77 at baseline). Subject4 failing a bar that was never
satisfiable for it given its own baseline is **not informative** about whether the mechanism
generalizes — this is why the task's own, weaker, direction-based falsifier (Sec.1) is the fair
test here, and why this doc reports both readings rather than picking the flattering one.

### 4.4 A near-miss against this repo's own historical 10% "real contributor" bar — a threshold artifact, not a sign the effect is unreal

Subject4's push-off relative-ratio-drop (9.1-9.3%, peak-level and waveform-level agree) sits just
under the ≥10% bar `docs/MECHANISM_FMAX_PCSA_VALIDATION.md` Sec.4 pre-registered and subject2 (17.1%)
/ subject3 (10.2%) both cleared. Per this repo's own "fixed-threshold-on-a-noisy-statistic-is-
regime-blind" convention: the smoothing-robustness check (9.10% smoothed vs 9.14% raw — a <0.1
percentage-point shift) shows the underlying effect size itself is stable to roughly ±1 percentage
point, so a rigid line at exactly 10% should not be read as "subject4 fails, subject3 passes" — both
are the **same real, replicated, order-of-magnitude effect** (9-10%), and subject4 landing a hair
under an arbitrary round-number cutoff is a threshold artifact, not evidence the mechanism is weaker
in kind for subject4 specifically (it IS weaker than subject2's 17.0%, which is the genuinely
informative comparison).

### 4.5 Hip specificity control replicates cleanly at both subjects

All 4 new hip-hump deltas are small (+0.4%, +2.1% for subject3; −1.0%, +2.3% for subject4), matching
`docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s already-established near-zero-net-effect pattern for hip
(mixed-direction per-muscle corrections cancel out). **A metric artifact flagged, not hidden**: the
generic `change_is_selective_not_uniform_25pct_rule` boolean reads `True` for all 4 hip series too
(mean_rel_drop is tiny — 1.22%/0.67% — so the spread/mean ratio blows up on a near-zero denominator,
the exact same degeneracy `docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md` Sec.4.2 already diagnosed for
subject2's hip). The substantively correct reading is the directly-computed
`hip_specificity_control.both_humps_within_5pct=True` for both subjects — no material change at
either hip hump, corroborating knee-crossing-specificity, not a second "selective" finding.

### 4.6 Geometric context (not re-derived here, cited from the already-verified source)

`docs/MECHANISM_FMAX_PCSA_VALIDATION.md` Sec.3 already established the geometric governor: Fmax is
frozen (never rescaled per-subject) while Handsfield-implied PCSA scales with each subject's own
height×mass, so the model/PCSA ratio is monotonic in H×M — subject3/4 (H×M≈107-105) sit further out
on that curve than subject2 (H×M=153.3). Measured directly here: the actual Fmax cut applied to each
subject's own 13 detected knee-crossing muscles is comparable in magnitude between subject3 and
subject4 (mean ratio 0.815 / 0.814 respectively, an ~18.5% cut for both — consistent with their near-identical anthropometrics,
1.69m/63.5kg vs 1.68m/62.6kg) — so the similar push-off gap-closed fractions (21.3%/23.8%) are not a
coincidence, they reflect a genuinely similar correction magnitude at a genuinely similar starting
mismatch. This doc does not re-derive why subject2's fractional gap-closed is so much larger; that
remains an open question (Sec.5).

## 5. Honest gaps

1. **The mechanism-level activation-redistribution check** (`docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md`
   Sec.4.4: gastrocnemius/rectus-femoris activation drops, soleus rises, at push-off) was **not
   re-run for subject3/subject4** in this doc — the already-committed `StaticOptimization_activation.sto`
   files exist for both (produced by the correction-test runs) and could support the identical check
   with zero new OpenSim runs, but this was scoped out to stay lean given the task's primary ask
   (hump/Z comparison) was already decisively answered without it. A real, disclosed scope limit, not
   an oversight.
2. **Why subject2's fractional gap-closed (53.8%) so far exceeds subject3/4's (21-24%) despite a
   comparable-or-larger knee-crossing Fmax cut** (Sec.4.6) is not resolved here — a candidate
   explanation (some other over-prediction driver, e.g. no-RRA/no-EMG-co-contraction, scales up
   faster than the Fmax-mismatch does as H×M shrinks) is plausible but untested.
3. **n=3 subjects total** (subject2 fully, subject3/4 this doc) — the same small, descriptive
   instance-space every cert in this family discloses. Sufficient to show the mechanism's DIRECTION
   generalizes and to quantify that its MAGNITUDE does not transfer 1:1; insufficient to fit a
   population model of the magnitude's own scaling law.
4. **Subject3's official hump1 boundary-clip artifact** (Sec.3.3, inherited from
   `docs/MECHANISM_WAVEFORM_CROSS_SUBJECT.md`) required the alt-reading workaround on BOTH baseline
   and corrected curves — re-verified fresh here (the corrected curve independently shows the same
   no-real-dip signature, not assumed to carry over), but this remains a per-subject hump-finder
   fragility this repo's `find_two_humps` was only ever validated against subject2's 6 curves.
5. **Hip's OrthoLoad ensemble is self-built** (raw-AKF, calibrated vs official knee only, r=0.928) —
   same inherited limitation as every prior doc in this family; affects only the specificity-control
   numbers, not the primary knee verdict.
6. **The R-toe-off textbook-band phase-alignment check fails for both subjects** (subject3: 70.1%,
   subject4: 69.6%, vs the 55-68% band) — inherited, unchanged from baseline (the mapping itself is
   identical between baseline/corrected by construction, confirmed via the 0.0 mapping-diff gate,
   Sec.3), already independently cross-validated sound via two decorrelated anchors in
   `docs/MECHANISM_WAVEFORM_CROSS_SUBJECT.md` Sec.3.2 — not re-litigated here.
7. **SO only, no CMC** for subject3/4 (no second-solve output exists for either) — same scope caveat
   `docs/MECHANISM_WAVEFORM_CROSS_SUBJECT.md` already discloses.
8. **Subject4's peak-level correction ratio-drop (9.3%) sits just under this repo's own historical
   10% "real contributor" bar** (Sec.4.4) — disclosed with the smoothing-robustness evidence that
   this is a threshold artifact, not treated as a clean pass or a clean fail.

## 6. Verification (machine-checked, not narrated)

- Section-0 premise check: `find` across `data/msk_models/` + `data/msk_smoketest/` for any
  subject4-fmax-corrected artifact returned zero hits before this session's own build — confirmed
  the task's premise was false for subject4, not assumed true.
- Reused-module integrity: `contact_waveform_analysis.py` md5 `9d54780f5146b8fea3821b3f92bac18d`
  (matches the value already on record in `docs/MECHANISM_WAVEFORM_CROSS_SUBJECT.md`) and
  `contact_waveform_cross_subject.py` md5 `b80ed220e0acce5513852ec372c651d8`, both asserted
  identical before/after this run (script would raise otherwise).
- Phase-mapping-vs-published-baseline cross-check: max absolute difference over
  `R_heel_strike_t`/`L_heel_strike_t`/`T_stride_estimated_s` = exactly `0.0` for both subjects
  (asserted in-script).
- Self-consistency: fresh-extracted corrected-model overall peaks match each subject's own
  already-published `jointreaction_peak_pct_bw` to **0.0000%** relative difference, all 4
  (subject3/4 × knee/hip) (asserted in-script, `rel_diff < 0.5%` gate, would raise otherwise).
- Convergence: `muscles_over_fmax_ratio_1p5={}` and `joint_reserve_leaning_hard={}` for both
  corrected runs (asserted, non-degenerate solve); `n_frames` matches each subject's own established
  baseline frame count exactly (132, 131).
- Hump-shape validity (`has_real_dip_between`): re-run fresh on the CORRECTED curves, not assumed
  to inherit the baseline curve's own diagnosis — subject3-knee shows the same boundary-artifact on
  both baseline and corrected curves (consistent, not cherry-picked).
- Smoothing-robustness: 5-point moving-average recomputation of all 4 knee hump values, both
  subjects, both baseline and corrected (all four push-off drops and both weight-acceptance drops
  reported both raw and smoothed, Sec.4.1-4.2).
- OrthoLoad knee raw-AKF-method calibration vs official ensemble: r=0.9278 (gate r>0.8, PASS),
  recomputed fresh in this run.
- All numbers in this document are pulled directly from
  `data/msk_smoketest/contact_waveform_fmax_corrected_cross_subject/contact_waveform_fmax_corrected_cross_subject_results.json`
  and `data/msk_smoketest/subject4_walking1/fmax_pcsa_correction_test/fmax_pcsa_correction_test_subject4_results.json`,
  never transcribed from console prose.

## 7. Confidence tier and falsifier (restated per report discipline)

**Confidence tier: in-vivo-anchored (OrthoLoad waveform target/verdict) × cadaveric/published-
plausibility (Handsfield-MRI-volume-via-Arnold-2010-specific-tension for the Fmax correction
itself)** — identical tier composition to `docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md`, extended to
n=3 subjects (still a small instance-space, Honest gap 3).

**Claim**: the Handsfield-MRI-PCSA Fmax correction's selective, push-off-concentrated repair
mechanism **generalizes in direction and pattern** to subject3 and subject4 (push-off ratio/Z move
toward OrthoLoad in both, robust to smoothing, self-consistent with independent peak-level
computations; hip specificity control replicates near-zero at both) but **generalizes weaker in
magnitude** than subject2 (21-24% vs 53.8% of the ratio-gap closed; both subjects remain
statistically significant push-off outliers post-correction) — and for subject4 specifically, the
correction does **not** meaningfully rescue the already-over weight-acceptance hump (a real but
practically-negligible 1.2-1.6% movement, smaller than subject2's own already-marginal hump1 change).

**Falsifier (what would have overturned "generalizes," none of which occurred)**: push-off
ratio-vs-OrthoLoad failing to move toward 1.0 at subject3-knee or subject4-knee; subject4's
weight-acceptance ratio failing to move toward 1.0 at all (it did, marginally); the hip specificity
control showing a comparably large shift at either subject (it did not, both <2.3%).

## 8. Files

- `scripts/msk/fmax_pcsa_correction_test_subject4.py` (new) — builds subject4's Fmax-corrected model
  (reusing `fmax_pcsa_correction_test.build_corrected_model` unedited) + fresh SO+JR re-solve +
  peak-level comparison vs subject4's cached baseline. Run with
  `source .venv-msk/bin/activate && python3 scripts/msk/fmax_pcsa_correction_test_subject4.py`.
- `scripts/msk/contact_waveform_fmax_corrected_cross_subject.py` (new) — the full waveform pipeline
  described above. Run with `/usr/bin/python3 scripts/msk/contact_waveform_fmax_corrected_cross_subject.py`
  (numpy+openpyxl, never imports OpenSim, NOT `.venv-msk`).
- `data/msk_models/subject4_scaled_handsfield_fmax_corrected.osim` (new) — subject4's Fmax-corrected
  model (80 muscles, geometry/Lopt/Lts/pennation/paths untouched, Fmax only rewritten).
- `data/msk_smoketest/subject4_walking1/fmax_pcsa_correction_test/{so,jr}/` (new) — subject4
  corrected-model SO+JR outputs.
- `data/msk_smoketest/subject4_walking1/fmax_pcsa_correction_test/fmax_pcsa_correction_test_subject4_results.json`
  (new) — subject4 peak-level correction-test full detail (Sec.2).
- `data/msk_smoketest/contact_waveform_fmax_corrected_cross_subject/contact_waveform_fmax_corrected_cross_subject_results.json`
  (new) — **the evidence JSON**: every number in this document (phase-mapping cross-checks,
  self-consistency gates, convergence gates, both subjects' corrected curves, before/after delta
  tables with gap-closed fractions, hump-shape diagnostics, smoothing-sensitivity values, hip
  specificity control, primary pre-registered verdict).
- Inputs read in place, never modified: `data/msk_smoketest/subject3_walking1/fmax_pcsa_correction_test/jr/`
  (already existed, built by `scripts/msk/fmax_pcsa_correction_test.py`), each subject's baseline
  `static_optimization/jr/` outputs, `data/msk_smoketest/contact_waveform_cross_subject/contact_waveform_cross_subject_results.json`
  (baseline waveform numbers for the delta table), `data/msk_smoketest/cross_subject_validation/subject{3,4}_walking1_cross_subject_results.json`
  (uncorrected-Fmax peak-level baselines), `data/external/orthoload/{knee,hip_gen1,hip_gen2}/...`
  (re-parsed fresh), each subject's GRF `.mot` (read-only external drive).
- **Not touched** (md5-verified, Sec.6): `scripts/msk/contact_waveform_analysis.py`,
  `scripts/msk/contact_waveform_cross_subject.py`, `scripts/msk/contact_waveform_fmax_corrected.py`,
  `scripts/msk/fmax_pcsa_correction_test.py`, `scripts/msk/fmax_pcsa_correction_test_subject2.py`,
  `scripts/msk/fmax_pcsa_validation.py`, `scripts/msk/subject_specific_scaling.py`,
  `scripts/msk/cross_subject_validation.py`, `scripts/msk/static_opt_knee.py`,
  `scripts/msk/validate_hip_force.py`, `scripts/msk/validate_joint_force.py`.
- Prior docs this builds on (context, not re-litigated): `docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md`,
  `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`, `docs/MECHANISM_WAVEFORM_CROSS_SUBJECT.md`,
  `docs/MECHANISM_CONTACT_WAVEFORM.md`, `docs/MECHANISM_CROSS_SUBJECT.md`.

Isolation respected throughout: bodytwin only; all inputs read in place (LabValidation external
drive, read-only); no OpenSim import in the waveform script (`.sto`/`.json`/`.xlsx`/`.akf` re-parsing
only); the model-building script ran under `.venv-msk` per this repo's own convention. No git commit,
no git add, no git push performed.
