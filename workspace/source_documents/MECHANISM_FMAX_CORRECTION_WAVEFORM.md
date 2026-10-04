# MECHANISM Fmax-CORRECTION WAVEFORM — does the MRI-PCSA fix SELECTIVELY repair push-off, or over-correct? (2026-07-21)

The sharpest available test of whether the Handsfield-MRI-PCSA Fmax correction
(`docs/MECHANISM_FMAX_PCSA_VALIDATION.md`) is a genuine, mechanistically-real fix or a coincidental/
over-tuned rescaling: run the **same full stance-phase waveform pipeline**
(`docs/MECHANISM_CONTACT_WAVEFORM.md`, `scripts/msk/contact_waveform_analysis.py`, reused
**unchanged**) on the **Fmax-corrected model** (`data/msk_models/subject2_scaled_handsfield_fmax_corrected.osim`)
and check BOTH humps, not just the peak.

**Headline: SELECTIVE, pre-registered falsifier NOT triggered.** The knee's push-off hump (hump2)
drops from 390.3 to 324.0 %BW — ratio-vs-OrthoLoad falls from **1.462x to 1.214x** (Z-score
2.67→1.23, no longer a stark statistical outlier) — a **17.0% drop**, closing **53.8%** of the gap
to a perfect ratio of 1.0. The weight-acceptance hump (hump1) barely moves: 245.4→240.6 %BW, ratio
**1.093x→1.071x** (Z 0.97→0.74) — a **2.0% drop**, staying comfortably inside "matches OrthoLoad"
territory and, if anything, getting marginally *better*, not worse. The relative-drop spread between
the two humps (15.0 percentage points around a 9.5% mean) is **8.7x larger at push-off than at
weight-acceptance** — far past this repo's own pre-registered 25%-relative-spread
selective/uniform threshold. **None of the three pre-registered falsifier conditions fired**
(hump1 Z did not reach 1.0; hump1's ratio did not fail to improve; hump2's ratio did not fail to
move toward 1.0). A parallel **specificity control at the hip** (whose Fmax correction was already
known to net to ~zero effect, `docs/MECHANISM_FMAX_PCSA_VALIDATION.md` Sec.4) shows both hip humps
changing by <1% — noise-floor, not a real shift — confirming the knee effect above is caused
specifically by the knee-crossing correction, not a generic "any re-solve moves everything"
artifact. **This is a genuine, phase-selective, mechanistically-explained correction** (Sec.4.4:
Static Optimization's minimum-cost recruitment measurably shifts load, at push-off only, from the
Fmax-cut biarticular knee-crossers onto a non-knee-crossing synergist) — **not** an over-correction
and **not** evidence the weight-acceptance match was luck. The correction does **not** fully close
the push-off gap (324.0 %BW is still 21.4% above OrthoLoad's 266.9 %BW) — a real, quantified,
partial fix, exactly as `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s own peak-level test already found
(1.255x, 25.5% over) — this waveform-level, independently-anchored re-derivation **corroborates**
that number (53.8% vs that test's own 50.4% gap-closed, converging within 3.4 points despite using
a completely different OrthoLoad-anchor convention, Sec.4.3) rather than contradicting it.

## Headline table

| joint | hump | %GC | baseline model (ratio, Z) | corrected model (ratio, Z) | Δ value | Δ ratio | Δ Z |
|---|---|---:|---:|---:|---:|---:|---:|
| **knee** | **hump1 (weight-accept.)** | 15% | 245.4 %BW (**1.093x**, Z=**0.97**) | 240.6 %BW (**1.071x**, Z=**0.74**) | −2.0% | −0.021 | −0.22 |
| **knee** | **hump2 (push-off)** | 44-45% | 390.3 %BW (**1.462x**, Z=**2.67**) | 324.0 %BW (**1.214x**, Z=**1.23**) | **−17.0%** | **−0.249** | **−1.44** |
| hip (control) | hump1 | 12-24% | 326.3 %BW (1.306x, Z=1.72) | 323.5 %BW (1.295x, Z=1.66) | −0.8% | −0.011 | −0.06 |
| hip (control) | hump2 | 41-48% | 386.9 %BW (1.666x, Z=4.56) | 388.2 %BW (1.672x, Z=4.60) | +0.3% | +0.005 | +0.04 |

OrthoLoad anchors (unchanged from `docs/MECHANISM_CONTACT_WAVEFORM.md`, re-verified fresh in this
run, not read from cache): knee hump1 224.6±21.6 %BW @15%GC, hump2 266.9±46.2 %BW @48%GC (official
`Dat_BW` ensemble, 8 subjects); hip hump1 249.9±44.4 %BW @12%GC, hump2 232.2±34.0 %BW @41%GC
(raw-AKF ensemble, 757 cycles, calibrated vs official knee at r=0.928).

**Pre-registered verdict: `H_SELECTIVE_CORROBORATED`** (knee, primary target). Falsifier
conditions checked and NOT triggered: Z1_corrected (0.744) < 1.0 ✓; ratio1_corrected (1.071) <
ratio1_baseline (1.093), i.e. did not fail to improve ✓; ratio2_corrected (1.214) < ratio2_baseline
(1.462), i.e. moved toward 1.0 ✓.

## 1. Pre-registration (stated before this script was run)

**H_selective:** hump2 (push-off) ratio-vs-OrthoLoad drops toward 1.0 relative to baseline AND
hump1 (weight-acceptance) Z-score stays < 1.0 (stays statistically indistinguishable from
OrthoLoad, as baseline already was at Z=0.97) ⇒ strong corroboration that push-off
plantarflexor/knee-extensor over-strength is the real, correctable mechanism.

**H_uniform:** hump1 is pushed down far enough that the corrected model falls BELOW the OrthoLoad
mean at that phase (ratio1 < 1.0, Z1 flips negative) by a drop comparable in relative size to
hump2's ⇒ the correction is not phase-specific; the baseline weight-acceptance "match" was partly
luck, or the correction is over-tuned.

**Falsifier for H_selective** (any one of these would have falsified it): Z1_corrected ≥ 1.0; OR
ratio1_corrected not lower than ratio1_baseline; OR hump2 fails to move toward 1.0 at all
(ratio2_corrected ≥ ratio2_baseline). **None fired** (Headline table).

**Decision statistic** (pre-registered, not chosen post-hoc): this repo's own established
25%-relative-spread convention (`docs/MECHANISM_CMC_SECOND_SOLVE.md`, reused verbatim by
`docs/MECHANISM_CONTACT_WAVEFORM.md`), applied here to the BEFORE/AFTER *change* itself
(`|Δratio1| − |Δratio2|` relative to their mean), not to the model/OrthoLoad snapshot ratio the
baseline doc's rule checks. Measured: knee spread/mean = 1.586 (159%) — far past the 0.25 bar.

## 2. Symmetric-QC — this is a prediction, not a fit

**The adversary:** "of course the corrected model looks better — you tuned Fmax to make OrthoLoad
agreement improve." **Forced, not asserted:** `grep -n -i orthoload` across
`scripts/msk/fmax_pcsa_correction_test.py`, `fmax_pcsa_correction_test_subject2.py`,
`fmax_pcsa_validation.py` shows every "orthoload" hit lives strictly inside `main()`'s
reporting/printout code, **after** `build_corrected_model()` (the actual Fmax-*setting* function)
has already run and returned. `build_corrected_model()` itself takes only `(model_file_in,
out_path, height_m, mass_kg)` and computes `Fmax_new = PCSA(Handsfield, height, mass) x
60 N/cm^2` — it never reads OrthoLoad data in any form. The correction was fixed entirely by an
independent MRI-volume anchor (Handsfield et al. 2014) and each subject's own anthropometrics;
OrthoLoad is consulted only afterward, to report a ratio, never to choose a value. **Whatever this
waveform comparison finds is therefore a prediction from an independent anchor, not a fit to the
OrthoLoad waveform target** — the same framing `docs/MECHANISM_FMAX_PCSA_VALIDATION.md` established
at the peak level, now re-verified structurally for the code path this task actually exercises.

## 3. Method

**Reuse, not reimplementation.** `scripts/msk/contact_waveform_fmax_corrected.py` (new) imports
`scripts/msk/contact_waveform_analysis.py` **as a module, zero edits** and calls its functions
directly: `model_pct_gc_mapping`/`verify_phase_alignment` (GRF-derived %GC mapping — model-
independent, so identical for baseline and corrected by construction), `model_curve_on_grid`
(generic `.sto`→%GC-curve extractor, pointed here at the corrected model's own JointReaction
output), `knee_official_ensemble`/`select_knee_level_walking_trials`/`select_hip_primary_trials`/
`raw_akf_ensemble` (OrthoLoad anchors, recomputed **fresh** in this run rather than read from the
cached baseline JSON — deliberately, to avoid the stale-cache trap
`scripts/msk/subject_specific_scaling.py`'s own `compute_baseline_fresh()` docstring already warns
about in this exact codebase), and `find_two_humps`/`compare` (hump-finder + Z-score/ratio/verdict
logic — byte-identical thresholds: `GRID_N=101`, 30%GC/62%GC hump split, 0.25 uniform/concentrated
band). This satisfies the task's own requirement to verify the corrected-model waveform uses the
**same** phase alignment and hump-finder as the baseline, by construction (same code object, not a
re-implementation that could silently drift).

**Controlled-experiment verification (checked, not assumed).** `diff` on the two JointReaction
setup XMLs (`data/msk_smoketest/subject2_walking1/{static_optimization,fmax_pcsa_correction_test}/jr/walking1_setup_jr_patched.xml`)
shows they differ **only** in `model_file` and `results_directory` — GRF, external loads, and time
window are identical. Direct comparison of both `.sto` time vectors: **158 frames, t=[0, 1.57s],
max abs difference 0.0** (byte-identical); column layout identical (same 28 columns, same order).
Body weight (766.88 N) and total mass (78.2 kg) identical. **Fmax is the only thing that changed
between these two SO+JR runs** — this is a real controlled experiment, not two independently-built
models being compared.

**Self-consistency gate** (forced before trusting anything else, same discipline
`contact_waveform_analysis.py` itself applies to the baseline): this script's own fresh
re-extraction of the corrected model's overall peak (knee 324.1093 %BW, hip 388.3556 %BW) matches
the already-published `fmax_pcsa_correction_test_subject2_results.json` `jointreaction_peak_pct_bw`
values to **0.0000% relative difference** (exact) on both joints — confirms no path/column/scaling
bug before anything downstream is trusted.

**Convergence / IK sanity at push-off** (task's explicit ask): the corrected model's own SO run
(`fmax_pcsa_correction_test_subject2.py`, already executed) reports `convergence_pass=true`,
`n_frames=158` (matches baseline), `muscles_over_fmax_ratio_1p5={}` (no muscle pinned at its
ceiling), `joint_reserve_leaning_hard={}` (no degenerate reserve-actuator dependency) — a real,
non-degenerate re-solve, re-surfaced and asserted in this script, not silently trusted.
Phase-alignment re-check (`verify_phase_alignment`) passes identically to baseline (both
independent toe-off checks in-band) because it depends only on the GRF file, which is unchanged.

## 4. Results in detail

### 4.1 Knee (primary target) — SELECTIVE, not uniform

Push-off (hump2) drops **17.0%** in absolute %BW (390.3→324.0), while weight-acceptance (hump1)
drops only **2.0%** (245.4→240.6) — an 8.7x differential. In ratio terms: hump2's ratio-vs-OrthoLoad
falls from 1.462x to 1.214x (closing **53.8%** of the gap to a perfect match), while hump1's ratio
falls only from 1.093x to 1.071x (closing 23.1% of its own, much smaller, gap). In Z-score terms
(model value vs the OrthoLoad ensemble's own inter-subject spread at that exact %GC): hump2's Z
falls from 2.67 (a real, significant outlier) to 1.23 (comfortably inside a typical ±2SD band, no
longer a stark statistical outlier); hump1's Z falls from 0.97 to 0.74 — **improves slightly**, not
worsens. The pre-registered decision statistic (relative-drop spread / mean) measures 1.586 (159%),
mechanically confirming visual "selective" reading against a fixed, pre-stated numeric bar (0.25),
not an eyeballed judgment.

### 4.2 Hip (specificity control) — near-zero at both humps, as predicted

Both hip humps move by <1% (hump1 −0.8%, hump2 +0.3%) — noise-floor, matching the direction and
magnitude `docs/MECHANISM_FMAX_PCSA_VALIDATION.md` Sec.4 already established at the peak level
(+0.3%, "no material effect") for the same reason: hip-crossing muscle Fmax corrections are a
genuine mix of up/down shifts (12/12 knee-crossing muscles were corrected uniformly down; the 25
hip-crossing muscles range 0.649–1.446x, roughly as many under 1.0 as over) that nets to
approximately zero, unlike the knee's uniformly-one-directional correction. This is a real
specificity control, not window-dressing: it demonstrates that re-solving SO on a model with
**many** muscles' Fmax simultaneously changed does **not**, by itself, produce a "push-off drops
more" pattern everywhere — only at the joint whose crossing muscles were actually, one-directionally
over-strong.

**A metric artifact caught and disclosed, not smoothed over:** the same generic
`change_is_selective_not_uniform` boolean this script computes for the knee also evaluates `true`
for the hip (spread/mean = 0.0117/0.00255 = 4.58, exceeding the 0.25 bar). This is **not** a real
selective effect — it is the formula's denominator (`mean_rel_drop` ≈ 0.00255, i.e. near zero)
blowing the ratio up when there is barely any change to divide by. The substantively correct
reading, from the separately-computed and more legible `hip_specificity_control` check
(`both_humps_within_5pct: true`, both raw drops <1%), is: **no material change at either hip
hump** — the selective/uniform dichotomy does not meaningfully apply when neither hump moved.
Flagged here explicitly so this near-zero-denominator degeneracy is not mistaken for a second
"selective" finding.

### 4.3 Anchor-convention cross-check — independent formulas converge (over-determination)

`docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s peak-level test used a **median-of-8-subjects'-own-peaks**
OrthoLoad anchor (258.22 %BW) and found the corrected/baseline ratio dropped 1.515x→1.255x
(**−17.1%** relative, closing 50.4% of the gap to 1.0). This waveform test uses a **completely
different** anchor convention — the ensemble **mean curve's own value at the fixed 48%GC push-off
location** (266.9 %BW) — and independently measures a **−17.0%** relative drop (closing 53.8% of
the gap to 1.0). These are different statistics (median-of-peaks-at-possibly-different-phases vs.
mean-curve-value-at-one-fixed-phase) computed by different scripts against different anchor
numbers, yet converge to within 0.13 percentage points on the relative drop and within 3.4 points on
the gap-closed fraction — a genuine over-determination, not a tautology (neither number was used to
compute the other; the waveform pipeline never reads the peak-test's ratio, and vice versa).

### 4.4 Mechanism check (secondary, bounded, geometric) — WHY it is selective, measured not narrated

Extracted per-muscle activation (already-committed `StaticOptimization_activation.sto`, both
baseline and corrected, no new OpenSim run) at each hump's %GC for the largest-corrected
knee-crossing muscles:

| muscle | role | baseline hump1 act. | baseline hump2 act. | push-off/weight-accept ratio |
|---|---|---:|---:|---:|
| gasmed_r | biarticular (knee+ankle) | 0.010 | 0.441 | **43.6x** |
| recfem_r | biarticular (knee+hip) | 0.010 | 0.227 | **22.6x** |
| gaslat_r | biarticular (knee+ankle) | 0.016 | 0.282 | **17.5x** |
| vaslat_r | knee-only (quad) | 0.128 | 0.010 | 0.08x |
| vasmed_r | knee-only (quad) | 0.064 | 0.010 | 0.16x |
| soleus_r | **mono-articular, ankle-only — does NOT cross the knee** | 0.042 | 0.072 | 1.7x |

The muscles that actually **drive** push-off in the baseline model (gastrocnemius, rectus femoris)
are exactly the biarticular knee-crossers whose Fmax got corrected downward; the quadriceps that
dominate weight-acceptance are already sitting at their activation floor (~0.01) during push-off, so
there is little further room for a knee-crossing Fmax cut to change push-off recruitment through
them. After correction, at push-off: gastrocnemius/rectus-femoris activation goes **down** (gasmed_r
0.441→0.408, recfem_r 0.227→0.170) — not up, ruling out a naive "same recruitment pattern, just
saturating harder against a lower ceiling" story — while soleus (mono-articular, contributes **zero**
knee-crossing force by anatomy) activation goes **up** (0.072→0.138, nearly doubling). This is
Static Optimization's minimum-activation-cost objective redistributing load, within the muscle
redundancy null-space, away from the now-more-expensive biarticular knee-crossers and onto a
synergist that satisfies the same required ankle moment **without** contributing to knee joint
reaction force at all — a geometric (moment-arm/redundancy-space) mechanism, directly measured from
the activation trajectories. At weight-acceptance, no such non-knee-crossing escape valve exists
among the dominant muscles (vasti and hamstrings are themselves knee-crossing), so quadriceps
activation rises slightly instead (vaslat_r 0.128→0.147) to partly compensate for their own reduced
capacity — the two effects (less force per muscle, more activation) partially cancel, leaving only
the small 2.0% net change measured at hump1.

### 4.5 Residual: a real, partial fix — not a full close

The corrected model's push-off hump (324.0 %BW) remains 21.4% above OrthoLoad (266.9 %BW,
Z=1.23) — the correction substantially shrinks but does not eliminate the over-prediction, exactly
matching `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s own disclosed finding at the peak level ("even
fully corrected, subject2's knee ratio remains a 25.5% over-prediction... something else remains
the dominant, unresolved driver"). This waveform-level test corroborates that conclusion rather than
revising it: the fix is real and phase-selective, but partial.

## 5. Honest gaps

1. **SO only — no CMC-corrected counterpart.** `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s causal
   correction test was built and run only through the Static-Optimization pipeline
   (`static_opt_knee.py`/`validate_hip_force.py`); no Fmax-corrected CMC model exists, so this
   analysis cannot say whether CMC's own concentrated-but-less-clean push-off signal
   (`docs/MECHANISM_CONTACT_WAVEFORM.md` Sec.3.2) would respond the same way. This was deliberate
   scope-matching (SO is "this repo's primary, most-used solve method," and the one showing the
   cleanest baseline signal this task's falsifier was built around), not an oversight, but it is a
   real, disclosed scope limit.
2. **n=1 subject/trial** (subject2/`walking1`) — same scope caveat inherited from every cert in
   this family; `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s own subject3 bonus check (peak-only) found
   a smaller but same-direction effect (−10.2%), suggesting the mechanism generalizes at least
   qualitatively, but that was never re-run at the waveform level here.
3. **Hump1's small 2.0% drop is directionally consistent but individually modest** — it sits well
   within this pipeline's own established smoothing-sensitivity band (3-9% per
   `docs/MECHANISM_CONTACT_WAVEFORM.md` Sec.3.2), so on its own it would not be distinguishable from
   noise; it is trusted here only in combination with the activation-level mechanism check
   (Sec.4.4), which shows the same small, physically-grounded direction (quad activation rises
   slightly to compensate), not because the 2.0% number alone is a strong signal.
4. **The `change_is_selective_not_uniform` generic flag is unreliable near-zero-mean-drop** (Sec.4.2)
   — disclosed and worked around with a directly-computed, more legible check for the hip, but the
   underlying formula (a ratio of two small deltas) should not be trusted blindly on future joints
   without checking the denominator magnitude first.
5. **Hip has no official OrthoLoad ensemble** (same inherited gap as
   `docs/MECHANISM_CONTACT_WAVEFORM.md` Sec.4.1 — the self-built raw-AKF ensemble, method-validated
   on knee at r=0.928) — affects only the specificity-control numbers (Sec.4.2), not the primary
   knee verdict.
6. **Hump location shifted by ≤1%GC between baseline and corrected** (knee hump2: 45%GC→44%GC) — a
   trivial, expected consequence of a genuinely different SO re-solve; the Z-scores above already
   account for this by reading the OrthoLoad ensemble's mean/SD at each curve's OWN hump location,
   not a shared fixed point.
7. **Mechanism check (Sec.4.4) is secondary and bounded** — 8 muscles at 2 time-points per model,
   not a full per-muscle/per-phase decomposition (a sibling in-progress analysis,
   `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`, looks adjacent but was not read or relied upon here,
   consistent with this repo's multi-instance convention) — sufficient to demonstrate the
   redistribution mechanism qualitatively and quantify its direction, not a formal claim about
   exact force contributions.

## 6. Verification (machine-checked, not narrated)

- Self-consistency: fresh-extracted corrected-model overall peaks (324.1093/388.3556 %BW) match
  already-published `jointreaction_peak_pct_bw` values to 0.0000% relative difference, both joints
  (script asserts `rel_diff < 0.5%`; a mismatch would raise, not silently pass).
- Controlled-experiment check: baseline vs corrected JointReaction `.sto` time vectors are
  byte-identical (158 frames, max abs diff 0.0); setup XMLs differ only in `model_file` and
  `results_directory` (verified via `diff`).
- Phase alignment: `verify_phase_alignment(...).all_pass == true` (both independent toe-off checks
  in-band), re-derived fresh in this script, not assumed from the baseline doc.
- OrthoLoad knee raw-AKF-method calibration vs official ensemble: r=0.9278 (gate r>0.8, PASS),
  recomputed fresh in this run (not read from cached baseline JSON).
- Corrected-model SO convergence gates re-surfaced and asserted:
  `convergence_pass=true`, `n_frames=158`, `muscles_over_fmax_ratio_1p5={}`,
  `joint_reserve_leaning_hard={}`.
- Symmetric-QC structural check: `grep -n -i orthoload` on `fmax_pcsa_correction_test*.py` shows
  zero references inside `build_corrected_model()` (the Fmax-setting function) — all references
  live in post-hoc reporting code.
- Pre-registered falsifier conditions for `H_selective` explicitly evaluated in-script
  (`z1_stays_below_1p0`, `ratio1_did_not_worsen`, `hump2_moved_toward_1`,
  `ratio1_fell_below_1p0_ortho_mean`, `z1_flipped_negative`) — all resolve in H_selective's favor,
  none of the falsifying conditions true.
- All numbers in this document are pulled directly from
  `data/msk_smoketest/subject2_walking1/contact_waveform_fmax_corrected/contact_waveform_fmax_corrected_results.json`,
  not transcribed from console prose.

## 7. Confidence tier and falsifier (restated per report discipline)

**Confidence tier: in-vivo-anchored (OrthoLoad, for the waveform target and verdict) combined with
cadaveric/published-plausibility (Handsfield-MRI-volume-via-Arnold-2010-specific-tension, for the
Fmax correction itself)** — inherits both anchor tiers already established in
`docs/MECHANISM_CONTACT_WAVEFORM.md` and `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`; no new anchor is
introduced by this document, only a new combination of the two on the same corrected model.

**Claim:** the Handsfield-MRI-PCSA Fmax correction SELECTIVELY reduces the knee push-off hump
(17.0% drop, closing 53.8% of the gap to OrthoLoad) while leaving the weight-acceptance hump's
already-good match intact (2.0% drop, Z stays at 0.74, still <1.0) — a real, phase-selective,
mechanistically-explained (Sec.4.4) partial fix, not a uniform rescaling and not evidence the
original weight-acceptance match was luck.

**Falsifier (what would have overturned this claim, none of which occurred):** the corrected
model's weight-acceptance Z-score reaching ≥1.0; the weight-acceptance ratio failing to improve
(or worsening) relative to baseline; or the push-off ratio failing to move toward 1.0 at all. A
secondary falsifier that also did not occur: the hip specificity control showing a comparably large
shift at either hump (which would have suggested a generic "any Fmax re-solve" artifact rather than
a knee-specific mechanism) — both hip humps moved by <1%.

## 8. Files

- `scripts/msk/contact_waveform_fmax_corrected.py` (new) — the full pipeline. Run with
  `/usr/bin/python3 scripts/msk/contact_waveform_fmax_corrected.py` (needs numpy+openpyxl, same
  convention as `contact_waveform_analysis.py`; never imports OpenSim). Imports
  `scripts/msk/contact_waveform_analysis.py` **unedited** for every shared function.
- `data/msk_smoketest/subject2_walking1/contact_waveform_fmax_corrected/contact_waveform_fmax_corrected_results.json`
  (new) — every number in this document: phase-alignment re-check, corrected-model curve
  metadata + self-consistency gate, fresh OrthoLoad ensemble diagnostics, both full comparison
  records (knee/hip, RMSE/Pearson-r/hump values/ratios/Z-scores), the before/after delta table with
  pre-registered verdicts, and the per-muscle mechanism-check activation table.
- Inputs read in place, never modified: `data/msk_smoketest/subject2_walking1/fmax_pcsa_correction_test/{jr,so}/`
  (corrected-model JointReaction + activation outputs, already committed by
  `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s own run), `data/msk_smoketest/subject2_walking1/static_optimization/{jr,so}/`
  (baseline outputs), `data/msk_smoketest/subject2_walking1/contact_waveform/contact_waveform_results.json`
  (baseline waveform comparison, for the delta table only),
  `data/msk_smoketest/subject2_walking1/fmax_pcsa_correction_test/fmax_pcsa_correction_test_subject2_results.json`
  (published corrected-model peak numbers, for the self-consistency gate),
  `data/external/orthoload/{knee,hip_gen1,hip_gen2}/...` (OrthoLoad source files, re-parsed fresh),
  `data/msk_models/subject2_scaled_handsfield_fmax_corrected.osim` (referenced only via its
  already-computed `.sto` outputs — this script never loads the `.osim` itself, no OpenSim import).
- **Not touched**: `scripts/msk/contact_waveform_analysis.py`, `scripts/msk/fmax_pcsa_correction_test.py`,
  `scripts/msk/fmax_pcsa_correction_test_subject2.py`, `scripts/msk/subject_specific_scaling.py`,
  `scripts/msk/validate_joint_force.py`, `scripts/msk/validate_hip_force.py`, `scripts/msk/static_opt_knee.py`.
- Prior docs this builds on (context, not re-litigated): `docs/MECHANISM_CONTACT_WAVEFORM.md`,
  `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`.

Isolation respected throughout: bodytwin only; all inputs read in place; no OpenSim import in this
script (`.sto`/`.json`/`.xlsx`/`.akf` re-parsing only, same environment convention as
`contact_waveform_analysis.py`). No git commit, no git add, no git push performed.
