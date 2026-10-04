# MECHANISM RESIDUAL DECOMPOSITION — what is the POST-Fmax-correction excess made of? (2026-07-21)

**Mid-task course correction, addressed and resolved (read this first):** this doc was originally
briefed to apply a JAM deformable-contact relief factor of "~0.7-0.9x" to the corrected knee number.
Partway through, the coordinating thread flagged that the JAM run backing that factor
(`docs/MECHANISM_JAM_CONTACT_DECORR.md`) had every muscle FROZEN at a flat 0.05 activation for the
whole trial — the real per-frame SO activation never reached the muscles. **This session
independently re-verified that claim** (read both h5 files directly, §3 below) and, further, **ran a
brand-new JAM re-solve on the Fmax-corrected model's own activation and found the identical frozen-
0.05 signature** — proving the bug is structural to the JAM invocation itself, not specific to one
input file. The contact-model leg of this decomposition is therefore reported as a **DIAGNOSED-GAP**,
not a number — exactly as directed. Everything else below (the Fmax recap, and a brand-new joint
Fmax+l_TS sweep on the corrected model) is a real, fresh measurement.

---

## 0. Pre-registration (stated before the l_TS-on-corrected sweep and the JAM re-run executed)

- **Knee falsifier:** PRIMARY-band (±10%) max–min swing of the push-off ratio (`ratio2`, vs
  OrthoLoad's push-off-hump anchor 266.9 %BW) computed **jointly** on the Fmax-corrected model >0.15
  absolute ⇒ l_TS remains FRAGILE even after Fmax-correction; ≤0.15 ⇒ ROBUST on this axis.
- **Closing criterion:** l_TS "explains" the residual only if a PRIMARY-band point moves `ratio2`
  materially toward 1.0 **without** hump1 (weight-acceptance) crossing below the OrthoLoad band
  (the same over-correction check `docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md` used), **and** the point
  is not an isolated knife-edge immediately flanked by a worse one.
- **No-fit rule:** the full pre-registered 9-point sweep (identical scales to
  `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md`: ±5/10/15/20%) is reported; no new scale point was
  searched post-hoc to find a better-looking answer.
- **Double-count guard:** l_TS is applied **on top of**, not instead of or additively alongside, the
  Fmax-corrected model — verified live at every sweep point (§2.2).
- **Symmetric-QC on the method itself:** if the joint sweep is fragile, that is reported as fragile,
  not smoothed into a point estimate; if the contact-model leg is broken, that is a diagnosed-gap, not
  silently patched with the old (now-known-bad) number.
- **Confidence tier:** method (attribution), anchored externally to OrthoLoad (knee/hip in-vivo
  medians) and to this subject's own real EMG (the hip's glmed1_r correction, inherited unchanged).

---

## 1. Headline

| joint | baseline | Fmax-correction | residual after Fmax | l_TS (joint, this session) | contact-model | net verdict |
|---|---:|---:|---:|---|---|---|
| **knee push-off** | 390.3 %BW (1.462×) | −17.0%, closes 53.8% of gap | **324.0 %BW (1.214×, +21.4%)** | **FRAGILE, knife-edge**: swing 0.573 (PRIMARY)/0.697 (FULL), ≥3.8× the 0.15 bar. Actual model sits near the LOCAL MAX of the sensitivity spike; ±5–20% l_TS moves the ratio anywhere from 0.616 (under by 38%) to 1.313 (over by 31%) | **DIAGNOSED-GAP** (JAM run invalid, §3) | **NOT a clean partial explanation** — the residual's apparent size is itself not well-determined once realistic l_TS uncertainty is admitted |
| **hip** | 386.77 %BW (1.412×) | ~0% (negligible, confirmed) | 1.086–1.194× (9–19% over) after glmed1_r EMG correction | measured on the PRE-correction model only (scope-disclosed): best l_TS point (scale 1.10) only reaches ratio 1.391 — **worse** than the already-EMG-corrected residual, i.e. l_TS cannot explain this residual | **no estimate exists at all** | **9–19% stands as the best-supported number**, not further reducible by any mechanism tested |

**Falsifier resolution:** neither pre-registered outcome ("fully explained" vs "clean irreducible
floor") is exactly what was found. The over-prediction is **not fully explained** (contact-model gap,
plus l_TS knife-edge for the knee) — but it is also **not a single, stable "floor" number** for the
knee, because the residual's magnitude is itself highly sensitive to a poorly-constrained nuisance
parameter (l_TS) at a scale comparable to the whole effect being measured. This is itself a **stronger
version of the Moissenet SO-then-subtract critique**: the method's own sensitivity to unmeasured Hill
parameters can swing the answer from "under-predicts" to "over-predicts" — an honest limitation of the
method's *resolution*, not just its central estimate. The hip's 9–19% residual, by contrast, **is** a
comparatively stable floor (l_TS tested and found unable to move it; no contact-model estimate exists
to test further).

---

## 2. Knee — full attribution chain

### 2.1 Fmax correction (done, recap — one correction to my own framing, disclosed)

`docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md`: push-off hump2 390.3 → 324.0 %BW, ratio 1.462× → 1.214×
(−17.0%, closes 53.8% of the gap to 1.0); weight-acceptance hump1 barely moves (245.4→240.6, ratio
1.093→1.071). **Scope correction, verified live this session** (not assumed from prior doc prose): the
Handsfield-PCSA Fmax correction (`fmax_pcsa_correction_test.py:build_corrected_model`) recalibrates
**all 80 muscles'** Fmax, not just the 12 knee-crossers — 80/80 muscles differ between
`subject2_scaled_handsfield_fmax_corrected.osim` and the original scaled model (checked directly via
the OpenSim API, both models loaded fresh). The correction's push-off-**selective** effect is a
Static-Optimization redistribution outcome (which muscles dominate at which phase), not a consequence
of a narrowly-scoped edit. This does not change any prior verdict — it corrects an imprecise
description of the correction's scope.

### 2.2 l_TS, applied JOINTLY on top of the Fmax-corrected model (new this session)

**Why this had to be re-run, not reused:** `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md` swept l_TS on
the **original** (Fmax-uncorrected) model and read a whole-trial argmax (which needed a post-hoc
window-fix for 5/9 points, `tendon_slack_sensitivity_diagnose.py`). Neither model nor extraction method
matches the corrected-model push-off headline. `scripts/msk/residual_decomposition_knee_lts.py` fixes
both: starts from `subject2_scaled_handsfield_fmax_corrected.osim`, and extracts push-off via the
%GC-hump-window method (`contact_waveform_analysis.find_two_humps`, hump2 = max over [30,62]%GC) —
phase-anchored by construction, immune to the boundary artifact from the start.

**Self-consistency gate (checked before trusting anything else):** the scale=1.00 control point (l_TS
literally unchanged, i.e. the corrected model re-solved from scratch) reproduces the published
324.0/1.214 to **0.015%/0.015% relative difference** — PASS.

**Full 9-point table** (same 27 prime-mover muscles, same ±10%/±20% bands as the original sweep,
scaled from the CORRECTED model's own l_TS values):

| scale | band | hump2 %BW | ratio2 | Z2 | hump1 %BW | ratio1 | Z1 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0.80 | STRESS | 198.08 | 0.742 | −1.49 | 242.14 | 1.078 | 0.81 |
| 0.85 | STRESS | 197.60 | 0.740 | −1.50 | 241.76 | 1.077 | 0.80 |
| 0.90 | PRIMARY | 197.48 | 0.740 | −1.50 | 236.18 | 1.052 | 0.54 |
| 0.95 | PRIMARY | 208.67 | 0.782 | −1.26 | 235.42 | 1.048 | 0.50 |
| **1.00 (actual)** | PRIMARY | **323.95** | **1.214** | **1.23** | 240.62 | 1.071 | 0.74 |
| 1.05 | PRIMARY | 350.37 | 1.313 | 1.81 | 238.67 | 1.063 | 0.65 |
| 1.10 | PRIMARY | 222.42 | 0.833 | −0.96 | 238.89 | 1.064 | 0.66 |
| 1.15 | STRESS | 169.53 | 0.635 | −2.11 | 237.05 | 1.056 | 0.58 |
| 1.20 | STRESS | 164.37 | 0.616 | −2.22 | 236.85 | 1.055 | 0.57 |

All 9 points: SO convergence PASS, `n_frames=158`, zero muscles >1.5×Fmax, zero reserve-actuator
degeneracy, **identical** pelvis residual (173.0 N / 45.1 N·m) at every point — a clean, non-degenerate
sweep, matching the original sweep's own validity gates.

**PRIMARY-band swing: 0.573. FULL-band swing: 0.697. Both ≥3.8× the pre-registered 0.15 fragility
bar ⇒ FRAGILE**, confirming and roughly matching in magnitude (same order, somewhat smaller in
absolute terms: 0.573 vs the original sweep's 0.809) the original l_TS fragility finding — **it
survives, essentially undiminished, jointly with Fmax-correction.**

**The knife-edge, stated plainly:** the model's own *actual* l_TS (scale=1.00, no perturbation at all)
sits almost exactly at the **local maximum** of a narrow sensitivity spike — ratio2=1.214, with an
*even worse* neighbor one step away (scale=1.05, ratio2=1.313). One step the OTHER way (scale=1.10,
still inside the ±10% primary band), the ratio drops to 0.833 (|Z2|=0.96, borderline-consistent with
OrthoLoad) — but scale=1.15/1.20 overshoot into a 37–38% *under*-prediction. **A physiologically
plausible l_TS uncertainty of ±5–10% is enough to swing the verdict from "worse than the current
21% over-prediction" through "statistically matches OrthoLoad" to "under-predicts by over a third" —
depending on a parameter this repo has already disclosed is poorly measured and not independently
re-verifiable this session** (`docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md` §2's confidence-tier note).
**No point here is picked as "the answer"** — picking scale=1.10 because it looks closest to 1.0 would
be exactly the forbidden tuning-to-the-anchor.

**Mechanism, confirmed not assumed (geometric, redundancy-null-space):** direct read of this sweep's
own SO activation output at push-off (t=0.51s) shows the SAME antagonist-pair discrete-switching
mechanism `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md` §5 diagnosed on the original model, now recurring
on the corrected one:

| scale | gasmed_r act. | recfem_r act. | soleus_r act. |
|---:|---:|---:|---:|
| 0.80 | 0.012 | 0.012 | 0.382 |
| 0.90 | 0.010 | 0.010 | 0.374 |
| 1.00 | **0.417** | **0.162** | 0.154 |
| 1.05 | **0.380** | **0.247** | 0.122 |
| 1.10 | 0.239 | 0.181 | 0.265 |
| 1.20 | 0.012 | 0.012 | 0.382 |

`gasmed_r`+`recfem_r` (a knee flexor/extensor antagonist pair) switch from near-silent (~0.01) to
dominant (0.16–0.42) inside a narrow l_TS window, with `soleus_r` (mono-articular, does not cross the
knee) taking the inverse pattern — Static Optimization's minimum-cost recruitment discretely
"discovering" a cheap co-contraction solution in a narrow parameter band. This is a real,
geometrically-grounded redundancy-space effect (not a solver artifact — all convergence/degeneracy
gates PASS at every point), and it is precisely the mechanism that makes the push-off contact force so
sensitive to l_TS in the first place.

### 2.3 Contact-model (rigid vs. deformable) — DIAGNOSED-GAP, not a number

**What was found (forced, not assumed):** the original JAM leg's own output
(`data/msk_smoketest/jam_contact_decorr/results/joint-mechanics/gait_driven.h5`) and a **brand-new**
JAM re-solve this session, driven by the Fmax-corrected model's own SO activation
(`gait_driven_fmax_corrected.h5`), were both read directly. In **both** files, every checked muscle —
the 13 real-activation-driven knee/pf-crossers (`gasmed_r`, `recfem_r`, `vasmed_r`, `vaslat_r`,
`gaslat_r`) and a non-crossing control (`addbrev_r`) — shows activation **frozen at exactly 0.0500**
for all 159 frames (`std` = 1.4e-17, i.e. bit-identical), even though the actually-fed
`actuator_input_file.sto` varies from 0.010 to 0.42+ over the trial (std 0.02–0.14 per muscle). The
ForsimTool never applied the real per-frame activation to any muscle, in either run.

**Why this matters for the interaction question the task asked about:** because both runs share the
identical frozen fallback, the fact that my fresh re-run's peak (330.29 %(lenhart-BW)) was close to
the original's (345.37 %(lenhart-BW), −4.4% relative) is **not** evidence that JAM is insensitive to
the Fmax-correction's activation shift — it is an artifact of both runs never having seen real
activation at all. **The previously-published JAM ratio (0.883× primary-norm / 0.700× alt-norm) is
retracted as a usable relief factor.** My own attempted joint measurement (which would have read
1.019× primary-norm / 0.808× alt-norm) is retracted for the same reason and is **not** reported as a
finding.

**Root cause: not yet identified**, bounded not chased (out of scope for this task, per the
coordinator's explicit direction). Ruled out (checked, not assumed): not the explicit
`<constant_muscle_control>0.02</constant_muscle_control>` fallback (frozen value is 0.05, not 0.02);
not an XML `<default_activation>`/`<activation>` tag inside `lenhart2015.osim` (zero grep matches). One
unconfirmed lead, disclosed as a lead only: `actuator_input_file.sto`'s columns are named
`<muscle>_activation` (e.g. `gasmed_r_activation`); if this JAM ForsimTool build expects bare actuator
names, the whole file could silently fail to match and fall back to an internal default — not verified.

**Consequence for the decomposition:** the contact-model leg contributes **zero** attributable %
to this residual — not because it was tested and found negligible, but because no valid measurement
of it exists yet. It must be re-run (with the frozen-activation bug fixed and independently confirmed
gone, e.g. `std(activation) > 0` and tracking the fed file) before it can enter any attribution.

### 2.4 Stacked attribution (knee push-off)

```
baseline over-prediction:            390.3 %BW / 1.462x  (+46.2%)
  − Fmax correction (done, measured): -17.0% abs, closes 53.8% of the gap-to-1.0
  = residual after Fmax:             324.0 %BW / 1.214x  (+21.4%)
      − l_TS (measured, JOINT):      NOT a stable subtraction -- ratio2 ranges 0.616-1.313
                                     across a plausible +/-20% band; the "true" l_TS is
                                     unknown, so no single %-closed number can be honestly stated
      − contact-model:                DIAGNOSED-GAP (0%, pending a JAM re-run with the
                                     frozen-activation bug fixed)
  = "irreducible structural floor":   NOT a clean number -- best honest statement is a RANGE
                                     (-38% to +31% around the naive 21.4%), reflecting the SO-
                                     then-subtract method's own sensitivity to a poorly-measured
                                     nuisance parameter, not a settled residual
```

---

## 3. Hip — what remains after the glute-med/EMG correction

### 3.1 Fmax correction: negligible, confirmed via JSON

`contact_waveform_fmax_corrected_results.json`: hip hump1 rel-drop +0.84%, hump2 rel-drop −0.33% (both
`<5%`, `both_humps_within_5pct: true`) — matches the already-published "nets to ~0" finding, re-checked
from the raw numbers here rather than trusted from prose.

### 3.2 glmed1_r EMG correction (already established, `docs/MECHANISM_HIP_STRUCTURAL_MECHANISM.md`)

Baseline 386.77 %BW / 1.412× (vs OrthoLoad hip median 273.93 %BW, n=162). `glmed1_r`'s SO tension
(585.31 N) vs this subject's own real EMG-implied tension (99.76 N, ratio 0.17) is over-recruited by
SO. A **moment-budget-preserving** reallocation (verified: `moment_eq_max_abs_resid` ~1e-15 N·m,
i.e. the exact same 3-DOF hip moment the original 25-muscle solution produced) closes:
- **LP-optimal (bang-bang):** 297.58 %BW, ratio **1.086** (79.0% of the gap to 1.0 closed)
- **Smooth QP (graded, physiologically smoother):** 327.04 %BW, ratio **1.194** (52.9% closed)

Both variants are moment-valid and clear the pre-registered 0.10-ratio materiality bar comfortably —
this is the one mechanism in the whole investigation that survives both a bang-bang **and** a smooth
adversary (§5 of that doc).

### 3.3 l_TS — scope-limited, does NOT explain the residual

**Disclosed scope limit:** unlike the knee, l_TS was **not** re-swept jointly with the glmed1_r EMG
correction this session (the EMG correction is an LP/QP tension-reallocation overlay, not a model-file
edit like Fmax was, so a joint re-sweep would require re-solving that LP at every l_TS scale point — a
larger lift, out of this task's scope; disclosed here, not hidden).

What **is** known (`docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md`'s original 9-point sweep, on the
PRE-correction baseline model): the hip ratio ranges from 1.391 (scale=1.10, the single best point in
the **entire** 9-point sweep) to 1.685 (scale=0.85) — l_TS perturbation **never** gets the hip anywhere
near 1.0, and its best case (1.391) is **worse** than the residual the glmed1_r correction alone already
achieves (1.086–1.194). **l_TS cannot plausibly explain the hip's residual**, even before accounting
for any interaction with the EMG correction.

### 3.4 Contact-model: no estimate exists

No JAM/deformable-contact leg has ever been attempted for the hip in this repo (the one JAM leg that
exists, `docs/MECHANISM_JAM_CONTACT_DECORR.md`, is knee-only: tibiofemoral + patellofemoral). This gap
is larger than the knee's (which at least has a — now invalidated — attempt).

### 3.5 Stacked attribution (hip)

```
baseline over-prediction:            386.77 %BW / 1.412x  (+41.2%)
  − Fmax correction (measured):      ~0% (negligible, confirmed)
  − glmed1_r EMG correction:         closes 53-79% of the gap (LP 1.086x / QP 1.194x)
  = residual:                        1.086x-1.194x  (+8.6% to +19.4%)
      − l_TS (measured, NOT joint):  cannot close it further (best l_TS point alone: 1.391x,
                                     WORSE than the already-corrected residual)
      − contact-model:                no estimate exists (bigger gap than knee)
  = best-supported residual:          9-19% over -- the most STABLE number in this whole
                                     decomposition (unlike the knee, nothing tested moves it)
```

---

## 4. Symmetric-QC checklist (what was forced, not assumed)

- **Double-counting guard (knee):** every l_TS-sweep point's Fmax was verified, live, against the
  CORRECTED model's own Fmax (not the original's) for all 80 muscles — `joint_application_check:
  all_match_corrected_fmax = True` at all 9 points. l_TS and Fmax act jointly, not additively, by
  construction.
- **No-fit guard:** the knee's 9-point sweep is reported in full; scale=1.10's closer-to-1.0 value is
  explicitly flagged as a knife-edge, not adopted as "the" answer.
- **JAM-factor guard:** the previously-instructed "apply 0.7–0.9× to the corrected number" step was
  **not executed as instructed** — it was independently re-verified first, found invalid on **both**
  the original run and a fresh joint re-run, and reported as a diagnosed-gap instead. This is the
  single most consequential correction in this document.
- **Self-consistency gates:** knee sweep's scale=1.00 control reproduces the published 324.0/1.214 to
  0.015%; the JAM re-run's mesh-pair agreement (tibia- vs femur-side casting) is 2.83%, consistent with
  the original run's own <10% gate.
- **Convergence gates:** all 9 knee sweep points, SO convergence PASS, zero saturated muscles, zero
  reserve-actuator degeneracy, bit-identical pelvis residual (a scale-independent invariant, confirming
  no drift in anything except the intended l_TS channel).

---

## 5. Honest gaps

1. **Contact-model leg is a diagnosed-gap for the knee, and non-existent for the hip** — the single
   largest open item. A real number requires the JAM ForsimTool's frozen-activation bug fixed and
   independently reconfirmed gone (activation `std>0`, tracking the fed file to a stated tolerance)
   before any deformable-contact relief factor can be trusted again.
2. **The knee's l_TS knife-edge means "the residual" is not a single number** — this document reports a
   range (0.616×–1.313× across the tested band) rather than a point estimate, which is itself the
   honest finding, not a failure to compute one.
3. **Hip l_TS × EMG-correction interaction untested** (§3.3) — a real, disclosed scope limit; the
   current reading (l_TS alone can't help) is very likely to hold jointly too (its best-case alone is
   already worse than the EMG-corrected residual), but this was not directly verified the way the
   knee's interaction was.
4. **Single subject (subject2), single trial (walking1)** — same scope caveat as every cert in this
   family; no claim of generality across subjects/trials.
5. **The JAM bug's root cause is unconfirmed** — one plausible, disclosed, UNVERIFIED lead
   (actuator-name column-naming convention) is offered for whoever fixes it next; this session did not
   chase it further, per explicit scope direction.
6. **Mechanism check (§2.2) covers 3 muscles at 6 of 9 scale points** — sufficient to confirm the
   antagonist-switching mechanism recurs and is not a solver artifact, not a full per-muscle/per-point
   decomposition.

---

## 6. Confidence tier and falsifier (restated per report discipline)

**Confidence tier: method (attribution), anchored externally to OrthoLoad** (knee/hip in-vivo medians
and push-off-hump ensemble anchors, re-verified live this session) **and to this subject's own real
EMG** (hip glmed1_r correction, inherited unchanged from `docs/MECHANISM_HIP_STRUCTURAL_MECHANISM.md`).
The contact-model leg carries no confidence tier at all this session — it is retracted, not merely
downgraded.

**Claim:** the knee's post-Fmax residual (+21.4%) is NOT a stable, fully-attributable number — l_TS
uncertainty alone (tested jointly with Fmax-correction, for the first time this session) can swing it
from a 38% under-prediction to a 31% over-prediction, and the contact-model leg that might have further
bounded it is currently unusable (diagnosed-gap, not a relief factor). The hip's post-correction
residual (+8.6% to +19.4%) is comparatively stable: l_TS (tested, if not jointly) cannot close it
further, and no contact-model estimate exists to test further — this is the best-supported "floor" in
the whole investigation.

**Falsifier (what would overturn this document):** a fixed JAM pipeline (activation verified
time-varying and tracking the fed file) producing a tf_contact ratio that, applied to the
Fmax-corrected model, falls clearly outside the knife-edge band already measured; or an independent
determination of subject2's TRUE l_TS (e.g. from imaging or a subject-specific cadaveric measurement)
that pins down which point in the 9-point sweep is actually correct — either would convert this
document's range into a point estimate.

---

## 7. Files

- `scripts/msk/residual_decomposition_knee_lts.py` — the joint l_TS-on-Fmax-corrected-model sweep
  (new this session). Reuses `tendon_slack_sensitivity.py` (`build_scaled_model`,
  `verify_model_edit_isolated`), `subject_specific_scaling.run_joint_forces_for_model`,
  `contact_waveform_analysis.py` (`model_pct_gc_mapping`, `model_curve_on_grid`, `find_two_humps`,
  `compare`, `knee_official_ensemble`) **all unchanged** — zero lines edited in any reused module.
- `scripts/msk/jam_contact_decorr_fmax_corrected.py` — the JAM re-run on the Fmax-corrected model's
  activation (new this session); reuses `jam_contact_decorr.py` unchanged via monkeypatched
  module-level path constants (same idiom `jam_contact_decorr_leftcheck.py` already established).
- `scripts/msk/residual_decomposition_consolidate.py` — pulls every number in this document live from
  its own source JSON (no hand-transcription).
- `data/msk_smoketest/subject2_walking1/residual_decomposition_knee_lts/residual_decomposition_knee_lts_results.json`
  — the full 9-point knee sweep: self-consistency gate, per-point convergence/joint-application checks,
  hump1/hump2 values/ratios/Z-scores, swing/verdict tables.
- `data/msk_smoketest/subject2_walking1/jam_contact_decorr_fmax_corrected/jam_contact_decorr_fmax_corrected_results.json`
  — the retracted JAM re-run, WITH the frozen-activation invalidation block and per-muscle
  activation-std diagnostic for both the original and corrected-activation runs (machine-checkable,
  `std≈1.4e-17` in both vs. the fed file's real `std=0.02–0.14`).
- **`data/msk_smoketest/subject2_walking1/residual_decomposition/residual_decomposition_evidence.json`**
  — the consolidated deliverable evidence file (pre-registration, knee table, hip table, both
  stacked attributions, all source pointers).
- Prior docs this builds on (context, not re-litigated): `docs/MECHANISM_FMAX_CORRECTION_WAVEFORM.md`,
  `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md`, `docs/MECHANISM_JAM_CONTACT_DECORR.md`,
  `docs/MECHANISM_HIP_STRUCTURAL_MECHANISM.md`, `docs/MECHANISM_OVERPREDICTION_DECOMP.md`.
- Inputs read in place, never modified: `data/msk_models/subject2_scaled_handsfield_fmax_corrected.osim`,
  the original (uncorrected) scaled model on the external drive, subject2's IK/GRF (external drive),
  `scratchpad/knee_lig/lenhart2015.osim` + Geometry (git-ignored JAM assets), the JAM `opensim-cmd`
  binary (`/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/opensim-build/opensim-cmd`,
  read-only per isolation).
- Additive, non-destructive environment change: `openpyxl` installed into `.venv-msk` (was absent;
  needed to reuse `contact_waveform_analysis.py`'s OrthoLoad-ensemble parsing in the same script as the
  OpenSim SO/JR pipeline) — same precedent as `hip_structural_mechanism.py`'s prior `scipy` install.

No git commit, no git add, no git push performed (isolation respected). All new files left untracked
for the coordinator to stage explicitly.
