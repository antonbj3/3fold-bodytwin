# MECHANISM — Fiber-length correction cross-subject: does the knee effect generalize beyond subject2?

**Question.** `docs/MECHANISM_FIBER_LENGTH_CORRECTION.md` found that replacing subject2's over-long
`optimal_fiber_length` with Ward, Eng, Smallwood & Lieber 2009 (Clin Orthop Relat Res, PMID 18972175,
cadaveric dissection Table 3) values for 6 muscles — pennation angle and `max_isometric_force` held at
baseline, machine-asserted unchanged — drops subject2's knee ratio-vs-OrthoLoad from **1.515× to
0.944×** (over-corrects, crosses *below* the in-vivo anchor) and *worsens* the hip (1.412×→1.735×).
Separately, `docs/MECHANISM_CROSS_SUBJECT.md` (no fiber-length correction) established subject2 is the
**mildest** of 3 subjects at baseline: subject3 knee 1.990×, subject4 knee 1.681×, both worse than
subject2's 1.515×. This doc repeats the *identical* fiber-length-only correction on subject3 and
subject4 and asks: does it also drop their knee toward/below OrthoLoad (over-long fibers = cohort-wide
knee driver), or does subject2's specific over-correction fail to reproduce (subject2 idiosyncratic)?

## Headline: DIRECTION generalizes robustly; subject2's specific OVER-correction does NOT

| subject | baseline knee ratio | corrected knee ratio (trustworthy reading) | relative drop | crosses ≤1.0? |
|---|---:|---:|---:|---:|
| subject2 (reference, already published) | 1.515× | 0.944× (self) / 0.960× (JR) | 37.6% | **YES — over-corrects** |
| subject3 | 1.990× | **1.504×** (self=JR, agree <0.01%) | 24.4% | NO |
| subject4 | 1.681× | **1.188×** (JR-consistent — see Sec.4, self-computed's naive 1.104× is a boundary-truncation *underestimate*) | 29.3% | NO |

**Pre-registered verdict: `C_GENERALIZES_DIRECTION` = TRUE, `C_GENERALIZES_OVERCORRECTION` = FALSE.**
Fiber-length correction is confirmed as a **real, quantified, cohort-wide knee-force driver** — both
new subjects clear the same ≥10%-relative-drop bar subject2 cleared, by a *larger* margin (24–29% vs
subject2's own qualifying threshold), in the same direction, with no degenerate solve. But subject2's
specific outcome — closing and then *overshooting* the OrthoLoad anchor — is subject2-specific, not
cohort-wide: subject3/4 start from a substantially worse baseline (1.99×/1.68× vs subject2's 1.52×)
and the correction narrows, but does not close, their gap.

## 1. Method (reuse discipline, zero edits to any existing cert script)

New script: `scripts/msk/fiber_length_correction_cross_subject.py`. Reused **unedited** (md5-verified
before/after, Sec.7): `cross_subject_validation.configure_for_subject` (path monkeypatch to
subject3/4's own data), `fiber_length_correction_test.build_fiber_length_corrected_model` (the
**identical** Ward2009 dict and pennation/Fmax/tendon-slack-length invariance asserts subject2's own
test used), `subject_specific_scaling.run_joint_forces_for_model` / `compute_joint_forces_from_so_output`
(fresh SO+JR resolve, `bw_n` always freshly recomputed from the model's own `BodySet`, never a stale
default), `static_opt_knee.check_so_convergence_and_sanity`, `ankle_spine_cross_subject.
SUBJECT_BW_N_EXPECTED` (trap-tripwire reference values, not re-transcribed). The only new code is the
cross-subject comparison/verdict logic and the compound-bonus orchestration.

### ★ Known trap (task-mandated, forced explicitly, not just avoided)

A def-time-default `bw_n` binding bug is real and documented elsewhere in this repo
(`contact_waveform_analysis.py`'s `joint_force_mag_pct_bw(..., bw_n=BW_N)`, `BW_N=766.88003`=subject2's
own body weight, bound once at function-definition time). **Verified this script's actual call chain
does not use that function** — `compute_joint_forces_from_so_output` recomputes
`total_mass = sum(model's own BodySet masses); bw_n = total_mass * G` fresh from the model object every
call (read directly from `scripts/msk/subject_specific_scaling.py` source). Enforced by two independent,
machine-checked layers anyway, never trusted on code-reading alone:

1. **TRAP TRIPWIRE** (`trap_tripwire()` in the new script, `raise`s on failure): every corrected run's
   own returned `body_weight_N` is asserted to equal that subject's independently-known value, not
   subject2's. **Result: subject3 622.7222750000001N == expected (exact); subject4 613.89629N ==
   expected (exact); both asserted NOT within 50N of subject2's 766.88003N.** Passed for all 4 runs
   (subject3/4 × primary/compound).
2. **SELF-CONSISTENCY GATE** (Part 0, run FIRST, gates the whole script — `RuntimeError` on failure):
   subject2's own already-published fiber-length-corrected knee number was reproduced through this
   script's `configure_for_subject` orchestration layer (never used by the original subject2 script,
   which relied on `vjf`/`sok`/`hvf`'s hardcoded subject2 defaults). **Result: published corrected knee
   ratio (self) = 0.9444272548177708; fresh, via `configure_for_subject` + a new SO/JR resolve on the
   *existing* corrected model = 0.9444272548177708. Relative diff = 0.0000% (bit-identical), well inside
   the pre-registered 1% tolerance. `bw_n` fresh=766.8800N == expected.** The orchestration layer
   introduces zero regression before being trusted on subject3/4.

## 2. Pre-registration (thresholds fixed before subject3/4 were run)

- **Q1 bar** (reused verbatim from `fiber_length_correction_test.py` / `MECHANISM_FMAX_PCSA_VALIDATION.md`,
  not loosened): ≥10% relative drop in self-computed knee ratio-vs-OrthoLoad = "real, quantified
  contributor" for that subject.
- **`C_GENERALIZES_DIRECTION`** (primary, the task's own falsifier): BOTH subject3 AND subject4 clear
  the Q1 bar, same direction as subject2. NOT-C: either subject drops <10% relative or moves the wrong
  direction (increases) — would falsify "cohort-wide driver", support "subject2-specific magnitude".
- **`C_GENERALIZES_OVERCORRECTION`** (secondary, stronger, subject2-mirroring): BOTH subjects' corrected
  ratio crosses to ≤1.0, matching subject2's specific 0.944× overshoot. Not required for the primary
  clause to hold.
- **Hip specificity control** (secondary, corroborating, NOT the primary falsifier — the task's own
  falsifier language is knee-specific): does the hip also worsen, matching subject2's direction?
  Predicted, before running, to be much smaller than the knee effect (only 1/25 hip-crossing prime
  movers, `bflh_r`, vs 5/12 knee-crossing prime movers, are touched by this correction).
- **Symmetric-QC forced adversaries**: (a) leaning-positive — a ratio drop that is really an
  SO-degeneracy artifact, forced via `muscles_over_fmax_ratio_1p5=={}` and
  `joint_reserve_leaning_hard=={}` on every corrected solve; (b) leaning-negative — a non-effect that is
  really a subject-specific pipeline defect, forced via the subject-agnostic solve-quality re-derivation
  (frames>100, zero NaN, zero out-of-bounds activations) since the raw `convergence_pass` gate fails *by
  construction* for any non-subject2 trial (`docs/MECHANISM_CROSS_SUBJECT.md` Sec.4); (c) the
  already-disclosed pelvis-residual-gate FAIL (no RRA) is inherited, restated, and machine-checked not to
  materially worsen under correction.
- **Confidence tier** (stated up front, per task): cadaveric (Ward2009) × in-vivo (OrthoLoad), same tier
  as subject2's own cert; solve-quality-caveated for subject3/4 (raw `convergence_pass` reads False by a
  subject2-specific 1.50s literal, pelvis-residual gate fails — both pre-existing, disclosed, inherited,
  not new).

## 3. Per-subject knee result (primary test: fiber-length only, Fmax/pennation at baseline)

Ward2009 targets are identical to subject2's own correction (verbatim, not re-transcribed):
`vaslat_r`→99.4mm, `bflh_r`→97.6mm, `bfsh_r`→110.3mm, `gasmed_r`→51.0mm, `gaslat_r`→58.8mm,
`soleus_r`→44.0mm. Pennation/Fmax/tendon-slack-length asserted bit-unchanged for all 6 muscles, both
subjects (machine-checked in `build_fiber_length_corrected_model`, not just claimed).

An honest per-muscle nuance (not present for subject2, whose 6 muscles were *uniformly* too long,
ratio 0.68–0.83×): subject3/4 show a **mixed** pattern — `bflh_r`/`bfsh_r`/`soleus_r` needed only small
(≤5%) *upward* adjustments (their pre-correction values were already close to, or slightly below, the
Ward2009 target), while `vaslat_r`/`gasmed_r`/`gaslat_r` needed large (12–16%) *downward* corrections,
same direction as subject2:

| muscle | subject3 old→new (ratio) | subject4 old→new (ratio) |
|---|---:|---:|
| bflh_r | 96.76→97.60mm (1.0087) | 92.95→97.60mm (1.0500) |
| bfsh_r | 109.38→110.30mm (1.0084) | 106.65→110.30mm (1.0342) |
| soleus_r | 42.88→44.00mm (1.0261) | 44.93→44.00mm (0.9793) |
| gaslat_r | 67.32→58.80mm (0.8735) | 70.17→58.80mm (0.8379) |
| gasmed_r | 57.57→51.00mm (0.8859) | 59.98→51.00mm (0.8503) |
| vaslat_r | 115.98→99.40mm (0.8570) | 113.29→99.40mm (0.8774) |

The aggregate knee-force effect is nonetheless large and real in both subjects — the high-force
quadriceps/gastrocnemius corrections dominate the net response even though 3/6 muscles barely move (or
move the other way):

| subject | baseline self / JR %BW (ratio) | corrected self / JR %BW (ratio) | rel. drop (self / JR) | self-vs-JR agreement (corrected) |
|---|---:|---:|---:|---:|
| subject3 | 513.75 / 513.75 %BW (1.9896×/1.9896×) | 388.29 / 388.32 %BW (1.5037×/1.5038×) | 24.42% / 24.41% | 0.007% |
| subject4 | 434.14 / 434.15 %BW (1.6813×/1.6813×) | 285.19 / 306.71 %BW (1.1044×/1.1878×) | 34.31% / 29.35% | 7.02% (naive) / 0.97% (whole-trial-consistent, Sec.4) |

Both clear the 10% Q1 bar by 2.4–3.4×; neither crosses ≤1.0. `convergence_pass` (raw) reads False for
both by the same pre-existing subject2-specific 1.50s literal already diagnosed in
`docs/MECHANISM_CROSS_SUBJECT.md`; the subject-agnostic solve-quality re-derivation (frames>100, zero
NaN, zero out-of-bounds activations) and the no-degenerate-solve gate
(`muscles_over_fmax_ratio_1p5=={}`, `joint_reserve_leaning_hard=={}`) both **PASS** for both subjects.
Pelvis residual force (the pre-existing, inherited, no-RRA limitation) is **bit-identical**
before/after correction for both subjects (subject3 566.999N, subject4 340.062N — a fiber-length-only
edit does not touch the reserve-actuator/residual system at all) — the correction is not silently
trading knee-force accuracy for dynamical-consistency.

## 4. Forced adversary: subject4's knee develops the SAME boundary-truncation confound already diagnosed at the hip

Subject4's corrected-knee self-vs-JR agreement (7.02%) broke the pattern subject2/subject3 both show
(<0.01%) and subject4's own *baseline* showed (0.0013%, already published). Per this repo's
"honest-negative is not a free pass" discipline, this was **forced through OODA**, not accepted or
hand-waved:

- **Observe**: pulled the raw per-frame self-computed contact-force signal for subject4's corrected
  knee near the trial boundary (t=1.30s, subject4's real, un-extrapolated trial end).
- **Orient**: found a smooth, **monotonic 15+-frame ramp** (156.4%BW at t=1.16s → 309.7%BW at t=1.30s,
  no discontinuity, still ascending at the last captured frame) — the model's own
  edge-exclusion convention (`Fbc_mag[e:-e]`, `e=5`, the same convention every cert in this family uses
  to avoid Savitzky-Golay differentiation edge effects) reports an earlier interior local peak
  (285.19%BW at t=0.15s) instead of this true, higher boundary value. This is **structurally identical**
  to the hip truncation confound `docs/MECHANISM_CROSS_SUBJECT.md` Sec.2.3 already diagnosed for
  subject3/4's *hip* (there, present at baseline; here, newly surfaced at the *knee*, caused by the
  correction shifting when in the gait cycle the knee force peaks). Cross-checked against JointReaction
  — a structurally independent computation path — evaluated at the SAME whole-trial convention: self
  (309.67%BW) vs JR (306.71%BW) agree to **0.97%**, vs 7.02% when self is read via its edge-excluded
  convention. Two independent pipelines agreeing tightly at the boundary is the same evidence
  `MECHANISM_CROSS_SUBJECT.md` used to conclude "a genuine kinematic/GRF event, not differentiation
  noise" — reproduced here independently, not assumed.
- **Decide**: the whole-trial/JR-consistent reading (306.71%BW, ratio 1.1878×) is the trustworthy
  number for subject4's corrected knee; the edge-excluded self-computed reading (285.19%BW, ratio
  1.1044×) **understates** the true peak and should not be reported as the primary number. Per the same
  logic `MECHANISM_CROSS_SUBJECT.md` already applied to the hip, the trial cuts off while the signal is
  **still ascending** — so even 1.1878× is itself potentially a lower bound on subject4's true corrected
  knee peak.
- **Act**: recomputed the cross-subject verdict using this reading. **Robustness result: unchanged.**
  `C_GENERALIZES_DIRECTION` stays TRUE (29.35% drop via JR-consistent reading still clears the 10% bar
  by 3×) and `C_GENERALIZES_OVERCORRECTION` stays FALSE (1.1878× is, if anything, *further* from ≤1.0
  than the naive 1.1044× reading suggested — forcing this adversary made the "does not over-correct"
  conclusion *more* solid, not less).

Systematically checked all 4 runs (subject3/4 × primary/compound) for the same effect: **subject3's
knee is clean in both** (edge-excluded == whole-trial, agreement 0.007–0.010%); **subject4's knee shows
the boundary shift in both** (primary 7.02%→0.97%, compound 7.66%→0.96%). The hip shows this same
convention mismatch in **all 4 runs**, both subjects — but that is the *already-disclosed, pre-existing*
`MECHANISM_CROSS_SUBJECT.md` finding, present already at baseline, not new here. Full diagnostic:
`data/msk_smoketest/fiber_length_correction_cross_subject/boundary_truncation_diagnostic.json`.

## 5. Hip secondary/corroborating check: subject2's "worsens the hip" does NOT generalize; the joint-differentiated (knee-specific) mechanism does

| subject | baseline hip ratio (self) | corrected hip ratio (self) | rel. drop, self-computed | rel. drop, JR-consistent (Sec.4 logic applied to hip) |
|---|---:|---:|---:|---:|
| subject2 (ref) | 1.412× | 1.735× | **−22.9% (worsens)** | — |
| subject3 | 1.664× | 1.629× | +2.09% (small improve) | +0.036% (flat) |
| subject4 | 1.529× | 1.528× | +0.07% (flat) | +0.027% (flat) |

Using the JR-consistent reading (the hip's boundary-truncation confound is present identically at
baseline *and* corrected for both subjects — an inherited, pre-existing property, not new — so the
whole-trial-native JR reading is the more apples-to-apples comparison here too, exactly Sec.4's logic):
**both subjects' hips are essentially flat (<0.04% change) under the fiber-length correction.**
Subject2's specific hip-*worsening* is subject2-specific, not cohort-wide. What DOES generalize is the
**joint-differentiated mechanism**: the knee moves by 24–34% while the hip moves by <0.05% in both new
subjects — exactly the anatomical prediction stated before running (5/12 knee-crossing vs 1/25
hip-crossing prime movers corrected). Subject2's hip is the outlier requiring its own explanation
(already attributed elsewhere, `docs/MECHANISM_HIP_STRUCTURAL_MECHANISM.md`'s glute-med recruitment
issue, not re-litigated here), not subject3/4.

## 6. Bonus (beyond the task's literal ask, cheap to add — reuses already-existing Fmax-corrected models)

Fiber-length correction stacked on each subject's already-published Handsfield-Fmax-corrected model
(`data/msk_models/subject{3,4}_scaled_handsfield_fmax_corrected.osim`, built by
`fmax_pcsa_correction_test.py`/`_subject4.py`, not by this script):

| subject | baseline | Fmax-only | fiber-length-only | compound (fiber-length+Fmax) |
|---|---:|---:|---:|---:|
| subject3 knee | 1.990× | 1.786× (10.2% drop) | 1.504× | **1.318×** |
| subject4 knee | 1.681× | 1.525× (9.3% drop, sub-threshold) | 1.188× (JR-consistent) | **~1.19×** (JR-consistent; naive self-computed 1.094× carries the identical Sec.4 boundary-shift artifact, confirmed in the systematic check) |

Subject3's Fmax-only correction independently clears its own 10% bar, so compounding adds a real,
further drop (1.504×→1.318×). Subject4's Fmax-only correction does *not* clear that bar (9.3%), and
correspondingly the compound result is statistically indistinguishable from fiber-length-alone
(~1.19× either way) — Fmax adds essentially nothing further for subject4's knee. Neither subject's
compound crosses ≤1.0. Both compound solves pass the same no-degenerate-solve gate
(`muscles_over_fmax_ratio_1p5=={}`, `joint_reserve_leaning_hard=={}`) and the same trap tripwire.

## 7. Verification (machine-checked, not narrated)

- Cert-script integrity: `validate_joint_force.py`/`static_opt_knee.py`/`validate_hip_force.py` md5
  identical before/after this run to `docs/MECHANISM_CROSS_SUBJECT.md`'s own recorded hashes
  (`e2bb15958cfe8d72f44a3e0ac3d95950` / `d00f962c1a8873606874fffd65e828c8` /
  `a6514b9f081e0b4767fb4b0f94e60200`) — confirmed **before this script ran at all**, meaning subject3/4's
  cached baseline (`cross_subject_validation`) JSON numbers used for comparison are guaranteed
  non-stale. `cert_scripts_all_unmodified: true` in the output JSON.
- Trap tripwire: PASS for all 4 runs (subject3/4 × primary/compound) — `body_weight_N` matched each
  subject's own independently-known value exactly, never subject2's.
- Self-consistency gate: PASS, 0.0000% relative diff reproducing subject2's published number through
  the new orchestration path.
- Read-only respect: only new files under `data/msk_models/`, `data/msk_smoketest/{subject3,subject4}
  _walking1/{fiber_length_correction_test,fiber_length_fmax_compound_test}/`,
  `data/msk_smoketest/subject2_walking1/fiber_length_correction_cross_subject_selfcheck/` (new, Part 0
  only), and `data/msk_smoketest/fiber_length_correction_cross_subject/` were written — verified absent
  before this run (Bash `ls` check per path, logged), so nothing already on disk (incl. any concurrent
  instance's files) was clobbered. No git operations.
- All numbers in this doc are pulled from
  `data/msk_smoketest/fiber_length_correction_cross_subject/fiber_length_correction_cross_subject_results.json`
  and its `boundary_truncation_diagnostic.json` sibling, not transcribed from console prose.

## 8. Confidence tier

Cadaveric (Ward2009, muscle-architecture literals) × in-vivo (OrthoLoad, instrumented-implant peak
contact force) — same tier as subject2's own cert. **Solve-quality-caveated**: raw `convergence_pass`
reads False for both subjects by the pre-existing subject2-specific 1.50s literal (not a real
convergence problem, `docs/MECHANISM_CROSS_SUBJECT.md` Sec.4); pelvis-residual gate fails the <75N band
in both subjects (566.99N/340.06N, no RRA — an inherited, disclosed, unaffected-by-this-correction
limitation). Subject-agnostic solve-quality and no-degenerate-solve gates **both PASS** for all 4 runs.

## 9. Honest gaps

1. **n=3 subjects total** (subject2 + 2 new) — descriptive, not a population estimate, same caveat
   `MECHANISM_CROSS_SUBJECT.md` already carries.
2. **Right side, walking1 only**, same inherited scope caveat as every cert in this family.
3. **Ward2009 absolute-mm values are applied identically regardless of subject height/limb-scale** — an
   inherited simplification from subject2's own original test (not new here); a cadaveric specimen's
   own fiber length need not transfer unscaled to a differently-sized living subject, and this was not
   corrected for in either subject2's or this test.
4. **Fascicle-decoupling (Q2b) and tendon-strain (Q2a) axes were intentionally out of scope** — no
   subject3/4 gastroc/quad-ham-fascicle-kinematics baseline exists to diff against, and the task's own
   ask was the knee contact-force ratio specifically. Not silently dropped — named here.
5. **Subject4's true corrected knee peak (and, by the same logic, both subjects' hip peaks) may still
   be under-captured** — the boundary-truncation diagnostic (Sec.4) shows the signal is still ascending
   at the trial's real end; the reported 1.188× for subject4's knee is, if anything, a lower bound.
6. **Root cause of the residual (post-correction) gap is not identified** — same open item
   `MECHANISM_CROSS_SUBJECT.md` Sec.6 already named (Fmax not subject-rescaled; subject3/4 are ~15–16kg
   lighter than subject2 but use the same absolute muscle-strength values, i.e. a *higher*
   strength-to-bodyweight ratio, yet still over-predict *more* — an unresolved, candidate, untested
   explanation for why subject2 achieves full closure and subject3/4 do not).
7. **Compound-bonus hip numbers were not separately re-derived with the JR-consistent reading** (Sec.6)
   — flagged as carrying the same Sec.4/5 caveat, not quantified, out of lean scope for a bonus
   measurement.

## 10. Files

- `scripts/msk/fiber_length_correction_cross_subject.py` (new) — orchestration + comparison script.
- `docs/MECHANISM_FIBER_LENGTH_CROSS_SUBJECT.md` — this doc.
- `data/msk_smoketest/fiber_length_correction_cross_subject/fiber_length_correction_cross_subject_results.json`
  — primary machine-readable results (self-consistency gate, per-subject trap tripwire/gates/comparison,
  compound bonus, cross-subject verdict, robustness-check-under-convention-choice).
- `data/msk_smoketest/fiber_length_correction_cross_subject/boundary_truncation_diagnostic.json` —
  supplementary Sec.4 diagnostic (edge-excluded vs whole-trial vs JR readings for all 4 runs, knee+hip).
- `data/msk_models/subject{3,4}_scaled_fiberlength_ward_corrected.osim` (new, primary test),
  `subject{3,4}_scaled_fiberlength_and_fmax_corrected.osim` (new, compound bonus test).
- `data/msk_smoketest/subject{3,4}_walking1/fiber_length_correction_test/{so,jr}/`,
  `.../fiber_length_fmax_compound_test/{so,jr}/` — fresh SO/JR outputs.
- `data/msk_smoketest/subject2_walking1/fiber_length_correction_cross_subject_selfcheck/{so,jr}/` — Part
  0 self-consistency-gate SO/JR output (new directory; the existing
  `subject2_scaled_fiberlength_ward_corrected.osim` model was read only, never rebuilt).
- **Not touched** (md5-verified, Sec.7): `scripts/msk/validate_joint_force.py`,
  `scripts/msk/static_opt_knee.py`, `scripts/msk/validate_hip_force.py`,
  `scripts/msk/subject_specific_scaling.py`, `scripts/msk/cross_subject_validation.py`,
  `scripts/msk/fiber_length_correction_test.py`.
- No git commit, no git push (isolation respected).
