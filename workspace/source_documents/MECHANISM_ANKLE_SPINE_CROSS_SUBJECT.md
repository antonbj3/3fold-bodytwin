# MECHANISM ANKLE + SPINE CROSS-SUBJECT — do the ankle/spine force certs generalize beyond subject2?

**Question.** `docs/MECHANISM_ANKLE_FORCE.md` and `docs/MECHANISM_SPINE_GAIT_VBR.md` were each validated on
exactly **one subject** (subject2/walking1) — the same single-subject gap `docs/MECHANISM_CROSS_SUBJECT.md`
already closed for knee/hip (finding a **1.412×–1.990×** over-prediction range vs. subject2's own 1.412×/1.515×).
This doc runs the **same** ankle and spine free-body machinery on **subject3** (63.5 kg/1.69 m/F) and
**subject4** (62.6 kg/1.68 m/F) — the identical two bodies the knee/hip cross-subject test already used — and
asks the task's own two falsifiers directly: does the ankle stay ~1.0× the published cadaveric anchor across
subjects, or diverge like knee/hip did? Does the spine's pure-kinematic-vs-muscle-subtraction split hold?

## Headline answer

**ANKLE generalizes, tightly — more tightly than subject2 alone suggested, and much more tightly than
knee/hip.** Self-computed/JointReaction ratio vs. the published cadaveric/literature midpoint (476.67 %BW):
subject2 **1.016×**, subject3 **1.001×**, subject4 **0.985×** — a **0.985×–1.016×** range, a 3.1-percentage-point
spread across two bodies that differ from subject2 by −15/−16 kg, −27/−28 cm, and opposite sex. This is the
same diverse instance-space the knee/hip test used, where the over-prediction ratio ranged 1.412×–1.990× (a
58-percentage-point spread) — the ankle's near-1.0× match is **categorically more stable across subjects**
than the knee/hip over-prediction magnitude is. Self-vs-JointReaction internal agreement stayed **<0.02%** in
every subject (subject2 0.0014%, subject3 0.0113%, subject4 0.0063%) — the same tight, structural agreement
the knee/hip test used as its own decisive anti-artifact check, reproduced here.

**SPINE: the pure-kinematic half of the split generalizes and is now independently cross-validated by
JointReaction in all 3 subjects (a check subject2's own cert could NOT do, because JointReaction was
numerically unusable on subject2's ligament-augmented substrate — Tier-1 on the plain base model does not
have that problem, verified). The muscle-subtraction half — the part that actually produces the over-prediction
finding subject2's spine cert reported (55.97 %BW → 210.22 %BW, a further ~3.75× jump) — is BLOCKED for
subject3/4, not silently skipped: the ligament+erector-spinae+trunk-flexor substrate
(`model_with_spine_ligaments.osim`) is built by grafting ~8 sequential, SUBJECT2-SPECIFIC model merges, and its
own build script hardcodes subject2's model path and subject2's own scaled local coordinates (§3). Building a
subject3/4 equivalent is a new model-SCALING task, not a "run the same machinery" task, and was correctly out
of this task's scope — diagnosed with direct, machine-checked evidence (§3), not assumed or hand-waved.**

**CONFIDENCE TIER**: ankle = **cadaveric/model-literature** (unchanged from subject2's own honest tier — no
in-vivo ankle program exists); spine Tier-1 = **in-vivo-VBR-anchored but muscle-subtraction-incomplete**; both
**solve-quality-caveated** for subject3/4 (§6) — same caveat class the task itself pre-registered, plus one
additional, newly-diagnosed caveat found and traced to its physical cause, not just flagged (§6.2).

---

## 1. Method (reuse discipline)

A **new** orchestration script, `scripts/msk/ankle_spine_cross_subject.py`, is the only new code. It imports
`validate_joint_force.py` (vjf), `static_opt_knee.py` (sok), `validate_ankle_force.py` (vaf),
`validate_spine_gait_force.py` (vsgf), and `cross_subject_validation.py` (csval) **unedited** — md5-identical
before and after every run (§5) — and calls their existing, parametrized functions directly
(`run_combined_frame_loop`, `reaction_and_contact_force`, `pure_reaction_peak_for_window`,
`extract_jr_force_pct_bw`, `classify_vs_published`, `check_so_convergence_and_sanity`, `torso_axis`,
`check_no_direct_external_force_on_above_cut`, `get_descendant_bodies`, `parse_mot`, `savgol_smooth_and_derivs`,
`configure_for_subject`). Two things this script does are genuinely new — the same class of "new work"
`cross_subject_validation.py`'s own docstring already licenses (orchestration + comparison logic) — because two
real, forced findings meant neither `vaf.main()` nor `vsgf.main()` could be called as a whole:

### 1.1 The task's pre-registered "def-time-default bw_n binding" trap — confirmed present, one level up

The exact trap (`joint_force_mag_pct_bw(reaction_sto, colprefix, bw_n=BW_N)`, a default bound once at
definition time) lives in `contact_waveform_analysis.py` (a knee/hip waveform tool), not in `vaf.py`/`vsgf.py`
directly — verified by reading: their own `bw_n` is always freshly recomputed from
`total_mass = sum(model's own body masses); bw_n = total_mass*G`, never a function default. But `vaf.py` has
the **same disease at a different site**: its module-level `MODEL_FILE`/`IK_MOT`/`GRF_MOT`/`SO_RESERVE_XML`
(its own lines ~127–139) are copied **once**, at `vaf`'s own import time, from `vjf`/`sok` into `vaf`'s
independent module dict (`MODEL_FILE = vjf.MODEL_FILE`) — repointing `vjf.MODEL_FILE`/`sok.MODEL_FILE`
afterward does **not** retroactively change `vaf.MODEL_FILE` (a bare name inside `vaf`'s own functions resolves
against `vaf`'s own `__dict__`). This script never relies on `vaf`'s copies: it builds the `osim.Model`/IK/GRF/
SO/JR objects itself from paths derived exactly as `cross_subject_validation.configure_for_subject` already
does, and passes them explicitly into `vaf`'s reusable functions (verified by reading: none of the six listed
above reads a bare subject-identifying module global — all take model/paths/data as explicit arguments).

The trap is additionally given a **live, decisive, machine-checked tripwire**, not just careful code: every
subject run asserts its own freshly-computed `body_weight_N` equals that subject's independently-known value
and is **not** close to subject2's —

| subject | computed body_weight_N (this run) | independently expected | tripwire |
|---|---:|---:|---|
| subject2 | 766.88003 | 766.88003 | n/a (reference) |
| subject3 | 622.7222750000001 | 622.7222750000001 | PASS — not within 50N of subject2's 766.88N |
| subject4 | 613.89629 | 613.89629 | PASS — not within 50N of subject2's 766.88N |

Had the trap fired, subject3/4's runs would have silently computed with subject2's 78.2 kg/766.88 N instead of
their own 63.5 kg/622.72 N or 62.6 kg/613.90 N — this assertion is the exact, direct falsifier of that failure
mode, and it never fired (the script raises `AssertionError` and aborts loudly if it does; it did not).

### 1.2 `vaf.main()`/`vsgf.main()` are not callable as a whole for subject3/4 — verified, not assumed

`vaf.main()`'s STEP 2 hard-stops (`return 1`, before computing **any** ankle number) unless
`sok.check_so_convergence_and_sanity(...)["convergence_pass"]` is `True`. That gate's own formula
(`static_opt_knee.py` L297–300) is `n_frames>100 AND covers_validated_peak_time_1p50s AND n_nan_act==0 AND
out_of_bounds==0`, where `covers_validated_peak_time_1p50s` is a hardcoded subject2-specific literal (does the
trial cover t=1.50s) that **fails by construction** for subject3 (trial ends 1.31s) and subject4 (ends 1.30s) —
the exact fact `docs/MECHANISM_CROSS_SUBJECT.md` already diagnosed for the knee/hip gate, and the reason that
script also never calls `vjf`/`sok`/`vhf`'s own `main()`s wholesale. This script follows the identical
precedent: call the lower-level functions directly, and replace the hard stop with the same "raw vs
subject-agnostic" dual reporting `cross_subject_validation.py` established
(`solve_quality_pass = n_frames>100 AND no NaN AND in-bounds`, dropping only the 1.50s literal — same formula,
same file, L318–321).

### 1.3 The spine's Tier-2 blocker — forced with direct evidence, not assumed

`add_spine_ligaments.py`'s own header hardcodes `UNIFIED_MODEL_PATH = data/msk_models/subject2_unified.osim`
and states (L137) *"Points are converted into LaiArnold subject2's own scaled pelvis/torso local
coordinates"* — a manual, subject-specific derivation, not a generic/parametrized rescaling. Verified live this
session: subject3's own base model
(`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` — the **same** model already used for
its existing knee/hip/ankle SO/JR) contains **zero** muscles matching
`erector|multifidus|iliocostalis|longissimus|oblique|rectus_abd` (`grep -ci` on the live `.osim` XML = 0).
`find data -iname "*subject3*" -o -iname "*subject4*"` (excluding walking1 trial dirs) turns up only
`subject{3,4}_scaled_handsfield_fmax_corrected.osim` — no unified/ligament/erector-spinae equivalent exists.
There is no subject3/4 analogue of `model_with_spine_ligaments.osim`, and building one (re-deriving ~8 chained,
subject-specific attachment-point conversions) is a new model-scaling **build** task, explicitly out of this
task's "run the same unmodified machinery" **reuse** scope. What *is* reusable unchanged is Tier-1 (pure
kinematic "back"-cut reaction, SO-independent by construction — `vsgf.py`'s own module docstring:
*"R = sum_i(m_i a_i) − gravity_above"*), which needs only the base model's mass distribution + that subject's
own IK + GRF — all already available for subject3/4. `run_spine_tier1()` re-derives `vsgf.py main()`'s own
inlined Tier-1 formula (L510–524 there, not factored into its own function) verbatim, reusing
`vsgf.torso_axis`/`vsgf.check_no_direct_external_force_on_above_cut`/`vjf.get_descendant_bodies`/
`vjf.parse_mot`/`vjf.savgol_smooth_and_derivs` unchanged.

**A cheaper Tier-2 proxy was considered and rejected, not silently skipped**: driving subject2's
already-built augmented model with subject3/4's own IK/GRF (rather than building a real subject3/4-scaled
substrate) would have been much cheaper, but was killed before running it — it would silently conflate
subject3/4's kinematics with subject2's anthropometry (mass distribution, segment lengths, muscle Fmax,
moment arms all remain subject2's), producing a number that is not attributable to subject3/4 in any
meaningful sense while superficially looking like a subject3/4 result. Rejected as invalid, not attempted.

## 2. Self-consistency gate (run FIRST, gates everything else — task's own explicit requirement)

Both new code paths were required to reproduce subject2's own **already-published, machine-recorded** numbers
(read from the on-disk results JSON, never from doc prose) before subject3/4 were trusted:

| quantity | mine (new code, subject2 inputs) | already-published (subject2) | rel. diff | verdict |
|---|---:|---:|---:|---|
| Ankle pure reaction | 107.5034 %BW | 107.5034 %BW | 0.00000% | PASS |
| Ankle self-computed | 484.4865 %BW | 484.4865 %BW | 0.00000% | PASS |
| Ankle official JointReaction | 484.4798 %BW | 484.4798 %BW | 0.00000% | PASS |
| Achilles scalar-sum | 320.2649 %BW | 320.2649 %BW | 0.00000% | PASS |
| Achilles vector-sum | 318.8285 %BW | 318.8285 %BW | 0.00000% | PASS |
| Spine Tier-1 (augmented model, the exact inputs the 55.97%BW number came from) | 55.96898 %BW @ t=0.68s | 55.96898 %BW @ t=0.68s | 0.00000% | PASS |

Bit-identical reproduction on all six quantities. **Understood as necessary but not sufficient alone**: a
subject-aliasing bug (always reading subject2's data) would trivially "pass" a subject2-vs-subject2 check too —
the decisive test for *that* specific adversary is §1.1's tripwire, which only fires on an actual subject3/4
run and did not.

## 3. ANKLE results (all subjects)

| subject | mass/height/sex | body weight (N) | self-computed peak (%BW) | JointReaction peak (%BW) | self-vs-JR agreement | ratio vs. published midpoint (476.67 %BW) |
|---|---|---:|---:|---:|---:|---:|
| subject2 | 78.2 kg / 1.96 m / M | 766.88 | 484.49 @ t=0.58s | 484.48 | 0.0014% | **1.016×** |
| subject3 | 63.5 kg / 1.69 m / F | 622.72 | 477.26 @ t=0.54s | 477.20 | 0.0113% | **1.001×** |
| subject4 | 62.6 kg / 1.68 m / F | 613.90 | 469.44 @ t=0.52s | 469.41 | 0.0063% | **0.985×** |

| subject | Achilles/triceps-surae scalar-sum peak (%BW) | ratio vs. Giddings 2000 (390 %BW) |
|---|---:|---:|
| subject2 | 320.26 @ t=0.59s | 0.821 |
| subject3 | 366.19 @ t=0.54s | 0.939 |
| subject4 | 340.87 @ t=0.53s | 0.874 |

Both candidates (self-computed, official JointReaction) land `PASS-plausible` inside the pre-registered
(0.30, 3.00) regime band in every subject (reused unchanged from `vaf.py`, never re-tuned for this task).
Anatomical crossing-muscle anchor: **11/11** expected ankle-crossing muscles detected, **100%** per-frame
coverage for the triceps-surae in every subject (0 missing, 0 dropout frames — same clean result as subject2).

**Free reproducibility gate** (this script's own re-derivation of the KNEE-r cut, pointed at each subject's OWN
already-published knee cert numbers from `cross_subject_validation.py`, not subject2's): **0.000% relative
difference on all three metrics (pure-reaction/self-computed/JointReaction), for both subject3 and subject4** —
an exact match, the strongest evidence this script's machinery (BFS, crossing-muscle detector, Newton
computation) is operating correctly on new subjects' data, fully decorrelated from whether the ankle-specific
literature anchor itself is right.

## 4. SPINE Tier-1 results (all subjects) — Tier-2 not computable, see §1.3

| subject | model | body weight (N) | Tier-1 interior peak (%BW) | JointReaction "back", SAME INSTANT as Tier-1 peak (%BW) | rel. diff | ratio vs. OrthoLoad level-walking median (53.37 %BW) |
|---|---|---:|---:|---:|---:|---:|
| subject2 | augmented (`model_with_spine_ligaments.osim`, = the already-published 55.97%BW) | 774.49 | 55.969 @ t=0.68s | n/a — JointReaction is numerically unusable on this substrate (subject2's own cert, §9: 61-66% NaN/overflow frames) | — | 1.049× |
| subject2 | **base model** (bonus, decorrelated comparison — what subject3/4 must use) | 766.88 | 55.283 @ t=0.68s | 55.679 | 0.716% | 1.036× |
| subject3 | base model | 622.72 | 63.557 @ t=0.61s | 60.015 | 5.572% | 1.191× |
| subject4 | base model | 613.90 | 58.158 @ t=1.12s | 57.005 | 1.982% | 1.090× |

The base-model detour incidentally produces something subject2's own spine cert never had: a **working**
JointReaction cross-check on Tier-1 (the augmented substrate's JointReaction was too numerically unstable to
use at all, §9 of `MECHANISM_SPINE_GAIT_VBR.md`) — here it agrees with the independently-computed, finite-
difference-based Tier-1 to **0.7–5.6%** in all three subjects, a genuine, decorrelated (different numerical
method: OpenSim's own internal reaction computation vs. this script's external Savitzky-Golay COM
differentiation) confirmation that Tier-1 is being computed correctly across subjects. The naive whole-trial
JointReaction peak (not shown as primary above) lands at a **different instant** than Tier-1's own interior
peak in subject3 (t=0.06s vs 0.61s) and subject4 (t=1.10s vs 1.12s) — checked directly rather than compared
naively, mirroring the same "compare at the same instant" discipline `docs/MECHANISM_CROSS_SUBJECT.md` §2.3
already established for the hip; the same-instant column above is the one load-bearing for the verdict.

**Tier-1 magnitude is consistent across subjects (55.3–63.6 %BW) and stays close to — modestly above, never
wildly above — the OrthoLoad level-walking median in every subject (1.04×–1.19×)**, the same near-1.0× pattern
subject2's own Tier-1 showed (1.05×). This is the **pure-kinematic half** of the task's falsifier, and it
generalizes cleanly. The **muscle-subtraction half** — subject2's own Tier-2b jumped a further **~3.75×** past
Tier-1 (55.97→210.22 %BW) once co-contraction was subtracted, which is what actually produced the "spine
over-predicts VBR by 3.9–5.4×" finding — is the part this test **cannot** reproduce for subject3/4 (§1.3): there
is no subject-specific augmented substrate to run it on. **Whether that further multiplication would also occur
in subject3/4 is open, not resolved either way** — reported honestly as a blocked, diagnosed gap, not silently
dropped and not guessed at.

## 5. Verification (machine-checked, not narrated)

- **md5 of all 5 reused files, before vs. after every run, identical**:
  `validate_joint_force.py=e2bb15958cfe8d72f44a3e0ac3d95950`,
  `static_opt_knee.py=d00f962c1a8873606874fffd65e828c8`,
  `validate_ankle_force.py=f381b439f032d9c9407044a10a17a4b4`,
  `validate_spine_gait_force.py=fd430376f9fa5385a8dfa01b449e02ef`,
  `cross_subject_validation.py=46c911a2a22cc81b6b636d951d612e76`. Zero bytes edited.
- All numbers in §3/§4 pulled directly from the on-disk JSON files (§7), not transcribed from console prose.
- Determinism: re-run after adding the same-instant JR comparison (a pure-addition code change) reproduced
  every pre-existing number bit-for-bit (self-consistency gate table unchanged, subject3/4 ankle numbers
  unchanged) — confirmed via diff of the two runs' logs, not assumed from "should be deterministic."

## 6. Solve-quality caveats — both the pre-registered one and one found by forcing a diagnosis

### 6.1 Pre-registered: raw SO convergence gate and pelvis residual (same as knee/hip cross-subject)

| gate | subject2 | subject3 | subject4 |
|---|---:|---:|---:|
| raw `convergence_pass` (incl. the 1.50s literal) | — | **False** | **False** |
| subject-agnostic solve-quality (frames>100, no NaN, in-bounds) | — | **True** | **True** |
| pelvis residual force peak (N), band <75N | 173.0 | **567.0** | **340.1** |

Identical numbers to `docs/MECHANISM_CROSS_SUBJECT.md` §4 (567.0N/340.1N) — confirms this script is reading the
exact same underlying SO solves, not a fresh/different run. Both are pre-existing, already-disclosed
limitations of this cert family (no RRA), not new to this test.

### 6.2 Newly found: a whole-body Newton residual gate reads elevated for subject3/4 — forced to a real diagnosis, not left as an unexplained "FAIL"

The ankle script's whole-body (all-bodies) Newton residual gate (<8 %BW, an internal-validity check on the
differentiation pipeline, separate from the ankle-specific self-vs-JR agreement) reads **17.87 %BW** (subject3,
vertical component) and **8.34 %BW** (subject4) — both above threshold, vs. subject2's own clean 4.10 %BW.
**Not accepted as an unexplained "solve-quality caveat" — forced to Orient**: printing the per-frame vertical
GRF sum (`R_ground_force_vy + L_ground_force_vy`) near the worst residual frame (t≈0.05s, in **all three**
subjects, including subject2) shows the measured total vertical ground force is implausibly **below** body
weight at the very start of the captured trial and climbs smoothly, monotonically toward the true
support-level over the next ~0.08s (subject3: 131 N → 482 N over t=0.01–0.09s, on a 623 N body) — the classic
signature of a **trial-onset force-plate-coverage transient** (the subject's foot(s) not yet fully registering
on the instrumented plates at the very first captured frames), not differentiation noise (noise would not decay
smoothly and monotonically) and not a bug in this script (the free reproducibility gate, §3, matched the
already-published knee numbers to 0.000%, which a differentiation bug would also have corrupted). Subject2's
own trial shows the identical-shaped transient at the identical instant (t=0.05s), just smaller in relative
terms (their trial happens to start with more of the true support force already registering) — the same
mechanism, not a new one specific to subject3/4, just more visible on their shorter trials. **This does not
affect the reported ankle numbers**: their peaks occur at t=0.52–0.58s, ~0.45–0.53s away from this onset
window, and the ankle-specific self-vs-JR agreement (which uses only the right-side GRF and only the
ankle-distal free body, not the whole-body sum) stayed <0.02% in every subject regardless.

## 7. Files

- `scripts/msk/ankle_spine_cross_subject.py` (new) — the orchestration script described in §1; the only new
  code this task wrote. Re-runnable; imports `validate_joint_force.py`/`static_opt_knee.py`/
  `validate_ankle_force.py`/`validate_spine_gait_force.py`/`cross_subject_validation.py` unmodified (md5-verified,
  §5).
- `data/msk_smoketest/subject3_walking1/ankle_force_validation/ankle_force_validation_results.json`,
  `data/msk_smoketest/subject4_walking1/ankle_force_validation/ankle_force_validation_results.json` — full
  ankle results, source for every ankle number in §3.
- `data/msk_smoketest/subject3_walking1/spine_gait_force_tier1/spine_gait_force_tier1_results.json`,
  `data/msk_smoketest/subject4_walking1/spine_gait_force_tier1/spine_gait_force_tier1_results.json` — full
  spine Tier-1 results, source for every spine number in §4.
- `data/msk_smoketest/ankle_spine_cross_subject/ankle_spine_cross_subject_results.json` — combined summary:
  self-consistency gate detail (§2), both subjects' full ankle+spine results, md5 before/after.
- Reused, read-only, unmodified: `scripts/msk/validate_joint_force.py`, `scripts/msk/static_opt_knee.py`,
  `scripts/msk/validate_ankle_force.py`, `scripts/msk/validate_spine_gait_force.py`,
  `scripts/msk/cross_subject_validation.py`; subject3/4's existing
  `data/msk_smoketest/subject{3,4}_walking1/static_optimization/{so,jr}/` (never re-run, only read); subject2's
  `data/msk_smoketest/subject2_walking1/ankle_force_validation/`, `.../spine_gait_force/`,
  `data/msk_smoketest/spine_ligaments/model_with_spine_ligaments.osim` (read-only, for the self-consistency
  gate only).
- Compared against: `docs/MECHANISM_ANKLE_FORCE.md`, `docs/MECHANISM_SPINE_GAIT_VBR.md` (the subject2-only
  certs this doc extends), `docs/MECHANISM_CROSS_SUBJECT.md` (the knee/hip precedent this doc's method and
  per-subject knee numbers reuse directly).
- No git commit, no git push performed (isolation respected).

## 8. Honest gaps

1. **Spine Tier-2 (the muscle+ligament-subtracted PRIMARY headline, the part that actually drives the
   over-prediction claim) is not tested for subject3/4** — a genuinely blocked, diagnosed gap (§1.3), not a
   silent skip. Building the missing subject-scaled ligament substrate is real future work, out of this task's
   reuse-only scope.
2. **n=3 subjects (ankle), n=3 subjects (spine Tier-1)** — descriptive, not a population estimate, same
   caveat the knee/hip cross-subject doc already stated for its own n=3.
3. **Right side, walking1 only** — same scope caveat inherited from every cert in this family.
4. **Whole-body Newton residual gate fails for subject3/4** (§6.2) — diagnosed to a real, physical,
   trial-onset force-plate-coverage transient rather than left as an unexplained caveat, but not eliminated
   (would require either a longer capture margin at trial start or an adaptive edge-exclusion window, neither
   attempted here — bounded scope).
5. **Pelvis residual force gate fails in all 3 subjects** (§6.1) — pre-existing, already-disclosed (no RRA
   run), not new here.
6. **Spine Tier-1 same-instant JointReaction agreement is weaker for subject3 (5.57%) than subject2/subject4
   (0.72%/1.98%)** — reported, not smoothed over; plausible contributors (shorter trial, different underlying
   OpenSim state-derivative behavior near a faster/different gait cycle) are not disentangled here.
7. **Static Optimization structural limitation** (inherited, unchanged): effort-minimizing, not measured EMG;
   real co-contraction can exceed it. Same for every joint in this cert family.
8. **Savitzky-Golay differentiation** (11-sample/110ms window) is a design choice, shown stable under a
   window sweep for the ankle (§3, reused unchanged from `vaf.py`) but not independently re-swept for spine
   Tier-1 in this task.
