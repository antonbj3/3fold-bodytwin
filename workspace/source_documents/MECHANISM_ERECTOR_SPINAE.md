# MECHANISM ERECTOR SPINAE — closing the lumbar-extensor muscle gap (2026-07-21)

Answers `docs/MECHANISM_MUSCLE_AUDIT.md` Anomaly #1: LaiArnoldModified2017 has a complete lower-limb
posterior chain (16/16 gastroc/soleus/hamstrings/glutes) but **erector spinae / multifidus /
longissimus / iliocostalis are ABSENT as muscles** — the trunk is driven only by 3 ideal 10 N·m
`CoordinateActuator`s (`lumbar_ext`, `lumbar_bend`, `lumbar_rot`). This is the operator's core-interest
gap: the posterior chain loaded from the ground up **includes the spinal erectors**, and a 10 N·m ideal
motor cannot supply a realistic trunk-extension moment. Executed via `scripts/msk/add_erector_spinae.py`,
run this session with `.venv-msk`'s OpenSim 4.6 (`opensim.GetVersion() == '4.6-2026-06-22-85aaf64'`).
Isolation respected: the working scaled model was **never modified** — a copy was made first
(`data/msk_models/`); both external NTFS drives were read-only inputs (mtime/hash-verified unchanged,
see §5); no git commit, no git push.

## 1. Verified source (checked live, not recalled)

**`Model_Pose2Sim_muscles_flex.osim`**, bundled with the Pose2Sim Python package already on disk at
`/media/anton/8838D60F38D5FBDE/mechanism_data/video2kin/envs/opensim_env/lib/python3.10/site-packages/Pose2Sim/OpenSim_Setup/`
(also present, byte-identical role, under the sibling `venv_pose2sim` copy). Its own `<credits>` tag,
grepped directly from the raw XML:

> "Full-body model adapted from Beaucage-Gauvreau et al., with Rajagopal et al. knee angles. Contact
> spheres and Muscles included."

This is the operator's own named candidate ("Actlab/Beaucage-Gauvreau lifting model"). Live citation
lookup this session (WebSearch was independently exhausted — a pre-existing constraint already logged in
`docs/MECHANISM_MSK_BUILD_PLAN.md` §6 — so `WebFetch` was used instead, and the DOI was independently
confirmed to resolve to a live publisher page, not just a search snippet):

> Beaucage-Gauvreau E, Robert-Lachaine X, Brandon SCE, Larivière C (2019). "Validation of an OpenSim
> full-body model with detailed lumbar spine for estimating lower lumbar spine loads during symmetric
> and asymmetric lifting tasks." *Computational Methods in Biomechanics and Biomedical Engineering*,
> 22(5):451-464. DOI [10.1080/10255842.2018.1564819](https://doi.org/10.1080/10255842.2018.1564819)
> (redirects live to tandfonline.com — verified reachable this session). PMID 30714401.

**Honest attribution gap:** the donor's own lumbar-spine architecture (5 vertebra bodies + sacrum, 6
intervertebral joints × 3 DOF = 18 DOF, `L1_T12`...`L5_S1` naming, `IL_`/`LTpL_`/`LTpT_`/`MF_` fascicle
prefixes) structurally matches the well-known Christophy, Faruk Senan, Lotz & Delp (2012) lumbar-spine
model lineage that Beaucage-Gauvreau's paper builds on — this specific sub-attribution is a **structural
inference from the on-disk model this session**, not independently re-verified against Christophy's own
paper/DOI. Flagged, not asserted as confirmed.

**No lumbar-muscle `.osim` was found anywhere else on either external drive** (395 `.osim` files total
searched; content-grepped for `erector|longissimus|multifidus|iliocostalis|spinae|semispinalis` — the
only other hit, `Neck3dof_point_constraint.osim`, is a **cervical**-spine test fixture with
iliocostalis-cervicis/longissimus-capitis neck fascicles, anatomically wrong region, not usable here).

## 2. The geometry problem (why this isn't a copy-paste) and its fix

The donor has 5 separate lumbar vertebra bodies (`lumbar1`..`lumbar5`) + `sacrum`. The target
LaiArnold model collapses the **entire** lumbar spine into **one** joint (`back`, pelvis→torso, 3 DOF:
`lumbar_extension`/`lumbar_bending`/`lumbar_rotation`). Real fascicle attachment points on `lumbar1-5`
have no body to live on in the target model.

**Fix — rigid kinematic-order reduction.** At each model's own default (neutral standing) pose, every
donor path point is re-expressed in the donor's **pelvis-local frame** via forward kinematics. This frame
is verified **identical** between donor and target (pelvis mass = 11.777 kg and mass-center =
`(-0.0707, 0, 0)` match to displayed precision in both — confirming shared Rajagopal-lineage skeletal
definition, checked live not assumed; a hard `assert` in the script aborts the whole build if this
cross-check ever fails on a different donor/target pair). Each point is then assigned to **either** the
target's `pelvis` **or** `torso` body by comparing its (unscaled) pelvis-frame height to the target's own
`back`-joint pelvis-offset height (0.0815 m) — i.e. exactly the same single-hinge simplification the
target model already uses for its 3 ideal actuators. Points assigned to torso are converted from donor
pelvis-frame into **target torso-frame** via the target's own verified pelvis→torso transform (not the
donor's torso frame — donor and target torso bodies do **not** share an origin convention; falsified
live: donor torso mass-center `(0.016, 0.180, 0)` vs target `(-0.03, 0.32, 0)`, clearly different origin
placement, so this conversion path is mandatory, not optional).

Muscles whose points **all** land on the same target body (short, single-level fascicles, e.g. a real
`lumbar4→lumbar5` multifidus fiber) have **zero moment arm** about the target's single hinge by
construction — excluded and logged, not silently dropped (this is genuine, disclosed information loss
from the reduction, the same class of thing as any coarser model omitting a finer DOF).

## 3. Two bugs found and fixed live (OODA, not a one-shot)

A first working version produced a plausible-looking 51.5 N·m ceiling — comfortably above the 10 N·m
placeholder, an easy place to stop. Per the "honest-negative/honest-result is not a free pass" rule, this
was forced further: the naive geometry-only ceiling (Σ Fmax·|momentArm|, ignoring force-length effects)
was 236 N·m, so the equilibrated result was only **22% of that ceiling** — too large a gap to accept
without a mechanism.

1. **Force-length mismatch (root-caused, not guessed):** per-muscle diagnostic showed fiber lengths
   badly off `optimal_fiber_length` (fiber/Lopt down to 0.44, one muscle — `LTpL_L5_r` — at literally
   **0%** active force). Cause: naively scaling `optimal_fiber_length`/`tendon_slack_length` by a
   per-muscle path-length ratio blows up for fascicles whose two real attachment points straddle the
   pelvis/torso split boundary (a genuinely short 1-level fascicle gets stretched across the model's much
   larger single-hinge gap). **Fix:** directly solve `tendon_slack_length` so `fiber_length ≈
   optimal_fiber_length` at the neutral reference pose (standard OpenSim tendon-slack-length tuning
   practice), with the length ratio itself clipped to the `[0.5, 2.0]` band already pre-registered in
   `MECHANISM_MUSCLE_AUDIT.md` §3 for this model family's own native muscles.
2. **Mixed-frame length bug (caught by the fix's own cross-check):** the fixed version still used a
   Euclidean distance computed directly between a pelvis-local point and a torso-local point — invalid,
   since those are two *different* frames (offset + rotated). Caught by comparing the manual calculation
   against the rebuilt model's own `muscle.getLength(state)` for a sample muscle (0.1454 m manual vs
   0.0576 m real — a 2.5× discrepancy). **Fix:** re-express every point in one common frame before
   computing any length, using the target *scaled* model's own pelvis↔torso transform.

After both fixes: **0/88** muscles need length-ratio clipping (confirming the earlier clipping was
compensating for the frame bug, not a real geometric distortion), force/Fmax fraction at the neutral pose
is **0.929–1.000 across all 88 muscles (mean 0.969)**, and the equilibrated total is **96.6% of the
geometry-only ceiling** — internally consistent, not just a bigger number.

## 4. What got added

| group | total in donor | kept (nonzero moment arm) | excluded (degenerate) |
|---|---:|---:|---:|
| Iliocostalis lumborum (`IL_`) | 24 | 22 | 2 |
| Longissimus thoracis pars lumborum (`LTpL_`) | 10 | 10 | 0 |
| Longissimus thoracis pars thoracis (`LTpT_`) | 42 | 28 | 14 |
| Multifidus (`MF_`) | 50 | 28 | 22 |
| **Total** | **126** | **88** | **38** |

Retention is lowest for multifidus (56%) — expected and anatomically sane: real multifidus is the
shortest, most segmentally-organized posterior muscle, so it's the most likely to be fully contained
within one side of a single-hinge reduction. Quadratus lumborum, psoas (multi-level), obliques, rectus
abdominis, and latissimus dorsi are all present in the donor too but **out of scope** (task was
explicitly "longissimus + iliocostalis; + multifidus if feasible").

All 88 added as `Millard2012EquilibriumMuscle` (per the task's request), Hill-type parameters from the
verified source: `max_isometric_force` and `pennation_angle_at_optimal` copied **unscaled** from the
donor (matching this model family's own measured, established convention —
`MECHANISM_MUSCLE_AUDIT.md` §3 found Fmax/pennation are *never* subject-rescaled for the native 80
muscles, bit-identical generic-vs-scaled); `optimal_fiber_length`/`tendon_slack_length` rescaled per the
solve in §3. Per-body scale factors used for geometry are **measured live from this model's own existing
muscles/mass properties**, never hardcoded: pelvis `(1.0429, 1.0429, 1.2780)` from 60+ existing
muscle-path-points (std ≈ 0, clean signal); torso `(0.9879, 0.9879, 0.9879-assumed)` from mass-center
ratio (torso has zero native muscle points to cross-check against — see honest gaps).

## 5. Verification

- **Source file untouched:** `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` mtime
  `2021-12-06 10:49:36` and md5 `0eb06c2a2b235b579b76855c8c8bff88` — unchanged (matches the date already
  logged in `docs/MECHANISM_MSK_ENV.md`).
- **Copy loads:** `data/msk_models/LaiArnoldModified2017_erector_spinae_subject2_scaled.osim` (1,365,569
  bytes) — `opensim.Model(path); model.initSystem()` → **exit 0**, reloaded fresh from disk (not just the
  in-memory build object).
- **Structural integrity of the original 80 muscles + 13 `CoordinateActuator`s (supplemented, not
  replaced — `lumbar_ext/bend/rot` still present at 10 N·m each):** 0/80 Fmax mismatches (bit-identical),
  13/13 `CoordinateActuator`s present unchanged. Muscle count 80→168 (+88), actuator count 93→181 (+88) —
  exact arithmetic match.
- **Moment-arm sign consistency (a real anatomical prediction, not tuned for):** all 88 ported muscles
  show a **positive** `computeMomentArm` sign about `lumbar_extension` at both poses tested (88/88 pos, 0
  neg, 0 zero) — i.e. every single one acts as a true extensor, exactly as real posterior-chain anatomy
  predicts, with zero sign outliers to explain away.

## 6. Max lumbar-extension moment (the headline number)

| pose | max available moment (full activation, equilibrated) | vs. 10 N·m placeholder |
|---|---:|---:|
| neutral standing (`lumbar_extension = 0°`) | **228.0 N·m** | ×22.8 |
| bent-over, 30° flexed (`lumbar_extension = -30°`) | **148.6 N·m** | ×14.9 |

This is a **theoretical mechanical ceiling** (every muscle at 100% activation, equilibrated fiber length,
no antagonist coactivation or voluntary-activation deficit) — not a claim about measured voluntary
effort. **External anchor (live-fetched this session, not recalled):** Arokoski et al. (2004, *Archives
of Physical Medicine and Rehabilitation*) measured maximal voluntary isometric trunk extension torque of
147.3±75.9 N·m (pre-intervention) / 170.1±72.3 N·m (post) in real subjects. A muscle-only theoretical
ceiling exceeding measured *voluntary* torque by roughly 1.3–1.5× is the expected, physiologically sane
relationship (voluntary activation is never 100%, and antagonist coactivation costs further reduce net
measured torque) — not a suspicious mismatch. The neutral-pose number sits just outside 1 SD above
Arokoski's post-intervention mean, consistent with a young/athletic OpenCap subject (subject2). The drop
from 228→149 N·m going into flexion reflects the tendon-slack-length calibration choice (fiber length =
optimal at the *neutral* reference pose, the standard convention) — it is a mechanical consequence of
that choice, not a separate physiological claim about the true optimum angle, and is disclosed as such.

## 7. Honest gaps

1. **Torso Z-axis scale factor is assumed, not measured.** Pelvis has 60+ existing muscle path-points to
   derive a clean per-axis scale from; torso has **zero** (all 80 native LaiArnold muscles are
   lower-limb, none cross the torso body). A marker-based check was tried and **rejected** — markers are
   individually re-registered during OpenSim scaling, not purely geometrically scaled (verified: marker
   ratios have std up to 0.53, vs 0.0 for the muscle-point-based pelvis derivation). Torso Z is therefore
   set equal to the measured torso XY factor (0.9879), an isotropic assumption. This mainly affects
   medial-lateral placement of torso-side attachment points, secondary to the flexion-extension moment
   arm reported on here (dominated by X/Y geometry).
2. **Single-hinge reduction cannot represent segmental (single-level) function**, most acutely for
   multifidus (44% of its fascicles excluded on exactly this basis, see §4) — its real clinical role as a
   segment-by-segment stabilizer is fundamentally not expressible in a 3-DOF lumped lumbar joint. This is
   a property of the *target model's* fidelity level (shared with its 3 original ideal actuators), not
   something a better port could fix.
3. **Christophy et al. 2012 lineage attribution is a structural inference** (18-DOF joint architecture +
   fascicle-naming convention match), not independently re-verified against that paper's own DOI this
   session — only the Beaucage-Gauvreau 2019 DOI was independently confirmed live.
3b. **The Beaucage-Gauvreau 2019 paper itself was verified by DOI/PMID/journal/year via WebFetch, not by
   reading its full text or methods this session** (WebSearch quota was independently exhausted before
   this task began — a pre-existing, disclosed constraint, not something this task caused).
4. **Activation-dynamics time constants use OpenSim's `Millard2012EquilibriumMuscle` class defaults**,
   not hand-matched to the native 80 muscles' own values — irrelevant to the static/isometric moment
   check performed here (no forward-dynamics excitation-to-activation lag is exercised), but worth
   checking before using these muscles in a dynamic (excitation-driven) simulation.
5. **Quadratus lumborum, multi-level psoas, obliques, rectus abdominis, latissimus dorsi** are present in
   the same donor model (218 more muscles) but were deliberately **not** ported — out of the task's
   explicit scope (erector spinae ± multifidus), not a capability gap; the same script's pattern extends
   to them directly if wanted later.
6. **The 228/149 N·m figures are ceilings from the newly-added muscles only**, evaluated independently of
   the pre-existing 10 N·m ideal actuators (which remain in the model, additive/supplementary, and would
   trivially add another 10 N·m of idealized capacity if summed — not included in the headline numbers
   above since they are the placeholder being superseded, not a real physiological contribution).

## 8. Next step

The lumbar-extensor muscle gap identified in `MECHANISM_MUSCLE_AUDIT.md` is closed on this copy: erector
spinae (iliocostalis lumborum + longissimus thoracis, 60/76 = 79% of donor fascicles retained) and a
meaningful share of multifidus (28/50, 56%) now exist as literature-parametrized
`Millard2012EquilibriumMuscle` actuators with a verified external anchor in a physiologically sane range,
superseding the 10 N·m placeholder's realism (while leaving it in place, harmless, for any existing
tooling that references it by name). Natural follow-ons, **not** executed here (scope discipline): (a)
re-run the existing IK/ID smoke-test pattern (`scripts/msk/verify_env.py`,
`scripts/msk/smoke_test_ik.py`) against this new copy to confirm it still tracks subject2's real
`walking1.trc` motion correctly with the added (currently zero-activation-by-default, so inert unless
driven) muscles present; (b) if a dynamic/Static-Optimization or Moco use case emerges, hand-tune the
activation/deactivation time constants (honest gap #4) and consider extending the same donor-port
pattern to quadratus lumborum + abdominals (honest gap #5) for a complete detailed-trunk model; (c) file
the same live-verification pattern used here (rather than recalled citations) as the default for any
future literature-model port in this project.

## Files

- `scripts/msk/add_erector_spinae.py` — the build script (idempotent: re-running it re-copies the
  pristine source and rebuilds from scratch).
- `scripts/msk/add_erector_spinae_evidence.json` — full machine-readable evidence (all 126 candidates,
  all 38 exclusions with reasons, per-muscle moment arms/forces at both poses, measured scale factors).
- `data/msk_models/LaiArnoldModified2017_erector_spinae_subject2_scaled.osim` — the upgraded model copy.
