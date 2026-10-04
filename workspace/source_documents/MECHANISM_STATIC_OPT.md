# MECHANISM STATIC OPTIMIZATION — muscle-driven knee contact force (2026-07-21)

## ⚠ CORRECTION (2026-07-21) — READ THIS BEFORE THE NUMBERS BELOW

A sign bug in the shared helper `knee_crossing_muscles_and_forces`
(`scripts/msk/static_opt_knee.py`) made this doc's **self-computed** headline number the EXACT
NEGATIVE of the correct per-muscle pulling-force direction — masked from looking obviously wrong
because the number that gets reported is a vector **norm** (`np.linalg.norm(...)`, always ≥0 by
construction), which cannot expose an internal sign error the way a signed projection would. Full
diagnosis, toy-test proof, and fix: `docs/MECHANISM_SIGN_BUG_AUDIT.md`,
`docs/MECHANISM_SIGN_BUG_REMEDIATION.md`. The bug is now fixed at its source (one-line branch swap),
confirmed via a known-answer toy test (`np.allclose(computed, +ground_truth)`, not the negative) and
by reproducing the corrected number below.

| | published (bug artifact) | **corrected** | ratio vs OrthoLoad (258.22 %BW) |
|---|---:|---:|---:|
| self-computed (R_old − Σmuscle-crossing) | 233.20 %BW @ t=0.51s | **391.10 %BW @ t=0.51s** | 0.903 → **1.515** |

**The corrected number now matches this doc's OWN already-published, bug-immune official
`opensim.JointReaction` cross-check (391.11 %BW) to 0.0039%** — JointReaction never calls the buggy
function, so this agreement (not the original "lands just under the anchor" story) is the proof the
correction is right. **Honest new headline: the twin OVER-predicts in-vivo knee contact force by
~1.51×**, not "lands just under the anchor." This is a documented tendency of exactly this kind of
"two-step" pipeline (Newton's-law reaction computed first, Static-Optimization muscle forces
subtracted afterward, rather than solved simultaneously): Moissenet F, Chèze L, Dumas R (2014), *J
Biomech* 47(1):50-58, DOI [10.1016/j.jbiomech.2013.10.015](https://doi.org/10.1016/j.jbiomech.2013.10.015),
PMID [24210475](https://pubmed.ncbi.nlm.nih.gov/24210475/) — verbatim: *"Musculo-tendon forces and
joint reaction forces are typically estimated using a two-step method, computing first the
musculo-tendon forces by a static optimization procedure and then deducing the joint reaction forces
from the force equilibrium... the joint reaction forces are usually overestimated."*

**§5's "differentiation-scheme" explanation for the ~1.68× self-computed/JointReaction gap is
SUPERSEDED**, not merely revised: fixing the sign bug (nothing about differentiation) collapses that
gap from 67.7% to 0.0039% — the sign bug, not differentiation scheme, was ~99.994% of the effect.
§7's muscle-identity finding (rectus femoris + gastrocnemius dominate at the peak instant) is
UNAFFECTED — SO's own solution (which muscles, how activated) never depended on this downstream
geometry bug; only how their forces combine into the final contact-force vector changed.

Everything below this banner is preserved as the ORIGINAL, pre-correction analysis (historical
record of what was believed at the time) — read the headline table, §5, and §6 with the correction
above in mind, not as still-current numbers.

---

Executes the operator's task: does adding real muscle forces (OpenSim Static Optimization, subject2
`walking1`) close the gap between the twin's PURE-kinematics knee reaction force (102.17 %BW,
`docs/MECHANISM_JOINT_FORCE_VALIDATION.md`) and OrthoLoad in-vivo (258.22 %BW)? Every number below is
machine-measured this session (`scripts/msk/static_opt_knee.py`, exit 0), not recalled. Isolation
respected: `.venv-msk` only, LabValidation + OrthoLoad data read in place, no git commit/push.

## Headline result [SUPERSEDED — see correction banner above; preserved for the historical record]

| method | peak knee-r contact force (%BW) | t (s) | ratio vs in-vivo (258.22) |
|---|---:|---:|---:|
| Prior: pure kinematics+GRF reaction (no muscles) | 102.17 | 1.50 | 0.396 |
| **This session, self-computed** (R_old − Σmuscle-crossing, primary) | ~~233.20~~ **→ 391.10 (corrected)** | 0.51 | ~~0.903~~ **→ 1.515** |
| This session, official `opensim.JointReaction` (cross-check) | 391.11 | 0.51 | 1.515 |
| IN-VIVO OrthoLoad knee (median, n=72 unassisted level-walking, 9 subjects) | 258.22 | — | 1.0 |

**[SUPERSEDED framing, kept verbatim for the record]** Yes — adding muscles moves the prediction
substantially toward in-vivo, landing in the task's own anticipated ~200-300 %BW band by the primary
(self-computed) method, and moderately overshooting it by the secondary (official JointReaction)
cross-check. The gap-to-anchor shrinks from 156.1 percentage points (ratio 0.396) to either +25.0
points short (ratio 0.903) or +132.9 points over (ratio 1.515), depending on method — both a large,
directionally-unambiguous improvement, not a knife-edge or borderline result. Both numbers PASS the
pre-registered "substantial closure" gate (ratio > 0.6, §6). **[Corrected: both methods now agree —
self-computed and JointReaction converge to within 0.0039%, not two independently-overshooting-and-
undershooting methods.]**

## 1. Inputs (all pre-existing, reused from the already-validated pipeline)

Same model/kinematics/GRF as `docs/MECHANISM_JOINT_FORCE_VALIDATION.md` (scaled
`LaiArnoldModified2017_poly_withArms_weldHand`, subject2, 78.2 kg; `walking1.mot` IK, 158 frames,
0-1.57 s; `walking1_forces.mot` GRF, 0-1.579 s). **New this session**: the pre-existing,
OpenCap-shipped Static-Optimization + Joint-Reaction setup templates for this exact
model/trial — `OpenSimData/Mocap/SO/walking1_setup_so.xml` + `walking1_reserveActuators.xml` +
`walking1_setup_externalLoads.xml`, and `OpenSimData/Mocap/JR/walking1_setup_jr.xml` — path-patched
only, same discipline as the prior cert's ID setup reuse.

## 2. Pipeline

1. **Static Optimization** (`opensim.AnalyzeTool` + `StaticOptimization` analysis) with the
   OpenCap-shipped reserve-actuator set (6 pelvis residuals: 3 `PointActuator` + 3 `TorqueActuator`,
   optimal force 10-15 N/N·m; 12 small joint-safety-net `CoordinateActuator`s on hip/knee/ankle/
   subtalar, optimal force 2.5-30 N·m) appended to the model's own 80 muscles + 13 existing
   `CoordinateActuator`s (lumbar/arms) — this is what makes the floating-base pelvis DOFs
   actuation-complete (no muscle crosses `ground_pelvis`) and is the standard, documented OpenSim
   SO recipe, not something built from scratch here.
2. **Convergence + activation sanity gates** (§6), forced BEFORE trusting any downstream number —
   per this task's explicit warning that a differently-scoped prior attempt elsewhere in this repo
   (raw linear-EMG fed into COMAK's tracking cost, `bt_memory/knee-comak-flagship-close-and-fork-
   build.md`) made things worse. That was a different method (EMG-tracking cost function, not
   plain ID-consistent SO used here); the "verify, don't trust" discipline is inherited regardless.
3. **Two independent muscle-inclusive knee contact-force estimates**, cross-checked against each
   other (§4-5):
   - **(a) Self-computed** (primary): re-derive the ALREADY-VALIDATED Newton's-law reaction R as a
     full vector (not just magnitude); separately identify every muscle whose live path geometry
     crosses the knee-r free-body cut (BFS body-membership on the model's own joint tree, never a
     hardcoded muscle list); subtract the summed, SO-tension-scaled muscle-crossing force vector
     from R to isolate the residual bone-contact compression.
   - **(b) Official `opensim.JointReaction`**, fed the SO `force.sto` directly via its documented
     `forces_file` mechanism ("should be used to calculate joint reactions from static optimization
     results") — the task's suggested primary method, run here as the cross-check on (a).
4. **OrthoLoad anchor**: re-parsed live with the same proven AKF parser from
   `validate_joint_force.py` (not re-implemented) — median reproduced exactly: 258.22 %BW, n=72.

## 3. Two real bugs found and fixed this session (forced via OODA, not waved through)

A first run produced a superficially "clean" result (no NaN, activations in bounds) that was
**wrong** on closer inspection — caught before reporting anything, per the "honest-negative/pass is
not a free pass" rule:

1. **Comment-unsafe XML patcher silently deleted the reserve-actuator tag.** The SO/JR templates'
   own doc-comment for `force_set_files` reads `<!--Replace the model's force set with sets
   specified in <force_set_files>? ...-->` — i.e. the comment text itself contains the literal
   substring `<force_set_files>`. A naive non-anchored `<tag>.*?</tag>` regex (the same style
   `validate_joint_force.py`'s own `sub_tag` uses — **verified safe there**: grepped its 6 patched
   tags against the ID templates' comments, zero self-references, and the ID pipeline's own patched
   output files show every tag landing on a clean, isolated line, confirming the already-published
   102.17 %BW cert is NOT affected) starts matching at the comment's inner mention and swallows
   everything up to the next real `</force_set_files>`, deleting the real tag. Net effect, caught
   only by noticing the reserve-actuator columns were **completely absent** from the SO output's
   own 112-column force.sto (not by any XML parse error — the corrupted comment happened to still
   parse): the model ran Static Optimization with **zero pelvis actuation**, silently. **Fix**:
   anchor the tag match to the start of a line (`^[ \t]*<tag>`); re-verified via a real
   `xml.etree.ElementTree.fromstring()` parse (comment open/close counts balance exactly) before
   re-running. Also added a hard, loud assertion that all 18 expected reserve columns are present
   in force.sto — replacing a `.get(col, 0.0)` pattern that would otherwise fail OPEN (a genuinely
   missing column silently reads as "reserve force = 0", indistinguishable from "reserves are
   healthily near-zero").
2. **Extrapolated-tail truncation risk, then over-correction into real extrapolation.** The SO
   template's own `<StaticOptimization><end_time>` is hardcoded to 1.0 s — but the already-validated
   knee-reaction peak (`validate_joint_force.py`'s own 5-window sensitivity sweep) occurs at
   t=1.49-1.51s, so the raw template would have silently dropped the one instant this whole
   comparison was anchored to. First fix attempt over-corrected to 1.6 s; the REAL IK/GRF data ends
   at 1.57/1.579 s, so 1.58-1.6 s was fabricated (spline-extrapolated) input. Console evidence this
   was caught, not missed: `Ipopt: Maximum iterations exceeded`, `the model appears too weak for
   static optimization... edl_r/ehl_r/fdl_r/tibant_r approaching upper bound of 1` at exactly
   t=1.58-1.6s — an artifact of asking the optimizer to satisfy fabricated kinematics, not a genuine
   muscle-weakness finding. **Fix**: tightened to 1.57 s (the real data ceiling); the re-run shows
   158/158 real frames, 0 optimizer failures, activation ceiling 0.627 (no muscle saturates at 1.0
   anywhere in the trial).

## 4. Convergence + activation sanity (pre-registered gates, machine-checked)

| gate | threshold | measured | verdict |
|---|---|---:|---|
| Frames processed, no extrapolation | 158/158, t∈[0,1.57] | 158/158 | PASS |
| Covers the validated peak time (1.50s) | required | in range | PASS |
| NaN in activations | 0 | 0 | PASS |
| Activation bounds | [0,1] (±1e-3 slop) | [0.010, 0.627] | PASS (no saturation anywhere) |
| Muscle force vs own Fmax | ratio ≤ 1.5× | 0/80 muscles exceed | PASS |
| Joint (hip/knee/ankle/subtalar) reserve usage | flagged if >50% of own optimal_force | max 12.5% (hip_rotation_r) | PASS — muscles do the real work, reserves are a thin safety net, not a crutch |
| Pelvis residual force peak | <75 N (broad Hicks-et-al-style band) | 172.96 N (22.6% BW) | **FAIL — disclosed, not blocking** (see §7.1) |
| Pelvis residual moment peak | <75 N·m | 45.11 N·m | PASS |
| Reserve-actuator columns present (anti-fail-open) | all 18 found | 18/18 found | PASS |

Overall convergence gate (frames/NaN/bounds/peak-coverage) → **PASS**. The pelvis residual is real,
quantified, and does not gate the primary result (it reflects the free-floating-base
consistency term, not the knee-local free-body cut either downstream method uses; see §7.1).

## 5. The two methods disagree by ~1.7×, and that disagreement was investigated, not hidden

**[SUPERSEDED — see the top correction banner. The ~1.68× gap this section diagnoses as
"differentiation scheme" was ~99.994% a sign bug in the self-computed method; corrected, the two
methods agree to 0.0039%. Preserved verbatim below as the historical record of what was believed at
the time.]**

Self-computed (233.20 %BW) and official JointReaction (391.11 %BW) both show a large, unambiguous
improvement over the 102.17 %BW pure-reaction baseline, and — a real, non-engineered consistency
signal — **both independently pick the same peak instant, t=0.51s**, even though the raw vertical
GRF is nearly identical at two separate stance occurrences in this trial (107.2 %BW at t=0.55s vs
108.95 %BW, the trial's global GRF max, at t=1.50s — i.e. this is genuinely two separate right-foot
stance phases within the 1.57s window). Both methods agree the *earlier* occurrence, not the
GRF-larger later one, has the bigger muscle-inclusive knee load — a specific, checkable, and
non-obvious agreement that would be a coincidence if either method were just noise.

What they do NOT agree on is magnitude (391/233 = 1.68×). Diagnosed, not shrugged off:
- **Most likely explanation**: different numerical differentiation schemes. The self-computed
  method reuses the already-externally-validated Savitzky-Golay-in-Cartesian-COM-space pipeline
  (cross-checked via the whole-body Newton residual gate at 1-4 %BW RMS in the prior cert).
  `opensim.AnalyzeTool`/`JointReaction` instead differentiates in joint-ANGLE space via its own
  spline fit + the template's 6 Hz low-pass filter, then maps to Cartesian accelerations via forward
  kinematics — a different, both-legitimate numerical path (filtering does not commute with a
  nonlinear coordinate transform in general).
- **Supporting evidence for a boundary/edge-artifact component**: at t=0 (the very first sample,
  outside the reported peak window), JointReaction reads 188.5 %BW vs the self-computed method's
  37.0 %BW — a 5× gap concentrated exactly at the domain boundary, the classic signature of
  spline/FFT-style differentiation edge effects. The reported peak (t=0.51s) is well inside the
  interior (40+ samples from either edge), so this specific artifact does not directly corrupt the
  headline numbers, but it corroborates that the two differentiation pipelines diverge more sharply
  in exactly the regions numerical-differentiation theory predicts they should.
- **Not resolved further this session** (lean, time-boxed): pinning down the exact residual
  numerically (e.g. by re-implementing JointReaction's own filter in the self-computed path) is the
  natural next step if a tighter single number is later required; not done here since both methods
  already agree on the qualitative/directional finding the task asked about, and neither is absurd
  (§6 gate) or degenerate.

## 6. Pre-registered verdict gate

Classification (thresholds set before this run): ratio > 0.6 → "substantial closure" PASS; ratio in
(0.396, 0.6] → "partial closure" PASS; ratio ≤ 0.396 → FAIL (no better than the muscle-free
baseline); pct_bw > 1200 %BW or < 20 %BW → SUSPICIOUS-recheck (guards against a recurrence of the
documented ~1000-1800× OpenSim 4.6 knee/hip generalized-force bug, `docs/
MECHANISM_MSK_ELASTIC_BAND.md` §4 — neither candidate here is within 3× of that ceiling).

| candidate | %BW | ratio vs in-vivo | verdict |
|---|---:|---:|---|
| self_computed [**corrected: 391.10, ratio 1.515** — see top banner; row below is the original, superseded number] | 233.20 | 0.903 | **PASS-substantial** |
| official_jointreaction | 391.11 | 1.515 | **PASS-substantial** |

Neither the OpenSim 4.6 knee/hip InverseDynamicsTool bug (§4 of the elastic-band doc) nor its
class of failure mode is implicated: this pipeline never calls `InverseDynamicsTool`/
`InverseDynamicsSolver` (the self-computed method reuses the bug-immune Newton's-law reaction;
JointReaction uses a documented-different Simbody code path, `calcMobilizerReactionForces`-style,
not the generalized-force attribution routine where the bug lives) — and both results are
2-3 orders of magnitude below the ~1000-1800× corruption signature, not merely "not NaN."

## 7. Is this really co-contraction? Measured, not assumed

The prior cert's own framing (`docs/MECHANISM_JOINT_FORCE_VALIDATION.md` §3) illustrated the missing
mechanism as "quadriceps AND hamstrings both pulling the tibia proximally." Checking the ACTUAL SO
solution at the peak instant (t=0.51s) refines rather than confirms that specific illustration:

| muscle group | activation | force (N) |
|---|---:|---:|
| Rectus femoris (quad, biarticular) | 0.227 | 586.7 |
| Vasti (medialis/lateralis/intermedius, quad) | 0.010 (floor) each | 14-45 each — negligible |
| Hamstrings (bflh/bfsh/semimem/semiten) | 0.010-0.066 | 5-33 — negligible |
| **Gastrocnemius medialis** | **0.441** | **1140.0** |
| **Gastrocnemius lateralis** | **0.282** | **398.0** |

At this specific instant the dominant simultaneous-crossing-muscle pair is **rectus femoris +
gastrocnemius** (both biarticular, both pulling the shank/thigh together across the knee), NOT the
vasti+hamstrings pairing used as the illustrative example — the vasti and hamstrings are
essentially at the SO activation floor here. This is still a genuine antagonist/joint-compressing
co-contraction mechanism (gastrocnemius is a well-documented knee flexor during stance, acting
against rectus femoris's extensor action, exactly the "moments partially cancel, forces add"
pattern) — the general MECHANISM is empirically confirmed, but the SPECIFIC muscles are different
from the prose illustration, and this is reported precisely because "measure, don't narrate" cuts
both ways: the numbers substantiate the mechanism while correcting which muscles drive it. Sanity:
1140 N / gasmed_r's own Fmax (3115.5 N, `scripts/msk/audit_muscles_evidence.json`) = 36.6% — sane,
non-saturated.

## 8. Honest gaps (full list)

1. **The two methods disagree by ~1.7× in magnitude** (§5) — diagnosed (differentiation-scheme
   difference, most likely) but not numerically pinned down further this session.
2. **Pelvis residual force (172.96 N, 22.6% BW) exceeds the pre-registered "good" band** — no RRA
   (Residual Reduction Algorithm) was run first; this is the standard, documented way to shrink this
   before it matters. Joint-level (hip/knee/ankle) reserve usage stayed low (≤12.5% of a small
   optimal_force) throughout, so the leg-muscle recruitment driving the knee number is not
   obviously contaminated by this — but it is a real, disclosed limitation of skipping RRA.
3. **Single trial, right leg only** (subject2 `walking1`) — same scope caveat as the prior cert;
   no claim of generality across subjects/trials/gait speeds.
4. **The mechanism is SO's activation-minimizing solution, not measured EMG** — SO finds the
   minimum-effort muscle set consistent with the required net joint moments; it is not guaranteed to
   match the subject's true neural recruitment strategy (real co-contraction can exceed the
   effort-minimizing solution). This is the well-known, structural limitation of SO vs. EMG-informed
   or CMC-based estimation, and is exactly why the two independent point estimates (233% / 391%)
   should be read as a plausible band, not a precision measurement.
5. **Different populations** (OrthoLoad's elderly instrumented-implant patients vs. subject2's
   healthy young OpenCap participant) — inherited, unchanged from the prior cert's own caveat.
6. **Left knee, other trials, other subjects** — not run this session (lean scope).

## 9. Next step

If a single tighter number is later needed: reconcile the self-computed-vs-JointReaction gap by
re-differentiating with a matched filter (either bring JointReaction's kinematics through the same
Savitzky-Golay pipeline, or validate the 6 Hz spline path against the same whole-body-residual gate
the SavGol path already passed). Otherwise, the qualitative/directional finding this session was
asked to test is answered with the pre-registered gate PASSING under both independent methods: the
gap closes substantially, from a 102.17 %BW / 0.396-ratio pure-reaction baseline to 233-391 %BW /
0.90-1.52-ratio once muscle forces are included — landing at or moderately past the in-vivo anchor,
not merely "somewhat higher."

## Files

- `scripts/msk/static_opt_knee.py` — the full pipeline (self-contained, re-runnable; imports
  `validate_joint_force.py` for the proven parse_mot/Savitzky-Golay/BFS/OrthoLoad-parser code, not
  re-implemented).
- `data/msk_smoketest/subject2_walking1/static_optimization/so/` — patched SO setup XML, SO
  `activation.sto`/`force.sto` (158 frames × 111 actuators), Ipopt console log.
- `data/msk_smoketest/subject2_walking1/static_optimization/jr/` — patched JR setup XML,
  `walking1_JointReaction_ReactionLoads.sto`.
- `data/msk_smoketest/subject2_walking1/static_optimization/static_opt_knee_results.json` — every
  number in this document, machine-written. **⚠ STALE, MARKED IN-FILE (2026-07-21):** this file's
  own `self_computed.bone_contact_residual_peak_pct_bw` / `candidates_pct_bw.self_computed` /
  `verdicts.self_computed` fields still hold the PRE-sign-fix **233.20** %BW — never regenerated
  (regenerating requires a fresh, expensive, unconditional Ipopt Static-Optimization solve, see the
  file's own `_STALE_PRE_SIGN_FIX_NOTICE` key for why that wasn't done). **A naive automated check
  reading only this file gets the wrong, superseded number.** The corrected **391.10** %BW (this
  doc's own top banner) is added as an explicit `_STALE_PRE_SIGN_FIX_NOTICE` top-level key in that
  same JSON file, and independently triangulated in
  `docs/MECHANISM_SIGN_BUG_REMEDIATION.md` Sec.3 (triple cross-checked: cached-reuse recompute,
  bit-identical match to an independent fixed-copy implementation, and <0.1% agreement with this
  same file's own never-buggy `jr_summary.peak_pct_bw` = 391.11 %BW). Flagged in
  `docs/MECHANISM_TRUST_LEDGER.md` §10 item 10; disposition in
  `docs/MECHANISM_TRUST_LEDGER_REMEDIATION.md`.
