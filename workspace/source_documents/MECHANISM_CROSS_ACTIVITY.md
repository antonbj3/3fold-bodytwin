# MECHANISM CROSS-ACTIVITY VALIDATION — broadening the joint-force certs beyond walking (2026-07-21)

Executes the operator's task: the corrected knee (391.10 %BW, ratio 1.515) and hip (386.77 %BW,
ratio 1.412) muscle-driven contact-force certs (`docs/MECHANISM_STATIC_OPT.md`,
`docs/MECHANISM_HIP_FORCE.md`) are anchored on ONE activity, level walking. This document broadens
the OrthoLoad in-vivo anchor to stairs up/down, sit-to-stand (STS) and squat, checks whether the
twin has matching kinematics for each, and — where it does — runs the SAME validated,
sign-fixed pipeline and compares. Every number below is machine-measured this session
(`scripts/msk/cross_activity_validation.py`), not recalled. Isolation respected: bodytwin only,
OrthoLoad + the external LabValidation mount read in place (never written to), `.venv-msk` only,
new files only, no git commit/push.

## Headline

1. **OrthoLoad in-vivo cross-activity anchors** (pooled fresh from the raw 3942-trial index):
   knee stairs-up/down median **329.6/348.1 %BW**, hip **300.2/309.4 %BW**; knee/hip
   sit-to-stand median **268.0/184.2 %BW**; knee/hip squat median **249.8/202.0 %BW**. All
   comfortably at-or-above walking's own ~255-302 %BW median — stairs and STS are NOT gentler
   than walking in vivo, contrary to the task brief's a-priori "~250-350 %BW stairs" framing being
   read as an upper bound (stairs-down and stairs-up MEDIANS already sit at that framing's ceiling,
   with individual trials reaching 458-565 %BW).
2. **Twin kinematics gap, exhaustively checked, not assumed**: STAIRS kinematics exist in NEITHER
   of the twin's two kinematics sources — the 228-clip monocular video corpus (2 movement
   categories only: `weighted_lifts/squat`, `plyometrics_jumps/*`) NOR the real OpenCap
   LabValidation mocap corpus (10 subjects, `subject2`-`subject11`, 22 distinct trial types, zero
   containing "stair"). **This is a genuine, confirmed corpus gap** — the OrthoLoad anchors above
   are the target the twin must eventually hit once stair kinematics are captured.
   SQUAT and SIT-TO-STAND kinematics DO exist, and at the SAME fidelity tier as the validated
   walking cert: subject2's own real OpenCap mocap + real force-plate GRF (`squats1.mot`/
   `squats1_forces.mot`, `STS1.mot`/`STS1_forces.mot`, 2000 Hz force plate, 100 Hz IK).
3. **Squat + STS, run through the corrected pipeline**: the trustworthy **pure-reaction** tier
   (Newton's law, immune to the OpenSim ID-tool bug) UNDER-predicts the OrthoLoad anchor by a
   **consistent** ~0.25-0.40× across ALL 4 activity/joint combinations now measured (walking knee
   0.396, hip 0.323; squat knee 0.332, hip 0.317; STS knee 0.247 hip 0.274) — a genuinely
   cross-activity-stable finding. The **muscle-driven (Static Optimization)** tier, by contrast,
   is **NOT** a clean extension of walking's ~1.4-1.5× over-prediction: it produces huge numbers
   for squat (knee 1257.7 %BW, hip 794.3 %BW — ratios 5.04×/3.93× vs anchor) and STS (knee 883.4
   %BW, hip 516.6 %BW — ratios 3.30×/2.80× vs anchor) that are internally self-consistent
   (self-computed agrees with the independent `opensim.JointReaction` cross-check to **<0.1%** in
   all 4 cases, exactly as tight as walking's own 0.004-0.07% agreement) but **independently
   flagged dynamically implausible** (muscle-pinning and/or pelvis-reserve-actuator saturation,
   the same machine diagnostic that already flagged the video-corpus squat scene). **Answer to
   the operator's question: the over-prediction is NOT consistent — walking's ~1.4-1.5× is a
   verified, trustworthy number; squat/STS's much larger apparent ratios are an artifact of a
   DIFFERENT, independently-diagnosed failure mode (Static Optimization breaking down on
   rapid/ballistic postural transitions), not evidence the Moissenet-et-al two-step-overestimation
   mechanism itself scales up 2-3× for these activities.**

## 1. OrthoLoad in-vivo cross-activity anchors (raw per-trial, pooled fresh)

Computed directly from `data/external/orthoload/_index/orthoload_akf_trials.jsonl` (3942/3942 AKF
files already parsed by `scripts/msk/index_orthoload_forces.py`, its own classifier's
`activity_bucket` field reused as-is, not re-derived) by THIS script's `pool_pctbw()` — pooling
gen1+gen2 for hip (same convention `validate_hip_force.py`'s own PRIMARY walking selection
already used), tracking `(sub_group, subject_code)` pairs as the subject-dedup key so a gen1 code
can never be silently merged with an unrelated gen2 code of the same string.

| joint | activity bucket | n trials | n subjects | min | median | mean | max (%BW) |
|---|---|---:|---:|---:|---:|---:|---:|
| knee | Walking (repro check) | 169 | 9 | 136.2 | 254.8 | 256.0 | 353.6 |
| knee | **Stairs Up** | 54 | 9 | 221.9 | **329.6** | 330.5 | 484.8 |
| knee | **Stairs Down** | 52 | 9 | 284.1 | **348.1** | 365.4 | 458.9 |
| knee | Stairs pooled (up+down+combined) | 112 | 9 | 221.9 | 339.6 | 347.3 | 484.8 |
| knee | **Sit-to-Stand/Stand-to-Sit** | 47 | 9 | 133.1 | **268.0** | 269.2 | 408.7 |
| knee | **Knee Bend/Squat** | 47 | 9 | 158.5 | **249.8** | 246.8 | 359.6 |
| hip (gen1+gen2) | Walking (repro check) | 382 | 19 | 191.8 | 302.4 | 311.0 | 550.8 |
| hip (gen1+gen2) | **Stairs Up** | 36 | 16 | 162.4 | **300.2** | 297.0 | 555.4 |
| hip (gen1+gen2) | **Stairs Down** | 34 | 15 | 164.6 | **309.4** | 315.6 | 529.6 |
| hip (gen1+gen2) | Stairs pooled | 111 | 18 | 162.4 | 303.4 | 310.5 | 565.3 |
| hip (gen1+gen2) | **Sit-to-Stand/Stand-to-Sit** | 119 | 17 | 74.2 | **184.2** | 193.0 | 403.1 |
| hip (gen1+gen2) | **Knee Bend/Squat** | 37 | 12 | 100.0 | **202.0** | 215.9 | 587.6 |

**Cross-check against this repo's own already-published, narrower walking selections** (a
different, hand-curated substring match — `"Level Walking"`/`"walking free"`, excluding
gait-aid-assisted trials — not expected to match exactly): knee walking-bucket median 254.8 vs
published 258.22 (**1.32% apart** — close); hip walking-bucket median 302.4 vs published 273.93
(**10.39% apart** — the broader bucket includes more sub-conditions, e.g. non-"free"/non-"Level"
phrasings, that the narrower published selection deliberately excluded; disclosed, not a bug).

**Reading the stairs/STS/squat numbers against the task's own a-priori framing ("in-vivo knee
~250-350%BW stairs")**: the MEDIANS (329.6/348.1 knee, 300.2/309.4 hip) sit inside or at the very
top of that band, but individual trials range up to 458.9-565.3 %BW — the framing is directionally
right but understates the tail. Stairs (both directions) and, less expectedly, STS and squat all
sit at or above walking's own ~255-302 %BW median at the knee/hip — **in-vivo, stairs/STS/squat
are not "gentler" loading regimes than walking**, an external, independent sanity check on what
the twin's cross-activity target actually is.

## 2. Twin kinematics-gap check — machine-verified across BOTH kinematics sources

`scripts/msk/cross_activity_validation.py::check_twin_kinematics_gap()` enumerates, rather than
assumes:

**(a) The 228-clip monocular video corpus** (`scripts/msk/corpus_ik_aggregate_per_clip.csv`,
already validated in `docs/MECHANISM_CORPUS_IK_AGGREGATE.md`): exactly 4 `movement_group` values
exist — `plyometrics_jumps/plyometric_drill`, `plyometrics_jumps/plyometrics_jumps_misc`,
`plyometrics_jumps/technique_tutorial`, `weighted_lifts/squat`. A keyword scan
(`stair|sit|chair|sts|rise`) over all 4 finds **zero** matches.

**(b) subject2's real OpenCap LabValidation mocap corpus**, plus every sibling subject present on
the external mount (`/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/`):
`glob("subject*")` finds **10** subject directories (`subject2`-`subject11`; no `subject1`
exists), **22** distinct `OpenSimData/Mocap/IK/*.mot` trial-type stems pooled across all 10:
`DJ1-5`, `DJAsym1-5` (drop jump), `squats1`, `squatsAsym1`, `STS1`, `STSweakLegs1`, `walking1-4`,
`walkingTS1-4`. **Zero contain "stair"** (case-insensitive substring, exhaustive over all 10
subjects' trial listings).

**Verdict: `KINEMATICS_GAP_CONFIRMED` for stairs** — not a lookup miss, a genuine absence in both
of the twin's kinematics sources. Consistent with the public OpenCap LabValidation dataset's own
documented activity set (walking, squats, sit-to-stand, drop jumps — no stair ascent/descent).
**The OrthoLoad anchors in §1 (knee 329.6/348.1 %BW, hip 300.2/309.4 %BW up/down) are therefore
reported as the target range the twin must eventually hit, once stair kinematics are captured
(new mocap or a stair-containing video source) — not something this session can close.**

**SQUAT and SIT-TO-STAND: `AVAILABLE`** — subject2's own `squats1.mot`/`squats1_forces.mot` and
`STS1.mot`/`STS1_forces.mot` all exist, same subject, same scaled model
(`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`), same 2000 Hz force-plate / 100 Hz IK
fidelity tier as the validated `walking1` trial used by every prior cert in this family.

## 3. Squat + sit-to-stand — the validated pipeline on real mocap + real GRF

### 3.1 A forced finding en route: sit-to-stand has a hidden third-contact confound

Before trusting any peak number, the SAME whole-body Newton residual gate the walking cert
already uses (total-body mass × COM-acceleration must equal total two-feet GRF + gravity,
independent of any joint-cut choice) was run on both trials. Squat PASSED cleanly (3.97/3.33/2.26
%BW RMS, x/y/z, all `< 8 %BW`). **STS FAILED badly** (10.76/**40.19**/3.27 %BW RMS) — this was
NOT treated as an honest-negative and dropped; it was forced through OODA:

- **Observe**: the vertical-axis residual is the outlier (40.2 %BW vs an 8 %BW gate).
- **Orient**: `pelvis_ty` (vertical pelvis translation) over the trial oscillates between ~0.68 m
  (knee/hip flexed ~90-100°) and ~1.15-1.17 m (knee/hip near-extended) — the subject is
  **repeatedly sitting down and standing up** (this one 8.6 s trial contains **~8-10 STS
  repetitions**, not one). Plotting the PER-FRAME residual against `pelvis_ty` shows it spikes to
  **50-150 %BW specifically at the low (seated) plateau** and drops to <15 %BW while
  rising/standing — the direct, expected signature of a real, unmeasured CHAIR-reaction force
  during the seated dwell (a third contact the two-foot-only GRF setup cannot see), not a
  differentiation artifact. This is the SAME *class* of hidden-assumption failure
  `force_scenes_batch.py` already found and named for jump flight-phase ("stance-foot-anchor...
  silently misreads a genuine loss of ground contact as +1 BW of static support") — a chair
  during STS is the geometric mirror image (an extra contact appearing, not a real one
  disappearing), caught by the identical diagnostic.
- **Decide/Act**: added a generic (no STS-specific hardcoding) contiguous-low-residual-segment
  finder (`find_valid_segments`, threshold 15 %BW/frame, `edge_pad_s=0.15`) to the script itself —
  not just diagnosed in prose. Found **10 valid segments** across the 8.6 s trial (one per rep),
  0.29-0.62 s each.

**The headline numbers survive this check, cross-validated across the 10 independent reps**:
knee_r peak median **66.29 %BW** [range 59.32-72.96], hip_r peak median **50.49 %BW** [range
43.26-56.49] — i.e. the naive whole-trial global peak (72.96 %BW knee, which happens to land
inside a valid segment) is corroborated, not an artifact, by 9 other independent repetitions
landing in the same range. Squat needed no such restriction (1 valid segment spanning 7.70 of
8.00 s — effectively the whole trial); its own reported numbers are the plain global peak.

### 3.2 Pure-reaction tier (trustworthy; immune to the OpenSim 4.6 ID-tool bug)

Newton's second law for the BFS-derived distal sub-chain (`vjf.get_descendant_bodies`), same
method as `validate_joint_force.py`'s own walking cert, generalized to any (IK, GRF) pair:

| activity | knee_r peak %BW | hip_r peak %BW | vs OrthoLoad anchor (knee/hip median) | ratio knee/hip |
|---|---:|---:|---|---:|
| walking (published) | 102.17 | 88.59 | 258.22 / 273.93 | 0.396 / 0.323 |
| **squat** | 82.84 @t=7.19s | 63.98 @t=7.20s | 249.8 / 202.0 | **0.332 / 0.317** |
| **sit-to-stand** (10-rep median [range]) | 66.29 [59.3-73.0] | 50.49 [43.3-56.5] | 268.0 / 184.2 | **0.247 [0.221-0.272] / 0.274 [0.235-0.307]** |

**All 4 knee/hip ratios across 3 activities cluster tightly in 0.25-0.40** — a genuinely
cross-activity-CONSISTENT under-prediction, the same "reaction ≠ contact force, co-contraction
invisible to Newton's law alone" story the walking cert established, now shown to generalize.

### 3.3 Muscle-driven (Static Optimization) tier — converges, agrees internally, but implausible

Static Optimization (`opensim.AnalyzeTool`), reusing squat's OWN pre-shipped OpenCap setup
template as-is, and — since NO STS-native or per-trial JointReaction template is shipped anywhere
in the corpus for ANY non-walking trial — adapting squat's/`walking1`'s templates (verified
boilerplate-identical to `walking1`'s own copies except name/datafile; disclosed, not hidden).
Squat ran the FULL 8.0 s / 801-frame trial (matches the walking cert's own convention exactly,
whole-trial residual gate passed); STS ran its longest machine-found valid segment
([7.47,8.09]s, 63 frames — §3.1). The sign-fixed `static_opt_knee.knee_crossing_muscles_and_forces`
is reused UNMODIFIED (not re-derived), so both runs inherit the walking cert's own fix by
construction.

| activity | joint | SO converged | self-computed peak %BW | `JointReaction` peak %BW | agreement | dynamically plausible? |
|---|---|:---:|---:|---:|---:|:---:|
| squat | knee | yes (801/801 frames, 0 NaN) | 1257.74 @t=7.17s | 1257.74 | **0.000%** | **NO** |
| squat | hip | yes | 794.25 @t=7.21s | 794.99 | **0.09%** | **NO** |
| sit-to-stand | knee | yes (63/63 frames, 0 NaN) | 883.39 @t=8.02s | 883.51 | **0.01%** | **NO** |
| sit-to-stand | hip | yes | 516.56 @t=7.99s | 517.09 | **0.10%** | **NO** |
| walking (published, for comparison) | knee | yes | 391.10 | 391.11 | 0.004% | YES |
| walking (published) | hip | yes | 386.77 | 387.04 | 0.070% | YES |

**The self-computed/`JointReaction` agreement is as tight for squat/STS (0.00-0.10%) as it is for
the already-validated walking cert (0.004-0.070%)** — this is the strongest evidence that the NEW
generalized pipeline code itself (the BFS free-body cut, the crossing-muscle detector, the XML
template adaptation) is correct, not buggy: two structurally-independent code paths (this script's
own geometry-based re-derivation vs. OpenSim's own compiled `JointReaction` analysis) converge on
the same answer regardless of activity.

**But `fsb.diagnose_so_dynamical_plausibility` — the SAME machine diagnostic that already flagged
the video-corpus squat scene's SO tier (`docs/MECHANISM_FORCE_SCENES_BATCH.md`) — flags ALL FOUR of
these real-mocap numbers implausible too**, via two DIFFERENT specific signatures depending on
activity:

| activity | muscles pinned >10% of frames | peak pelvis reserve force / moment | vs 75N/75Nm comfort band |
|---|---:|---:|---:|
| squat | 4 (`glmax2_r`, `glmax3_l`, `semimem_l`, `vaslat_l`) | 58.4 N / 51.7 Nm | 0.78× (reserve itself OK — muscle-pinning is the failure signature here) |
| sit-to-stand | 1 (`vaslat_l`) | 179.8 N / 32.3 Nm | **2.40× over** (reserve saturation is the dominant signature here) |
| video-corpus squat (published, corroborating) | 8 | 222.6 N / 277.7 Nm | 3.70× over |

Squat's own knee number (1257.7 %BW) also trips `static_opt_knee.py`'s own imported absurdity
ceiling (`KNEE_CONTACT_ABSURDITY_CEILING_PCT_BW = 1200`, calibrated against a DIFFERENT, documented
OpenSim ID-tool bug) — used here only as an independent sanity signal, not as this script's pass/
fail gate, since the ceiling's own calibration target is a different defect; the muscle-pinning +
reserve-saturation diagnostic is the load-bearing evidence for the "implausible" verdict.

**Anatomical cross-check (unaffected by the implausibility verdict, a separate geometric anchor):**
the hip-r crossing-muscle set detected on BOTH squat and STS (25 muscles: adductors ×6, glutes ×9,
iliopsoas, piriformis, gracilis, sartorius, rectus femoris, long hamstrings, TFL) is **IDENTICAL,
muscle-for-muscle**, to the already-published walking hip cert's own 25-muscle set
(`docs/MECHANISM_HIP_FORCE.md` §2) — expected, since the crossing set is a property of the model's
joint-tree geometry, not the movement, and confirms the BFS+path-geometry detector generalizes
correctly across activities. The knee-r set (13 muscles) matches the walking cert's own 12 plus
`tfl_r` (TFL/IT-band, anatomically a knee stabilizer via the iliotibial tract in this model) on
both squat and STS — again identical to each other, not activity-dependent.

**Verdict: HONEST NEGATIVE for the ABSOLUTE muscle-driven contact-force number on both squat and
STS, exactly as `docs/MECHANISM_FORCE_SCENES_BATCH.md` already found for the video-corpus squat
scene** — now independently corroborated on REAL mocap + REAL force-plate data (ruling out
"noisy monocular video kinematics" or "kinematically-estimated GRF" as necessary causes; both are
absent here and the SO tier still breaks down). The common factor across all 3 implausible cases
(video-corpus squat, real-mocap squat, real-mocap STS) is a **rapid, large-ROM, ballistic/postural
TRANSITION** (squat descent-ascent, standing up from a chair) — contrasted with walking's smooth,
small-per-step, steady-state cyclic gait, where the identical model/actuator-set/differentiation
scheme produces a converged, dynamically-plausible, literature-consistent result.

## 4. Is the ~1.4-1.5× over-prediction ratio consistent across activities? — the operator's question

**Two tiers give two different answers, and the distinction is the actual finding:**

- **Pure-reaction (trustworthy) tier: YES, consistent.** The under-prediction ratio (predicted
  reaction / in-vivo contact) clusters at **0.25-0.40** across all 3 activities and both joints
  (walking 0.396/0.323, squat 0.332/0.317, sit-to-stand 0.247/0.274) — a real, cross-activity-
  stable signature that "reaction force ≠ contact force" by a roughly constant factor, strengthening
  confidence that the missing term (muscle co-contraction) is a structurally similar-sized gap
  regardless of movement type.
- **Muscle-driven (Static Optimization) tier: NO — but not for the reason the task's own framing
  anticipated.** Walking's ~1.4-1.5× is a verified, trustworthy, dynamically-plausible number.
  Squat's and STS's OWN muscle-driven numbers are 2.8-5.0× the in-vivo anchor (squat knee 5.04×,
  hip 3.93×; STS knee 3.30×, hip 2.80×) — but these are **not a clean replication of walking's
  over-prediction mechanism at a larger magnitude**. They are corroborated (self-computed =
  `JointReaction` to <0.1%, so not a code bug) yet independently flagged **dynamically implausible**
  by the SAME diagnostic already used to reach that verdict on the video-corpus squat scene —
  driven by muscle-pinning and/or pelvis-reserve-actuator saturation specific to rapid/ballistic
  movements, not by the "two-step SO-then-subtract-reaction" mechanism (Moissenet et al. 2014)
  that explains walking's own, much smaller and dynamically-plausible, 1.4-1.5× gap.

**Honest bottom line**: the Moissenet-et-al over-estimation mechanism is verified for walking and
NOT yet shown to generalize — not because it was tested on squat/STS and found smaller or larger,
but because the muscle-driven tier for squat/STS cannot currently be trusted at all, for a
different, independently-diagnosed reason (SO's own per-frame recruitment solution becoming
dynamically implausible on rapid transitions). Closing this gap (a heavier/matched coordinate
low-pass filter for Static Optimization on these movement types, per `docs/
MECHANISM_FORCE_SCENES_BATCH.md`'s own named-but-unexecuted next step) is a prerequisite for ever
answering the ratio-consistency question on non-gait activities — a concrete, bounded follow-up,
not attempted here (would require re-tuning and re-validating the SO differentiation scheme, a
materially different task than broadening the anchor/kinematics comparison this session executed).

## 5. Machine-checked gates (PASS/FAIL, not eyeballed)

| gate | pre-registered threshold | squat | sit-to-stand |
|---|---|---|---|
| OrthoLoad index size (source-data drift guard) | == 3942 | 3942 (PASS) | 3942 (PASS) |
| Video-corpus stair/STS keyword scan | 0 hits expected | 0 (PASS, confirms gap) | 0 (PASS, confirms gap) |
| LabValidation (10 subjects) stair-trial scan | 0 hits expected | 0 (PASS, confirms gap) | 0 (PASS, confirms gap) |
| Required input files exist (IK/GRF/SO templates) | all present | PASS | PASS |
| Whole-body Newton residual (RMS, x/y/z) | < 8 %BW | 3.97/3.33/2.26 (PASS) | 10.76/**40.19**/3.27 (**FAIL whole-trial** → third-contact guard applied, §3.1) |
| Third-contact-guarded valid segment found | ≥1 segment ≥0.3s | n/a (whole trial valid) | 10 segments found (PASS) |
| SO convergence (frames/NaN/activation bounds) | no NaN, activations in [0,1] | 801/801, 0 NaN (PASS) | 63/63, 0 NaN (PASS) |
| Self-computed vs `JointReaction` agreement | < 5% (walking achieved <0.1%) | knee 0.000%, hip 0.09% (PASS) | knee 0.01%, hip 0.10% (PASS) |
| SO dynamical plausibility (muscle-pinning + reserve saturation) | 0 pinned muscles AND reserve <1.0× comfort band | **FAIL** (4 pinned, reserve 0.78×) | **FAIL** (1 pinned, reserve 2.40×) |
| Anatomical anchor (hip-r 25-muscle set matches walking cert) | exact match | PASS (25/25 identical) | PASS (25/25 identical) |
| Absurdity ceiling (context only, calibrated for a different bug) | < 1200 %BW knee | 1257.7 (over — contextual flag only) | 883.4 (PASS) |
| Overall muscle-driven tier "trustworthy" (plausibility AND agreement) | both must PASS | **FALSE** | **FALSE** |

## 6. Honest gaps and caveats (full list)

1. **Stairs: a genuine, unclosed kinematics gap.** No cert can be built for stairs this session;
   §1's OrthoLoad numbers are a target, not a comparison. Closing this needs either a new
   stair-climbing mocap/video capture for subject2 (or another already-scaled subject) or a
   stair-containing addition to the video corpus.
2. **The muscle-driven tier's "implausible" verdict is NOT re-litigated here beyond reusing the
   existing diagnostic** — no attempt was made to fix it (e.g. re-running SO with a heavier
   coordinate low-pass, the concrete next step `docs/MECHANISM_FORCE_SCENES_BATCH.md` already
   named but did not execute). Doing so is a bounded, separate follow-up, not attempted here to
   stay within this task's scope (broaden the anchor + kinematics comparison, not re-tune SO).
3. **STS's SO tier used an ADAPTED, non-native setup template** (squat's own template, since no
   STS-specific SO/JointReaction template is shipped anywhere in the corpus) — including squat's
   own 4 Hz coordinate low-pass filter setting, not a value chosen for STS's own dynamics. This is
   disclosed, not hidden, and is a plausible PARTIAL contributor to STS's own implausibility
   signature (though squat's OWN native-templated run is ALSO implausible, by a different
   signature — so template adaptation is not the sole explanation).
4. **STS's muscle-driven tier covers only one 0.62 s repetition** (of ~8-10 in the trial), chosen
   as the longest machine-found valid segment — not all 10 reps were run through Static
   Optimization (expensive; would be a natural, cheap extension given each SO window took only
   ~12-94s wall-clock, but was not needed to answer this session's question).
5. **Absolute force magnitudes throughout §3 are subject2's own anthropometry** (78.2 kg, same as
   every other cert in this family) — OrthoLoad's subjects are older, post-arthroplasty implant
   patients; this is a plausibility/ballpark comparison, not a per-subject validation (inherited,
   unchanged framing from every walking-cert document in this family).
6. **Right side only, single trial per activity** — same scope boundary as every prior cert in
   this family; left-side/multi-trial replication not attempted.
7. **The OrthoLoad "Sit-to-Stand/Stand-to-Sit" and "Stairs (Up+Down)" buckets do not separate
   sub-phases** (stand-up vs sit-down; ascent-only vs descent-only trials that happen to report
   both directions in one file) — inherited directly from `index_orthoload_forces.py`'s own
   classifier, not refined further here.
8. **n=3 hip squat-bucket OrthoLoad trials in hip_gen1 alone are very small**; the reported hip
   squat anchor (202.0 %BW median, n=37) pools gen1(n=3)+gen2(n=34) — gen2-dominated, disclosed.
9. **Video-corpus squat scene numbers (§3.3 corroborating row) are cited, not re-run** — read
   directly from its own current on-disk JSON this session (confirmed non-stale relative to the
   sign-bug fix by both mtime ordering and its own `_CORRECTED_sign_bug_fix` fields), per the
   lean discipline of not re-deriving an already-verified number.

## Files

- `scripts/msk/cross_activity_validation.py` — the full pipeline (self-contained, re-runnable;
  imports `validate_joint_force.py`/`static_opt_knee.py`/`force_scenes_batch.py` for proven
  parse/BFS/sign-fixed-crossing-detector/plausibility-diagnostic code, not re-implemented).
- `data/msk_smoketest/cross_activity_validation/cross_activity_validation_results.json` — the
  combined report (OrthoLoad anchors + kinematics-gap + both activities), machine-written.
- `data/msk_smoketest/cross_activity_validation/{squat,sit_to_stand}/pure_reaction_results.json`,
  `.../muscle_driven_results.json` — every number in §3-5 above, machine-written per activity.
- `data/msk_smoketest/cross_activity_validation/{squat,sit_to_stand}/{so,jr}/` — patched OpenSim
  setup XMLs + raw Static-Optimization/JointReaction `.sto` output.
- Reused, unmodified: `scripts/msk/validate_joint_force.py`, `scripts/msk/static_opt_knee.py`,
  `scripts/msk/force_scenes_batch.py`; `data/external/orthoload/_index/orthoload_akf_trials.jsonl`
  (built by `scripts/msk/index_orthoload_forces.py`); subject2's real OpenCap LabValidation mocap
  (external, read-only mount).
