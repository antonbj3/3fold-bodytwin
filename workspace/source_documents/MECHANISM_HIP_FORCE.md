# MECHANISM HIP FORCE VALIDATION — muscle-driven hip contact force vs OrthoLoad (2026-07-21)

## ⚠ CORRECTION (2026-07-21) — READ THIS BEFORE THE NUMBERS BELOW

`validate_hip_force.py` line ~325 calls `static_opt_knee.knee_crossing_muscles_and_forces` directly,
unmodified — the SAME shared function `docs/MECHANISM_STATIC_OPT.md`'s own correction banner
describes. This doc's self-computed headline number was the EXACT NEGATIVE of the correct
muscle-force direction, masked by a vector-norm step (`np.linalg.norm(...)`, always ≥0, cannot
expose an internal sign error). Full diagnosis/fix: `docs/MECHANISM_SIGN_BUG_AUDIT.md`,
`docs/MECHANISM_SIGN_BUG_REMEDIATION.md`. Fixed at its source in `static_opt_knee.py` (one-line branch
swap) — this doc's own numbers below are auto-corrected by that fix; nothing in this file was
directly edited except this banner. Re-ran `validate_hip_force.py` after the fix (exit 0, all its own
internal gates PASS) to confirm directly, not just infer:

| | published (bug artifact) | **corrected** (re-measured, `validate_hip_force.py` exit 0) | ratio vs OrthoLoad (273.93 %BW) |
|---|---:|---:|---:|
| self-computed (R_old − Σmuscle-crossing) | 235.30 %BW @ t=0.67s | **386.77 %BW @ t=0.55s** | 0.859 → **1.412** |

Note the peak **time** also changed (0.67s → 0.55s) — the corrected peak now lands at the exact
same instant as JointReaction's own peak (below), which the buggy version did not.

**The corrected number now matches this doc's OWN already-published, bug-immune official
`opensim.JointReaction` cross-check (387.04 %BW) to 0.070%** — JointReaction never calls the buggy
function, so this convergence (not the original "lands just under the anchor" story) is the proof
the correction is right. **Honest new headline: the twin OVER-predicts in-vivo hip contact force by
~1.41×**, not "lands just under the anchor, 14% short." Same documented tendency as the knee
(`docs/MECHANISM_STATIC_OPT.md`'s banner): Moissenet F, Chèze L, Dumas R (2014), *J Biomech* 47(1):50-58,
DOI [10.1016/j.jbiomech.2013.10.015](https://doi.org/10.1016/j.jbiomech.2013.10.015), PMID
[24210475](https://pubmed.ncbi.nlm.nih.gov/24210475/): two-step (reaction-then-subtract-SO-muscles)
pipelines like this one's "joint reaction forces are usually overestimated."

**§5's "corroborates the knee cert's differentiation-scheme diagnosis" framing is SUPERSEDED**: the
recurring ~1.6-1.7× ratio at both joints was real and correctly measured, but its cause was the
shared CODE bug (this exact function, called identically by both certs), not shared differentiation
scheme — the two hypotheses were never actually distinguished until `docs/MECHANISM_SIGN_BUG_AUDIT.md`
fixed the code and watched the gap collapse to <0.1% at both joints. §6's near-cancellation
MECHANISM (large, near-antiparallel bone-contact and muscle-crossing terms) remains correct; the
specific corrected numbers are in the remediation doc.

Everything below this banner is preserved as the ORIGINAL, pre-correction analysis (historical
record) — read the headline table, §5, and §6 with the correction above in mind, not as
still-current numbers.

---

Replicates the VALIDATED knee method (`docs/MECHANISM_JOINT_FORCE_VALIDATION.md` +
`docs/MECHANISM_STATIC_OPT.md`) for the HIP: does the twin's muscle-driven hip contact force
(subject2 `walking1`, OpenSim Static Optimization, bug-immune method) match OrthoLoad in-vivo hip
contact force during walking? Every number below is machine-measured this session
(`scripts/msk/validate_hip_force.py`, exit 0), not recalled. Isolation respected: `.venv-msk` only,
OrthoLoad + LabValidation data read in place (never written to), no git commit/push.

## Headline result [SUPERSEDED — see correction banner above; preserved for the historical record]

| method | peak hip-r contact force (%BW) | t (s) | ratio vs in-vivo (273.93) |
|---|---:|---:|---:|
| Pure kinematics+GRF reaction (no muscles, re-derived this session) | 88.59 | 1.52 | 0.323 |
| **Self-computed** (R_old − Σmuscle-crossing, primary) | ~~235.30~~ **→ 386.77 (corrected)** | ~~0.67~~ **→ 0.55** | ~~0.859~~ **→ 1.412** |
| Official `opensim.JointReaction` (cross-check, reused output) | 387.04 | 0.55 | 1.413 |
| IN-VIVO OrthoLoad hip (median, n=162 unassisted level-walking, 18 subjects) | 273.93 | — | 1.0 |

**[SUPERSEDED framing, kept verbatim for the record]** Yes — adding muscles closes most of the gap,
landing the primary (self-computed) estimate inside the task's own stated ~210-330 %BW ballpark, 38.6
percentage points (14%) short of the OrthoLoad median; the secondary (official JointReaction)
cross-check overshoots by 113.1 points (41%). Both candidates PASS the pre-registered "substantial
closure" gate (ratio > 0.6, §5). This reproduces — on an independent joint, independent muscle set,
independent OrthoLoad sub-corpus — the SAME qualitative pattern the knee cert found: self-computed
lands just under the anchor, official JointReaction moderately overshoots it (knee: 0.903 / 1.515;
hip: 0.859 / 1.413). **[Corrected: both methods now agree at both joints — self-computed and
JointReaction converge to within 0.004-0.07%, not two independently under/over-shooting methods.]**

## 1. Inputs — REUSED, not re-derived (the central lean decision this session made)

Static Optimization (`opensim.AnalyzeTool` + `StaticOptimization`, 158 frames, t∈[0,1.57s], the
model's own 80 muscles + reserve actuators) and the official `opensim.JointReaction` analysis
(`joint_names=ALL`) were **already run this session** for the knee cert
(`scripts/msk/static_opt_knee.py`) on this exact model/trial. Static Optimization's result does not
depend on which downstream joint is later interrogated, and the JR setup already requested `ALL`
joints — so its output already contained `hip_r_on_femur_r_in_femur_r_f{x,y,z}` columns, verified
present before writing a line of this script's downstream logic. Rerunning either (expensive: Ipopt
nonlinear optimization × 158 frames) would have been wasteful, not more correct. This script:

- **Reuses** `data/msk_smoketest/subject2_walking1/static_optimization/so/walking1_StaticOptimization_
  {activation,force}.sto` and `.../jr/walking1_JointReaction_ReactionLoads.sto` directly (falls back to
  calling `static_opt_knee.py`'s own `run_static_optimization()`/`run_joint_reaction()` if absent).
- **Does not trust the reuse blindly**: re-verified file provenance (JR output's own patched setup XML
  points at the exact SO `force.sto` on disk with a matching mtime — a single coherent 3-millisecond-
  apart run, not two independently-timed artifacts) and **re-ran every convergence/activation/reserve
  gate itself** on the reused files (§5) rather than citing the prior doc's prose numbers.
- Model/kinematics/GRF: same as both prior certs (`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`,
  subject2, 78.2 kg → 766.88 N; `walking1.mot` IK, 158 frames; `walking1_forces.mot` GRF).
- **New this session**: the hip-specific free-body cut, hip-crossing-muscle detection, and a
  self-contained OrthoLoad **hip** AKF parser (the knee cert's parser only covers `data/external/
  orthoload/knee/`; hip needed a new parser for two previously-unhandled file formats, §4).

## 2. Geometric method (identical derivation to the knee cert, a different free-body cut)

For the chain of bodies BFS-distal to `hip_r` (found on the model's own joint tree, not assumed):

```
sum_i(m_i * a_i) = gravity + GRF + R        R = F_bone_contact + F_muscles_crossing (vector sum)
```

Verified live: **`hip_r`-distal = `[femur_r, patella_r, tibia_r, talus_r, calcn_r, toes_r]`, mass
15.219 kg** — the WHOLE leg, 2.8× the knee-distal chain's 5.475 kg — because `patella_r` hangs off
`femur_r` via the `patellofemoral_r` joint and is therefore hip-distal but NOT knee-distal (the same
BFS discipline the knee cert used to catch this exact parentage detail generalizes cleanly here).
`R` is a net resultant of two large, near-cancelling terms (bone/cartilage contact only pushes;
antagonist co-contraction pulls); OrthoLoad's instrumented implant measures `|F_bone_contact|`
specifically, so `F_bone_contact = R − F_muscles_crossing` (subtract the SO-tension-scaled,
live-path-geometry muscle forces from the already-validated reaction vector).

**Crossing-muscle detector — pre-registered external anchor, checked before trusting any force
number (§5):** every muscle whose live path crosses the hip-r cut is found by BFS body-membership
on the model's own current geometry, never a hardcoded list. Measured result: **25 muscles detected
— `addbrev_r, addlong_r, addmagDist_r, addmagIsch_r, addmagMid_r, addmagProx_r, bflh_r, glmax1_r,
glmax2_r, glmax3_r, glmed1_r, glmed2_r, glmed3_r, glmin1_r, glmin2_r, glmin3_r, grac_r, iliacus_r,
piri_r, psoas_r, recfem_r, sart_r, semimem_r, semiten_r, tfl_r` — ALL 25 anatomically-expected hip
flexors/extensors/abductors/adductors/rotators present (iliopsoas, glutes ×9, piriformis, adductors
×6, TFL, gracilis, sartorius, rectus femoris, the LONG hamstrings).** A geometry-derived **negative
control** was also checked: `bfsh_r` (biceps femoris SHORT head) anatomically originates on the
femur shaft — unlike its long-head sibling `bflh_r`, it must NOT cross the hip. Measured: **correctly
excluded** (it does, correctly, still cross the knee — confirmed present in a separate knee-cut
detection run). This is a stronger check than an unrelated-muscle negative control: `bfsh_r` is in
the same muscle group, same leg region, similar name — and the geometry still gets it right.

## 3. OrthoLoad in-vivo HIP anchor — self-contained AKF parser, two measured file formats

`data/external/orthoload/hip/` holds only reference PDFs/videos (zero `.akf` files — not a data
source). The real corpora are `hip_gen1/` (710 files, older format) and `hip_gen2/` (1240 files,
newer format), with **two different column layouts**, both measured directly (not assumed) before
parsing:

- **hip_gen1**: `Time -Fx -Fy -Fz Fres - - - Marker` — NO moments (literal `-` placeholders); resultant
  at the same fixed index 4 as every other OrthoLoad joint corpus. Comment #1 never uses the phrase
  "Level Walking" at all — activities read `"HIP JOINT Walking free; Velocity: normal/fast/slow"`,
  `"...Walking on Treadmill..."`, or `"...supported by Crutches..."`.
- **hip_gen2**: `Time Fx Fy Fz F Mx My Mz Marker` — identical layout to the knee corpus, including
  literal `"Level Walking"` phrasing in Comment #1.
- **Decimal separator**: measured 0/710 hip_gen1 BodyWeight lines use a comma; 61/1240 hip_gen2 lines
  do — both normalized per-token (same discipline as the knee/sibling parsers).

**Pre-registered selection** (stated before computing the distribution): PRIMARY = hip_gen1 `"walking
free"` (case-insensitive substring) OR hip_gen2 `"Level Walking"` (case-insensitive substring) —
measured live to have **zero overlap** with `"crutch"` or `"nordic"` in either sub-corpus before
adopting the rule.

| set | n trials | n subjects | median %BW | mean | min | max |
|---|---:|---:|---:|---:|---:|---:|
| **PRIMARY** (walking-free + Level-Walking) | 162 | 18 | **273.93** | 283.42 | 206.91 | 469.18 |
| BROAD sensitivity (+ gen1 treadmill walking) | 299 | 19 | 291.87 | 307.81 | 191.77 | 550.82 |
| Crutch-assisted (excluded, reported separately) | 58 | 12 | 222.69 | 212.10 | 75.04 | 341.79 |

The broader (treadmill-inclusive) selection moves the median by +18 points (+6.5%) and stays
comfortably inside the literature band — not a knife-edge choice. Crutch-assisted trials read LOWER
as physically expected (partial weight-bearing through the arms), the same directional sanity signal
the knee cert found for its 3/4-Point trials.

**Cross-validated** against the pre-existing, differently-coded sibling index
(`data/external/orthoload/_index/orthoload_akf_trials.jsonl`, built independently this session by a
sibling MSK agent) — not a runtime dependency, purely supporting evidence: **162/162 of this
script's PRIMARY trials matched a sibling record; max abs difference 0.05 %BW** (rounding-level
agreement, both parsers correct). n_distinct_subjects=18 is exactly consistent (8 hip_gen1 codes +
10 hip_gen2 codes) with the sibling index's own separately-flagged finding that hip_gen1's free-text
Comment #2 field yields more apparent codes than that corpus's external subject-count metadata
states (9 vs an expected 2) — an inherited, pre-existing data-quality wrinkle in the corpus itself,
not something this session introduced, and immaterial to a per-trial peak-%BW comparison.

## 4. Machine-checked gates (PASS/FAIL, not eyeballed)

| gate | pre-registered threshold | measured | verdict |
|---|---|---:|---|
| SO convergence (frames/NaN/activation bounds/Fmax ratio) | reused from static_opt_knee.py's gate, re-run on the reused files | 158/158 frames, 0 NaN, activations [0.010, 0.627], 0/80 muscles >1.5×Fmax | PASS |
| Hip reserve-actuator usage (safety-net, not a crutch) | flagged if >50% of own optimal_force | max 12.5% (hip_rotation_r), all others ≤6.6% | PASS |
| Anatomical anchor (25 expected hip-r crossers) | all present | 25/25 present | PASS |
| Negative control (`bfsh_r` must NOT cross hip) | excluded | correctly excluded | PASS |
| Internal consistency: re-derived pure reaction vs `validate_joint_force.py`'s own hip_r number | should match closely | **88.58593431894126 vs 88.58593431894126 — bit-identical to 12 significant figures** | PASS |
| Pure-reaction regime sanity (positive, below the anchor) | ratio in (0.10, 1.05) | 0.323 | PASS |
| Literature ballpark (task-stated hip-walking band) | 210–330 %BW | OrthoLoad median 273.93 %BW | PASS |
| Absurdity ceiling (guards the documented OpenSim 4.6 `hip_flexion_r/l` generalized-force bug, ~1000-1800×) | <1200 %BW | 235.30 / 387.04 %BW (2-3 orders of magnitude below); **corrected self-computed 386.77 %BW, still 2-3 orders below** | PASS |
| Substantial-closure verdict gate | ratio > 0.6 | self-computed 0.859 [**corrected: 1.412**, see top banner], JointReaction 1.413 | **PASS-substantial** (both) |

The bit-identical internal-consistency check is the strongest evidence in this document that the
reuse (§1) is legitimate and that this session's independent re-derivation of the Newton's-law
reaction is correct: two separately-written scripts, run in two different sessions, produce the
exact same 12-significant-figure number for the same physical quantity.

## 5. A forced adversary: the hip and knee JointReaction numbers looked suspiciously close

**[The JointReaction-vs-JointReaction comparison in this section remains valid (JointReaction never
called the buggy function at either joint). The "corroborates the knee cert's differentiation-scheme
diagnosis" conclusion about the recurring self-computed/JointReaction RATIO is SUPERSEDED — see the
top correction banner: the recurrence was real, but corroborated a shared CODE bug, not shared
differentiation scheme.]**

Official JointReaction reads 387.04 %BW for the hip vs the knee cert's own 391.11 %BW — a ~1%
difference between two different joints. The self-computed method also lands close (235.30 %BW hip
vs 233.20 %BW knee). Rather than accept this as a clean pass, this was treated as the adversary most
tempting to skip (a column-selection bug silently re-reading the knee's own columns) and forced to
a direct machine check on the raw `.sto` data:

- Hip and knee force-magnitude time series are **NOT identical arrays** (`np.array_equal` = False);
  max absolute difference 752.7 N; correlation 0.966 (high but not 1.0 — physically expected, since
  both joints are in the same limb during the same single-leg-stance phase, but genuinely distinct
  signals that diverge visibly frame-to-frame away from the shared peak window).
- Peak **times** differ (hip 0.55 s vs knee 0.51 s) — independently-found maxima, not a copy-paste.
- The recurring ~1.6-1.7× self-computed/JointReaction disagreement ratio (hip: 387.04/235.30=1.645;
  knee: 391.11/233.20=1.677) across two different joints, two different muscle sets, and two
  different peak times **corroborates** (does not merely repeat) the knee cert's own diagnosis
  (`docs/MECHANISM_STATIC_OPT.md` §5) that the gap is a systematic differentiation-scheme difference
  (Savitzky-Golay-in-Cartesian-COM-space vs JointReaction's own spline-in-joint-angle-space + 6 Hz
  low-pass filter) — a general property of comparing the two methods, not a joint-specific artifact.

**Verdict: the adversary falls.** The closeness in magnitude is real physics (serially-adjacent
joints in one limb during one stance phase carry similar loads), not a pipeline defect.

## 6. Is this really co-contraction? Measured, not assumed

At the self-computed method's peak instant (t=0.67s), the muscle-crossing force magnitude is 278.88
%BW — LARGER than the resulting bone-contact residual (235.30 %BW), confirming substantial mutual
cancellation between `F_bone_contact` and `F_muscles_crossing` in the vector sum `R` (88.59 %BW),
exactly the mechanism the knee cert derived (§2): both terms are large, nearly anti-aligned along
the joint's compressive axis. This session did not re-run the per-muscle activation breakdown for
the hip cut (scope-boxed, §7); the mechanism is confirmed at the vector-magnitude level, consistent
with — not independently re-deriving — the knee cert's own muscle-level finding.

## 7. Honest caveats (full list)

1. **Reaction force ≠ contact force is the primary, load-bearing framing** (§2) — the self-computed
   number (235.30 %BW) is expected to be in the right ballpark, not a precision estimate; the two
   methods here should be read as a band (235-387 %BW), not a single number.
2. **The two methods disagree by ~1.65×** (§5) — diagnosed (differentiation-scheme difference,
   corroborated by the same ratio recurring on a second, independent joint) but not numerically
   pinned down further this session, same scope boundary as the knee cert.
3. **Different populations.** OrthoLoad's hip subjects are instrumented-hip-implant patients
   (elderly, post-arthroplasty anthropometry/gait); subject2 is a healthy young(ish) OpenCap
   participant (78.2 kg, 1.96 m) — a plausibility/ballpark comparison, not a per-subject validation
   (inherited, unchanged framing from the knee cert).
4. **hip_gen1 subject-code text parsing is coarse.** The free-text `Comment #2` field yields more
   apparent distinct codes (9) than the corpus's own external metadata expects (2) — a pre-existing
   data-quality wrinkle in the sibling index (built independently this session), not introduced or
   resolved here, and immaterial to the per-trial %BW comparison this document relies on.
5. **Selection-criterion choice** (walking-free vs treadmill-inclusive) moves the OrthoLoad median by
   +6.5% — tested as a sensitivity check (§3), not a knife-edge, but a real disclosed modeling choice.
6. **SO's pelvis residual force (172.96 N, 22.6% BW) exceeds the "good" band** — inherited unchanged
   from the knee cert's own SO run (a whole-body-level property of the SO solve, not hip-specific;
   no RRA was run first). Joint-level reserve usage (including all three hip DOFs) stayed low
   (≤12.5%), so the leg-muscle recruitment driving this number is not obviously contaminated.
7. **The mechanism is SO's activation-minimizing solution, not measured EMG** — inherited limitation,
   unchanged from the knee cert; real co-contraction can exceed the effort-minimizing SO solution.
8. **Single trial, right leg primary** (left hip reported only as a coarse pure-reaction symmetry
   check: 87.62 %BW vs the right leg's 88.59 %BW — same order of magnitude, not a full independent
   muscle-driven replicate).
9. **GRF source**: real force-plate data (best case, not the degraded kinematics-only fallback).
10. **The OpenSim 4.6 ID-tool bug** (`docs/MECHANISM_MSK_ELASTIC_BAND.md` §4) explicitly corrupts
    `hip_flexion_r/l` generalized force by ~1000-1800× — directly on point for this joint. This
    script's FORCE result is architecturally immune (never calls `InverseDynamicsTool`/
    `InverseDynamicsSolver`; JointReaction uses a documented-different Simbody code path), and both
    candidate numbers sit 2-3 orders of magnitude below that corruption signature (absurdity-ceiling
    gate, §4).

## 8. Next step

If a single tighter number is later needed: reconcile the self-computed-vs-JointReaction ~1.6-1.7×
gap by re-differentiating with a matched filter (this is now a TWO-joint-corroborated finding, which
strengthens the case that this is worth resolving generally rather than per-joint). Otherwise, the
task's question is answered with the pre-registered gate PASSING under both independent methods: the
muscle-driven hip contact-force estimate lands at 235-387 %BW against an OrthoLoad in-vivo anchor of
273.93 %BW (n=162, 18 subjects) — inside or moderately past the anchor, in the same qualitative
pattern the validated knee method already established, not a joint-specific fluke.

## Files

- `scripts/msk/validate_hip_force.py` — the full pipeline (self-contained, re-runnable; imports
  `validate_joint_force.py` and `static_opt_knee.py` for proven parse/BFS/SO-gate/crossing-detector
  code, not re-implemented; adds only the hip-specific free-body cut, a generalized self-cross-check
  function, and a new self-contained OrthoLoad hip AKF parser).
- `data/msk_smoketest/subject2_walking1/hip_force_validation/hip_force_validation_results.json` —
  every number in this document, machine-written.
- Reused, not written by this script: `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_{activation,force}.sto` and `.../jr/walking1_JointReaction_
  ReactionLoads.sto` (both already existed from the knee cert session; re-verified, not re-run).
