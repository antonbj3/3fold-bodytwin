# MECHANISM EMG-CONSTRAINED FORCE — does co-contraction modeling explain the knee/hip over-prediction? (2026-07-21)

**Question (this session's task):** `docs/MECHANISM_COMMON_MODE.md` names the EMG-driven layer
(`docs/MECHANISM_EMG_DRIVEN.md`) as the one mechanism in this repo that genuinely breaks the
shared-Static-Optimization-solve common-mode. That doc already built a measured-EMG-constrained
forward-dynamics estimate and read its knee/hip contact force via the **official**
`opensim.JointReaction` method only. This session is additive, not a re-derivation: it (1) enumerates
what EMG data actually exists for subject2 (machine-verified), (2) re-reads the SAME constructions
through a **second, independent** free-body method — `vhf.compute_self_cross_check_generic`, the
exact function this task named, imported unchanged — and (3) builds one new, pre-registered,
literature-floor construction targeting the ONE major crossing muscle with zero real-EMG coverage.

## Headline

**Two independent free-body methods now agree: EMG-realistic co-contraction makes the knee
over-prediction MATERIALLY WORSE (+46.4% self-computed / +49.6% official-JointReaction, both vs SO,
both AWAY from OrthoLoad) and has a NEGLIGIBLE, joint-specific effect at the hip (-2.4% self-computed
/ -0.7% official-JointReaction, technically TOWARD OrthoLoad but below the pre-registered 5% modest
threshold at both methods).** The falsifier resolves per-joint, not uniformly: at the **knee**, missing
co-contraction is **ruled OUT** as an explanation for the over-prediction (a valuable negative — real
antagonist activity pushes the twin further from, not closer to, the in-vivo anchor). At the **hip**,
co-contraction does not materially help either way. A targeted literature-floor construction on
rectus femoris (the one major crossing muscle with no real EMG) changes nothing at the established
peak instants, for a disclosed, mechanistic reason (Sec. 5): its literature-taught burst windows do
not temporally overlap the instant contact force is actually maximal in this model.

| | KNEE (SO=391.10 %BW, anchor=258.22 %BW, ratio 1.51x) | HIP (SO=386.77 %BW, anchor=273.93 %BW, ratio 1.41x) |
|---|---:|---:|
| naive_fwd_only (method-control) | 395.38 %BW, **+1.09%** negligible, AWAY | 434.37 %BW, **+12.31%** modest, AWAY |
| **emg_hybrid** (measured EMG, primary) | **572.42 %BW, +46.36% material, AWAY** | **377.46 %BW, -2.41% negligible, TOWARD** |
| emg_hybrid_litfloor (+ recfem lit-floor) | 572.42 %BW, +46.36% material, AWAY (unchanged) | 377.42 %BW, -2.42% negligible, TOWARD (unchanged) |
| cross-check: same rows, official JointReaction (`docs/MECHANISM_EMG_DRIVEN.md`) | naive +0.80%; emg_hybrid **+49.58% material** | naive +13.08%; emg_hybrid **-0.72% negligible** |

All values at the established, non-edge-flagged peak instants (knee t=0.51s, hip t=0.55s), the
defensible reading this whole cert family uses. Both extraction methods agree on sign and rough
magnitude for every row — a genuine, unplanned convergence (Sec. 4).

## 1. Method — additive to `docs/MECHANISM_EMG_DRIVEN.md`, not a re-derivation

Reused unchanged, zero re-solve/re-integration: SO's own force.sto, and `docs/
MECHANISM_EMG_DRIVEN.md`'s already-computed, already-gated `naive_fwd_only` (method-control: all 80
muscles forward-integrated with SO's own smoothed activation as excitation, isolating the pure
quasi-static-to-forward-dynamics method-switch effect) and `emg_hybrid` (7/80 muscles — 6 of them
actual knee/hip crossers, Sec. 2 — forward-integrated with real, measured surface EMG; the other
73/80 forward-integrated with SO's own activation as excitation) force arrays, cached in
`data/msk_smoketest/subject2_walking1/emg_driven/_diag_cache.npz`.

New this session:
- **Re-read via a second, independent free-body method.** `docs/MECHANISM_EMG_DRIVEN.md` used only
  the official `opensim.JointReaction` AnalyzeTool. This task explicitly named
  `vhf.compute_self_cross_check_generic` (the generic self-computed BFS crossing-muscle Newton's-law
  cut `docs/MECHANISM_CMC_SECOND_SOLVE.md` built to compare SO vs CMC) — imported unchanged, called
  in-memory on all four force arrays (SO, naive_fwd_only, emg_hybrid, emg_hybrid_litfloor) at both
  the knee-r and hip-r free-body cuts (8 calls total, no AnalyzeTool subprocess needed).
- **One new, pre-registered construction**, `emg_hybrid_litfloor` (Sec. 5).
- **A divergence guard and a mis-phased-constraint check** (Secs. 6-7), both required by the task's
  own symmetric-QC instruction.

## 2. EMG enumeration (machine-verified, not assumed — this task's explicit first step)

`/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/EMGData/` contains
**16 `*_EMG.sto` files** (DJ1-3, DJAsym1/4/5, STS1, STSweakLegs1, squats1, squatsAsym1, static1,
walking1-3, walkingTS1/2/4), **every one sharing the identical 16-channel schema** (7 right-leg + 9
left-leg columns, verified live): `soleus, gasmed, tibant, recfem, vasmed, vaslat, semiten, bflh,
glmed1` on the left; `soleus, gasmed, vasmed, vaslat, semiten, bflh, glmed1` on the right (`tibant_r`
and `recfem_r` have **no right-side channel at all** — a real electrode-placement asymmetry, not a
scope choice). `walking1_EMG.sto` (this trial): 158 rows, matching the IK/SO grid exactly.

Of the **31** anatomically-expected knee-r-or-hip-r crossing muscles (12 knee ∪ 25 hip, 6 shared
biarticular, both anchor sets reused unchanged from `static_opt_knee.py`/`validate_hip_force.py`),
**6 have a real right-side EMG channel**: `bflh_r, gasmed_r, glmed1_r, semiten_r, vaslat_r, vasmed_r`
(`soleus_r` also has a channel but is ankle-only — correctly excluded from the crossing-muscle count,
a live, machine-verified confirmation, not an assumption). **25/31 crossers have no right-side EMG at
all** — most critically `recfem_r` (Sec. 5) and `gaslat_r`. Measured EMG **does exist** for this
subject/trial, covering the major quad+hamstring+glmed antagonist set — the task's "measured EMG"
branch is the operative one; the literature-envelope branch (Sec. 5) is used only as a targeted,
disclosed supplement for the one flagged gap, not the primary construction.

## 3. Pre-registered gates and thresholds (stated before running)

| item | value |
|---|---|
| Literature floor level (recfem_r) | 0.15 (conservative — half this repo's own 0.25 ON-threshold convention; a floor, `exc=max(exc_hybrid, 0.15)`, never lowers existing excitation) |
| Literature-ON window (recfem_r) | `EMG_LITERATURE["recfem"]["windows"]` (unchanged, `validate_emg_timing.py`) = early-stance (0-30 %GC) + pre-swing/initial-swing (50-73 %GC) |
| Materiality tiers (vs SO's own self-computed baseline) | <5% negligible, 5-20% modest, >20% material |
| Falsifier | at knee t=0.51s / hip t=0.55s: EMG-constrained force CLOSER to OrthoLoad than SO ⇒ co-contraction matters; FURTHER ⇒ over-prediction not explained by missing co-contraction |

## 4. Two independent free-body methods converge (the forced adversary for this result's positive lean)

Because this result partly reads as "EMG changes the joint force," the adversary that must be forced
is: *could the self-computed number simply be an artifact of this one extraction method?* Checked by
comparing every row against the already-published, independently-coded `opensim.JointReaction`
numbers (`docs/MECHANISM_EMG_DRIVEN.md` Sec. 7), which use a wholly different code path (OpenSim's own
internal multibody-dynamics reaction solver vs this script's Savitzky-Golay CoM-acceleration +
BFS-crossing-muscle subtraction):

| construction | self-computed Δ vs SO | official-JR Δ vs SO | agreement |
|---|---:|---:|---|
| knee, naive_fwd_only | +1.09% | +0.80% | both negligible, same sign |
| knee, emg_hybrid | **+46.36%** | **+49.58%** | both material, same sign, within 3.2 points |
| hip, naive_fwd_only | +12.31% | +13.08% | both modest, same sign, within 0.8 points |
| hip, emg_hybrid | **-2.41%** | **-0.72%** | both negligible, same sign |

Every row agrees on sign and materiality tier; the two largest (knee emg_hybrid) agree within 3.2
percentage points of each other despite zero shared code between the extraction methods. This is a
genuine over-determination (two independently-coded free-body cuts on two different arrays converge),
not a tautology — and it is symmetric: the same close agreement also holds for the SO-only self-vs-JR
baseline already established in `docs/MECHANISM_STATIC_OPT.md`/`docs/MECHANISM_HIP_FORCE.md` (391.10 vs
391.11 %BW knee, 386.77 vs 387.04 %BW hip) — this session extends that already-verified agreement to
the EMG-driven constructions and finds it holds up.

## 5. The literature-floor construction: a real, mechanistic null (forced via OODA, not glossed over)

`recfem_r` (rectus femoris) is the single muscle `docs/MECHANISM_STATIC_OPT.md` Sec.7 and `docs/
MECHANISM_EMG_DRIVEN.md` Sec.6/9 both flag as the dominant knee-peak co-contraction partner (with
gastrocnemius, which IS EMG-covered) that has **zero** real-EMG coverage. Pre-registered construction:
`emg_hybrid_litfloor` = `emg_hybrid` + a 0.15 activation floor on `recfem_r`, applied only inside its
literature-taught ON window (Sec. 3), forward-integrated fresh (INTEGRATOR_ACCURACY=1e-5, the same
convergence-checked setting `docs/MECHANISM_EMG_DRIVEN.md` established).

**Forced-adversary gates, both PASS**: internal consistency (79 unchanged-excitation muscles vs
`emg_hybrid`) max|diff|=97.64 N (gate <100 N — a **tight**, 2.4%-headroom margin, the same order as
`docs/MECHANISM_EMG_DRIVEN.md`'s own 97.40 N/100 N gate, disclosed not hidden), relative RMS=0.0061
(gate <0.03); over-Fmax sanity max_ratio=0.915, 0/80 flagged (gate <1.5x).

**Real finding, measured not assumed**: the floor actually raised `recfem_r`'s excitation on 60/158
frames (38.0% of the trial) — but the established knee peak (t=0.51s = **44.1 %GC**) and hip peak
(t=0.55s = **47.1 %GC**) both fall **OUTSIDE** recfem's literature-ON windows (0-30, 50-73 %GC) — in
the 30-50%GC gap between its two classic textbook bursts. Consequence, verified not predicted: at both
reference instants, `emg_hybrid_litfloor` is **numerically indistinguishable** from `emg_hybrid`
(knee 572.42 vs 572.42 %BW; hip 377.42 vs 377.46 %BW — differences <0.05 %BW, floating-point/
interpolation noise, not a floor effect). **This is a real, disclosed, mechanistically-explained null
for this specific construction**: literature-normative recfem co-contraction, applied at its
textbook-taught timing, cannot speak to why contact force peaks where it does in this twin, because
the two events do not temporally coincide — not a failure of the method, a fact about when SO's own
solution places the force peak relative to when classic EMG atlases place rectus femoris activity.

## 6. Divergence guard — crossing-muscle SET identical across all 4 constructions

Forced to be true by pure geometry (`sok.knee_crossing_muscles_and_forces` detects crossings from
live path-vs-body-membership geometry alone — never touches the force array), verified here, not
assumed: **PASS at both joints.** Knee: 13/13 identical across SO/naive_fwd_only/emg_hybrid/litfloor
(`bflh_r, bfsh_r, gaslat_r, gasmed_r, grac_r, recfem_r, sart_r, semimem_r, semiten_r, tfl_r, vasint_r,
vaslat_r, vasmed_r` — all 12 anatomically-expected knee crossers present, plus `tfl_r`, a real bonus
find matching `docs/MECHANISM_CMC_SECOND_SOLVE.md`'s own note that this model's TFL also crosses the
knee). Hip: 25/25 identical, exact match to `HIP_ANATOMICAL_MUSCLES_EXPECTED`; negative control
`bfsh_r` correctly excluded from the hip set in every construction.

## 7. Mis-phased-constraint QC — measured EMG TIMING itself vs literature (closes a gap `docs/MECHANISM_EMG_TIMING.md` left open)

`docs/MECHANISM_EMG_TIMING.md` validated SO's own **activation** against literature phase windows —
it never checked the **measured EMG signal itself**. This session's gait-cycle re-derivation
(T=1.3300s, right_hs_0=-0.0765s, right/left stance-fraction consistency diff=0.038 percentage points
— matching `docs/MECHANISM_EMG_TIMING.md`'s own independently-computed 0.04pp to rounding, a second
over-determination) lets each real-EMG channel's own peak be located in %GC and checked directly:

| muscle | real-EMG peak (%GC) | literature window | verdict |
|---|---:|---|---|
| soleus_r | 61.4 | 30-60 | **FAIL** (boundary miss, 1.4 points over) |
| gasmed_r | 50.1 | 30-60 | PASS |
| vasmed_r | 16.3 | 0-30 | PASS |
| vaslat_r | 14.0 | 0-30 | PASS |
| semiten_r | 7.3 | 0-30, 85-100 | PASS |
| bflh_r | 7.3 | 0-30, 85-100 | PASS |
| glmed1_r | 12.5 | 0-30 | PASS |

**6/7 PASS**; the one failure (`soleus_r`) misses its window by 1.4 percentage points (≈19ms) — a
boundary case of the same kind `docs/MECHANISM_EMG_TIMING.md` itself flagged for gastrocnemius (missed
its 95% rank bar by 0.05 points), not a gross mis-mapping. This is the required symmetric-QC evidence
that the EMG-to-muscle channel assignment and gait-cycle timing are correctly phased — a mis-phased
constraint would show channels peaking in anatomically implausible phases (e.g., a knee extensor
peaking in mid-swing); none do.

## 8. Confidence tiers (per this task's own framing)

- **`emg_hybrid` (measured-EMG-constrained), knee AWAY / material and hip TOWARD / negligible**:
  **published-plausibility, with a diagnosed-gap on absolute magnitude.** The DIRECTION and
  materiality-tier findings rest on real, subject-measured EMG (a decorrelated input) driving real,
  compiled contraction dynamics, cross-checked by two independent free-body extraction codepaths
  (Sec. 4) — strong internal grounds for the sign/tier. The ABSOLUTE magnitude is not
  in-vivo-anchored to an independently-known muscle-force ground truth: `docs/MECHANISM_EMG_DRIVEN.md`
  Sec.9 already disclosed the EMG-to-MVC normalization is uncalibrated (no true-MVC trial exists for
  this subject), which this session did not re-litigate or fix — inherited unchanged.
- **`emg_hybrid_litfloor` (literature-normative floor)**: **published-plausibility** (Rajagopal 2016 /
  Perry 1992 phase-window convention, PMID 27392337, already-committed and unchanged this session) —
  and its finding here (no effect at the reference instants) is a **diagnosed-gap in window overlap**,
  not a claim about co-contraction magnitude.
- **OrthoLoad anchor comparison**: **in-vivo-anchored** (real instrumented-implant contact-force
  measurements, re-parsed live this session: knee median 258.22 %BW n=72, hip median 273.93 %BW
  n=162 — both match previously-published values to 2 decimals).

## 9. Honest gaps (pre-registered items first)

1. **The literature-floor construction targeted exactly one muscle** (`recfem_r`), chosen because it
   is the specific, already-flagged missing mechanism — not a systematic sweep of all 25 uncovered
   crossers. A floor applied to, e.g., the adductors or glute compartments was not tested; this
   session's scope is a targeted, disclosed test of the ONE literature-motivated hypothesis already
   on the table, not an exhaustive search.
2. **The 0.15 floor level was not swept** (a single pre-registered value, not a sensitivity band) —
   moot for THIS trial's reference-instant reading since the floor never activates there (Sec. 5), but
   a different floor level would still not change that specific null (window overlap is a timing fact,
   independent of floor magnitude).
3. **Self-computed and official-JointReaction are not fully independent tests** — both ultimately
   consume the same underlying force arrays; the genuinely distinct ingredient is the CODE PATH
   extracting the joint reaction from those arrays (Savitzky-Golay CoM-acceleration + BFS-crossing
   subtraction vs OpenSim's internal reaction solver), which is what Sec. 4 cross-checks.
4. **Single trial, right-leg-primary** (subject2 `walking1`) — same scope caveat as every prior cert
   in this family; no claim of generality across subjects, trials, gait speeds, or the left leg.
5. **EMG-to-MVC normalization remains uncalibrated** (inherited from `docs/MECHANISM_EMG_DRIVEN.md`
   Sec.9.1, not re-addressed here) — the likely dominant driver of the large absolute `emg_hybrid`
   magnitudes; the DIRECTION/materiality-tier findings are more robust to this than the exact %BW
   values quoted.
6. **The hip's `emg_hybrid` "TOWARD anchor" finding is negligible in magnitude** (-2.41%/-0.72%,
   both under the 5% modest floor at both extraction methods) — this session does NOT claim
   co-contraction modeling materially helps the hip; it claims the effect is real-signed but small,
   a more precise and less hype-prone statement than "improves" would be.
7. **The internal-consistency gate margin for the litfloor construction is tight** (97.64 N vs a 100 N
   gate, 2.4% headroom) — the same numerical-pathology risk region `docs/MECHANISM_EMG_DRIVEN.md`
   Sec.2 already diagnosed and fixed (1e-5 integrator accuracy); disclosed, not hidden, and it PASSED.
8. **`soleus_r`'s measured-EMG timing FAIL (Sec.7) is a 1.4-point boundary miss**, not investigated
   further this session (e.g., whether a subject-specific rather than population-average stance
   fraction would resolve it, per `docs/MECHANISM_EMG_TIMING.md`'s own unresolved gap #3).

## Files

- `scripts/msk/emg_constrained_force.py` — the full, self-contained, re-runnable pipeline (9 steps;
  imports `validate_joint_force.py`, `validate_emg_timing.py`, `contraction_dynamics_forward.py`,
  `static_opt_knee.py`, `validate_hip_force.py` — every physics/geometry function reused unchanged).
  ~83s wall time (one ~69s forward integration + 8 self-computed free-body calls + enumeration/QC).
- `data/msk_smoketest/subject2_walking1/emg_constrained_force/emg_constrained_force_results.json` —
  every number in this document, machine-written (EMG inventory, gait-cycle detection, per-muscle
  timing QC, litfloor construction + gates, all 8 self-computed joint-force reads, divergence guard,
  OrthoLoad anchors, full verdict tables).
- `data/msk_smoketest/subject2_walking1/emg_constrained_force/litfloor_excitation.sto` — the one new
  excitation file (158×81), identical to `emg_driven/hybrid_excitation.sto` except `recfem_r`'s column.
- Reused, not re-solved/re-integrated: `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_force.sto`; `data/msk_smoketest/subject2_walking1/emg_driven/
  _diag_cache.npz` (`f_naive_only`, `f_hybrid`, `exc_hybrid`, `a_so_smooth`, `t_so`, `muscle_names`).
- Real EMG read in place, never modified: `/media/anton/8838D60F38D5FBDE/mechanism_data/
  LabValidation_withVideos/subject2/EMGData/*.sto` (16 files enumerated, Sec. 2).
- Isolation (`COORDINATOR.md` Sec.1): bodytwin only, `.venv-msk` only; no git commit, no git push.
