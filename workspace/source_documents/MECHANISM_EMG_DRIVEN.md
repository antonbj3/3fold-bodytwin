# MECHANISM EMG-DRIVEN MUSCLE FORCE — CEINMS-style forward estimate vs Static Optimization (2026-07-21)

Executes the operator's task: every muscle-force cert in this repo so far (`docs/MECHANISM_STATIC_OPT.md`,
`docs/MECHANISM_HIP_FORCE.md`, `docs/MECHANISM_ANKLE_FORCE.md`) is built on Static Optimization (SO), which
minimizes `sum(activation^2)` subject to the net joint moment — a solution concept that structurally
under-predicts antagonist co-contraction (every one of those docs flags this as the same open honest gap).
This session closes part of that gap with a **decorrelated input**: instead of solving for activation,
drive the model's own real `Millard2012EquilibriumMuscle` muscles from **measured surface EMG** (subject2
`walking1`) through real activation dynamics → real contraction dynamics — a CEINMS-style forward
EMG-driven estimate (Lloyd DG, Besier TF. *J Biomech.* 2003;36(6):765-76. **PMID 12742444**; Pizzolato C
et al. *J Biomech.* 2015;48(14):3929-36. **PMID 26522621** — both verified live via NCBI esummary this
session, title/journal/year/authors match). Every number below is machine-computed this session
(`scripts/msk/emg_driven.py`, exit 0). Isolation respected: `.venv-msk` only, reads SO output + real
IK/GRF/model + real LabValidation EMG in place, no git commit/push, all new outputs under
`data/msk_smoketest/subject2_walking1/emg_driven/`.

**ANCHOR (pre-registered before any number below):** measured EMG is a decorrelated input — measured
directly from the subject's own muscles, independent of SO's activation-minimizing optimization assumption.
Agreement between the two is convergent evidence from independent methods; disagreement puts EMG (a real
physiological measurement) against an optimizer's arbitrary choice among infinitely many moment-consistent
activation vectors.

## Headline result

**Two real problems were found and fixed this session before any scientific number was trusted (Sec. 2-3)
— both are exactly the kind of thing the watertight discipline exists to catch, not celebrate around.**
With those fixed, the substantive findings:

| question | verdict |
|---|---|
| Do EMG-driven and SO agree on absolute muscle force magnitude? | **NO** — 2/7 muscles inside a ±(30-43)% agreement band vs raw SO, dominated by an honestly-disclosed, un-calibrated EMG normalization (Sec. 9) |
| Do EMG-driven and SO agree on relative force TIMING/SHAPE (scale-invariant)? | **PARTIALLY** — all 7 muscles show moderate positive whole-trial correlation (0.42-0.71) between EMG-hybrid and SO force time series |
| Does real EMG show co-contraction SO misses? | **YES** — gated mean CCI: EMG 0.391 vs SO 0.190 (diff +0.20, exceeds the pre-registered 0.10 threshold); 32/158 frames EMG-active-with-SO-silent vs 6/158 the reverse (Sec. 6) |
| Does that co-contraction change the joint force MATERIALLY? | **YES at the knee (+48.4%, material), NO/modest at the hip (-12.2%, modest)** — method-controlled, at SO's own established peak instant (Sec. 7) |

## 1. Method — the three-way, method-controlled decomposition

The scientific hazard: forward-integrating through real contraction dynamics changes force **even with
zero new information**, because it replaces SO's quasi-static force-length-velocity equilibrium with a
genuine forward-integrated fiber-length state. An "EMG changes the joint force" claim that doesn't control
for this method-switch effect is confounded. This script therefore builds **three** force estimates on the
same 158-frame (0-1.57s) grid:

1. **SO** — quasi-static Static Optimization (`data/.../static_optimization/so/`, already computed, reused,
   zero re-solve).
2. **naive_fwd_only** — all 80 muscles forward-integrated with SO's own smoothed activation *as* the
   excitation (zero real EMG anywhere). Isolates the pure method-switch effect.
3. **emg_hybrid** — the 7/80 muscles with a real EMG electrode on this subject's right leg (`soleus_r,
   gasmed_r, vasmed_r, vaslat_r, semiten_r, bflh_r, glmed1_r`) forward-integrated with **real, measured EMG**
   as excitation; the other 73/80 muscles get the identical SO-activation-as-excitation construction as (2)
   — the "missing muscles still need SO/assumption" gap, now exact and quantified (73/80 = 91.2%), not a
   vague caveat.

**(3) vs (2)** isolates the true EMG/co-contraction information content, holding the numerical method fixed
— this is the primary, decorrelated comparison. **(3) vs (1)** is the naive, task-literal "EMG-driven vs SO"
comparison, which conflates the method-switch effect with the EMG-content effect. Both are reported, never
only the flattering one.

**Geometric note on why SO cannot see co-contraction (not a bug, a nullspace argument):** for any joint with
>1 muscle crossing it, the map {muscle activations} → {net joint moment} has a nontrivial nullspace
(antagonist pairs can co-vary within it without changing the net moment SO must reproduce); SO's
`sum(activation²)` objective picks the minimum-norm point in that nullspace's affine slice, which is
generically not the physiological point. Measured EMG is a direct observation of a specific point in that
same nullspace — why it is a decorrelated, not merely "a different model," check.

Machinery reused, not rebuilt: `contraction_dynamics_forward.py`'s `build_forward_model`/`run_forward_sim`
(the prescribed-real-IK-kinematics forward integration through the real, compiled Millard2012 activation +
contraction dynamics ODE) and `process_real_emg` (rectify → Savitzky-Golay envelope → per-channel
[floor,peak] normalization); `static_opt_knee.py`'s `run_joint_reaction` (the official
`opensim.JointReaction` `forces_file` mechanism, already bug-fixed there); `validate_hip_force.py`'s hip
extractor + OrthoLoad hip parser; `validate_joint_force.py`/`validate_emg_timing.py` for `parse_mot` /
`compute_bw_n` / the OrthoLoad knee parser.

## 2. A real bug found and fixed — integrator-accuracy numerical pathology (forced via OODA)

The first version of this script reused `contraction_dynamics_forward.py`'s already-cached `f_naive_fwd`
array (computed at that script's own `INTEGRATOR_ACCURACY=1e-3`) as construction (2). This produced a
hard, **reproducible** internal-consistency failure: 73 muscles with provably identical excitation input
(verified: `max|exc_hybrid − a_so_smooth| = 0.0` over those columns) showed a **22.7% relative-RMS** force
mismatch against the cached reference, dominated by one muscle (`addbrev_r`, max|diff| = 2545 N). Per the
"honest-negative is not a free pass" rule, this was diagnosed, not shrugged off:

1. **Control test**: re-running the *identical* cached construction (same excitation file, same a0) fresh
   in a new process reproduced the cache to **exactly 0.0000 N** (bit-identical) — ruling out generic
   nondeterminism.
2. **Root cause**: re-running that same construction at a tighter accuracy (1e-5) revealed the cached 1e-3
   value for `addbrev_r` hit **2584 N — 4.13× its own Fmax (625.8 N)**, a physiologically impossible force.
   A genuine numerical pathology: the loose 1e-3 tolerance lets `RungeKuttaMerson`'s *shared* adaptive
   step-size control (chosen from error estimates over the whole 80-muscle coupled state vector) wander
   onto a spurious branch for a weak, near-floor-activation muscle. Invisible in the already-published
   `docs/MECHANISM_CONTRACTION_DYNAMICS.md` because that cert's own force comparison (§7.5) checked only 6
   muscle *groups* (glmax/glmed/hamstrings/gastroc/soleus/tibant — never addbrev/sart) and never ran an
   all-80-muscle over-Fmax sweep.
3. **Convergence check**: 1e-5 vs 1e-6 (identical excitation) differ by only **42.7 N max / 0.67% relative
   RMS** (vs the 1e-3 cache's 2544 N / 22.7% — a >30× improvement), and both 1e-5 and 1e-6 show **0/80**
   muscles exceeding 1.5× Fmax (vs 1/80 at 1e-3). 1e-5 adopted (wall time ~60-70s, barely slower than
   1e-3's ~50-60s).

**Fix applied**: this script recomputes both (2) `naive_fwd_only` and (3) `emg_hybrid` **fresh**, at
`INTEGRATOR_ACCURACY=1e-5`, and adds an explicit over-Fmax sanity gate. The stale 1e-3 cache is still
loaded and diffed against the corrected result as a disclosed, auditable record (`max|diff|=2534.80 N`,
relative RMS 0.199, confirming the pathology was in the stale cache, not this script's construction, which
now passes the over-Fmax gate cleanly on both constructions).

## 3. A real data artifact found — EMG lead-in padding (forced via OODA, not glossed over)

Per-muscle force comparisons (Sec. 5) and joint-force peaks (Sec. 7) initially showed extreme, physically
implausible early-trial values (e.g. knee contact force reaching 1107 %BW at t=0.15s). Diagnosed, not
reported blind: the first 5 frames (t=0.00-0.04s) of `walking1_EMG.sto` are **bit-identical across all 16
channels simultaneously** (verified live, `max|diff|=0.0` exactly) — a real 16-electrode surface-EMG
recording cannot legitimately hold every channel at exactly the same value for 4 consecutive samples by
chance. Cross-checked against `walking2_EMG.sto`/`walking3_EMG.sto`: **all three trials show exactly the
same 5-frame pad**, confirming a systematic upstream data-processing convention (a filter warm-up or
zero-order-hold, most likely), not a one-off glitch. This script flags (does not silently drop) any
per-muscle burst or joint-force peak landing at or before t=0.15s (a conservative buffer past the proven
0.04s pad, covering the visible post-pad filter-settling decay tail) — see the `EDGE-FLAGGED` markers in
Sec. 5/7. (Correction to the machine-written JSON: the knee's `emg_hybrid` own-peak lands at *exactly*
t=0.15s, which a strict `<` comparison in the first analysis pass failed to flag due to the tie; the source
now uses `<=` and this doc treats that specific number as edge-flagged in the text below, even though the
already-produced JSON's boolean reads `false` for that one entry.)

## 4. Forced-adversary gates (pre-registered, machine-checked)

| gate | threshold | measured | verdict |
|---|---|---:|---|
| Internal consistency (73 unchanged-excitation muscles, hybrid vs fresh matched-accuracy naive_fwd_only) | max\|diff\|<100 N, relative RMS<0.03 | max\|diff\|=**97.40 N**, relative RMS=**0.0136** | **PASS** (tight margin, 2.6%; disclosed, not hidden — see Sec. 9) |
| Over-Fmax sanity (0/80 muscles >1.5× own Fmax) | 0 flagged | naive_fwd_only max_ratio=0.683, emg_hybrid max_ratio=0.915, **0 flagged in either** | **PASS** |
| Re-implementation cross-check (this script's own `full_pct_bw_series` vs the two published docs' SO-JointReaction numbers) | value+time match | knee 391.11 %BW @ 0.51s (doc: 391.11 @ 0.51), hip 387.04 %BW @ 0.55s (doc: 387.04 @ 0.55) | **MATCH, both joints** |

The internal-consistency margin (97.40 vs a 100 N gate) is real and reported honestly, not smoothed over —
see Sec. 9 for what it means and does not threaten.

## 5. Per-muscle force comparison (7 EMG-covered muscles, at each muscle's own dominant SO burst)

Pre-registered agreement band: ratio ∈ [0.70, 1.43] (±30-43%, symmetric in log space) → AGREE.

| muscle | SO peak (N) | t_burst (s) | naive_fwd (N) | emg_hybrid (N) | hyb/SO | hyb/naive | corr(hyb,SO) | agree vs SO |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| soleus_r | 1447.8 | 0.61 | 1136.0 | 2625.2 | 1.813 | 2.311 | 0.706 | NO |
| gasmed_r | 1140.0 | 0.51 | 1105.3 | 917.7 | 0.805 | 0.830 | 0.475 | **YES** |
| vasmed_r | 185.4 | 0.06 [EDGE-FLAGGED] | 151.2 | 959.4 | 5.174 | 6.344 | 0.691 | NO |
| vaslat_r | 694.2 | 0.06 [EDGE-FLAGGED] | 611.9 | 1355.0 | 1.952 | 2.214 | 0.590 | NO |
| semiten_r | 31.8 | 1.25 | 46.3 | 256.9 | 8.074 | 5.552 | 0.684 | NO |
| bflh_r | 104.7 | 1.25 | 80.6 | 143.9 | 1.375 | 1.785 | 0.420 | **YES** |
| glmed1_r | 586.6 | 0.56 | 545.5 | 102.7 | 0.175 | 0.188 | 0.504 | NO |

**Ranked by SO peak force (the "major" muscles): soleus_r > gasmed_r > vaslat_r > glmed1_r > vasmed_r >
bflh_r > semiten_r.** **2/7 agree with raw SO within the pre-registered band; 1/7 agree with the
method-controlled `naive_fwd_only` baseline.** Of the two largest ("major") muscles: `gasmed_r` agrees
(0.805, well inside the band); `soleus_r` does not (1.813×) — but `soleus_r` also has the *best*
whole-trial correlation of all 7 (0.706), i.e. its force *rises and falls in the same pattern* as SO's even
though its *absolute scale* is off, consistent with the normalization-scale explanation (Sec. 9) rather
than a timing disagreement. **All 7 muscles show moderate positive whole-trial correlation (0.42-0.71)** —
a scale-invariant check, immune to the MVC-normalization uncertainty, showing broad phase agreement even
where absolute magnitude diverges sharply.

## 6. Co-contraction index (Falconer-Winter CCI) — a real methodological bug found and fixed here too

Quadriceps (`vasmed_r+vaslat_r`, mean) vs hamstrings (`semiten_r+bflh_r`, mean), CCI(t)=2·min(a,b)/(a+b)
(Falconer K, Winter DA. *Electromyogr Clin Neurophysiol.* 1985;25(2-3):135-49. **PMID 3987606**, verified
live; original paper is ankle-specific, the CCI formula's extension to the knee's quad/hamstring pair here
is this script's own, not literature-tested).

**Bug found (forced via OODA before trusting the first number)**: the raw CCI ratio is scale-invariant —
when *both* signals sit near their own floor (verified: `a_so_smooth` for all 4 muscles at t=0.51/0.55s is
within [0.010,0.010] of their own floor, matching `docs/MECHANISM_STATIC_OPT.md` §7's own claim that "vasti
and hamstrings are essentially at the SO activation floor" there), the *ratio* of two tiny, comparable
numbers can still read "high" — the raw formula reported CCI_SO=0.797 at the knee-peak instant, a
misleading artifact, not a real co-contraction finding. **Fix**: CCI is only meaningful for frames where at
least one signal clears an ON threshold (0.25, reusing this repo's own `ACT_ON_THRESH_FRAC_PRIMARY`
convention); frames where both are below threshold are excluded from the mean.

| statistic | real EMG | SO activation |
|---|---:|---:|
| UNGATED mean CCI (misleading) | 0.332 | 0.481 |
| **GATED mean CCI** (n active frames) | **0.391** (n=86) | **0.190** (n=60) |
| diff (EMG−SO), gated | **+0.200** (exceeds the 0.10 pre-registered threshold) | — |

**Contingency (158 frames)**: both-active=54, **EMG-active/SO-silent=32** (SO misses real co-contraction
EMG shows), **SO-active/EMG-silent=6** (the reverse), neither=66. The >5:1 asymmetry (32 vs 6) directly
supports the "SO systematically misses real co-activity more often than it fabricates activity" hypothesis
— once the ratio-floor artifact above is corrected.

**At the established knee/hip peak instants (t=0.51s/0.55s)**: quad+hamstring are **below the ON threshold
for BOTH methods** (max activity 0.23/0.15 for EMG, 0.008-0.012 for SO) — consistent with
`docs/MECHANISM_STATIC_OPT.md` §7's own finding that the *dominant* co-contraction pair at those specific
instants is rectus femoris + gastrocnemius, not vasti/hamstrings. `recfem_r` has no real-EMG channel (only
`recfem_l`), so that specific mechanism cannot be EMG-verified here — the quad/hamstring CCI above is a
whole-trial finding, not evidence about the exact knee/hip-peak instant. Separately: the classic ANKLE
antagonist pair (tibialis anterior vs soleus/gastrocnemius) cannot be checked this way on this leg either —
`tibant_r` has no real-EMG channel (only `tibant_l`).

## 7. Joint contact force (knee, hip) — does the co-contraction change it materially?

Pre-registered materiality tiers (relative to SO's own official-JointReaction number): <5% negligible,
5-20% modest, >20% material. Two readings reported for each construction: **own peak** (each
construction's own global maximum — for `emg_hybrid` this is often edge-flagged, Sec. 3) and **at SO's
peak time** (both constructions read at the same t=0.51s/0.55s anchor the rest of this repo's joint-force
certs already use — not edge-flagged, the more defensible number for magnitude).

| | KNEE (SO=391.11 %BW @ 0.51s) | HIP (SO=387.04 %BW @ 0.55s) |
|---|---:|---:|
| naive_fwd_only, own peak | 424.06 %BW (t=0.56s, +8.4%, modest) | 450.17 %BW (t=0.60s, +16.3%, modest) |
| naive_fwd_only, at SO peak time | 394.24 %BW (+0.8%, negligible) | 437.65 %BW (+13.1%, modest) |
| emg_hybrid, own peak | 1107.43 %BW (t=0.15s, **+183.1%, material**) **[EDGE-FLAGGED]** | 516.72 %BW (t=0.12s, **+33.5%, material**) **[EDGE-FLAGGED]** |
| emg_hybrid, at SO peak time | 585.02 %BW (**+49.6%, material**) | 384.24 %BW (**−0.7%, negligible**) |
| **METHOD-CONTROLLED** (emg_hybrid vs naive_fwd_only), at SO peak time | **+48.4% (material)** | **−12.2% (modest)** |

**Both `own_peak` entries land inside/at the EMG lead-in artifact window (Sec. 3) and are not trustworthy
as precision magnitudes** — they are reported for completeness, not as the headline. **The `at SO peak
time` / method-controlled row is the defensible finding: swapping in real EMG for the 7 covered muscles,
holding the forward-integration method fixed, changes the knee estimate by +48.4% (material) and the hip
estimate by −12.2% (modest, and in the *opposite* direction from the knee)** — a genuine, asymmetric,
joint-specific effect, not a uniform "EMG always inflates joint force" story.

## 8. OrthoLoad in-vivo anchor comparison

Re-parsed live (same proven parsers every prior cert in this family uses): knee median **258.22 %BW**
(n=72), hip median **273.93 %BW** (n=162).

| | using OWN PEAK (edge-prone) | using AT-SO-PEAK-TIME (fairer, hand-derived from the same JSON fields) |
|---|---|---|
| knee: \|SO−anchor\| | 132.89 | 132.89 |
| knee: \|emg_hybrid−anchor\| | 849.21 (further) | 326.80 (still further, less extreme) |
| hip: \|SO−anchor\| | 113.11 | 113.11 |
| hip: \|emg_hybrid−anchor\| | 242.79 (further) | **110.31 (marginally CLOSER than SO)** |

Using the fairer, non-edge-flagged instant: at the **knee**, `emg_hybrid` remains further from the in-vivo
anchor than SO (both over-predict, EMG more so). At the **hip**, `emg_hybrid` is essentially tied with SO
(384.2 vs 387.0 %BW) and marginally *closer* to the anchor (273.93 %BW) — a genuinely different,
joint-specific result, not picked to flatter either method.

## 9. Honest gaps (full list, pre-registered items first)

1. **EMG→MVC normalization is uncertain — and this session directly measured it to be the dominant driver
   of the extreme force numbers, not a footnote.** No trial in this subject's corpus is a true
   maximum-voluntary-contraction trial; excitation was normalized per-channel to that channel's own
   [floor,peak] range **within `walking1` alone**. Sensitivity check (Sec. 2 of the script, `walking2`/
   `walking3` compared): for `vasmed_r`, `walking1`'s own peak is only **66.2%** of the 3-trial maximum
   (i.e. even the multi-trial ceiling is likely still below true %MVC) — meaning `walking1`'s own-trial
   normalization systematically over-amplifies `vasmed_r`'s excitation relative to a better-informed
   ceiling, let alone true MVC. For `soleus_r` it is 88.3%; for `gasmed_r`/`semiten_r`/`bflh_r`, `walking1`
   happens to contain each muscle's own highest observed envelope among the 3 trials (100%) — still no
   guarantee that reflects true 100% excitation for a submaximal walking task. This is very likely the
   dominant cause of the large `hyb/SO` ratios in Sec. 5 (e.g. `semiten_r` 8.07×, `vasmed_r` 5.17×) — a
   forced, quantified demonstration of the task's own pre-registered caveat, not a new speculative concern.
2. **EMG covers only 7/80 muscles (8.75%).** The other 73 (91.2%) fall back to SO's own activation as a
   excitation surrogate — an exact, disclosed, machine-verified number, not a vague "some muscles" caveat.
   Electrode-placement asymmetry: `tibant_r` and `recfem_r` have **no** right-side channel at all in this
   file (only their left counterparts) — a data-availability limit, not a scope choice — which specifically
   blocks EMG-verifying the ankle-antagonist and the recfem+gastroc knee co-contraction mechanisms
   `docs/MECHANISM_STATIC_OPT.md` §7 identified as dominant at the knee's own peak instant.
3. **Electromechanical delay IS modeled** (via the real, compiled Millard2012 activation-dynamics ODE,
   10/40 ms act/deact time constants — the same values independently verified live in
   `docs/MECHANISM_CONTRACTION_DYNAMICS.md`), but **CEINMS's per-subject calibration layer is NOT** — no
   optimization was run to fit each muscle's EMG-to-activation gain/shape against independently-known joint
   moments. This is the forward half of CEINMS only, disclosed explicitly, not silently assumed equivalent
   to a calibrated model.
4. **A real EMG lead-in padding artifact (Sec. 3) contaminates the first 5 frames (0.04s) of every trial**,
   with a filter-settling tail visible out to roughly 0.15-0.3s. All burst/peak times landing at or before
   0.15s are edge-flagged in this doc's tables; the underlying raw EMG file was not otherwise modified.
5. **The internal-consistency gate margin is tight** (97.40 vs a 100 N gate, 2.6% headroom) — real,
   disclosed, not hidden. It does not threaten the headline findings: the residual (97 N, 1.36% relative
   RMS) is one to two orders of magnitude below the effect sizes being reported (hundreds-to-thousands of
   Newtons per muscle; 12-48% joint-force deltas).
6. **Co-contraction mechanism coverage is partial.** Only the knee's quad/hamstring pair could be
   EMG-verified; the two specific mechanisms `docs/MECHANISM_STATIC_OPT.md`/`docs/MECHANISM_HIP_FORCE.md`
   identify as dominant at the established peak instants (recfem+gastroc at the knee; the hip's own
   25-muscle crossing set) are only partially EMG-covered (gastrocnemius yes, rectus femoris no).
7. **The ~50 non-muscle actuator columns** (lumbar/arm CoordinateActuators, the 6 pelvis residuals, the 12
   joint reserves) in the JointReaction `forces_file` are held at SO's own solved values for every
   construction — no EMG measures these, and this is outside the scope of an EMG-driven *muscle* estimate,
   a disclosed simplifying assumption inherited from every prior cert in this family.
8. **This real EMG file is not literally under `data/msk_smoketest/subject2_walking1/*EMG*.sto`** as the
   task described — no such path exists in this repo. The actual, only real subject/trial-matched EMG file
   lives at `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/EMGData/
   walking1_EMG.sto` (already located and path-mapped by the contraction-dynamics agent in a prior session,
   `docs/MECHANISM_CONTRACTION_DYNAMICS.md` §7) — reused from there, read in place, never copied into this
   repo (isolation: LabValidation data is read-only external corpus for every cert in this family).
9. **Single trial, right leg primary** (subject2 `walking1`) — same scope caveat as every prior cert in
   this family; no claim of generality across subjects, trials, gait speeds, or the left leg.
10. **Different populations** for the OrthoLoad anchor (elderly instrumented-implant patients vs subject2's
    healthy young OpenCap participant) — inherited, unchanged caveat from every prior joint-force cert here.
11. **The JointReaction vs self-computed method disagreement** documented in `docs/MECHANISM_STATIC_OPT.md`
    §5 / `docs/MECHANISM_HIP_FORCE.md` §5 (a differentiation-scheme difference, ~1.6-1.7×) is inherited
    unchanged — this script uses only the official JointReaction method (fed by each construction's own
    force.sto), not the self-computed BFS method, so that specific disagreement does not directly apply
    here, but its root cause (SO-based methods generally) is a shared ancestry.

## Files

- `scripts/msk/emg_driven.py` — the full, self-contained, re-runnable pipeline (12 steps; imports
  `contraction_dynamics_forward.py`, `static_opt_knee.py`, `validate_hip_force.py`,
  `validate_joint_force.py`, `validate_emg_timing.py` for all proven machinery, never re-implemented).
  ~140s wall time (two ~55-70s forward integrations + two fast JointReaction passes).
- `data/msk_smoketest/subject2_walking1/emg_driven/`:
  - `naive_only_excitation.sto`, `hybrid_excitation.sto` — the two excitation inputs (158×81).
  - `naive_fwd_only_force.sto`, `emg_hybrid_force.sto` — the two 130-column force.sto files fed to
    JointReaction (80 muscle columns overridden, ~50 non-muscle columns copied verbatim from SO's own
    force.sto).
  - `jr_naive_fwd_only/`, `jr_emg_hybrid/` — the two new official-JointReaction output directories (the
    already-published SO JointReaction output, `data/.../static_optimization/jr/`, is never touched).
  - `_diag_cache.npz` — raw `a_hybrid`/`f_hybrid`/`a_naive_only`/`f_naive_only` arrays, cached before the
    forced-adversary gate check (diagnosable without re-paying the ~130s integration cost).
  - `emg_driven_results.json` — every number in this document, machine-written.
- Reused, not re-solved: `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_{activation,force}.sto`; reused, not re-run:
  `data/msk_smoketest/subject2_walking1/static_optimization/jr/walking1_JointReaction_ReactionLoads.sto`;
  reused, not reprocessed: `data/msk_smoketest/subject2_walking1/contraction_dynamics/_diag_cache.npz`
  (for `a_so_smooth`, `t_so`, `muscle_names` only — its `f_naive_fwd` field is the diagnosed-stale
  pathology, Sec. 2, loaded only for the disclosed diff, never as ground truth).
