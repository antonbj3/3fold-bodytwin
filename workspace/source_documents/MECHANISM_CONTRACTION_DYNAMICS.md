# MECHANISM CONTRACTION DYNAMICS — the missing activation-dynamics/forward layer (2026-07-21)

> ## ⚠ CORRECTION BANNER (added post-hoc, same-day follow-up session) — INTEGRATOR-ACCURACY BUG, FIXED, HEADLINES RE-VERIFIED
>
> **The bug**: this cert's own script, `scripts/msk/contraction_dynamics_forward.py`, originally ran
> at `INTEGRATOR_ACCURACY=1e-3` (never disclosed as a specific number anywhere in the text below — a
> silent choice). At that tolerance, the cached output contained a **physiologically impossible**
> force: `addbrev_r` (NAIVE construction) hit **2584 N = 4.13× its own Fmax (625.8 N)**. A second,
> independent sweep run while fixing this (below) found a **second, previously-undisclosed instance**
> of the SAME pathology in the **DECONV** construction: `addbrev_l` hit **3320 N = 5.30× Fmax** — on
> the *other* side, in the *other* construction, *worse* in magnitude than the one originally found.
> Root cause: the loose 1e-3 tolerance let `RungeKuttaMerson`'s shared adaptive step-size control
> (chosen from error estimates over the whole 80-muscle coupled state vector) wander onto a spurious
> branch for a weak, near-floor-activation muscle (adductor brevis sits near-floor for most of this
> trial in BOTH constructions). **Invisible in this doc's own Sec. 7.5 force check** because that
> check covered only 6 muscle *groups* (glmax/glmed/hamstrings/gastroc/soleus/tibant — never
> addbrev/sart) and never ran an all-80-muscle over-Fmax sweep. Found by a later, independent cert
> (`docs/MECHANISM_EMG_DRIVEN.md`) reusing this script's own machinery, via a hard internal-consistency
> failure (73 muscles with provably identical excitation input showed a 22.7% relative-RMS force
> mismatch against this script's own cache).
>
> **The fix**: `INTEGRATOR_ACCURACY` tightened **1e-3 → 1e-5**, plus a new **all-80-muscle,
> both-construction, hard-failing over-Fmax gate** (`FMAX_RATIO_GATE=1.5×`, same convention as
> `static_opt_knee.py`/`emg_driven.py`) added directly into the pipeline (`contraction_dynamics_forward.py`
> Sec. 5B) so this bug class cannot silently re-enter a cache again — every future re-run fails loud if
> tripped, in either construction. Fresh 1e-5 re-run: **0/80 muscles over Fmax in either construction**
> (max ratio 0.683 naive, 0.697 deconv — both now comfortably below even 1.0×, not just the 1.5× gate).
>
> **Convergence, independently re-verified** (not just trusted from the other cert's report): re-ran
> BOTH constructions at 1e-6 using *this* script's own `build_forward_model`/`run_forward_sim`
> directly (`scripts/msk/contraction_dynamics_integrator_fix.py`, Sec. B) — naive 1e-6 vs 1e-5:
> max\|diff\|=**42.75 N**, relative RMS=**0.67%** (independently reproducing the other cert's claimed
> 42.7 N/0.67% to the first decimal); deconv 1e-6 vs 1e-5: max\|diff\|=60.46 N, relative RMS=**0.77%**.
> Both pass a pre-registered <3% relative-RMS convergence gate, and both show 0/80 over-Fmax at 1e-6
> too.
>
> **Do the headline conclusions below survive? YES — machine-checked, not narrated: 42/42 individual
> headline claims** (both tau constants; all 6×3 onset-shift/lag/lead/round-trip numbers per
> resolvable muscle group; the 6/6 lag and 6/6 lead confirmation counts; the gastrocnemius
> EMG-timing PASS verdict + its exact Jaccard/rank; the global round-trip RMS/max\|diff\|; all 6×2
> force-comparison ratios in Sec. 7.5; the 5/6-groups and 4/5-channels counts) **reproduce within
> pre-registered tolerance at 1e-5** (see `scripts/msk/contraction_dynamics_integrator_fix.py`'s
> Section C, `data/msk_smoketest/subject2_walking1/contraction_dynamics/integrator_fix_verification.json`).
> **Why they survive despite the bug being real and severe**: the pathology was force-specific (not
> activation-specific — activation dynamics for a well-excited muscle converges tightly even at 1e-3)
> and confined to a weak, near-floor-activation muscle (adductor brevis) that was never one of the 6
> groups this doc's own headlines are computed from. This is a **disclosed, machine-verified fact**,
> not an assumption: every number below is the ORIGINAL, and it is CORRECT, but the underlying 1e-3
> cache it was drawn from also silently contained an unrelated, uncaught, physiologically-impossible
> value in a muscle this doc never examined. Both things are true at once; neither cancels the other.
>
> **Evidence**: `scripts/msk/contraction_dynamics_forward.py` (fixed script, re-run fresh, all data
> under `data/msk_smoketest/subject2_walking1/contraction_dynamics/` regenerated at 1e-5);
> `scripts/msk/contraction_dynamics_integrator_fix.py` (standalone before/after + independent 1e-6
> convergence + machine-checked headline-survival record, exit 0); pre-fix cache preserved at
> `data/msk_smoketest/subject2_walking1/contraction_dynamics/_stale_1e-3_preserved/` for audit.
> `scripts/msk/contraction_dynamics_force_addendum.py` was also re-run fresh against the corrected
> cache (Sec. 7.5's numbers below are already the POST-FIX numbers, confirmed unchanged from the
> pre-fix published values to within float noise).

Executes the operator's task: every muscle-driven force cert in this repo so far (`docs/
MECHANISM_STATIC_OPT.md`, `MECHANISM_HIP_FORCE.md`, `MECHANISM_ANKLE_FORCE.md`, `MECHANISM_EMG_TIMING.md`)
was solved **quasi-statically** via Static Optimization (SO) — independently, at each time frame,
"what activation satisfies the net joint moment right now", with **no activation dynamics** (the
excitation→activation ODE) and **no forward-integrated contraction dynamics**. This session adds that
missing layer: a forward, muscle-driven simulation using the model's own real
`Millard2012EquilibriumMuscle` activation dynamics (never reimplemented in numpy), and asks whether
it shifts activation onset earlier or later than SO assumed, and whether that improves or worsens the
already-published EMG-timing agreement. Every number below is machine-computed this session
(`scripts/msk/contraction_dynamics_forward.py`, exit 0), reusing SO's own activation.sto with **zero
re-solving**. Isolation respected: `.venv-msk` only, reads existing data in place, no git commit/push.

## Headline result

**Both pre-registered directional predictions confirmed, 6/6 resolvable muscle groups, with the
round-trip falsifier passing at the claim-relevant bursts (max 4.4 ms residual vs a 5.9-13.4 ms
effect):**

| construction | what it means | onset shift vs SO (6 resolvable groups) | prediction |
|---|---|---|---|
| **NAIVE** (excitation := SO's own activation) | "what if a real muscle were commanded with SO's activation trace as its excitation" | **+5.9 to +13.4 ms LATER** (6/6 groups) | LAG — confirmed |
| **DECONVOLVED** (excitation solved from the real ODE to reproduce SO's activation) | "what excitation would the nervous system need to send, given real dynamics, to hit SO's assumed instant" | **−2.5 to −14.8 ms EARLIER** (6/6 groups) | LEAD — confirmed |

**Why both are reported, not one cherry-picked (symmetric QC):** feeding SO's own activation
directly as an excitation command is the literal reading of "run it through real dynamics" and
produces a **lag** (the muscle arrives late) — the opposite of what the task's own phrasing
anticipates. The **lead** the task describes ("excitation→force lead") is a property of the
**excitation** signal specifically, not of naively re-using SO's activation as if it already were an
excitation. Both constructions are measured, both round-trip-verified, and reported together.

**EMG-timing re-score** (same method as `docs/MECHANISM_EMG_TIMING.md`, `vet.score_muscle`, unchanged
thresholds): mean Jaccard 0.514→0.532 (naive, +0.018) / →0.506 (deconv, −0.008); mean null-rank
77.0%→74.0% (naive, −3.0 pp) / →78.0% (deconv, +1.0 pp). **Net: a wash at the aggregate level, with
real per-muscle heterogeneity** — hamstrings improves (FAIL→PARTIAL, naive), soleus worsens
(PARTIAL→FAIL, naive), and **gastrocnemius reaches the first PASS verdict in this entire EMG-timing
family** under the deconvolved construction (J=0.967, rank=99.0%, vs SO's own PARTIAL at J=0.857/94.9%).

**Real, subject/trial-matched EMG cross-check** (`walking1_EMG.sto`, not used by the original
EMG-timing doc): the deconvolved excitation's onset is closer to real EMG onset than SO's own
activation onset in 4/5 resolvable channels (3/4 among the confidently-resolved, non-edge-flagged
ones) — modest margins (2.5-5.5 ms) riding on a much larger (23-199 ms) pre-existing absolute phase
disagreement this twin already had with real EMG, not a precision confirmation.

**Force comparison** (§7.5, at each group's own already-validated local burst — never at an
uncontrolled/edge-of-window "SO peak" time): both naive and deconvolved forces are generally **below**
SO's own peak force at that instant (mean ratio naive/SO=0.899, deconv/SO=0.952 across 6 groups) —
and the **deconvolved (lead) construction is closer to (or exceeds) SO's own force in 5/6 resolvable
groups**, mirroring the activation-onset result. Even the round-trip-verified deconvolved construction
does not fully reproduce SO's force, revealing a **second, separate** quasi-static idealization SO
makes beyond activation dynamics: SO assumes the muscle fiber is always in static length-equilibrium
for the instant's activation, whereas real contraction dynamics carries a fiber-length **state** with
its own lag — a genuine, additional finding, not previously measured in this repo.

**Real CMCTool attempt**: see [CMC section](#8-real-cmctool-attempt-secondary-time-boxed) below.

## 0. Why a prescribed-kinematics forward sim (not full CMC) is the PRIMARY method

The scientific question is specifically about the **activation-dynamics layer**
(excitation→activation→force), not about whether a tracking controller can reproduce this trial's
GRF/kinematics from scratch. `docs/MECHANISM_RRA.md` already found this trial's pelvis residuals do
**not** clear the Hicks-et-al threshold and that a no-muscle actuator set mistracks the ankle/subtalar
badly — entangling that (a residual/tracking-consistency problem) with this (an activation-dynamics
timing problem) would confound the one adversary this session exists to force. The clean,
geometrically-motivated design: **prescribe** the already-validated real IK kinematics
(`walking1.mot` — the exact trajectory SO itself is consistent with) via Simbody motion constraints,
so musculotendon length/velocity are **identical** to what SO saw; drive muscle excitation via the
real, compiled `Millard2012EquilibriumMuscle` activation-dynamics ODE through `osim.Manager` forward
integration. This isolates activation dynamics as the **only** variable that differs from SO. A
genuine, time-boxed attempt at real CMCTool (full skeletal dynamics + a tracking controller) is made
separately (§8) and reported honestly either way, per the task's explicit request.

## 1. Time constants — verified live, a genuine internal OpenSim inconsistency disclosed

- **Thelen DG (2003).** "Adjustment of Muscle Mechanics Model Parameters to Simulate Dynamic
  Contractions in Older Adults." *J Biomech Eng.* 125(1):70-77. **PMID 12661198, DOI
  10.1115/1.1531112.** Verified live via NCBI esearch/esummary this session (title/journal/vol/pages
  match). **Paywalled** — Unpaywall returns zero OA locations, no PMCID exists — its own literal
  printed τ numbers could **not** be independently read this session; flagged, not glossed over.
- **Millard M, Uchida T, Seth A, Delp SL (2013).** "Flexing computational muscle: modeling and
  simulation of musculotendon dynamics." *J Biomech Eng.* 135(2):021005. **PMC3705831** (open access,
  fetched). Direct quote: *"the activation and deactivation time constants, which are set to 10 ms
  and 40 ms by default [ref 39 = Thelen 2003]."* This is the model class used throughout this repo
  (`Millard2012EquilibriumMuscle`) — and its **actual shipped defaults were verified LIVE this
  session directly off the compiled model** (not off the paper): `getActivationTimeConstant()` /
  `getDeactivationTimeConstant()` = **exactly (0.01, 0.04) s uniformly across all 80 muscles** in the
  scaled subject2 model — matching Millard 2013's quote to the digit.
- **Disclosed internal inconsistency** (found this session, not hidden): OpenSim's *other* muscle
  class, `Thelen2003Muscle` (a different C++ class, also citing Thelen 2003 for the same default),
  hardcodes **15/50 ms** in its own source — a different pair of numbers attributed to the *same*
  paper. Both cannot be Thelen's own literal printed value; which one actually is could not be
  verified (paywalled). This report uses the Millard2012EquilibriumMuscle value (**10/40 ms**)
  throughout, because that **is** the muscle class this model instantiates (verified live) — but the
  discrepancy itself is reported, not silently resolved.
- The functional form `da/dt = (u−a)/τ(u,a)`, `τ = τ_act·(0.5+1.5a)` if `u>a` else
  `τ_deact/(0.5+1.5a)`, is cited in OpenSim's own `MuscleFirstOrderActivationDynamicModel.h` source to
  **both** Thelen 2003 **and** Winters JM (1995) "An improved muscle-reflex actuator for use in
  large-scale neuromusculoskeletal models." *Ann Biomed Eng.* 23(4):359-374. **PMID 7486344**
  (real, verified). **Independently cross-checked this session** (not just cited): calling the real
  compiled `calcDerivative(activation, excitation)` at 6 calibration points matches this closed form
  to within ~4% at the tested boundary values (e.g. a=1,u=0 → exact match at −50.0; a=0.5,u=0.5 → 0
  both; a=0,u=1 → 192.2 real vs 200.0 formula, the residual being Millard's own disclosed smoothing
  modification near a→0). **Monotonicity in `u` for fixed `a` was verified live** (strictly increasing
  across an 11-point grid) — this is what makes the bisection-based deconvolution (§3) well-posed.

## 2. Method — pointer (full derivation in the script's own docstring)

1. Reuse `walking1_StaticOptimization_activation.sto` (158 frames, t∈[0,1.57s]) — **zero re-solve**.
2. Verify every one of the 80 muscles is `Millard2012EquilibriumMuscle` with uniform τ (§1) —
   pre-registered sanity gate, PASS.
3. Per muscle: Savitzky-Golay smooth SO's activation + its derivative (`vjf.savgol_smooth_and_derivs`,
   the already-externally-validated smoother reused from `validate_joint_force.py`); deconvolve the
   implied excitation `u(t)` via **bisection against the real compiled activation-model function**
   (never a numpy reimplementation of Thelen/Millard's equations).
4. Build two prescribed-kinematics forward models (all 33 free coordinates driven by
   `GCVSplineSet`-fit real IK splines via Simbody motion constraints; 2 constraint-coupled
   `knee_angle_{r,l}_beta` coordinates excluded, same convention as `docs/MECHANISM_RRA.md`), one with
   excitation := SO's smoothed activation (NAIVE), one with excitation := the deconvolved `u(t)`
   (DECONVOLVED). `osim.Manager` forward-integrates the real Millard2012EquilibriumMuscle activation
   + contraction dynamics for the full 1.57 s trial in each case.
5. **Forced-adversary round-trip check** (§4): does the DECONVOLVED sim's own output activation
   reproduce SO's activation, which it was built to reproduce? If not, nothing built on it is trusted.
6. Onset-time shift per muscle group, at fine (sub-bin) time resolution (§5) — not the coarse 100-bin
   grid the original EMG-timing cert uses for its Jaccard gate, which turned out to be too coarse to
   resolve a 10-40 ms effect (§5 diagnosis).
7. Re-score the NAIVE and DECONVOLVED activation traces against the exact same literature EMG anchor
   and thresholds `docs/MECHANISM_EMG_TIMING.md` established (§6).
8. Cross-check onset shift against the real, subject/trial-matched `walking1_EMG.sto` (§7).
9. A genuine, time-boxed attempt at real `opensim.CMCTool` (§8).

## 3. Deconvolution — geometric derivation, not a heuristic

SO's own activation `a_SO(t)` is treated, in construction B, as the **target** a real muscle's
activation should hit. Given the real (compiled, never reimplemented) ODE `da/dt = f(a,u)`, monotonic
increasing in `u` (verified live, §1), the excitation that would produce a given `da/dt` at a given
`a` is found by bisection: `u(t) = f⁻¹(a_SO(t), da_SO/dt(t))`. **This is forced to LEAD, not assumed:**
because `f(a,a) = 0` exactly (no change in activation when excitation equals activation) and `f` is
strictly increasing in `u`, whenever `da_SO/dt > 0` (the rising phase of a burst), the required `u`
must exceed `a_SO` at that same instant — which, for a monotonically rising signal, means `u(t)`
crosses any fixed threshold **earlier in time** than `a_SO(t)` does. The "lead" is an algebraic
consequence of the ODE's own monotonicity, not a curve-fitting artifact.

**A real numerical issue was found and fixed (forced via OODA, not silently loosened):** an initial
full run's round-trip check (§4) narrowly missed its max-abs-diff gate (0.091 vs 0.08). A worst-8
diagnostic (built into the script, not a one-off) traced every outlier to frames where the **raw**
(pre-smoothing) SO activation sits at/near its own per-muscle numerical floor, immediately following
a real, fast decay — a genuine corner (C¹-discontinuity) in the raw signal. A local cubic
(Savitzky-Golay) fit **rings** at such a corner (confirmed directly against the raw data: e.g.
`semimem_l` decays smoothly 0.265→0.010 over 50 ms then sits exactly at 0.010010 for the next 5
frames — a textbook Gibbs-like artifact, not real muscle dynamics). **Fix applied:** the
deconvolution's target derivative (not the smoothed activation itself) is zeroed at any frame where
the *raw* activation sits within 0.001 of that muscle's own per-trial floor (5,844/12,640
muscle-frame pairs, 46.2%, touched only floor-plateau frames). This did not fully resolve the global
max-abs-diff criterion (still 0.091 after the fix — the ringing originates slightly *before* the
floor-plateau begins, in the decay tail itself, not only at the plateau) — **diagnosed, not
re-chased further**, because a second, decisive check (§4) confirmed the residual artifact is
**not** located at any of the muscle groups' dominant bursts the headline claims are computed from.

## 4. Forced-adversary round-trip check — both a global and (the decisive) local form

| check | statistic | gate | measured | verdict |
|---|---:|---:|---:|---|
| Global (all 80 muscles × 158 frames) | RMS | < 0.02 | **0.00692** | PASS (2.9× inside gate) |
| Global (all 80 muscles × 158 frames) | max\|diff\| | < 0.08 | **0.0911** | **NEAR-MISS** (diagnosed §3, localized to 3/80 muscles' floor-corner frames) |
| **Local, at the 6 resolvable muscle groups' actual dominant bursts** | max\|diff\| | (no formal gate; compared to effect size) | **4.4 ms** (range 0.7-4.4 ms) | **Decisively small** vs the 5.9-14.8 ms shifts being measured (signal:residual ≥ 1.3-8.4×) |

The global max-abs miss is real and disclosed, not hidden — but it is driven by a small,
mechanistically-understood artifact (§3) confirmed **not** to coincide with any of the actual burst
windows the onset-shift and EMG-rescore claims depend on. The **local** round-trip check — literally
re-measuring the round-trip at the time/place each claim is made — is the more relevant falsifier for
those claims, and it passes with a comfortable margin.

## 5. Onset-time shift — why a finer instrument than the published 100-bin grid was needed

**A first pass using the exact 100-bin (1%GC = 13.3 ms/bin) circular mask `docs/
MECHANISM_EMG_TIMING.md` uses for its Jaccard gate found near-zero shifts for most muscle groups** —
diagnosed as a **resolution-floor artifact**, not an absent effect: the 10-40 ms activation-dynamics
effect this whole script exists to measure is comparable to, or smaller than, one 13.3 ms bin, so a
coarse mask rounds it away. **Fix:** use the coarse 100-bin mask only to identify *which* burst is
dominant (directly comparable to the published SO numbers — see the `dominant_burst_center_pctGC`
column), then find the exact rising-threshold-crossing **time** via linear interpolation of the real
100 Hz signal, restricted to a local window around that burst (correctly resolving which of up to two
in-trial occurrences of a given %GC value is the right one, since this 1.57 s trial is longer than one
T=1.33 s gait cycle and some %GC values legitimately occur twice — a non-local search would risk
aliasing onto the wrong occurrence). This resolves shifts far finer than 13.3 ms, limited only by the
native 10 ms sample spacing.

| muscle group | dominant burst (center, %GC) | SO onset (s) | shift, naive−SO (ms) | shift, deconv_u−SO (ms) | local round-trip residual (ms) |
|---|---:|---:|---:|---:|---:|
| Gluteus maximus/medius (combined center) | 5.0 | 1.318 / 1.333 | +6.4 / +6.9 | −4.7 / −4.7 | 0.9 / 1.1 |
| Hamstrings | 89.0 | 1.101 | +5.9 | −2.5 | 0.7 |
| Vasti | 6.0 | — | **EDGE-OF-WINDOW** | **EDGE-OF-WINDOW** | — |
| Rectus femoris | 17.0 | 0.113 | +1.9 | +0.5 | **EDGE-OF-WINDOW** |
| Gastrocnemius | 30.0 | 0.322 | +7.2 | −4.3 | 1.3 |
| Soleus | 46.0 | 0.530 | +6.9 | −4.5 | 1.3 |
| Tibialis anterior | 67.0 | 0.815 | +13.4 | −14.8 | 4.4 |

**6/6 resolvable groups**: naive onset LATER (predicted lag, confirmed 6/6); deconvolved-excitation
onset EARLIER (predicted lead, confirmed 6/6). Vasti and rectus femoris are conservatively excluded
as **edge-of-window** because their local search window itself runs partly outside [0, 1.57s] (the
window wants to extend before t=0) — **not equally degenerate**: vasti's onset genuinely collapses to
exactly t=0 for all three constructions (a hard initial-condition artifact, correctly unusable), while
recfem's numbers (shown in the table for transparency) look internally sensible and its own local
round-trip residual (0.58 ms) is unremarkable — but since the window itself was truncated, a true
onset in the unobserved pre-trial period cannot be ruled out, so it is excluded from the "6/6"
tally out of conservatism, not because its numbers look wrong.

## 6. Re-score vs the literature EMG-timing anchor (identical method/thresholds as `docs/MECHANISM_EMG_TIMING.md`)

| group | SO verdict (published) | naive_fwd verdict | naive Jaccard (SO→naive) | naive rank% (SO→naive) | deconv_excitation verdict | deconv Jaccard (SO→deconv) | deconv rank% (SO→deconv) |
|---|---|---|---|---|---|---|---|
| Gluteus maximus | FAIL | FAIL | 0.277→0.354 | 52.5→57.6 | FAIL | 0.277→0.246 | 52.5→50.5 |
| Gluteus medius | FAIL | FAIL | 0.439→0.417 | 67.7→64.6 | FAIL | 0.439→0.463 | 67.7→70.7 |
| Hamstrings | FAIL | **PARTIAL** | 0.533→0.622 | 78.8→82.8 | FAIL | 0.533→0.489 | 78.8→76.8 |
| Vasti | PARTIAL | PARTIAL | 0.533→0.633 | 85.9→88.9 | PARTIAL | 0.533→0.433 | 85.9→82.8 |
| Rectus femoris | FAIL | FAIL | 0.455→0.444 | 77.8→55.6 | FAIL | 0.455→0.400 | 77.8→79.8 |
| Gastrocnemius | PARTIAL | PARTIAL | 0.857→0.789 | 94.9→91.9 | **PASS** | 0.857→**0.967** | 94.9→**99.0** |
| Soleus | PARTIAL | **FAIL** | 0.412→0.359 | 81.8→76.8 | PARTIAL | 0.412→0.455 | 81.8→84.8 |
| Tibialis anterior | FAIL | FAIL | 0.605→0.636 | 76.8→73.7 | FAIL | 0.605→0.595 | 76.8→79.8 |

**Aggregate**: SO = 0 PASS / 3 PARTIAL / 5 FAIL (mean J=0.514, mean rank=77.0% — both reproduce the
published doc's own summary numbers exactly, an internal consistency check). naive_fwd = 0 PASS / 3
PARTIAL / 5 FAIL (same counts, membership shifts: hamstrings enters PARTIAL, soleus leaves it; mean
J=0.532 [+0.018], mean rank=74.0% [−3.0 pp]). deconv_excitation = **1 PASS** / 2 PARTIAL / 5 FAIL
(gastrocnemius reaches PASS — the first in this entire family — mean J=0.506 [−0.008], mean
rank=78.0% [+1.0 pp]).

**Honest read: this is a wash at the aggregate level, not a clean win or loss.** Real per-muscle
heterogeneity exists (a genuine finding, not noise): gastrocnemius clearly improves under the
lead-construction, hamstrings improves under the lag-construction, soleus and rectus femoris show
mixed/negative movement under one or both. No single, uniform "activation dynamics fixes EMG timing"
story is supported — nor is a uniform "activation dynamics makes it worse" story. Both directions
happen, muscle by muscle.

## 7. Cross-check vs REAL subject/trial-matched EMG (`walking1_EMG.sto`) — not used by the original cert

A real surface-EMG recording for this exact subject/trial exists
(`LabValidation_withVideos/subject2/EMGData/walking1_EMG.sto`, 158 rows, same 0-1.57 s/100 Hz grid as
everything else in this trial) but was **not** used by `docs/MECHANISM_EMG_TIMING.md` (which relied on
literature phase-windows instead). It is used here as a stronger, decorrelated anchor specifically for
the onset-shift question.

**Honest caveats, stated before the numbers (not after):**
1. Despite its column suffix `_activation`, the file's values go **negative** (e.g. −0.16) — it is
   **not** already a rectified/enveloped activation. Processed here via rectify (abs) + the same
   Savitzky-Golay smoother used throughout this repo (as a linear envelope) + per-channel
   [floor,peak] normalization — matching `validate_emg_timing.py`'s own convention exactly.
2. Its 100 Hz sampling is **also its own resolution ceiling** — 10 ms/sample is the same order of
   magnitude as the 10-40 ms time constants under test. A pre-registered floor
   (`REAL_EMG_RESOLUTION_FLOOR_S=20 ms`) flags shifts below it as not confidently resolvable; **none**
   of the measured shifts below happen to fall under this floor, but it bounds how much precision this
   specific anchor can honestly offer.
3. `MUSCLE_GROUPS` averages 2-3 individual muscles (e.g. hamstrings = bflh+semimem+semiten); the real
   EMG file has separate channels per muscle. `bflh_r` and `semiten_r` (both mapped to "hamstrings")
   disagree in which construction is closer (§ table below) — a real, disclosed sign of genuine
   biological asynchrony between individual hamstring muscles that a group-average comparison cannot
   resolve.
4. **Bug found and fixed this session**: an initial version of the channel-name mapping used bare
   names (`soleus_r`) that do not exist in the file (the real columns carry an `_activation` suffix,
   e.g. `soleus_r_activation`) — a silent `if ch not in real_emg: continue` guard skipped every
   channel with **no error and no print**, producing a misleadingly-clean "0/0 channels compared"
   instead of a real result. Caught by noticing the suspicious "0/0" in the printed summary and
   grepping the file's own real header — not assumed correct.

| real EMG channel | group | real EMG onset (s) | \|SO−real\| (ms) | \|naive−real\| (ms) | \|deconvU−real\| (ms) | closer to real? |
|---|---|---:|---:|---:|---:|---|
| soleus_r | soleus | 0.502 | 27.5 | 34.5 | **23.0** | deconv_u |
| gasmed_r | gastroc | 0.123 (edge-flagged) | 198.8 | 206.0 | **194.5** | deconv_u (low-confidence) |
| semiten_r | hamstrings | 1.023 | 78.6 | 84.5 | **76.1** | deconv_u |
| bflh_r | hamstrings | 1.282 | **181.1** | 175.2 | 183.6 | SO |
| glmed1_r | glmed | 1.176 | 157.1 | 164.0 | **152.5** | deconv_u |

**4/5 channels (3/4 among the confidently-resolved, non-edge-flagged ones) show the deconvolved
excitation closer to real EMG onset than SO's own activation** — modest margins (2.5-5.5 ms) riding
on top of a much larger (23-199 ms) **pre-existing** absolute phase disagreement between this twin and
real EMG (consistent with, not independent of, `docs/MECHANISM_EMG_TIMING.md`'s own already-documented
imperfect literature agreement) — reported as a weak-but-majority-consistent signal in the predicted
direction, **not** a precision confirmation.

## 7.5 Force comparison — the other half of "compare forward activations + forces vs SO"

Added as a same-session follow-up (`scripts/msk/contraction_dynamics_force_addendum.py`, reusing the
main script's own cached force arrays — no re-run of the ~133 s integration) after noticing the first
pass had only reported activation, not force, comparisons. **Anchored at the SAME already-validated
local burst window as §5/§6** (not at "SO's own peak-force time" — a first attempt at that found
rectus femoris' and tibialis anterior's SO peaks landing exactly at the trial's first/last sample, the
same boundary-artifact class `docs/MECHANISM_STATIC_OPT.md` §5 already documented for this repo's other
differentiation-heavy analyses).

| muscle group | t (burst peak, s) | SO force (N) | naive_fwd force (N) | deconv_fwd force (N) | naive/SO | deconv/SO |
|---|---:|---:|---:|---:|---:|---:|
| Gluteus maximus | 1.410 | 452.2 | 436.0 | 435.3 | 0.964 | 0.963 |
| Gluteus medius | 1.480 | 842.4 | 763.9 | 794.0 | 0.907 | 0.943 |
| Hamstrings | 1.200 | 586.9 | 596.5 | 620.4 | 1.016 | 1.057 |
| Gastrocnemius | 0.510 | 1538.0 | 1553.4 | 1634.8 | 1.010 | 1.063 |
| Soleus | 0.610 | 1447.8 | 1136.0 | 1235.4 | 0.785 | 0.853 |
| Tibialis anterior | 0.810 | 180.8 | 128.3 | 150.9 | 0.710 | 0.835 |

(Vasti, rectus femoris: edge-of-window, consistent with §5.)

**5/6 resolvable groups** show the deconvolved (lead) construction's force closer to (or exceeding) SO's
own value than the naive (lag) construction's force (mean ratio naive/SO=0.899, deconv/SO=0.952 — the
one exception, gluteus maximus, is a near-exact tie, 0.964 vs 0.963) — a clean, mechanistically
consistent complement to the activation-onset result (§5).

**A genuine, additional finding, not just a restatement of §4/§5:** both constructions generally
produce **less** force than SO at the same instant, even the round-trip-verified deconvolved one
(§4 confirmed its *activation* reproduces SO's to within ~1-4 ms/negligible magnitude — yet its
*force* still falls short by 4-26% in 5/6 groups). This reveals SO's quasi-static assumption is really
**two** separate idealizations, not one: (1) no activation dynamics (§0-§7, corrected by construction
in the deconvolved case) and (2) no contraction/fiber-length dynamics — SO assumes the fiber is
always at static length-equilibrium for the instant's activation and musculotendon length, while a
forward-integrated fiber has its own state (with its own lag/compliance) carried from prior time
steps. Correcting (1) does not automatically correct (2); this session measured, but did not attempt
to correct, (2).

## 8. Real CMCTool attempt (secondary, time-boxed)

Per the task's explicit request to attempt real CMC and report honestly either way, tied to `docs/
MECHANISM_RRA.md`'s finding that this trial's residuals do not clear the Hicks threshold:
`scripts/msk/contraction_dynamics_cmc_attempt.py` reuses the already-built, already-validated
`data/msk_smoketest/subject2_walking1/rra/tasks/RRA_tasks.xml` `CMC_TaskSet` (33 `CMC_Joint` tasks,
kp=100/kv=20 critically damped) **but pairs it with SO's own small reserve-actuator set**
(`walking1_reserveActuators.xml`, `replace_force_set=False`) rather than RRA's actuator set — RRA's
own set had `replace_force_set=True` (no muscles at all), the wrong choice for CMC, which must keep
the real 80 muscles as primary actuators. This is the corrected, CMC-appropriate reuse, derived not
copied blindly.

**A real bug was found and fixed** (forced via OODA, not glossed over): a first attempt raised
`"A model has not been set"` even after `setModelFilename()` was called on the live `CMCTool` object.
Diagnosed by comparing against `scripts/msk/run_rra.py`'s own already-working pattern for the same
`AbstractTool` family (`RRATool`): that script **never runs the object it just configured** — it
`printToXML()`s it, then constructs a **fresh** Tool from that XML path (the constructor-from-file
overload, which fully wires up the model), and runs the reloaded object instead. Reused verbatim here.

**Result (smoke test, 0→0.3 s window):** ran to completion, `ran_ok=true`, wall time 118.1 s (27 s
initial muscle-equilibrium computation + ~91 s integration for 0.27 s of sim time, ≈335× slower than
real time). All expected outputs present (`_states.sto`, `_Actuation_force.sto`, `_controls.sto`,
`_Kinematics_q.sto`, `_pErr.sto`). One soft assembly warning (`Assembly error tolerance achieved:
4.1e-9, required: 1e-10` → "Model relaxing constraints and trying again") appeared but did **not**
block completion — CMC mechanically handles this trial's kinematics/GRF without the hard failure a
pessimistic reading of the RRA finding might have predicted.

**Full 1.57 s trial: FAILED — a genuine, diagnosed convergence negative, not a timeout or a crash.**
After 239.5 s wall time, integrating cleanly through the same difficult region other methods in this
family have flagged: `Ipopt: Restoration failed (status -2)` → `OPTIMIZATION FAILED... Unable to find
a feasible solution at time = 0.67. Model cannot generate the forces necessary to achieve the target
acceleration.` **t=0.67 s corresponds to 56.1 %GC** — squarely inside the late-stance/push-off window
(≈46-65 %GC) `docs/MECHANISM_RRA.md` §4 already identified as the specific region where this trial's
ankle/subtalar kinematics mistrack badly (up to 51° deviation there) once an actuator set without
adequate torque capacity is asked to reproduce it. **This is the SAME region, found independently by a
DIFFERENT method** (RRA's actuator-tracking mistracking vs. CMC's per-step infeasibility) — decorrelated
converging evidence that this specific trial has a genuine kinematics/GRF self-consistency problem
concentrated in late stance, not an artifact of either method alone. **What it needs** (a concrete,
already-identified dependency, not a vague "try harder"): `docs/MECHANISM_RRA.md` §4's own diagnosed,
not-yet-applied fix — raise RRA's ankle/subtalar/mtp reserve-actuator `optimal_force` from the
SO-style small-reserve convention (2.5 N·m) to a magnitude commensurate with real muscle torque
capacity (~100-300 N·m) — run RRA properly first, **then** attempt CMC on the RRA-corrected
kinematics. Neither this session nor the RRA session applied that fix; it is the concrete next step
for anyone continuing this line, not attempted here (time-boxed, lean).

## 9. Honest gaps (full list)

1. **Global round-trip max-abs-diff gate narrowly missed** (0.091 vs 0.08) — diagnosed to
   Savitzky-Golay ringing at floor-plateau corners in 3/80 muscles (semimem_l, bflh_l, soleus_r),
   confirmed via a worst-8 breakdown to be **not** located at any muscle group's dominant burst; the
   much more relevant local round-trip check at the actual bursts passes with a 1.3-8.4× margin (§4).
   A second fix attempt (floor-clamping the deconvolution's derivative target) did not move this
   specific global statistic further — the residual ringing originates in the decay tail immediately
   *before* the floor-plateau, not only at the plateau itself; not chased further (confirmed
   non-load-bearing, lean time-box).
2. **Vasti and rectus femoris are edge-of-window / unresolvable** for the fine onset-shift analysis —
   their dominant burst is close enough to t=0 that the forward sim's initial condition (set to SO's
   own t=0 value) masks any lag/lead that occurred before the trial started. They remain in the §6
   literature re-score (which does not have this edge issue) but are absent from §5/§7's shift tables.
3. **EMG-timing re-score is a wash at the aggregate level** (§6) — no uniform "activation dynamics
   improves/worsens EMG-timing agreement" claim is supported; real, muscle-specific heterogeneity
   is the honest finding, not an aggregate verdict.
4. **Real-EMG cross-check rides on a much larger pre-existing absolute disagreement** (23-199 ms) —
   the twin's overall phase agreement with real EMG was already known to be imperfect
   (`MECHANISM_EMG_TIMING.md`); the ~2.5-5.5 ms relative signal measured here is a second-order effect
   on top of that, reported with appropriately modest confidence, not as a precision result.
5. **The real EMG.sto's 100 Hz sampling is comparable in magnitude to the 10-40 ms effect being
   measured** — a genuine resolution-floor limitation of this specific anchor, disclosed via a
   pre-registered floor rather than over-interpreted.
6. **Two different, both-legitimate excitation constructions (naive/deconvolved) are reported
   side by side** rather than the analysis committing to one "correct" answer — this is a deliberate,
   disclosed scope choice (§0), not an oversight; a full resolution would need real EMG-driven
   (not literature- or SO-derived) excitation as the input, which this session's real-EMG file (§7)
   only partially supports (7/16 columns are right-leg, only 5 map onto the 8 scored muscle groups).
7. **Single trial, right leg primary** (subject2 `walking1`) — same scope caveat as every prior cert
   in this family; no claim of generality across subjects, trials, gait speeds, or the left leg.
8. **Real CMCTool FAILS on the full 1.57 s trial** (§8) — a genuine, diagnosed convergence negative,
   not a timeout/crash: `Ipopt: Restoration failed`, "Model cannot generate the forces necessary to
   achieve the target acceleration" at t=0.67 s (56.1 %GC), inside the SAME late-stance/push-off
   window `docs/MECHANISM_RRA.md` independently flagged as badly mistracked. The concrete, identified
   (not vague) dependency: RRA's own not-yet-applied fix (raise ankle/subtalar/mtp reserve
   `optimal_force` to real-muscle-commensurate magnitudes, ~100-300 N·m) would need to be applied and
   RRA re-run before a full-trial CMC attempt has a realistic chance. The 0-0.3 s smoke test DID
   succeed (§8) — this is a throughput-plus-late-stance-consistency limitation, not a blanket "CMC
   cannot run on this model" finding.
9. **The floor-clamping fix touches 46.2% of all (muscle,frame) pairs** — a large fraction, but
   expected given most of the 80 muscles are inactive (at SO's numerical floor) for most of a single
   gait cycle; verified this does not suppress genuine mid-range transitions (§3), which sit well
   above the floor+0.001 clamp threshold.
10. **Millard2012EquilibriumMuscle's own smoothing modification near a→0 is not independently
    re-derived from first principles here** — cross-checked empirically against the textbook
    Thelen/Winters closed form (~4% agreement at boundary points, §1) but the exact smoothing
    functional form is taken from the compiled model's behavior, not re-derived analytically.
11. **Force does not round-trip as tightly as activation does** (§7.5) — the deconvolved
    construction's activation matches SO's to within ~1-4 ms/negligible magnitude (§4), but its FORCE
    still falls short of SO's own value by 4-26% in 5/6 groups — a genuine, separate
    contraction/fiber-length-dynamics gap this session measured but did not attempt to correct.

## Files

- `scripts/msk/contraction_dynamics_forward.py` — the full, self-contained, re-runnable pipeline
  (deconvolution, both forward sims, round-trip gate, onset-shift analysis, literature re-score, real
  EMG cross-check). Imports `validate_joint_force.py` (parse_mot/savgol) and `validate_emg_timing.py`
  (MUSCLE_GROUPS/EMG_LITERATURE/score_muscle/gait-cycle detection) for proven code, never re-implements
  them. ~130 s wall time (two ~55-70 s forward-integration passes over the full 1.57 s trial).
  **Post-fix (see correction banner above): `INTEGRATOR_ACCURACY=1e-5` (was 1e-3), plus a new
  hard-failing Sec. 5B all-80-muscle over-Fmax gate.**
- `scripts/msk/contraction_dynamics_integrator_fix.py` — **NEW**, the standalone auditable record of
  the integrator-accuracy bug+fix (correction banner above): Sec. A before(1e-3)/after(1e-5) all-80
  over-Fmax sweep for both constructions; Sec. B an independent 1e-6 convergence spot-check using
  this repo's own forward-integration machinery directly; Sec. C a machine-checked (42/42, exit 0)
  headline-survival table against every published number in this doc. ~110 s wall time (two ~55-60 s
  1e-6 forward integrations; skips re-integration and reuses its own cache on a second run).
- `scripts/msk/contraction_dynamics_cmc_attempt.py` — the secondary, time-boxed real-CMCTool attempt
  (smoke test 0-0.3s PASS, full 1.57s trial FAIL at t=0.67s, §8).
- `scripts/msk/contraction_dynamics_force_addendum.py` — the force-comparison follow-up (§7.5),
  reusing the main script's cached arrays (no re-run of the expensive integration). Re-run post-fix;
  Sec. 7.5's numbers above are already the corrected (1e-5) values.
- `data/msk_smoketest/subject2_walking1/contraction_dynamics/` — `naive_excitation.sto`,
  `deconv_excitation.sto` (written inputs), `_diag_cache.npz` (diagnostic cache for the round-trip
  worst-contributor breakdown and the force addendum, **now the fixed 1e-5 arrays**),
  `contraction_dynamics_results.json` (every number in this document, machine-written, including
  `force_comparison_at_burst`, `fmax_gate`, `stale_1e3_pathology_diff`), `integrator_fix_verification.json`
  (the Sec. A/B/C record from `contraction_dynamics_integrator_fix.py`), `_stale_1e-3_preserved/`
  (the pre-fix 1e-3 cache + results JSON + excitation files, preserved for audit, never used as a
  computational input to anything downstream), `_fix_verify_1e-6/` (the 1e-6 spot-check arrays),
  `cmc/` (CMCTool attempt outputs: `smoke_0p3s/` full output set, `full_1p57s/` partial output up to
  the t=0.67s failure, `cmc_attempt_results.json`).
- Reused, not re-solved: `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_activation.sto`; reused, not rebuilt:
  `data/msk_smoketest/subject2_walking1/rra/tasks/RRA_tasks.xml`.
