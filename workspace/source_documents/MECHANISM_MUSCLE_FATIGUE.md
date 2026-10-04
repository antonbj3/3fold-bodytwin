# MECHANISM MUSCLE FATIGUE — the history-state layer Hill/Millard cannot have (2026-07-21)

Adds **muscle fatigue** to the twin: a validated motor-unit fatigue model coupled to the
already-computed Static Optimization (SO) muscle activations for subject2/walking1 — **no
re-solve**. Script: `scripts/msk/muscle_fatigue.py`. Evidence:
`data/msk_smoketest/subject2_walking1/muscle_fatigue/muscle_fatigue_results.json`.

## 0. Why this, why now

`docs/MECHANISM_CROSSBRIDGE_MODEL.md` (this repo, same day) verified **live**, on the real
`Millard2012EquilibriumMuscle` curves, that Hill-type muscle force is a **pure algebraic function**
of instantaneous (length, velocity, activation): *"RFE_Hill = 0.0000000000%... exactly zero to
floating-point precision, and (mathematically) the relaxation time constant is exactly zero too,
since there is no ODE state to relax... a structural property of the equation's form, regardless of
parameter values."* Every muscle in this twin (168 instances, verified there) is exactly this model.
Fatigue is the same class of gap one level up: `Muscle.getActivation()` has no memory of how long or
how hard it has already worked. This script adds that missing **history state**.

## 1. Model — Xia & Frey Law (2008) three-compartment motor-unit model

The population of motor units in one muscle is partitioned into three pools — **Resting (MR)**,
**Active (MA)**, **Fatigued (MF)** — with `MR + MA + MF = 1` at all times (a 2-simplex conservation
law, not a soft constraint: the RK4 integration holds it to **4e-14** absolute error over the
longest run, see §4). First-order flow between pools:

```
dMR/dt = -C(t) + R·MF
dMA/dt =  C(t) - F·MA
dMF/dt =  F·MA - R·MF
```

driven by a controller `C(t)` that tries to track a target load `TL(t)` (target activation, %max):

```
C(t) = L·(TL-MA),   if MA < TL and MR > (TL-MA)     [enough resting pool -- recruit toward TL]
C(t) = L·MR,          if MA < TL and MR <= (TL-MA)     [NOT enough resting pool -- recruit everything left]
C(t) = L·(TL-MA),   if MA >= TL                        [de-recruit back toward TL]
```

`F` (fatigue rate) and `R` (recovery rate) are joint-specific rate constants (1/s); `L=10` is an
"arbitrary constant tracking factor" (both papers' own phrase). Cases 1 and 3 are algebraically
identical (both `L·(TL-MA)`) — only case 2 (insufficient resting pool) differs, so the whole
controller reduces to one `np.where`. This governing system was **not** reconstructed from memory:
it is reprinted verbatim in ref [2]'s open-access Appendix (fetched live, quoted above verbatim),
and this project's own prior finding is a measured **~62% drift rate** on recalled citation details
— everything below was fetched live this session, not recalled.

**Intermittent-task correction (ref [4]):** gait is cyclic (activation troughs), not a sustained
isometric hold — the regime the base model above was built/validated for. Ref [4]'s own finding:
applying the base model (no rest-boost) to intermittent tasks **"over-predicted fatigue by 23.3%
torque decline on average."** Their fix: during rest (`TL≈0`), multiply `R` by a generic factor
`r=15` (ankle/knee/elbow optimum). This script applies exactly that correction for the
sustained/repeated-gait protocol (§5), using `TL < 0.05` as its own (disclosed, not
paper-specified) operational "at rest" gate, and reports the `r=1` vs `r=15` sensitivity explicitly
rather than silently using the out-of-regime base model.

## 2. Citations — every PMID/DOI verified LIVE (NCBI eutils + PMC full text), not recalled

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Xia T, Frey Law LA (2008). "A theoretical approach for modeling peripheral muscle fatigue and recovery." *J Biomech* 41(14):3046-52. | **18789445**, `10.1016/j.jbiomech.2008.07.013` | **Primary model** (3-state ODE + controller). No PMC copy (verified via esummary) — paywalled; equations used here are ref [2]'s verbatim reprint. |
| 2 | Frey-Law LA, Looft JM, Heitsman J (2012). "A three-compartment muscle fatigue model accurately predicts joint-specific maximum endurance times for sustained isometric tasks." *J Biomech* 45(10):1803-8. | **22579269**, `10.1016/j.jbiomech.2012.04.018`, PMC3397684 (OA, full text fetched) | Verbatim Appendix equations; joint-specific F/R (Table 1); the model's own validation-error magnitude (Table 2). |
| 3 | Frey Law LA, Avin KG (2010). "Endurance time is joint-specific: a modelling and meta-analysis investigation." *Ergonomics* 53(1):109-29. | **20069487**, `10.1080/00140130903389068`, PMC2891087 (OA, full text fetched) | **External anchor**: directly-fit power-law `ET = b0·intensity^b1` per joint (Table 2). |
| 4 | Looft JM, Herkert N, Frey-Law L (2018). "Modification of a three-compartment muscle fatigue model to predict peak torque decline during intermittent tasks." *J Biomech* 77:16-25. | **29960732**, `10.1016/j.jbiomech.2018.06.005`, PMC6092960 (OA, full text fetched) | Rest-recovery multiplier `r` for intermittent/cyclic tasks (gait-relevant); quantifies base-model over-prediction. |
| 5 | Liu JZ, Brown RW, Yue GH (2002). "A dynamical model of muscle activation, fatigue, and recovery." *Biophys J* 82(5):2344-59. | **11964225**, `10.1016/S0006-3495(02)75580-X`, PMC1302027 (OA) | Task's named **alternative** model. Verified, cited, **not implemented** — a real citation, not a fabricated fallback. |

Fetch trail: NCBI eutils `esearch`→`esummary` (`version=2.0`, for DOI/PMC ids) resolved all 6 PMIDs
(a 3-candidate `esearch` initially returned 18789445/22579269/29960732 together — disambiguated by
title in the same `esummary` call, not assumed from ranking). PMC full text fetched directly (not
just abstracts) for refs 2, 3, 4 to pull verbatim equations/tables.

## 3. Coupling — reuses the existing SO output, no re-solve

Reads `data/msk_smoketest/subject2_walking1/static_optimization/so/walking1_StaticOptimization_activation.sto`
(158 frames, t=0–1.57s, convergence-gated PASS per `static_opt_knee_results.json`) and treats each
muscle's own SO activation as `TL(t)` — the literal reading of "target task intensity" in refs
[1][2]. `opensim.AnalyzeTool`/`StaticOptimization` is **never called**. The only OpenSim use is a
**static parameter read** — `Muscle.getMaxIsometricForce()` on the subject's already-scaled model
(`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`) — to convert dimensionless `MA(t)` into
Newtons; no dynamic simulation runs on the model.

**40 unique muscle stems → joint-parameter group** (disclosed simplification; ref [2] Table 1 has no
dedicated hip row): Ankle (soleus, gastrocs, tib ant/post, etc., 11 muscles) · Knee (vasti,
hamstrings, rectus femoris, sartorius, gracilis, 10 muscles — biarticular hip-knee muscles classed
by knee action) · General (glutes, iliopsoas, adductors, TFL, piriformis — 19 muscles, the fallback
for anything primarily hip/trunk-crossing). Filtering muscle columns by this mapping automatically
excludes the SO file's other 31 non-muscle columns (lumbar/shoulder/elbow ideal torque actuators,
pelvis residuals, joint reserves) — verified: exactly 80 real `Muscle` columns recovered, asserted
in code.

## 4. Implementation self-check (before trusting any downstream number)

| check | method | result |
|---|---|---|
| dt convergence | Tlim at dt=0.05s vs dt=0.01s, General/TL=0.5 | 120.60s vs 120.58s (**0.017%** diff, gate <1%) |
| Reduced-form cross-check | independent closed-form `Tlim = -(1/R)·ln(1 - R(1-TL)/(F·TL))`, derived from a **singular-perturbation** argument (L=10 ≫ F,R~1e-3, so recruitment is a fast servo nested inside slow fatigue/recovery — see docstring for the full derivation) vs. the RK4 numeric integration | ratio **1.114** (gate [0.4, 2.5]) |
| Conservation | max\|MR+MA+MF−1\| over the longest (60-min, 362k-step) run | **4.3e-14** (gate <1e-6) |

The reduced-form check is a genuinely independent method (closed-form algebra vs. numerical ODE
integration) — not just re-running the same code twice.

## 5. Protocol C — validation against the external anchor (the real test)

Constant-`TL` "sustained isometric" runs (the regime refs [1][2][3] themselves validate against),
swept over **all 7** joint groups × 5 intensities (20/35/50/65/80% MVC) — deliberately including the
3 groups outside this trial's own muscle set (Trunk/Shoulder/Grip) as a diverse-instance-space test
of the model's own fidelity, decorrelated from this script's muscle-mapping choices. Compared
against ref [3]'s **directly-fit** power-law endurance-time curve (a genuinely different functional
form — empirical regression vs. mechanistic-ODE reduction — though both papers share an author and
likely overlapping meta-analysis source data, disclosed, not oversold as fully independent).

| joint | TL=0.20 | TL=0.35 | TL=0.50 | TL=0.65 | TL=0.80 |
|---|---|---|---|---|---|
| General (R²=0.81) | 535.8s / 530.6s | 224.1s / 175.2s | 120.6s / 86.5s | 68.6s / 51.4s | 37.4s / 34.1s |
| Ankle (R²=0.88) | 894.9s / 955.7s | 370.9s / 301.8s | 199.1s / 144.7s | 113.1s / 84.3s | 61.6s / 55.0s |
| Knee (R²=0.79) | 352.8s / 399.4s | 146.0s / 139.5s | 78.4s / 71.3s | 44.6s / 43.6s | 24.3s / 29.5s |

(format: **Tlim_model / Tlim_anchor**, seconds; full 7×5 = 35-point table in the JSON evidence)

**Result: ratio range 0.77–1.57 across all 35 (joint, intensity) pairs, monotonically decreasing in
intensity for both curves at every joint (verified, not eyeballed). Gate (pre-registered [0.2, 5.0]):
PASS.** Context that calibrates how tight this gate is fair to set: ref [2]'s **own** model-vs-351-
real-data-point validation (same paper, same fit) had RMS errors of 2.7–28.2s (roughly −20%..+16%
relative, their Table 2) — this script's cross-paper spread (0.77–1.57×) is the **same order of
magnitude**, not an artificially tighter (or looser) bar than the source literature's own standard.

## 6. Protocol A — single real gait cycle (158 frames, 1.57s) — pre-registered TINY

Pre-registered (before integrating): `ΔMF ≈ F·mean(activation)·duration`. Measured: **max ΔMF =
0.34%** (`glmed1_r`, the trial's highest-mean-activation muscle at 25.1%), matching the
linearization to **15.3%** relative error (gate <30%, itself an appropriately loose bound for a
linearization). **Gate (MF < 2% for all 80 muscles): PASS — single-cycle fatigue is genuinely tiny,
exactly as pre-registered, not narrated.**

**A forced adversary caught a real confound before it became a false claim.** A naive
last-frame/peak `MA/TL` ratio initially showed up to **66% "force decline"** within one cycle —
alarming, and *wrong* to attribute to fatigue. OODA: **Observe** the surprising number → **Orient**:
this trial's own activation keeps *rising* for the first 0.57s (peak at frame 57/158) while the
simulation starts from a cold `MA(0)=0` — a recruitment-servo (finite bandwidth, τ≈1/L=0.1s) chasing
a still-rising target has nothing to do with fatigue → **Decide**: run a **zero-fatigue control**
(`F=0` for all 80 muscles, identical recruitment dynamics, same real trace) → **Act**: measure.
Result: the **F=0 control reproduces 99.8%** of the summed real-run "decline" — proving it is a
cold-start/bandwidth artifact. The genuinely fatigue-attributable component (real minus control) is
**max 0.12%** (`sart_r`) — tiny, confirming the pre-registered claim for the right reason instead of
an accidental one. **Gate (fatigue-attributable decline < 0.5%): PASS.**

## 7. Protocol B — sustained/repeated gait (5–60 min continuous walking)

The same measured 158-frame stride tiled to cover up to 60 minutes (2293 repeats), at the native
0.01s sampling resolution — this is where fatigue becomes genuinely visible, exactly because
single-cycle fatigue (§6) is too small to be interesting on its own.

| time | MF mean (r=15) | MF max (r=15) | MF max (r=1, no rest-boost) |
|---|---|---|---|
| 1 cycle | 0.08% | 0.34% | 0.34% |
| 5 min | 7.0% | 38.6% | 62.4% |
| 10 min | 8.2% | 47.6% | 84.7% |
| 30 min | 8.7% | 49.9% | 88.5% |
| 60 min | 8.8% | **50.0%** | 89.9% |

**r=1 vs r=15 sensitivity (forced adversary variant, matching ref [4]'s own finding):** at 60 min,
`MF_r1/MF_r15 = 1.80×` — the base (sustained-task-calibrated) model over-predicts fatigue on this
cyclic gait input relative to the intermittent-corrected variant, **directionally consistent** with
ref [4]'s own measured 23% over-prediction on intermittent tasks (gate: PASS).

**A second, deconfounded reporting decision, learned directly from §6's finding:** force decline is
reported as **peak-per-cycle capacity vs. a "settled" baseline** (mean peak-MA of cycles 1–4, which
have already passed the cold-start transient but have accrued <1% fatigue) — **not** an
instantaneous single-frame ratio, for the identical reason §6 forced.

**Geometric insight (derived from the MR+MA+MF=1 simplex, not asserted, and machine-verified):**
visible peak-force decline is **dissociable** from accumulated fatigue (MF) — it only appears once a
muscle's *peak* demand collides with its own MF via the conservation ceiling `peak_TL + MF → 1`.
Concretely: `glmed1_r` (peak activation 0.582, the highest in this muscle set) is this run's worst
visible-decline muscle (**7.96%** at 60 min) with `peak_TL+MF = 0.582+0.500 = 1.082` — over the
ceiling. By contrast, `tibant_r` accumulates **more** hidden fatigue debt (MF=0.479, second-highest
of the 80) yet shows **0% visible peak-decline**, because its peak demand is modest
(`peak_TL+MF = 0.477+0.479 = 0.956`, comfortably under 1.0). A muscle can be substantially fatigued
"underneath" (MR depleting) while still hitting every peak it's asked for — until, nonlinearly, it
can't. This is a genuine structural property of the model (a direct consequence of its own
conservation law), not a tuning artifact — and it argues for reporting **both** MF (the complete,
monotonic fatigue-state trajectory) and peak-decline (the functionally-relevant, threshold-like
consequence) rather than either alone.

## 8. Force decline in Newtons — 8 representative major leg muscles (60-min sustained, r=15)

Real, subject-scaled `Fmax` (from `Muscle.getMaxIsometricForce()`), baseline vs. final (60 min)
peak-capacity:

| muscle | group | Fmax (N) | MF @60min | baseline peak capacity (N) | final peak capacity (N) | peak deficit |
|---|---|---:|---:|---:|---:|---:|
| soleus_r | Ankle | 6194.8 | 8.6% | 1908.2 | 1908.2 | 0.0% |
| gasmed_r | Ankle | 3115.5 | 9.8% | 1227.8 | 1227.8 | 0.0% |
| tibant_r | Ankle | 1227.5 | **47.9%** | 302.8 | 302.8 | 0.0% |
| vaslat_r | Knee | 5148.8 | 2.9% | 595.6 | 595.6 | 0.0% |
| bflh_r | Knee | 1313.2 | 1.8% | 96.8 | 96.8 | 0.0% |
| glmax2_r | General | 1406.0 | 6.6% | 311.9 | 311.9 | 0.0% |
| psoas_r | General | 1426.8 | 16.1% | 514.1 | 514.1 | 0.0% |
| iliacus_r | General | 1021.1 | 25.3% | 528.0 | 528.0 | 0.0% |

**Honest reading (§7's insight applied directly):** none of these 8 hand-picked "major" muscles
happens to cross the peak-demand ceiling within 60 min at these generic parameters — `tibant_r`
accumulates the most hidden fatigue (47.9% MF, the highest of these 8, because it stays active with
little true rest through gait: only 13% of the cycle is below the 5%-of-max "rest" threshold) but
its peak demand (47.7%) stays under its own MF ceiling. The muscle that **does** cross it
(`glmed1_r`, §7) simply wasn't on this pre-registered representative list — named here explicitly
rather than swapped in after the fact.

## 9. Pre-registered gates — all 12 PASS (machine-verified, `overall_pass: true`)

```
dt_convergence_spotcheck, conservation_spotcheck, conservation_sweep_35case,
conservation_protocolA, conservation_protocolB,
external_anchor_ratio_all_7joints_x_5intensities, monotonic_decreasing_tlim_vs_intensity,
single_cycle_MF_is_tiny, single_cycle_force_decline_is_negligible,
single_cycle_matches_preregistered_linearization, sustained_protocol_shows_visible_fatigue,
r_modification_directionally_consistent_with_ref4
```

## 10. Honest gaps (disclosed, not hidden)

- **Generic (published, joint-group-level) F/R parameters** — not fit to this subject; explicit
  first-step scope, per the task brief.
- **Biarticular hip-knee muscles** (hamstrings, rectus femoris, sartorius, gracilis) assigned to
  "Knee"; ref [2]'s Table 1 has no hip row, so all primarily hip/trunk-crossing muscles fall back to
  "General" — a disclosed simplification, not derived from the source papers.
- **The "at rest" threshold (TL<0.05)** gating ref [4]'s recovery multiplier is this script's own
  operational choice — the source paper uses exact `TL=0` between discrete contractions, which
  gait's continuous, never-quite-zero (SO floor ≈0.01) activation never literally satisfies.
- **Anchor independence is partial**: ref [3] (the anchor) and ref [2] (the parameters) share an
  author and likely overlapping meta-analysis source studies. The cross-check uses genuinely
  different functional forms (mechanistic ODE-reduction vs. direct empirical power-law regression),
  but is not a fully independent-lab replication.
- **Gait is cyclic/intermittent**; the base model (refs [1][2]) was built/validated for sustained
  isometric holds. This script applies ref [4]'s own correction (r=15) and reports the r=1-vs-r=15
  sensitivity explicitly (§7) rather than silently using the out-of-regime base model — but this is
  still a first-step application of a generic correction to a new task class, not a gait-specific
  re-validation.
- **SO's own known bias** (minimum-effort activation, no antagonist co-contraction — this repo's own
  `metabolic_cost.py` finding, externally verified against Koelewijn et al. 2019) propagates into
  `TL(t)` here: if true muscle activation is higher than SO's solution, true fatigue accrual is
  **under-estimated** by this whole pipeline.
- **No central fatigue, metabolic/pH, or motor-unit size-principle recruitment order** is modeled —
  this is the peripheral, lumped-compartment phenomenological layer only, exactly what refs [1][2]
  themselves are.
- **The repeated-gait protocol tiles the same measured stride** — it does not model real
  within-session gait adaptation (altered kinematics/synergies as an actual walker compensates for
  fatigue).
- **Distinct from `MSK-MUSCLE-FATIGUE-EMG`** (this repo's existing cell, per `docs/MECHANISM_STATE.md`
  §7): that capability detects a fatigue *signature* in surface EMG (median-frequency-down /
  RMS-up spectral shift). This script is a *predictive dynamical model* (motor-unit compartment
  ODE) coupled to SO activation — a complementary, not overlapping, capability.

## 11. Files

- `scripts/msk/muscle_fatigue.py` — the model + all 3 protocols + self-checks + gates. Run with
  `.venv-msk/bin/python3 scripts/msk/muscle_fatigue.py` (~38s wall time).
- `data/msk_smoketest/subject2_walking1/muscle_fatigue/muscle_fatigue_results.json` — full evidence
  (citations, per-joint anchor table, per-muscle single-cycle and sustained trajectories, gates).
