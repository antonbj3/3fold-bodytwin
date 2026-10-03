# CLOUD-W6-STRENGTH-ID — Results

**Verdict: FAIL.** Neither part of the preregistered gate holds.
- **Held-out RMSE:** the proposed mechanistic curve does not reduce held-out 90 deg/s RMSE by ≥10% in both directions (extension −2.6%, flexion +9.4%, i.e. worse).
- **Identifiability:** no individual *muscle* parameter (Fmax, activation, moment arm, Lopt, Lslack, vmax) has full-rank sensitivity from these measurements.

The only quantity that is well determined and stable across strokes is a joint-level one: the angle of peak moment at 60 deg/s, θ_opt ≈ 53–57° knee flexion for extension and 52–56° for flexion. This is **not** a validated individual muscle model.

## Data, provenance, conventions
- **Source:** Grand Challenge Competition to Predict In Vivo Knee Loads, subject SC, bundled Biodex CSVs. DOI 10.1002/jor.22023; https://pmc.ncbi.nlm.nih.gov/articles/PMC4067494/. This is one implanted knee, so none of it is a population estimate.
- **SHA-256 checks:** `SC_isokin60_biodex.csv` d743d150…c757 and `SC_isokin90_biodex.csv` 1252b575…5b3c. Both match `SOURCE_SHA256.txt`. `PREREG.md` e1e23987… matches `PREREG.sha256`.
- **Sign convention** (declared before scoring):
  - velocity > 0 means the Biodex position is increasing and the knee flexion angle is decreasing. This is an **extension stroke**, and torque > 0 is an extension moment.
  - velocity < 0 is a **flexion stroke**, and torque < 0 is a flexion moment.
  - Torques are used as logged. No gravity correction is applied.
- **Measured quantity:** the external dynamometer moment in Nm. It is not muscle force. The moment arm, force–length and fascicle data are absent. The CSVs do contain a "Filtered vaslat EMG" column, which conflicts with the brief's statement that EMG is not in the bundle. It has no MVC normalisation or units, so I did not use it or fit it.
- **Re-derivation from raw CSVs:**
  - Qualified peaks are 53.625 Nm (60 deg/s) and 40.512 Nm (90 deg/s), matching the prior report.
  - My 10° bin medians match the prior `curve_60.csv` to within 0.0005 Nm. The prior file's `increasing_angle` label means increasing Biodex position, i.e. extension.

## Method (`strength_id.py`)
- **Valid samples:** |v| ≥ 0.8 × nominal speed. Direction comes from the sign of v. A stroke is a contiguous valid run of at least 10 samples. That gives 60 deg/s: 3 extension and 3 flexion strokes (493 and 568 samples), and 90 deg/s: 4 extension and 5 flexion strokes (435 and 521 samples).
- **Preregistered baseline:** per direction, medians in 10° bins (edges at multiples of 10°, n ≥ 5) from the 60 deg/s data, linearly interpolated in angle, with scale = 1. No public, subject-independent speed prior was allowed, so the "speed-scaled baseline" with scale = 1 is identical to the angle-only baseline. They are one baseline here.
- **Mechanistic joint curve:** τ(θ) = A·exp(−((θ−θ_opt)/w)²) + c, fit by least squares on the 60 deg/s samples only. Here A lumps Fmax·a·r·f_v(60)·(unit factor).
- **Scoring:** held-out RMSE (Nm) on all valid 90 deg/s samples inside the 60 deg/s bin-centre range. Both models are scored on the same samples, with no extrapolation. That leaves 419 of 435 extension samples and 495 of 521 flexion samples; the rest are UNKNOWN. Scoring on 10° bins is secondary.
- **Uncertainty:**
  - leave-one-60-stroke-out refits;
  - a cluster bootstrap over 90 deg/s strokes (2000 reps, seed 20260924);
  - Jacobian singular values.
- **Structural identifiability:** the numeric Jacobian SVD of a Hill-type single equivalent muscle, at 5 random linearisation points. Parameters are Fmax, a, r, Lopt, Lslack, L0_mt, vmax factor and Hill curvature k, with a rigid tendon and a constant moment arm.

## Primary results (threshold 0.8, 10° bins)
| Direction | Baseline RMSE 90 (Nm) | Mechanistic RMSE 90 (Nm) | Reduction | Bootstrap 95% CI | Bin-level reduction | ≥10%? |
|---|---|---|---|---|---|---|
| Extension | 23.90 | 23.27 | +2.6% | [+0.9%, +6.2%] | +9.3% | no |
| Flexion | 5.49 | 6.00 | −9.4% | [−12.9%, −7.0%] | −7.4% | no |

**60 deg/s-only fit (training RMSE):**
- Extension: 7.63 Nm mechanistic vs 7.35 Nm baseline.
- Flexion: 3.56 Nm vs 3.50 Nm.

**Fitted parameters:**
- **Extension:** A = 48.7 Nm, θ_opt = 54.7°, w = 41.7°, c = 0.35 Nm.
- **Flexion:** A = −60.8 Nm, θ_opt = 54.2°, w = 73.4°, c = +40 Nm. The offset c sits at its bound.

**Repeated strokes (leave-one-60-stroke-out), which show practical non-identifiability of amplitude vs offset:**
- **Extension:** A ranges 41.9–87.5 Nm and c ranges −40 to +10.7 Nm, but θ_opt only 53.0–57.3°.
- **Flexion:** c stays at the bound in every refit, A ranges −58.5 to −62.6 Nm, and θ_opt 52.0–56.1°.
- The Jacobian condition numbers (32 extension, 297 flexion) are finite. So A, θ_opt, w and c are locally full rank, but A and c are not practically identifiable.

**Held-out RMSE per 90 deg/s extension stroke:** 35.9, 23.6, 10.6 and 24.2 Nm. The error is dominated by the counterexamples below, not by curve shape.

**Diagnostic only (uses 90 labels, not part of the gate):**
- The least-squares 90/60 scale is 0.50 for extension (RMSE would drop to 14.2 Nm) and 0.68 for flexion (3.07 Nm).
- The ratios of qualified peaks are 0.755 (extension) and 0.717 (flexion).
- So the dominant held-out error is the untrained speed effect. Training at one speed cannot identify it.

## Sensitivity to angle bins and the valid-velocity threshold
18 cases, all completed in 0.6 s (see `results.json` → `sensitivity_sweeps`): threshold 0.7, 0.8 or 0.9; bin width 5°, 10° or 15°; bin offset 0 or half a bin.
- **Extension:** the sample-level reduction ranges from +0.1% to +4.3%.
- **Flexion:** it ranges from −3.0% to −13.7%.
- **Bin-level reduction:** extension reaches ≥10% in 4 of 18 cases (at most 10.5%), and flexion is negative in all 18.

No case passes in both directions, so the FAIL does not depend on the bin choice.

## Identifiability and gauges (structural, Hill-type muscle)
- **Training at 60 deg/s only:** the Jacobian rank is 3 of 8 at all 5 points, leaving **5 null directions**. In log space these span:
  1. **Fmax ↔ activation gauge.** Only Fmax·a enters, so these two are exactly degenerate.
  2. **Length-scale gauge.** Moment arm r, Lopt and (L0−Lslack) scale by c while Fmax·a scales by 1/c. This leaves every fl and fv argument, and Fmax·a·r, unchanged.
  3. **L0_mt ↔ Lslack.** Only fibre length lf = L0 − rθ − Lslack is observable. **Lopt and Lslack therefore cannot be inferred uniquely.** From the angle of peak moment you get L0 − Lslack − rθ_opt = Lopt: one equation in four unknowns (L0, Lslack, Lopt, r).
  4. and 5. **vmax factor and Hill curvature k** both collapse into the amplitude at a single speed.
- **With 60 and 90 deg/s both used:** the rank is 4, leaving **4 null directions**. Gauges 1–3 remain, plus one force–velocity direction (vmax ↔ k). Two speeds give two amplitudes, which is not enough for three unknowns (isometric amplitude, vmax, k).
- **What is identifiable:** at most the joint-level quantities A₆₀ (Nm), θ_opt (deg) and w (deg). Of these, A is practically confounded with the offset c.

## Concrete next measurement
**Measurement:** isometric maximal voluntary moment at about 55° knee flexion (the stable θ_opt), in both directions, at least 3 repetitions. Together with 60 and 90 deg/s this gives three speeds for the three force–velocity unknowns. It removes the vmax ↔ k null direction, but not gauges 1–3.

**Predicted signal:** this is a design calculation using the 90 deg/s peak ratio and a joint-space Hill curve, not a validation (`next_measurement.py`). Isometric/60 deg/s torque would be:
- **Extension:** 1.85, 1.96 or 2.11 for k = 0.40, 0.25 or 0.15. On the 51.6 Nm bin median that is roughly 95, 101 or 109 Nm.
- **Flexion:** 2.06, 2.25 or 2.49. On the 20.1 Nm bin median that is roughly 41, 45 or 50 Nm.

The ~14 Nm (extension) and ~9 Nm (flexion) spread between curvature hypotheses is about 2× the per-sample residual SD at 60 deg/s (7.2 and 3.5 Nm). A few repetitions should therefore separate them.

**Other gauges:**
- The Fmax ↔ a gauge needs an interpolated-twitch activation measurement. I cannot derive a magnitude for it from this bundle, so it is UNKNOWN.
- The r/Lopt scale gauge needs an imaged moment arm or ultrasound fascicle length. Its magnitude is also UNKNOWN from the bundle.

## Counterexamples and limitations
- **Non-isokinetic, spiky 90 deg/s extension strokes.** They contain large negative torques in the middle of the stroke (e.g. −18.4 Nm at 52.7°, −26.6 Nm at 23°), and the early velocity lags at 62–86 deg/s. These look like impact or oscillation artefacts, and no smooth angle model can predict them. Stroke 1 has an RMSE of 36 Nm.
- **Offset at bound.** The flexion offset sits at the ±40 Nm bound. A broad Gaussian plus an offset is nearly a quadratic, so A and c trade off.
- **Scope.** This is one subject and 3 training strokes per direction, so stroke-level uncertainty is coarse. There is no gravity correction, and the EMG column was not used.
- **What model fit shows.** Model fits alone do not validate muscle parameters. Everything here is measured external joint moment.

## Preregistered criteria
| Criterion | Status |
|---|---|
| Split by speed; 60 train / 90 held-out | PASS (done) |
| Directions analysed separately; sign convention declared first | PASS (done) |
| ≥10% held-out RMSE reduction, extension | **FAIL** (+2.6%) |
| ≥10% held-out RMSE reduction, flexion | **FAIL** (−9.4%) |
| Full-rank sensitivity for claimed individual muscle parameters | **FAIL** (5 null directions at 60 deg/s only; 4 with 60 + 90) |
| 60-only fit, stroke uncertainty and gauges reported | PASS (done) |
| 90 deg/s angles outside the training range | UNKNOWN (16 extension and 26 flexion samples excluded) |
| Independent measurements (moment arm, fascicle length, activation) | UNKNOWN (missing) |

**Command:** `python3 strength_id.py && python3 next_measurement.py` (Python 3, numpy 2.4.6, scipy 1.17.1; runtime < 1 s; all sweeps COMPLETE)
