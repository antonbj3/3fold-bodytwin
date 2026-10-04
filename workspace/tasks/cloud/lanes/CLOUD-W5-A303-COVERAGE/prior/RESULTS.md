# CLOUD-W4-A303-BOUND: Emulator error bounds at active-set switches

**Verdict: FAIL.** The preregistered criterion is not met. Coverage part (b) passed. Max-error part (a) failed: the worst relative error was 22.8%, against a limit of 3%.


## What was done
- **Frozen plan.** `ANALYSIS_PLAN.md` was committed together with the code before the evaluation run (commit `333d335`). PREREG.md named no strata, so the plan fixed them before any results were seen.
- **Families** (3 instances each, deterministic seeds; θ ∈ [0,1]², normalised units):
  - **F1 `muscle`:** a synthetic planar-arm static-optimisation QP. 8 muscles; minimise Σ(fᵢ/Nᵢ(q))² subject to R(q)ᵀf = τ(q, load) and 0 ≤ f ≤ N(q). The output is y = Σf (N). It has 5–8 distinct active sets on a 21×21 grid.
  - **F2 `generic`:** a random strictly convex pQP with n = 6 and 11 inequalities. Q(θ), c(θ) and h(θ) are nonlinear in θ. The output is y = wᵀx (arbitrary units). It has 11–34 active sets.
- **Exact reference.** quadprog (Goldfarb–Idnani), then an active-set KKT polish, then KKT verification: primal and dual feasibility, and stationarity. All 11,400 evaluated points passed.
- **Emulators.** All use the same 11×11 anchor grid (spacing 0.1) and the nearest anchor.
  - **E1 (primary), one-sided derivative:** y(θ₀) + h·D_u y(θ₀). D_u is the exact one-sided directional derivative from the directional-derivative QP: strongly active constraints become equalities and weakly active ones become inequalities.
  - **B1, global smooth:** thin-plate-spline RBF interpolation over the anchors.
  - **B2, set-aware:** a KKT re-solve at θ with the anchor's active set frozen.
- **90% intervals.** Mondrian split conformal, calibrated on 6,000 uniform points kept separate from the test sets. Groups are family × predicted-switch flag. The score is the relative error divided by h² (flag 0) or by h (flag 1).
- **Strata.** True distance d from the query θ to the nearest active-set change: the minimum over 32 rays of the exit from the frozen-set KKT region, found by a 0.001 scan followed by bisection. S1 [0, .005), S2 [.005, .02), S3 [.02, .05), S4 ≥ .05.
- **Test sets.** 2,400 held-out stratified points (100 per instance per stratum) for the criterion. A secondary set of 3,000 uniform points was also evaluated.

## Results (relative error |ŷ−y|/|y|, dimensionless)

| Emulator | Stratified test: max | p99 | median | share >3% | Uniform test: max |
|---|---|---|---|---|---|
| E1 one-sided | **0.2283** | 0.0830 | 0.00104 | 4.75% | 0.2273 |
| B1 global smooth | 0.1640 | 0.0659 | 0.00171 | 4.29% | 0.1366 |
| B2 set-aware | 0.2124 | 0.0795 | 0.0 (exact within the region) | 4.54% | 0.2269 |

E1 90% interval coverage on the stratified test (n = 600 per stratum; Wilson 95% CI; range across the 6 instances):

| Stratum | Coverage | 95% CI | Muscle / generic | Instance range | Max rel. err |
|---|---|---|---|---|---|
| S1 [0, .005) | 0.902 | 0.875–0.923 | 0.903 / 0.900 | 0.81–0.98 | 0.183 |
| S2 [.005, .02) | 0.892 | 0.864–0.914 | 0.887 / 0.897 | 0.81–0.99 | 0.210 |
| S3 [.02, .05) | 0.870 | 0.841–0.895 | 0.877 / 0.863 | 0.72–0.97 | 0.228 |
| S4 ≥ .05 | 0.920 | 0.896–0.939 | 0.910 / 0.930 | 0.83–1.00 | 0.127 |

The baselines on the same intervals:
- **B1** under-covers near switches: 0.81 in S1, 0.80 in S2 and 0.85 in S3. It over-covers in S4 (0.97).
- **B2** over-covers away from switches: 0.95 in S2, 0.97 in S3 and 0.998 in S4.

Neither baseline would pass (b).

## Preregistered criterion
| Part | Result | Pass? |
|---|---|---|
| (a) E1 max relative error ≤ 3% on ≥ 1000 held-out points | 22.83% max over 2,400 points; 114 points exceed 3% | **FAIL** |
| (b) 90% intervals cover 85–95% in every stratum (pooled over families) | 0.870–0.920 | PASS |
| (b, stricter reading, per family × stratum; secondary) | 0.863–0.930 | PASS |
| **Overall** | | **FAIL** |

## Independent checks
- **SLSQP (scipy) as a second solver**, on 300 stratified test points:
  - The maximum relative difference in y from the reference solution was 8.8e-8.
  - The median relative difference was 6e-14.
  - SLSQP reported success on only 257 of the 300 points. The other 43 did not report success but still agreed to within that bound.
- **Finite-difference check of one-sided derivatives** at 180 anchor/direction pairs: the maximum relative mismatch was 3.7e-5, which is consistent with the O(ε) error of the finite difference.
- **KKT verification** passed at every reference point.

## Counterexamples and failure analysis
- **Worst case.** On generic0 at θ = (0.848, 0.467), y = 1.503 and E1 gives ŷ = 1.160, a relative error of 22.8%. The query is d = 0.023 from the nearest switch, and the step from the anchor is h = 0.058. The other four worst cases are also on generic0, with errors of 15–21%.
- **Crossing a switch is the main cause.** This split was exploratory and not preregistered; a switch counts as crossed when B2 is not exact. Of the 114 points above 3%, 110 cross a switch. The E1 max is 22.8% with a crossing (n = 761) and 5.7% without one (n = 1,639). So even with no switch crossed, curvature at h ≤ 0.071 already breaks the 3% limit on the generic family.
- **The flag catches most crossings.** The linearised switch flag caught 96.7% of crossings. It fired on 1.8% of points that did not cross.
- **Distance from the query is not the right stratifier for error.** The largest errors are in S3 and S4, where the query is far from any switch but the anchor-to-query segment crosses one. Coverage by query-distance stratum can look fine while the maximum error stays large.
- **Conditional coverage is fragile.** On the secondary uniform set, pooled E1 coverage in S3 was 0.829 (CI 0.798–0.856), outside 85–95%. It fell as low as 0.65 on one instance. The per-instance range in the stratified test reached 0.72. Pooled coverage passing (b) does not mean every instance is covered.

## Uncertainty and limitations
- **Coverage.** Wilson CIs are shown above. With 600 points per stratum the half-width is about ±2.5 pp. S3 (0.870, lower CI 0.841) could be below 0.85 in truth.
- **The max error is a sample extreme.** It depends on the anchor spacing (0.1), the domain and the family design. A finer grid would shrink both the curvature error and how often switches are crossed. This was not tested, because that would be a post-hoc change.
- **Distance to switch is approximate.** It is a 32-ray approximation (overestimates by at most about 0.5% for straight boundaries). Slivers thinner than the 0.001 scan step could be missed.
- **The model is synthetic.** The muscle family is only a stylised stand-in for musculoskeletal static optimisation, with invented geometry and strengths.
- **Stratified sampling changes the distribution.** The test set oversamples near-switch points compared with the uniform calibration distribution, so marginal coverage on it is not the deployment coverage.
- **Relative error depends on output scale.** Generic instances have y as low as about 1.26.

## Sources
- **No external datasets or papers were accessed or checked in this session.**
- **Software** (installed from PyPI): numpy 2.4.6, scipy 1.17.1, quadprog 0.1.13 (https://pypi.org/project/quadprog/).
- **Standard methods used without checking a citation here:** Goldfarb–Idnani dual active-set QP, parametric-QP directional sensitivity, and split/Mondrian conformal prediction. No DOIs are given because none were checked.

## Reproduce
```
pip install numpy scipy quadprog
python3 emu_bound.py              # ~9.5 min on 4 cores; writes raw_points.json (gitignored) + results.json
python3 emu_bound.py --summarize  # recompute results.json from raw_points.json
```
