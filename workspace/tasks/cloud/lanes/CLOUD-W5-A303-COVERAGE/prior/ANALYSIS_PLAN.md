# Analysis plan (written before any evaluation run; supplements, does not modify, PREREG.md)

PREREG.md froze the criterion but did not enumerate the switch strata or the QP
families. The following operational choices were fixed *before* the first
evaluation run and committed in the same commit as the code (see git log).

## Families (original synthetic, no external data)
- F1 `muscle`: planar 2-link arm static-optimisation QP. 8 muscles, variables f (N);
  min sum (f_i/N_i(q))^2 s.t. R(q)^T f = tau(q, load), 0 <= f <= N_i(q).
  theta1 -> elbow angle q = 0.3 + 1.8*theta1 rad, theta2 -> hand load 0..m_max.
  m_max = 0.9 x (min over a 41-point elbow grid of the LP-maximal feasible load).
  Output y = sum_i f_i (N).
- F2 `generic`: random strictly convex pQP, n=6, x>=0, sum x >= 1+0.5 theta1,
  4 random general constraints with guaranteed-feasible nonlinear RHS, Q(theta), c(theta)
  nonlinear. Output y = w^T x (w in [0.5,1.5]). Seeds rejected if <4 distinct active
  sets on a 21x21 grid.
- 3 instances per family (deterministic seeds), theta in [0,1]^2 (normalised units).

## Exact solver and checks
quadprog (Goldfarb-Idnani) -> active set -> KKT polish -> KKT verification
(primal/dual feasibility, stationarity). Independent check: scipy SLSQP on a
random subset of test points; finite-difference check of one-sided derivatives.

## Emulators (all use the same 11x11 anchor grid per instance, spacing 0.1)
- E1 one-sided derivative emulator (primary): nearest anchor theta0,
  y_hat = y(theta0) + h * D_u y(theta0), D_u = exact one-sided directional derivative
  from the directional-derivative QP (strongly active -> equalities, weakly active -> inequalities).
- B1 global smooth: thin-plate-spline RBF interpolation of y over anchors.
- B2 set-aware: KKT re-solve at theta with the anchor's active set frozen.

## Intervals (90%)
Mondrian split conformal on relative score |y_hat - y|/|y_hat| divided by s,
groups = family x predicted-switch flag (flag = linearised slack or multiplier of the
anchor changes sign within the step); s = h^2 (flag 0) or h (flag 1).
Calibration: 1000 uniform points per instance (separate seed).

## Switch strata (true distance d to nearest active-set change, normalised theta units)
d = min over 32 ray directions of the first exit from the frozen-active-set KKT region
(scan step 0.001 to 0.06, then 40-step bisection).
S1 [0, 0.005), S2 [0.005, 0.02), S3 [0.02, 0.05), S4 [0.05, inf).
Test set: stratified acceptance sampling, 100 points per instance per stratum
(1200 per family, 2400 total, all held out). Secondary: 500 uniform test points/instance.

## Criterion evaluation (frozen in PREREG)
PASS iff (a) max relative error of E1 over the 2400 stratified test points <= 3% AND
(b) E1 90% interval coverage in [85%, 95%] in each of S1..S4 (pooled over families).
Per-family x stratum coverage is reported as a secondary, stricter reading.
Uncertainty: Wilson 95% CIs on coverage; instance-level ranges.
