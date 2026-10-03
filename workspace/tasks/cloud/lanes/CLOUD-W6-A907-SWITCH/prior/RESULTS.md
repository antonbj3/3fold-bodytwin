# CLOUD-W4-A907-PARITY — Implicit KKT derivative parity (synthetic method check)

**Verdict against PREREG: PASS on synthetic cases, with one post-hoc correction to the test harness (see below). Parity with the private A907 model: UNKNOWN.**


## Exact command
```
pip install numpy          # tested with numpy 2.4.6, Python 3 (Linux)
python3 kkt_parity.py      # writes results.json (deterministic; sha256 10474ea5…5673558 on this machine)
python3 kkt_parity.py --raw   # also writes full per-case record raw_results.json (~110 KB, not committed)
```
Runtime is about 1.5 s.

## Derivation (done independently)
Force solver: `f*(p) = argmin ½ fᵀWf  s.t.  A(p) f = b(p)`, with W diagonal and positive definite. That makes the problem strongly convex. A (6×10) is full row rank.
KKT system: `K z = r`, with `K = [[W, Aᵀ],[A, 0]]`, `z = [f; λ]`, `r = [0; b]`.
Differentiate `Wf + Aᵀλ = 0` and `Af = b` with respect to a geometry parameter p_j:
```
K [∂f/∂p_j ; ∂λ/∂p_j] = [ −(∂A/∂p_j)ᵀ λ ;  ∂b/∂p_j − (∂A/∂p_j) f ]
```
K is nonsingular exactly when A has full row rank (W ≻ 0). If A loses rank, the derivative does not exist in general. In the bounded mode (`nn`, muscle f ≥ 0), the same formula applies on a fixed active set with strict complementarity. At an active-set switch the solution map has a kink.
∂A/∂p and ∂b/∂p are taken by complex step (h = 1e-30i) on the analytic geometry. This is exact to rounding, and the solver itself is never differentiated.

## Synthetic model (seeded)
- **Setup:** a planar 2-body chain (ground–body 1–body 2), with 2 revolute joints under static equilibrium.
- **Equations:** 6 Newton–Euler equations, 3 per body.
- **Unknowns (10):** 6 muscle forces (2 mono-articular at each joint, 2 bi-articular) plus 4 joint-reaction components.
- **Weights:** W = diag(1/s_i², s_i ∈ [300,1500] N; reactions 1e-2).
- **Varied by seed:** segment lengths, masses m1 ∈ [1.5,3] kg and m2 ∈ [1,2] kg, external wrench at the distal tip (F ∈ [−30,30] N per component, M ∈ [−3,3] N·m), muscle attachments and posture.
- **Geometry parameters p:** q1 and q2 (rad), and d (m), the insertion offset of muscle 1. The gravity and wrench moments depend on p, so both A and b change.
- **Units:** forces in N; Jacobian in N/rad and N/m.
- **Relative error:** ‖(J_fd − J_imp)·diag(1,1,0.1)‖_F / ‖J_imp·diag(1,1,0.1)‖_F. This is the Frobenius norm, with columns scaled to the step scale 1 rad, 1 rad, 0.1 m.
- **Finite differences:** central differences of the complete solver at 5 steps, h ∈ {1e-2, 1e-3, 1e-4, 1e-5, 1e-6} × scale.
- **Stable step:** picked without using the implicit result. It is the finer member of the adjacent pair (h, h/10) whose FD Jacobians differ least, skipping any step whose stencil changes the active set.

## Results — 20 regular seeded cases (seeds 0–19, equality QP, the preregistered target)
| quantity | value |
|---|---|
| cases with relerr < 1e-4 at the identified stable step | **20 / 20** |
| relerr at stable step: min / median / max | 1.88e-10 / 6.22e-10 / 1.24e-9 |
| identified stable step (all cases) | 1e-6 (scaled) |
| independent check: implicit vs complex step through the whole solver, max | 2.0e-14 |
| cond(A) range; cond(K) range | 38.5–66.3; 1.15e4–1.30e6 |
| flags raised on regular cases | none |

### Error / cancellation budget (maximum over the 20 cases)
| h | relerr max (median) | truncation est. ‖D(h)−D(h/10)‖ | rounding est. ε·max‖f‖/(h·rms‖J‖) |
|---|---|---|---|
| 1e-2 | 1.45e-2 (1.27e-4) | 1.43e-2 | 9.2e-14 |
| 1e-3 | 1.47e-4 (1.27e-6) | 1.45e-4 | 9.2e-13 |
| 1e-4 | 1.47e-6 (1.27e-8) | 1.45e-6 | 9.2e-12 |
| 1e-5 | 1.47e-8 (1.72e-10) | 1.48e-8 | 9.2e-11 |
| 1e-6 | 1.24e-9 (6.22e-10) | ~1.5e-10 (h² extrapolation) | 9.2e-10 |

From 1e-2 down to 1e-5 the error falls by a factor of 100 per decade, which is the O(h²) truncation expected of central differences. At 1e-6, rounding and cancellation dominate. The ~1e-9 error there matches the rounding estimate: about 6 of 16 digits are lost in f(p+h) − f(p−h). The implicit and complex-step Jacobians agree to about 1e-14, so the residual FD error is finite-difference error, not a defect in the derivative. The 5-step grid does not reach deep into the rounding-dominated region, so the "optimal" step may lie between 1e-5 and 1e-6. Either choice satisfies the criterion by more than 4 orders of magnitude.

**Uncertainty:** the FD reference is itself uncertain at about 1e-9 relative (from the budget above). The implicit derivative is therefore verified to about 1e-9 by FD and to about 1e-14 by complex step. These figures cover this synthetic model family only. The worst case at the coarse step h = 1e-2 is 1.4e-2 relative error, which does not pass. A careless single-step FD would fail the criterion.

## Failure cases / counterexamples (all flagged)
| case | construction | flag(s) | behaviour |
|---|---|---|---|
| seed 100 | every joint-2 muscle line passes through the joint-2 centre, chain straight | `SINGULAR_A` (cond A = 1.6e17) | not differentiable; no derivative reported |
| seed 101 | the same with 1e-9 m attachment offsets | `ILL_CONDITIONED_A` (cond A = 3.1e9, cond K = 2.7e12) | FD disagrees with implicit by about 100% at every step, while implicit vs complex step is 1.6e-8. This is a counterexample: in this case the FD "stable step" is meaningless. |
| seeds 0, 2, 6, 8 (nn mode) | q1 bisected onto an active-set switch | `WEAK_COMPLEMENTARITY` (margin ≤ 8e-13 N), `ACTIVE_SET_CHANGE` at all 5 steps, `NO_STABLE_STEP` | FD vs one-sided implicit differ by 49–79% at every step (a kink) |

Exploratory `nn` mode on the 20 regular seeds: 9 are feasible and all 9 pass (relerr ≤ 7.8e-10). The other 11 are flagged `SOLVE_FAILED: no feasible active set`. I confirmed that this infeasibility is real, not a solver bug. In this generator every muscle has a positive moment arm at joint 2, so nonnegative muscle forces can balance only one sign of joint-2 load. In the equality mode, muscle forces change sign within the stencil in 2 cases at h = 1e-2 and in 1 case at h = 1e-3 and 1e-4. This does not matter for the equality QP. It would be an active-set change if bounds were present.

## Preregistered criterion
- "20 of 20 seeded regular cases have relative error below 1e-4 at an identified stable step": **PASSED** (20/20, max 1.24e-9).
- "all singular/active-set changes are flagged": **PASSED** in the final run (6 of 6 constructed cases are flagged). This depends on the correction below.
- Private A907 numerical parity (A907 reports 2.6e-8 on its private model; I used that only as context): **UNKNOWN**. I did not reproduce or compare it.

## Post-hoc corrections (recorded separately; the criterion itself was not changed)
1. **First run: overall `PASS=false`.** The criterion part was 20/20, but one constructed active-set test case (seed 102) could not be built. I had tried to place it at a kink by perturbing M_ext, and that made the nonnegative problem infeasible. The resulting record had no flag, so "all flagged" evaluated false. This was a failure to build the test case, not a missed kink.
2. **Second attempt (base seeds 0 and 2, still tuning M_ext).** This landed 7.5e-4 N from the switch (flagged at h = 1e-2 only), and seed 2 again became infeasible.
3. **Final construction.** Bisect q1 onto an actual active-set switch (60 iterations) on nn-feasible base seeds 0, 2, 6 and 8. The flag thresholds (cond A > 1e8, σ_min/σ_max < 1e-10, margin < 1e-9 N) and the regular-case definition were fixed before the first run and did not change.

## Limitations
- This is a synthetic, planar, static, 2-body model with 10 unknowns. It is not a validated anatomical model, and cond(K) is at most 1.3e6. Large 3-D models with many muscles and worse conditioning may lose more digits.
- ∂A/∂p comes from complex step on my own geometry code. It checks the KKT formula, not a separate hand-coded analytic ∂A/∂p.
- The flag thresholds are heuristic. A near-singularity just below cond 1e8, or an active-set switch farther than 1e-2 × scale from p0, would not be flagged. The second of these does not affect the FD comparison.
- The `nn` mode uses exhaustive active-set enumeration (2⁶ subsets). This is exact but only practical for small models.
- **Sources:** I could not check any external source. The sandbox proxy denied outbound requests to doi.org and arxiv.org (HTTP 403 on CONNECT, checked 2026-09-24). No value in this report depends on an external source. The derivation is standard implicit differentiation of KKT conditions (related background, *not checked here*: Amos & Kolter 2017, OptNet, arXiv:1703.00443; Squire & Trapp 1998, complex-step, doi:10.1137/S003614459631241X). The only package fetched was numpy 2.4.6 from pypi.org.
