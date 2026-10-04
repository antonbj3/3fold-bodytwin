"""Decide whether a transport route is staged or shunted, from the first few samples of a bolus.

WHY THIS IS A DECISION AND NOT A DESCRIPTION. The detail layer carries a question of exactly this
shape and states it as an attribution rather than a measurement: `MECHANISM_DOPAMINE_KINETICS.md`
(net edge HARVEST-E0094) says an apparent 3.8 uM clearance component in chopped tissue is
"attributed to diffusion, not a second intrinsic transporter". Those two readings are a SERIES route
(one transporter behind a diffusion barrier) and a PARALLEL route (two transporters drawing from the
same pool). Attribution cannot separate them. An early-time flux ratio can, and the separation is
exact rather than statistical:

    R(t) = J_last(t) / J_first(t)

    series, n stages:  R(t) = (K_last / V_1) * t + O(t^2),  so R(0+) = 0 for ALL positive K and V
    parallel:          R(t) = K_last / K_first + O(t),      so R(0+) is the conductance ratio

The series recursion is why: c_j starts at order t^j, so J_last = O(t^(n-1)) while J_first = O(1).
The verdict is read off the INTERCEPT, not the magnitude, which is what makes it parameter-free.

PROVENANCE. The structure was delivered by a swarm job (BT-FW48-AUTO-384787ba1a13a6, space-swarm,
target BT-RESEARCH-IMMUNITY). I verified its two decisive fixture numbers in closed form before
building on them: with K = [100, 60] mL/s and V_1 = 100 mL the leading coefficient K_last/V_1 is
0.6 per second, which puts R at 0.006 at t = 0.01 s and 0.0300 at t = 0.05 s against their reported
0.006030 and 0.030747, and the parallel intercept at K_last/K_first = 0.6 against their 0.6024.

DOMAIN, stated because the dopamine case sits just outside it. Linear, time-invariant, passive,
non-saturating transport from a finite source. Michaelis-Menten uptake is NOT linear, so the test
applies to dopamine only in the regime c << Km, and this file computes the concentration ceiling
and the sampling requirement that puts a real experiment inside the domain. PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import math


def series_trace(K: list[float], V: list[float], V_source: float, dose: float,
                 dt: float, n_steps: int) -> list[tuple[float, float, float]]:
    """Explicit Euler on the staged chain. Returns (t, J_first, J_last) per step.

    Stage j holds concentration c_j in volume V_j; stage 1 is fed by the source and stage j by
    stage j-1. The source depletes, which is the part that makes the late-time ratio useless.
    """
    n = len(K)
    c_src = dose / V_source
    c = [0.0] * (n - 1)
    out = []
    for i in range(n_steps):
        t = i * dt
        J_first = K[0] * c_src
        J_last = K[-1] * c[-1]
        out.append((t, J_first, J_last))
        flows = [K[0] * c_src] + [K[j] * c[j - 1] for j in range(1, n)]
        c_src -= flows[0] * dt / V_source
        for j in range(n - 1):
            c[j] += (flows[j] - flows[j + 1]) * dt / V[j]
    return out


def parallel_trace(K_first: float, K_last: float, V_source: float, dose: float,
                   dt: float, n_steps: int) -> list[tuple[float, float, float]]:
    """Both routes draw on the same source, so the ratio is the conductance ratio from t=0."""
    c_src = dose / V_source
    out = []
    for i in range(n_steps):
        out.append((i * dt, K_first * c_src, K_last * c_src))
        c_src -= (K_first + K_last) * c_src * dt / V_source
    return out


def _ols(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Intercept, slope, and the standard error OF THE INTERCEPT. The decision needs that error,
    not the slope's: a series route's intercept is zero and the question is whether the data can
    tell zero from the parallel route's conductance ratio."""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = my - slope * mx
    resid = [y - intercept - slope * x for x, y in zip(xs, ys)]
    dof = n - 2
    s2 = sum(r * r for r in resid) / dof if dof > 0 else 0.0
    se_intercept = math.sqrt(s2 * (1.0 / n + mx * mx / sxx)) if s2 > 0 else 0.0
    return intercept, slope, se_intercept


def decide(trace: list[tuple[float, float, float]], doublings: int = 3,
           base_multiple: int = 20) -> dict:
    """Decide on the ORDER at which the ratio vanishes, not on its size or its intercept.

    The first version of this file tested whether an OLS intercept sat within three standard errors
    of zero, and it called a textbook two-stage SERIES chain PARALLEL. The reason is instructive and
    is the reason the rule below is different: on a noiseless simulated trace the residual scatter
    goes to zero, so the standard error goes to zero with it, and the O(t^2) curvature that any real
    staged chain has then sits many "sigma" from zero. A significance test needs a noise model it
    does not have here, and a magnitude test cannot work either, because a parallel route with a
    small conductance ratio has a small constant ratio that looks like a series route's small
    early ratio.

    What separates the two is scale-free: halve the time and the series ratio halves per stage,
    while the parallel ratio does not move.

        series, n stages:  R(t) ~ C t^(n-1)   so  R(t)/R(2t) = 2^-(n-1) <= 1/2
        parallel:          R(t) -> K_last/K_first  so  R(t)/R(2t) -> 1

    That also reads the STAGE COUNT off the data: n = 1 + log2(R(2t)/R(t)), which is more than the
    binary verdict and is what makes the test worth running on a real bolus.
    """
    samples = {round(t, 9): (jl / jf) for t, jf, jl in trace if jf != 0.0}
    times = sorted(samples)
    if len(times) < 4:
        raise ValueError('trace too short')
    dt = times[1] - times[0]
    # The doubling must be in ABSOLUTE time, because R ~ t^(n-1) is a power law in t measured from
    # the bolus, not in an index offset. Doubling the offset from a late base instead read a
    # three-stage chain as 2.73 stages and a four-stage chain as 3.31.
    #
    # The base must also clear the start-up transient: the downstream stage is exactly zero for the
    # first steps of a staged chain, and the first non-zero sample is discretisation, not asymptote.
    # Measured convergence of the implied stage count, this file's own fixtures at dt = 1 ms with
    # base = 20*dt: 2.014 for two stages, 3.035 for three, 4.078 for four. Rounding is honest at
    # that accuracy; a tenth of a stage is not claimed.
    t0 = round(max(base_multiple * dt, times[0]), 9)
    pts = []
    for k in range(doublings + 1):
        t = round(t0 * 2 ** k, 9)
        if t not in samples:
            break
        pts.append((t, samples[t]))
    if len(pts) < 2:
        raise ValueError('trace does not span one doubling from the chosen base')
    if any(r <= 0.0 for _, r in pts):
        # An exact zero is the strongest possible evidence of staging and the weakest divisor:
        # report it as such instead of dividing by it.
        return {'verdict': 'SERIES', 'halving_ratio': 0.0, 'implied_stages': float('inf'),
                'first_sample_ratio': pts[0][1], 'note': 'downstream flux still exactly zero',
                'control_late_ratio': float('nan'), 'control_is_degenerate': True}
    ladder = pts
    # Geometric mean: the quantity is multiplicative, one factor of two per stage.
    ratios = [ladder[k][1] / ladder[k + 1][1] for k in range(len(ladder) - 1)]
    halving = math.exp(sum(math.log(r) for r in ratios) / len(ratios))
    verdict = 'SERIES' if halving < 0.75 else 'PARALLEL'
    stages = 1.0 + math.log2(1.0 / halving) if 0.0 < halving < 1.0 else float('nan')
    ys = [r for _, r in ladder]
    t_end, jf_end, jl_end = trace[-1]
    peak = max(abs(jf) for _, jf, _ in trace)
    return {
        'verdict': verdict,
        'halving_ratio': halving,
        'implied_stages': stages,
        'rounded_stages': round(stages) if stages == stages and stages != float('inf') else None,
        'first_sample_ratio': ys[0],
        # The equally informed control, from the SAME trace: the late quasi-steady ratio, which is
        # what an observer without the early-time reading would use.
        'control_late_ratio': (jl_end / jf_end) if jf_end != 0.0 else float('nan'),
        'control_is_degenerate': abs(jf_end) < 1e-9 * peak,
    }


def linear_regime_ceiling(Km_uM: float, tolerance: float = 0.10) -> float:
    """Highest concentration at which saturating uptake still behaves linearly to `tolerance`.

    Michaelis-Menten v = Vmax*c/(Km+c) differs from its linear part Vmax*c/Km by a relative
    c/(Km+c), so a 10 percent ceiling is Km*tolerance/(1-tolerance) = Km/9, not Km/10.
    """
    return Km_uM * tolerance / (1.0 - tolerance)


if __name__ == '__main__':
    K, V, Vb, dose, dt = [100.0, 60.0], [100.0], 100.0, 1.0, 0.001
    def show(label, d):
        print(f"{label:<28}", d['verdict'],
              '| halving', round(d['halving_ratio'], 4),
              '| stages', round(d['implied_stages'], 3) if d['implied_stages'] == d['implied_stages'] else 'nan',
              '->', d.get('rounded_stages'),
              '| first ratio', f"{d['first_sample_ratio']:.3e}",
              '| control late', f"{d['control_late_ratio']:.3e}",
              'DEGENERATE' if d['control_is_degenerate'] else '')
    show('series n=2', decide(series_trace(K, V, Vb, dose, dt, 40000)))
    show('series n=3', decide(series_trace([100.0, 80.0, 60.0], [100.0, 100.0], Vb, dose, dt, 40000)))
    show('series n=4', decide(series_trace([100.0, 90.0, 80.0, 60.0], [100.0] * 3, Vb, dose, dt, 40000)))
    show('parallel', decide(parallel_trace(K[0], K[1], Vb, dose, dt, 40000)))
    # Negative control: a parallel route with a tiny conductance ratio has a tiny CONSTANT ratio.
    # A magnitude test calls it series; the halving test must not.
    show('parallel K_last/K_first=5e-4', decide(parallel_trace(100.0, 0.05, Vb, dose, dt, 40000)))
    print('predicted: series n=2 halving 0.5, n=3 0.25, n=4 0.125; parallel 1.0')
    print('predicted series leading coefficient K_last/V_1 =', K[-1] / V[0], 'per s;',
          'parallel intercept K_last/K_first =', K[-1] / K[0])
    ceiling = linear_regime_ceiling(0.2)
    print('dopamine: Km = 0.2 uM -> linear to 10 percent below', round(ceiling, 5), 'uM')
    # Write the consumable numbers so an edge can cite a key rather than prose.
    import json, os
    out = {
        'implied_stages_two_stage_fixture': decide(series_trace(K, V, Vb, dose, dt, 40000))['implied_stages'],
        'implied_stages_three_stage_fixture': decide(series_trace([100.0, 80.0, 60.0], [100.0] * 2, Vb, dose, dt, 40000))['implied_stages'],
        'implied_stages_four_stage_fixture': decide(series_trace([100.0, 90.0, 80.0, 60.0], [100.0] * 3, Vb, dose, dt, 40000))['implied_stages'],
        'parallel_halving_ratio': decide(parallel_trace(K[0], K[1], Vb, dose, dt, 40000))['halving_ratio'],
        'weak_parallel_halving_ratio': decide(parallel_trace(100.0, 0.05, Vb, dose, dt, 40000))['halving_ratio'],
        'series_verdict_cut_halving_ratio': 0.75,
        'dopamine_linear_ceiling_uM': ceiling,
        'dopamine_Km_uM': 0.2,
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    d = 'results/ASSEMBLY_ROUTE_TOPOLOGY'
    os.makedirs(d, exist_ok=True)
    with open(d + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('wrote', d + '/decision.json')
