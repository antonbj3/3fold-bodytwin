"""Does any decision we built tonight get WORSE when one of its inputs gets better?

WHY. Three unrelated chains showed the same failure on 2026-10-04: a favourable number that holds
only while an internal cancellation is retained. The toric system margin is certified at
rho = -0.796 and degrades 1.65x if a stage is repaired; a swarm error budget understated its spread
by 1.90x to 2.74x through opposite-sign superposition; and the laser lane's sharp endpoint bounds
exist "only after spatial cancellation is retained", with the aggregate strictly understating in
14 412 of 15 413 instances. Standing gate A17 therefore requires a repair response on any favourable
aggregate. This file is that gate applied to the decisions outside the eye.

THE TEST, one at a time per input: scale the input toward its ideal (smaller error, or a value
closer to the no-effect point) and recompute the headline. A decision whose headline improves
monotonically as each input improves is cancellation-free on that input. A decision whose headline
WORSENS somewhere is cancellation-dependent, and its number may not be quoted as a bound that
survives improvement.

WHAT A PASS MEANS, and it is narrow: it means the headline is monotone in each input separately over
the swept range. It does not mean the decision is right, that the inputs are correct, or that joint
moves are safe. The toric margin is included as the POSITIVE CONTROL precisely so a reader can see
the screen firing on a case we already know fails. PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import math
import os


def monotone_report(name: str, headline, inputs: dict, recompute, better_is_lower: bool) -> dict:
    """Sweep each input toward its ideal and record whether the headline ever moves the wrong way."""
    base = headline
    rows = {}
    for key in inputs:
        worsened = []
        for frac in (0.9, 0.75, 0.5, 0.25, 0.1, 0.0):
            trial = dict(inputs)
            trial[key] = inputs[key] * frac
            try:
                val = recompute(trial)
            except Exception as exc:
                rows[key] = {'error': f'{type(exc).__name__}: {exc}'}
                break
            wrong_way = (val > base) if better_is_lower else (val < base)
            if wrong_way:
                worsened.append({'input_scaled_to': frac, 'headline': val,
                                 'factor_vs_base': val / base if base else None})
        else:
            rows[key] = {'worsens_when_improved': bool(worsened), 'cases': worsened[:3]}
    return {'decision': name, 'baseline_headline': base, 'better_is_lower': better_is_lower,
            'per_input': rows,
            'cancellation_dependent': any(r.get('worsens_when_improved') for r in rows.values())}


def toric_positive_control() -> dict:
    """Known failure: two stages at rho = -0.796. Repairing either stage raises the combined sd."""
    s1, s2, rho = 0.45562, 0.59698, -0.7961

    def combined(d):
        a, b = d['stage_posterior_sd'], d['stage_residual_sd']
        return math.sqrt(a * a + b * b + 2 * rho * a * b)

    base = combined({'stage_posterior_sd': s1, 'stage_residual_sd': s2})
    return monotone_report('toric_system_margin_POSITIVE_CONTROL', base,
                           {'stage_posterior_sd': s1, 'stage_residual_sd': s2},
                           combined, better_is_lower=True)


def renal_secretion() -> dict:
    """Headline: blocked secretion fraction. Inputs: the clearance ratio and the filtration terms."""
    CL = 30.0 * 1000 / 60
    inputs = {'clearance_ratio': 0.5867, 'gfr_shift': 129.8 - 129.0}

    def blocked(d):
        gfr_drug = 129.8 - d['gfr_shift']
        return 1.0 - (d['clearance_ratio'] * CL - gfr_drug) / (CL - 129.8)

    base = blocked(inputs)
    # A SMALLER clearance ratio means a bigger block, so "better" here is a larger headline.
    return monotone_report('renal_blocked_secretion_fraction', base, inputs, blocked,
                           better_is_lower=False)


def complement_margin() -> dict:
    """Headline: critical deposition efficiency. Smaller kcat or lifetime should RAISE the ceiling."""
    inputs = {'kcat_per_s': 1.78, 't_half_s': 90.0}

    def ceiling(d):
        tau = d['t_half_s'] / math.log(2.0)
        return 1.0 / (d['kcat_per_s'] * tau)

    base = ceiling(inputs)
    return monotone_report('complement_critical_deposition_efficiency', base, inputs, ceiling,
                           better_is_lower=False)


def load_redistribution() -> dict:
    """Headline: fraction of redistribution a uniform scaling explains, from the four measured rows."""
    rows = [(3.15, 3.0876624781292343), (4.29, 3.9896), (3.58, 3.5134), (4.88, 4.5172)]
    inputs = {'residual_scale': 1.0}

    def explained(d):
        errs = [abs((p - t) * d['residual_scale']) for t, p in rows]
        change = [abs(t - p) for t, p in rows]
        return 1.0 - sum(errs) / sum(change)

    base = explained(inputs)
    return monotone_report('load_redistribution_explained_fraction', base, inputs, explained,
                           better_is_lower=False)


def main() -> None:
    out = [toric_positive_control(), renal_secretion(), complement_margin(), load_redistribution()]
    for r in out:
        flag = 'CANCELLATION-DEPENDENT' if r['cancellation_dependent'] else 'monotone, no cancellation'
        print(f"{r['decision']:<48} {flag}")
        for k, v in r['per_input'].items():
            if v.get('worsens_when_improved'):
                worst = max(v['cases'], key=lambda c: abs(c['factor_vs_base'] - 1))
                print(f"    {k}: worsens to {worst['headline']:.5f} "
                      f"({worst['factor_vs_base']:.4f}x base) when scaled to {worst['input_scaled_to']}")
    d = 'results/ASSEMBLY_CANCELLATION_SCREEN'
    os.makedirs(d, exist_ok=True)
    with open(d + '/screen.json', 'w') as fh:
        json.dump({'screen': out, 'standing_gate': 'A17',
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'}, fh, indent=2)
    print('wrote', d + '/screen.json')


if __name__ == '__main__':
    main()
