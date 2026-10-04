"""Does any decision we built tonight get WORSE when one of its inputs gets better?

READ THIS FIRST, 2026-10-04 afternoon: an proof_lane round on the cancellation question found a THIRD
mechanism that this screen cannot see, and it limits what a hit here means. On the 69 toric rows,
replacing only the stored intermediate predictor leaves EVERY total error exactly unchanged while
changing the measured stage-2 repair effect. So part of the anticorrelation is created by the
ADDITIVE DECOMPOSITION into stages, not by the system, and it survives arbitrarily large samples.
A hit in this screen is therefore a property of how the error was SPLIT as much as of the chain.

proof_lane also corrected the direction of my own reading, and the correction is the useful part:

  exact K-stage condition   dV = d^2 v_j - 2 d g_j,  with v_j = Sigma_jj and
                            g_j = (Sigma 1)_j = Cov(e_j, S)
  a finite repair harms iff g_j < d v_j / 2;  first-order harm needs g_j < 0
  interior optimum          sigma_j* = max(0, -sum_{k != j} rho_jk sigma_k)

For the toric pair that gives sigma_2* = 0.36272 against the current 0.59698, a reduction of
39.24 percent -- so stage 2 IMPROVES the system up to a 39 percent repair and only harms beyond it.
My earlier reading, that repairing stage 2 is 1.2594x worse, was the FULL-removal case stated as if
it were the repair case. Partial repair helps. I verified the 39.24 percent and the variance drop of
0.0507528 D^2 on halving against proof_lane's 0.0507480 independently.

And the metric matters: halving stage 2 reduces variance by 0.0507480 D^2 while INCREASING mean
squared error by 0.00570006 D^2, because component means do not vanish. For MSE the covariance has
to be replaced by E[e e^T].

Finally, the honest limit: no finite observational test universally separates physical compensation
from an identical constructed split. A conditional test exists and is specified in
results/PROOF_LANE_CANCELLATION/. So this screen stays useful as a WARNING that a margin depends on a
cancellation, and must not be read as evidence that the cancellation is physical.


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
import pathlib
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
    """Headline: the fraction of the spine redistribution a uniform load scaling explains.

    The first version of this arm was WRONG and its baseline said so: it scaled the model's own
    residual and computed 1 - sum|residual|/sum|change|, which is identically 0 at scale 1 because
    those are the same numbers. The headline being screened is the SURROGATE's error against the
    measured change, so the surrogate has to be rebuilt from the muscle data rather than imitated
    from a residual. A screen that reports a baseline of 0.0 for a decision whose headline is 0.702
    is screening the wrong quantity, and the baseline is the cheapest place to catch that.
    """
    import statistics
    base_dir = pathlib.Path('source_repository/data/msk_smoketest/subject2_spine_stoop_lift')
    NON = ('reserve', 'residual', 'box', 'FX', 'FY', 'FZ', 'MX', 'MY', 'MZ')

    def read(arm):
        lines = (base_dir / arm / 'walking1_StaticOptimization_force.sto').read_text().splitlines()
        h = next(i for i, l in enumerate(lines) if l.strip().lower() == 'endheader')
        cols = lines[h + 1].split('\t')
        rows = [[float(x) for x in l.split('\t')] for l in lines[h + 2:] if l.strip()]
        return {c: statistics.mean(abs(r[j]) for r in rows)
                for j, c in enumerate(cols)
                if j and not any(k.lower() in c.lower() for k in NON)}

    mn, mb = read('so_no_box'), read('so_with_box')
    shared = sorted(set(mn) & set(mb))
    inputs = {'surrogate_error_scale': 1.0}

    def explained(d):
        scale = sum(mb[m] for m in shared) / sum(mn[m] for m in shared)
        errs = [abs(mn[m] * scale - mb[m]) * d['surrogate_error_scale'] for m in shared]
        change = [abs(mb[m] - mn[m]) for m in shared]
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
