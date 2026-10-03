"""Can the power decision be made BEFORE the operation? The control it was compared against could.

The problem, found while looking for further quantities to strengthen. The chain predicts refraction
from `CASIA_POR_AQD post` -- the aqueous depth measured AFTER the lens is in the eye. The implant-power
decision was then scored against the power the surgeon actually implanted, and reported as closer to
the hindsight-correct power in 11 of 20 eyes against 2. But a surgeon choosing a lens cannot know where
that lens will sit afterwards. So the twin held information the control did not, and the comparison was
not between equally informed decisions. That does not make the result wrong; it makes it not a
prospective capability, which is what a decision has to be to be usable.

What this cell does. It replaces the postoperative depth with a PREDICTION from preoperative data only,
the way practice must, and reruns the same decision and the same scoring. The predictor is a
leave-one-out linear fit on quantities available before surgery -- preoperative anterior chamber depth,
axial length, lens thickness and corneal thickness -- so no eye contributes to its own prediction.

Three numbers come out of it, and all three matter: how well the postoperative depth can be predicted
at all, how much the recommendation moves when the prediction is used instead, and whether the decision
still beats the implanted power when both sides know only what was knowable beforehand. The third is
the one that says whether there is a usable capability here.
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

import numpy as np
import openpyxl

W = Path('.')
R15 = W / 'results/LANE_EYE_OPTICAL_TWIN/r15'
OUT = W / 'results/ASSEMBLY_PROSPECTIVE_POWER'
TARGET_D = 0.0
TOL_D = 0.25
GRID_D = 0.5
FEATURES = ['IOLM_ACD', 'IOLM_AL', 'IOLM_LT', 'IOLM_CCT']

sys.path.insert(0, str(R15))
sys.path.insert(0, str(R15.parent))
from clinical_anchor import predict                                          # noqa: E402


def choose_power(row: dict, aqd: float, lo: float = 5.0, hi: float = 40.0) -> float:
    """Bisection on the chain for the power that lands the eye on target, with a fixed lens position."""
    r = dict(row)
    r['CASIA_POR_AQD post'] = aqd
    f = lambda p: (r.update(IOLP=p), predict(r, 'thick')[0]['refraction_infinity_D'])[1] - TARGET_D  # noqa
    a, b = f(lo), f(hi)
    if a * b > 0:
        return float('nan')
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if f(lo) * f(m) <= 0:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def main() -> int:
    rows = list(openpyxl.load_workbook(next(R15.glob('*.xlsx')), data_only=True).active.values)
    h = [str(x) for x in rows[1]]
    need = FEATURES + ['IOLM_R1', 'IOLM_R2', 'IOLM_rpmean', 'IOLP', 'CASIA_POR_AQD post', 'sph', 'cyl']
    recs = []
    for n, raw in enumerate(rows[2:], 3):
        d = dict(zip(h, raw))
        if not all(isinstance(d.get(k), (int, float)) for k in need):
            continue
        d['sheet_row'] = n
        d['measured_refraction_D'] = d['sph'] + 0.5 * d['cyl']
        recs.append(d)

    X = np.array([[1.0] + [r[f] for f in FEATURES] for r in recs])
    y = np.array([r['CASIA_POR_AQD post'] for r in recs])
    n = len(recs)

    # Leave-one-out prediction of the postoperative lens position from preoperative data only.
    pred = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        beta, *_ = np.linalg.lstsq(X[m], y[m], rcond=None)
        pred[i] = X[i] @ beta
    resid = pred - y

    out = []
    for i, r in enumerate(recs):
        # sensitivity of the chain's refraction to the lens position, at this eye
        base = predict(r, 'thick')[0]['refraction_infinity_D']
        bumped = dict(r)
        bumped['CASIA_POR_AQD post'] = r['CASIA_POR_AQD post'] + 0.1
        dR_dAQD = (predict(bumped, 'thick')[0]['refraction_infinity_D'] - base) / 0.1

        p_post = choose_power(r, r['CASIA_POR_AQD post'])
        p_pre = choose_power(r, float(pred[i]))
        # the power that would in hindsight have been right, from the measured refraction
        bump_p = dict(r)
        bump_p['IOLP'] = r['IOLP'] + 1.0
        dR_dP = predict(bump_p, 'thick')[0]['refraction_infinity_D'] - base
        hindsight = r['IOLP'] - (r['measured_refraction_D'] - TARGET_D) / dR_dP
        grid = lambda p: round(p / GRID_D) * GRID_D                                        # noqa: E731
        out.append({
            'sheet_row': r['sheet_row'],
            'implanted_power_D': r['IOLP'],
            'aqd_measured_mm': round(r['CASIA_POR_AQD post'], 4),
            'aqd_predicted_mm': round(float(pred[i]), 4),
            'aqd_prediction_error_mm': round(float(resid[i]), 4),
            'dR_dAQD_D_per_mm': round(float(dR_dAQD), 4),
            'power_with_measured_position_D': round(grid(p_post), 2),
            'power_with_predicted_position_D': round(grid(p_pre), 2),
            'hindsight_correct_power_D': round(float(hindsight), 4),
        })

    def miss(key):
        return [abs(r[key] - r['hindsight_correct_power_D']) for r in out]

    m_post, m_pre = miss('power_with_measured_position_D'), miss('power_with_predicted_position_D')
    m_impl = miss('implanted_power_D')
    shift = [abs(r['power_with_predicted_position_D'] - r['power_with_measured_position_D'])
             for r in out]

    summary = {
        'eyes': n,
        'problem': ('the chain consumed the POSTOPERATIVE aqueous depth, which the surgeon it was '
                    'compared against could not know; the earlier 11-of-20 result is therefore not a '
                    'prospective capability'),
        'lens_position_prediction': {
            'features': FEATURES,
            'method': 'leave-one-out linear fit, no eye in its own fit',
            'mean_abs_error_mm': round(float(np.mean(np.abs(resid))), 4),
            'sd_of_error_mm': round(float(np.std(resid, ddof=1)), 4),
            'max_abs_error_mm': round(float(np.max(np.abs(resid))), 4),
            'spread_of_the_quantity_itself_mm': round(float(np.std(y, ddof=1)), 4),
        },
        'chain_sensitivity_to_lens_position': {
            'median_D_per_mm': round(float(np.median([r['dR_dAQD_D_per_mm'] for r in out])), 4),
            'implied_median_refraction_error_D': round(float(np.median(
                [abs(r['dR_dAQD_D_per_mm'] * r['aqd_prediction_error_mm']) for r in out])), 4),
        },
        'power_recommendation_shift_D': {
            'mean': round(statistics.mean(shift), 4), 'max': round(max(shift), 2),
            'eyes_where_the_grid_choice_changes': sum(1 for s in shift if s > 1e-9)},
        'scored_against_hindsight_correct_power': {
            'mean_miss_with_measured_position_D': round(statistics.mean(m_post), 4),
            'mean_miss_with_predicted_position_D': round(statistics.mean(m_pre), 4),
            'mean_miss_of_implanted_power_D': round(statistics.mean(m_impl), 4),
            'prospective_closer_than_implanted': sum(1 for a, b in zip(m_pre, m_impl) if a < b),
            'implanted_closer_than_prospective': sum(1 for a, b in zip(m_pre, m_impl) if b < a),
            'within_tolerance_prospective': sum(1 for x in m_pre if x <= TOL_D / abs(1.0)),
            'of': n},
        'claim_type': 'capability',
        'control': ('the power actually implanted, now compared against a twin that knows only what was '
                    'knowable before the operation, which is the equally informed comparison the earlier '
                    'result did not make'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'rows': out,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'PROSPECTIVE_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    lp, sc = summary['lens_position_prediction'], summary['scored_against_hindsight_correct_power']
    print(f"  eyes {n}")
    print(f"  lens position: predicted to {lp['mean_abs_error_mm']} mm MAE "
          f"(sd {lp['sd_of_error_mm']}, max {lp['max_abs_error_mm']}), "
          f"quantity's own spread {lp['spread_of_the_quantity_itself_mm']} mm")
    print(f"  chain sensitivity {summary['chain_sensitivity_to_lens_position']['median_D_per_mm']} D/mm "
          f"-> median {summary['chain_sensitivity_to_lens_position']['implied_median_refraction_error_D']} D")
    print(f"  recommendation shifts by {summary['power_recommendation_shift_D']['mean']} D mean, "
          f"max {summary['power_recommendation_shift_D']['max']}; grid choice changes in "
          f"{summary['power_recommendation_shift_D']['eyes_where_the_grid_choice_changes']} of {n}")
    print(f"  miss vs hindsight: measured-position {sc['mean_miss_with_measured_position_D']}, "
          f"prospective {sc['mean_miss_with_predicted_position_D']}, "
          f"implanted {sc['mean_miss_of_implanted_power_D']} D")
    print(f"  prospective closer than implanted in {sc['prospective_closer_than_implanted']} of {n}, "
          f"implanted closer in {sc['implanted_closer_than_prospective']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
