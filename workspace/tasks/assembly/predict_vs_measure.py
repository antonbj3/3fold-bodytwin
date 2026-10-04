"""When is a quantity's own measurement noisier than predicting it from the others?

The thing that started this, and it is strange enough to chase. Feeding the chain a PREDICTED lens
position beat feeding it the position actually measured after surgery: 0.5345 D against 0.5690 D on 89
eyes. A guess outperformed a measurement of the same quantity. The easy explanation is that the
leave-one-out fit smooths measurement noise, and if that is what happened then it is not a curiosity at
all -- it is a general rule with reach across the whole twin.

The rule, stated so it can be wrong: for any quantity we measure, compare how much its own measurement
disagrees with itself against how well it can be predicted from the OTHER quantities. Where the
prediction is tighter than the measurement's own disagreement, measuring that quantity harder is the
wrong investment, and the chain should take the prediction. Where it is not, the measurement is carrying
information no regularity can supply and must be kept.

What makes this testable here rather than philosophical: for four biometric quantities we hold two
independent instruments on the same 89 eyes, so the measurement's disagreement with itself is observed
rather than assumed. The prediction is leave-one-out from the other three plus axial length, so no eye
is in its own fit, and the comparison is like for like -- both numbers are in the quantity's own unit.

The control is the quantity's own mean: a prediction that cannot beat predicting the cohort average is
not a prediction. That matters because three of these quantities are tightly distributed, and a model
can look accurate purely by returning something near the middle.
"""
from __future__ import annotations

import json
import statistics
from pathlib import Path

import numpy as np
import openpyxl

W = Path('')
R15 = W / 'results/LANE_EYE_OPTICAL_TWIN/r15'
OUT = W / 'results/ASSEMBLY_PREDICT_VS_MEASURE'

PAIRS = [('CCT', 'IOLM_CCT', 'CASIA_PR_CCT', 'um'),
         ('ACD', 'IOLM_ACD', 'CASIA_PR_ACD', 'mm'),
         ('LT', 'IOLM_LT', 'CASIA_PR_LT', 'mm'),
         ('WTW', 'IOLM_WTW', 'CASIA_PR_WTW', 'mm')]
EXTRA = ['IOLM_AL', 'IOLM_K1', 'IOLM_K2']


def main() -> int:
    rows = list(openpyxl.load_workbook(next(R15.glob('*.xlsx')), data_only=True).active.values)
    h = [str(x) for x in rows[1]]
    need = [c for _, a, b, _ in PAIRS for c in (a, b)] + EXTRA
    recs = [dict(zip(h, r)) for r in rows[2:]]
    recs = [r for r in recs if all(isinstance(r.get(k), (int, float)) for k in need)]
    n = len(recs)

    results = []
    for name, a, b, unit in PAIRS:
        # what the quantity's own measurement does: two instruments on the same eye
        disagreement = [abs(r[a] - r[b]) for r in recs]
        target = np.array([r[a] for r in recs])

        # predict it from the OTHER quantities, leave-one-out
        feat_cols = [c for nm, ca, cb, _ in PAIRS if nm != name for c in (ca,)] + EXTRA
        X = np.array([[1.0] + [r[c] for c in feat_cols] for r in recs])
        pred = np.empty(n)
        for i in range(n):
            m = np.ones(n, bool)
            m[i] = False
            beta, *_ = np.linalg.lstsq(X[m], target[m], rcond=None)
            pred[i] = X[i] @ beta
        pred_err = np.abs(pred - target)

        # control: predict the cohort mean, also leave-one-out
        mean_err = np.array([abs(target[i] - target[np.arange(n) != i].mean()) for i in range(n)])

        results.append({
            'quantity': name, 'unit': unit, 'eyes': n,
            'measurement_disagrees_with_itself_mean': round(statistics.mean(disagreement), 5),
            'prediction_error_from_other_quantities_mean': round(float(pred_err.mean()), 5),
            'cohort_mean_control_error': round(float(mean_err.mean()), 5),
            'prediction_beats_the_cohort_mean': bool(pred_err.mean() < mean_err.mean()),
            'prediction_tighter_than_the_measurement': bool(pred_err.mean()
                                                            < statistics.mean(disagreement)),
            'ratio_prediction_over_disagreement': round(
                float(pred_err.mean() / statistics.mean(disagreement)), 4),
        })

    summary = {
        'question': ('for which quantities is the prediction from other quantities tighter than the '
                     'quantity own measurement disagreement'),
        'why': ('a predicted lens position beat the measured one in the decision, 0.5345 against 0.5690 '
                'D, and if that is noise smoothing then it generalises: measuring a quantity harder is '
                'the wrong investment wherever its regularity is tighter than its instrument'),
        'measurement_disagreement_is_observed': ('two independent instruments on the same 89 eyes, not '
                                                 'an assumed repeatability'),
        'control': ('the cohort mean, leave-one-out. A prediction that cannot beat the average is not a '
                    'prediction, and three of these quantities are tightly distributed enough that a '
                    'model can look accurate by returning the middle'),
        'results': results,
        'claim_type': 'capability',
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'PREDICT_VS_MEASURE_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    print(f"  {n} eyes with all quantities\n")
    print(f"  {'quantity':6s} {'unit':>5s} {'measurement disagreement':>14s} {'prediction':>11s} "
          f"{'cohort mean':>12s}  tighter?")
    for r in results:
        print(f"  {r['quantity']:6s} {r['unit']:>5s} "
              f"{r['measurement_disagrees_with_itself_mean']:>14.5f} "
              f"{r['prediction_error_from_other_quantities_mean']:>11.5f} "
              f"{r['cohort_mean_control_error']:>12.5f}  "
              f"{'YES' if r['prediction_tighter_than_the_measurement'] else 'no'}"
              f"{'' if r['prediction_beats_the_cohort_mean'] else '  (does not beat cohort mean)'}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
