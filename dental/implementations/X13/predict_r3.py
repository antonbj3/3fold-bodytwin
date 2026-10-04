from pathlib import Path
import numpy as np, json, time, datetime
from predict_r1 import check, sha, dump
P = Path(__file__).resolve().parent

def run():
    check('R3')
    started = time.perf_counter()
    rows = json.loads((P / 'measurements.json').read_text())
    pred = []
    models = []
    for study in ['PMC10246932', 'PMC10721348']:
        q = sorted([d for d in rows if d['study'] == study and d['region'] == 'marginal'], key=lambda d: d['internal_spacer_um'])
        anchors = q[:2]
        s = np.array([d['internal_spacer_um'] for d in anchors])
        y = np.array([d['measured_mean_um'] for d in anchors])
        A = np.stack([np.ones(2), 1 / s], axis=1)
        beta = np.linalg.solve(A, y)
        model = dict(study=study, b_um=float(beta[0]), a_um_squared=float(beta[1]), anchor_rows=[d['row_id'] for d in anchors], anchor_spacers_um=s.tolist(), anchor_means_um=y.tolist(), anchor_sd_um=[d['reported_sd_um'] for d in anchors], matrix_rank=int(np.linalg.matrix_rank(A)), interpolation_domain_um=s.tolist())
        models.append(model)
        for d in q[2:]:
            st = d['internal_spacer_um']
            v = np.array([1, 1 / st])
            mu = float(v @ beta)
            w = v @ np.linalg.inv(A)
            radius = float(np.abs(w) @ np.array(model['anchor_sd_um']))
            lam = (1 / st - 1 / s[0]) / (1 / s[1] - 1 / s[0])
            ctrl = float((1 - lam) * y[0] + lam * y[1])
            aff = float(y[0] + (st - s[0]) / (s[1] - s[0]) * (y[1] - y[0]))
            pred.append(dict(row_id=d['row_id'], study=study, source_family=study, region='marginal', spacer_um=st, prediction_um=mu, reciprocal_control_um=ctrl, affine_control_um=aff, constant_control_um=float(y[1]), explicit_nominal_um=d['marginal_spacer_um'], anchor_sd_sensitivity_radius_um=radius, extrapolation=True))
    record = dict(round='R3', created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(P / 'PREREG_R3.json'), source_measurements_sha256=sha(P / 'measurements.json'), models=models, predictions=pred, seconds=time.perf_counter() - started, prospective_physical_measurement=False)
    f = P / 'FROZEN_PREDICTIONS_R3.json'
    if f.exists():
        assert sha(f) == f.with_suffix('.sha256').read_text().strip()
        old = json.loads(f.read_text())
        assert max((abs(a['prediction_um'] - b['prediction_um']) for (a, b) in zip(old['predictions'], pred))) < 1e-07
    else:
        dump(f, record)
        f.with_suffix('.sha256').write_text(sha(f) + '\n')
    print('Frozen R3', len(pred), 'higher-dose predictions')
if __name__ == '__main__':
    run()
