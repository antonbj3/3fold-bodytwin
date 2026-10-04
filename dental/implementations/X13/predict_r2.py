from pathlib import Path
import numpy as np
import json, datetime, time, hashlib, resource
from calibrate_gap import feature, REGIONS
from predict_r1 import sha, dump, check
P = Path(__file__).resolve().parent

def anchor_key(d):
    return (d['internal_spacer_um'], d.get('layer_um', 0), d['arm'])

def run():
    check('R2')
    started = time.perf_counter()
    rows = json.loads((P / 'measurements.json').read_text())
    fr = json.loads((P / 'FROZEN_PREDICTIONS_R1.json').read_text())
    assert sha(P / 'measurements.json') == fr['source_measurements_sha256']
    lookup = {d['row_id']: d for d in fr['predictions']}
    out = []
    anchors = []
    for region in REGIONS:
        for study in sorted(set((d['study'] for d in rows if d['region'] == region))):
            groups = sorted([d for d in rows if d['region'] == region and d['study'] == study], key=anchor_key)
            a = groups[0]
            fm = next((d['model'] for d in fr['folds'] if d['region'] == region and d['held_family'] == a['source_family']))
            b = np.array(fm['beta'])
            C = np.array(fm['beta_cov'])
            xa = np.array(feature(a))
            mua = xa @ b
            vara = fm['sigma_um'] ** 2 + fm['tau_um'] ** 2 + xa @ C @ xa
            anchors.append(dict(study=study, region=region, row_id=a['row_id'], arm=a['arm'], observed_group_mean_um=a['measured_mean_um'], n_specimens=a['n_specimens'], has_evaluation=len(groups) > 1))
            for d in groups[1:]:
                x = np.array(feature(d))
                base = lookup[d['row_id']]
                cross = fm['tau_um'] ** 2 + x @ C @ xa
                var = fm['sigma_um'] ** 2 + fm['tau_um'] ** 2 + x @ C @ x
                mu = x @ b + cross / vara * (a['measured_mean_um'] - mua)
                condvar = var - cross * cross / vara
                assert condvar >= 0
                cov = np.array([[vara, cross], [cross, var]])
                ctrl = float(x @ b + cov[1, 0] * np.linalg.solve(cov[:1, :1], np.array([a['measured_mean_um'] - mua]))[0])
                out.append(dict(row_id=d['row_id'], study=study, source_family=d['source_family'], region=region, arm=d['arm'], prediction_um=float(mu), lower_um=float(mu - 1.96 * np.sqrt(condvar)), upper_um=float(mu + 1.96 * np.sqrt(condvar)), R1_prediction_um=base['prediction_um'], gaussian_control_um=ctrl, constant_anchor_control_um=a['measured_mean_um'], additive_control_um=float(x @ b + a['measured_mean_um'] - mua), anchor_id=a['row_id']))
    record = dict(round='R2', created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(P / 'PREREG_R2.json'), source_measurements_sha256=sha(P / 'measurements.json'), r1_prediction_sha256=sha(P / 'FROZEN_PREDICTIONS_R1.json'), predictions=out, anchors=anchors, seconds=time.perf_counter() - started, peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, prospective_physical_measurement=False)
    f = P / 'FROZEN_PREDICTIONS_R2.json'
    if f.exists():
        assert sha(f) == f.with_suffix('.sha256').read_text().strip()
        old = json.loads(f.read_text())
        assert max((abs(a['prediction_um'] - b['prediction_um']) for (a, b) in zip(old['predictions'], out))) < 1e-07
    else:
        dump(f, record)
        f.with_suffix('.sha256').write_text(sha(f) + '\n')
    print('Frozen R2:', len(out), 'evaluation rows;', len(anchors), 'regional anchor values')
if __name__ == '__main__':
    run()
