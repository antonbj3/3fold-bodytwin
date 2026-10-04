from pathlib import Path
import json, hashlib, datetime, time, resource
from calibrate_gap import fit, predict, gls_control, simple_controls, REGIONS
P = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, x):
    Path(p).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def check(tag):
    assert sha(P / f'PREREG_{tag}.json') == (P / f'PREREG_{tag}.sha256').read_text().strip()

def public_input(d):
    return {k: v for (k, v) in d.items() if k not in ['measured_mean_um', 'reported_sd_um', 'source_cell']}

def run():
    check('R1')
    check('R1_LINEAGE')
    started = time.perf_counter()
    rows = json.loads((P / 'measurements.json').read_text())
    out = []
    folds = []
    for region in REGIONS:
        subset = [d for d in rows if d['region'] == region]
        for family in sorted(set((d['source_family'] for d in subset))):
            train = [d for d in subset if d['source_family'] != family]
            test = [public_input(d) for d in subset if d['source_family'] == family]
            t = time.perf_counter()
            model = fit(train)
            fitsec = time.perf_counter() - t
            (mu, sd) = predict(model, test)
            t = time.perf_counter()
            gls = gls_control(train, model, test)
            ctrlsec = time.perf_counter() - t
            ctl = simple_controls(train, test)
            folds.append(dict(region=region, held_family=family, model=model, fit_seconds=fitsec, control_whitening_seconds=ctrlsec, train_families=sorted(set((d['source_family'] for d in train))), n_test=len(test)))
            for (i, d) in enumerate(test):
                nominal = d['marginal_spacer_um'] if region == 'marginal' else d['internal_spacer_um']
                out.append(dict(row_id=d['row_id'], study=d['study'], source_family=family, region=region, arm=d['arm'], internal_spacer_um=d['internal_spacer_um'], prediction_um=float(mu[i]), lower_um=float(mu[i] - 1.96 * sd[i]), upper_um=float(mu[i] + 1.96 * sd[i]), gls_um=float(gls[i]), controls_um={k: float(v[i]) for (k, v) in ctl.items()}, practice_um=nominal, internal_identity_um=d['internal_spacer_um'], nominal_definition='explicit_regional_spacer' if nominal is not None else 'unknown_marginal_nominal'))
    record = dict(round='R1', created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), source_measurements_sha256=sha(P / 'measurements.json'), prereg_sha256=sha(P / 'PREREG_R1.json'), lineage_prereg_sha256=sha(P / 'PREREG_R1_LINEAGE.json'), predictions=out, folds=folds, total_seconds=time.perf_counter() - started, peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, prospective_physical_measurement=False)
    f = P / 'FROZEN_PREDICTIONS_R1.json'
    if f.exists():
        assert sha(f) == f.with_suffix('.sha256').read_text().strip()
        old = json.loads(f.read_text())
        assert old['source_measurements_sha256'] == record['source_measurements_sha256']
        assert max((abs(a['prediction_um'] - b['prediction_um']) for (a, b) in zip(old['predictions'], out))) < 1e-07
        dump(P / 'R1_REPLAY_CHECK.json', dict(replay_agrees=True, seconds=record['total_seconds'], max_error_um=max((abs(a['prediction_um'] - b['prediction_um']) for (a, b) in zip(old['predictions'], out)))))
    else:
        dump(f, record)
        f.with_suffix('.sha256').write_text(sha(f) + '\n')
    dump(P / 'CURRENT_WORK_STATE.json', dict(lane='X13-cement-gap', phase='R1_PREDICTIONS_FROZEN', latest_gate='UNSCORED', next_operation='score frozen study-family predictions; preserve failure then change information', updated_utc=record['created_utc']))
    print('Frozen', len(out), 'predictions, seconds', record['total_seconds'])
if __name__ == '__main__':
    run()
