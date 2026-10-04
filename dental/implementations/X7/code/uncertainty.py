import json, time, hashlib, math, datetime
from pathlib import Path
import numpy as np
from benchmark import P, labels, write, sha, state
from parser import FIELDS

def quantile(scores, alpha=0.1):
    k = math.ceil((len(scores) + 1) * (1 - alpha))
    return float(np.partition(np.asarray(scores), min(k, len(scores)) - 1)[min(k, len(scores)) - 1]) if k <= len(scores) else 1.0

def scores_for(rows, ys, field, all_reports):
    vals = []
    for r in rows:
        if field not in r['fields']:
            continue
        dd = r['fields'][field]
        prob = dict(zip(dd['classes'], dd['probabilities']))
        yy = [o['labels'][field] for o in ys[r['case_id']] if o['labels'][field] is not None]
        if not all_reports:
            yy = [ys[r['case_id']][0]['labels'][field]] if ys[r['case_id']] and ys[r['case_id']][0]['labels'][field] is not None else []
        if yy:
            vals.append(dict(case_id=r['case_id'], score=max((1 - prob.get(y, 0) for y in yy)), labels=sorted(set(yy))))
    return vals

def run():
    start = time.perf_counter()
    r2 = json.load(open(P / 'PREREG_R2.json'))
    fp = json.load(open(P / 'FROZEN_PREDICTIONS.json'))
    assert sha(P / fp['prediction_path']) == fp['prediction_sha256']
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    pred = json.load(open(P / fp['prediction_path']))
    cal = [r for r in pred if r['split'] == 'calibration']
    test = [r for r in pred if r['split'] == 'test']
    yscal = labels(m['calibration'])
    write('raw/CALIBRATION_LABELS.json', yscal)
    thresholds = {}
    calscores = {}
    for f in FIELDS:
        full = scores_for(cal, yscal, f, True)
        first = scores_for(cal, yscal, f, False)
        q = quantile([v['score'] for v in full], r2['alpha'])
        q1 = quantile([v['score'] for v in first], r2['alpha'])
        k = math.ceil((len(full) + 1) * (1 - r2['alpha']))
        ctrl = sorted((v['score'] for v in full))[k - 1] if k <= len(full) else 1.0
        assert abs(q - ctrl) < 1e-12
        thresholds[f] = dict(patient_max_q=q, first_report_q=q1, strong_control_q=ctrl, cal_n=len(full))
        calscores[f] = dict(patient_max=full, first_report=first)
    write('raw/CALIBRATION_SCORES_R2.json', calscores)
    write('raw/THRESHOLDS_R2.json', thresholds)
    sets = []
    for r in test:
        row = dict(case_id=r['case_id'], fields={})
        for (f, d) in r['fields'].items():
            scores = 1 - np.asarray(d['probabilities'])
            cats = d['classes']
            tt = thresholds[f]
            row['fields'][f] = {meth: [y for (y, s) in zip(cats, scores) if s <= tt[key] + 1e-12] for (meth, key) in [('candidate', 'patient_max_q'), ('first_reference_conformal', 'first_report_q'), ('strong_conformal_control', 'strong_control_q')]}
            row['fields'][f]['point'] = [d['point']]
        sets.append(row)
    write('raw/PREDICTION_SETS_R2.json', sets)
    frozen = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), set_sha256=sha(P / 'raw/PREDICTION_SETS_R2.json'), threshold_sha256=sha(P / 'raw/THRESHOLDS_R2.json'), prereg_sha256=sha(P / 'PREREG_R2.json'), phase='After calibration labels, before this R2 evaluator opens test labels. R1 test results may be known; R2 method/gates froze before R1 evaluation.')
    if not (P / 'FROZEN_PREDICTIONS_R2.json').exists():
        write('FROZEN_PREDICTIONS_R2.json', frozen)
    else:
        assert json.load(open(P / 'FROZEN_PREDICTIONS_R2.json'))['set_sha256'] == frozen['set_sha256']
    ystest = labels(m['test'])
    all_y = dict(yscal)
    all_y.update(ystest)
    train = json.load(open(P / 'raw/TRAIN_LABELS.json'))
    all_y.update(train)
    write('raw/ALL_STRUCTURED_REPORTS.json', all_y)
    out = {}
    raw = []
    for f in FIELDS:
        ii = [r for r in sets if f in r['fields'] and any((o['labels'][f] is not None for o in ystest[r['case_id']]))]
        summary = {}
        for meth in ['candidate', 'first_reference_conformal', 'strong_conformal_control', 'point']:
            cov = []
            sz = []
            reportcov = []
            discord = []
            singletonwrong = []
            for r in ii:
                yy = [o['labels'][f] for o in ystest[r['case_id']] if o['labels'][f] is not None]
                S = r['fields'][f][meth]
                cov.append(set(yy).issubset(S))
                sz.append(len(S))
                reportcov.append(np.mean([y in S for y in yy]))
                discord.append(len(set(yy)) > 1)
                if len(S) == 1:
                    singletonwrong.append(not set(yy).issubset(S))
                if meth == 'candidate':
                    raw.append(dict(case_id=r['case_id'], field=f, reference_set=sorted(set(yy)), predicted_set=S, all_reports_covered=set(yy).issubset(S)))
            n = len(cov)
            p = float(np.mean(cov)) if n else 0
            z = 1.96
            den = 1 + z * z / n if n else 1
            mid = (p + z * z / (2 * n)) / den if n else 0
            rad = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den if n else 0
            summary[meth] = dict(n=n, all_report_patient_coverage=p, coverage_wilson95=[mid - rad, mid + rad], mean_size=float(np.mean(sz)), singleton_fraction=float(np.mean(np.array(sz) == 1)), mean_within_patient_report_coverage=float(np.mean(reportcov)), discordant_patient_fraction=float(np.mean(discord)), singleton_all_report_error=float(np.mean(singletonwrong)) if singletonwrong else None)
        summary['matched_control_max_score_error'] = 0.0
        out[f] = summary
    macro = {k: float(np.mean([out[f]['candidate'][k] for f in FIELDS])) for k in ['all_report_patient_coverage', 'mean_size', 'singleton_fraction']}
    parser_pass = json.load(open(P / 'raw/PARSER_REVIEW_CORRECTED.json'))['summary']['ours']['pass_gate']
    gate = bool(macro['all_report_patient_coverage'] >= 0.85 and macro['mean_size'] <= 2 and (macro['singleton_fraction'] >= 0.2) and all((out[f]['candidate']['n'] >= 40 for f in FIELDS)) and parser_pass)
    injected_sets = []
    injcovered = []
    for row in test:
        for (f, d) in row['fields'].items():
            yy = [r['labels'][f] for r in ystest[row['case_id']] if r['labels'][f] is not None]
            if not yy:
                continue
            badS = [y for (y, p) in zip(d['classes'], d['probabilities']) if 1 - p <= -0.01]
            injected_sets.append(badS)
            injcovered.append(set(yy).issubset(badS))
    injcov = float(np.mean(injcovered))
    assert injcov < 0.85 and all((len(s) == 0 for s in injected_sets))
    write('raw/SET_COMPARISONS_R2.json', raw)
    write('raw/RESULTS_R2.json', dict(metrics=out, macro=macro, primary_gate='PASS' if gate else 'FAIL', method_gain='TIE_WITH_STANDARD_PATIENT_MAX_CONFORMAL', external_referent=r2['external_referent'], injected_bad_threshold=dict(value=-0.01, result='All prediction sets empty', coverage=injcov, n_scored_patient_fields=len(injcovered), n_empty_sets=sum((len(s) == 0 for s in injected_sets)), gate='FAIL', refutes=True), wall_seconds=time.perf_counter() - start))
    disagreement = {}
    for f in FIELDS:
        pairs = []
        for (c, rs) in all_y.items():
            yy = [r['labels'][f] for r in rs if r['labels'][f] is not None]
            if len(yy) >= 2:
                pairs.append(dict(case_id=c, first=yy[0], second=yy[1], n_readable=len(yy), all_values=yy))
        from sklearn.metrics import cohen_kappa_score
        a = [p['first'] for p in pairs]
        b = [p['second'] for p in pairs]
        disagreement[f] = dict(n_patients=len(pairs), first_second_agreement=float(np.mean(np.array(a) == np.array(b))) if pairs else None, kappa=float(cohen_kappa_score(a, b)) if pairs else None, all_report_discordance=float(np.mean([len(set(v['all_values'])) > 1 for v in pairs])) if pairs else None, pairs=pairs)
    write('raw/REPORT_DISAGREEMENT.json', disagreement)
    if not (P / 'raw/RESULTS_R2_INITIAL.json').exists():
        write('raw/RESULTS_R2_INITIAL.json', json.load(open(P / 'raw/RESULTS_R2.json')))
    (P / 'HANDOFF_R2.md').write_text('R2 primary ' + ('PASS' if gate else 'FAIL') + '; ' + str(macro) + '. Strong matched patient-maximum conformalTIE. Changes calibration observation from a chosen single reference to all readable reports per patient. Numeric landmarks remainUNKNOWN. Next construction: independent tooth landmark port and full vector-frame retraining with a fresh holdout; no retuned thresholds here.\n')
    state('R2_DECIDED', 'R2 ' + ('PASS' if gate else 'FAIL') + '; matched conformal TIE', 'R3 continuous vector frame; tooth identity remains missing')
    print('R2', gate, macro, flush=True)
if __name__ == '__main__':
    run()
