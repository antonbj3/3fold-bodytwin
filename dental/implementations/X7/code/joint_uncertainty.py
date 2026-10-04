"""Joint patient-level report-set calibration; frozen before its calibration computation."""
import json, datetime, hashlib, time, math
import numpy as np
from benchmark import P, write, sha
from parser import FIELDS

def preregister():
    r = dict(round='R4', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), capability='A single uncertainty statement covering all readable bite quantities and reports of one patient simultaneously', obstacle='Independent field sets passed macro coverage but jointly covered only0.4596 of198testpatients; per-field uncertainty is inadequate for whole-bite assessment', changed_operation='Calibration state is a patient vector; max nonconformity across fields and all reports, rather than9independent marginal scores', equation='s_i=max_{field,report}(1-p_field(y_report|x_i)); q=ceil((n_cal+1)*0.9)-th order statistic; use q in every field set', alpha=0.1, primary_gate=dict(joint_all_field_all_report_patient_coverage_min=0.85, macro_mean_set_size_max=2.0, macro_singleton_fraction_min=0.2, min_test_n=40, tolerance=1e-12, decision='All simultaneous criteria required; broad abstention is a useful failure, not diagnostic success'), strong_control='Standard multi-output split conformal with patient-maximum scores over identical labels/probabilities, independently sorted score quantile; expected TIE', practice='R2 marginal sets and forced point categories on same stored probabilities; joint coverage recomputed', falsifier='q=1/uninformative sets, mean size>2, singleton<.2, or joint coverage<.85. Actual injected q=-.01 must fail coverage.', external_referent=json.load(open(P / 'PREREG_R2.json'))['external_referent'], leaf_status=dict(patient_max='DERIVED_UNDER_ASSUMPTIONS', categorical_judgments='EXTERNALLY_MEASURED', probability_model='CONSTITUTIVE_CLOSURE', exchangeability='UNKNOWN', numeric_landmark_accuracy='UNKNOWN'), stopping_argument='Maximum score composes marginal inclusion into a vector event; it cannot improve class probabilities or identify missing teeth. More detailed model/physical claims need independent landmark observations.', cost=dict(preparation='Reuse complete R1/R2 upstream scans/models/reports, charged in COSTS', fit='One calibration order statistic, timed', discovery='UNKNOWN', validation='Read same held-out source reports, timed', queries='0 new physical queries', fallback='Independent per-tooth landmark acquisition UNKNOWN'), evaluation_scope='Sequential construction after R1/R2 test inspection; same fixed test population, no hyperparameter tuning; descriptive replay, not fresh independent confirmation', consumer='Whole-bite lab report with explicit review/abstention')
    p = P / 'PREREG_R4.json'
    if not p.exists():
        p.write_text(json.dumps(r, indent=2) + '\n')
        p.with_suffix('.json.sha256').write_text(sha(p) + '\n')

def run():
    preregister()
    r = json.load(open(P / 'PREREG_R4.json'))
    start = time.perf_counter()
    pred = json.load(open(P / 'raw/PREDICTIONS_R1.json'))
    yscal = json.load(open(P / 'raw/CALIBRATION_LABELS.json'))
    ysall = json.load(open(P / 'raw/ALL_STRUCTURED_REPORTS.json'))
    cal = [x for x in pred if x['split'] == 'calibration']
    test = [x for x in pred if x['split'] == 'test']
    scores = []
    for row in cal:
        ss = []
        for (f, v) in row['fields'].items():
            pp = dict(zip(v['classes'], v['probabilities']))
            yy = [z['labels'][f] for z in yscal[row['case_id']] if z['labels'][f] is not None]
            ss.extend((1 - pp.get(y, 0) for y in yy))
        if ss:
            scores.append(dict(case_id=row['case_id'], score=max(ss)))
    n = len(scores)
    k = math.ceil((n + 1) * 0.9)
    q = float(np.partition([v['score'] for v in scores], k - 1)[k - 1]) if k <= n else 1.0
    control = sorted((v['score'] for v in scores))[k - 1] if k <= n else 1.0
    assert abs(q - control) < 1e-12
    output = []
    for row in test:
        ff = {f: [y for (y, p) in zip(v['classes'], v['probabilities']) if 1 - p <= q + 1e-12] for (f, v) in row['fields'].items()}
        output.append(dict(case_id=row['case_id'], fields=ff))
    write('raw/PREDICTION_SETS_R4.json', output)
    frozen = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prediction_sha256=sha(P / 'raw/PREDICTION_SETS_R4.json'), prereg_sha256=sha(P / 'PREREG_R4.json'), phase='Joint calibration computed before this evaluator rereads test reference values; prior test labels already seen, descriptive sequential replay')
    if not (P / 'FROZEN_PREDICTIONS_R4.json').exists():
        write('FROZEN_PREDICTIONS_R4.json', frozen)
    else:
        assert json.load(open(P / 'FROZEN_PREDICTIONS_R4.json'))['prediction_sha256'] == frozen['prediction_sha256']
    cov = []
    sizes = []
    single = []
    rows = []
    badcov = []
    for (row, pr) in zip(output, test):
        ok = []
        bad = []
        for (f, S) in row['fields'].items():
            yy = {z['labels'][f] for z in ysall[row['case_id']] if z['labels'][f] is not None}
            if not yy:
                continue
            ok.append(yy.issubset(S))
            sizes.append(len(S))
            single.append(len(S) == 1)
            pp = pr['fields'][f]
            badS = [y for (y, p) in zip(pp['classes'], pp['probabilities']) if 1 - p <= -0.01]
            bad.append(yy.issubset(badS))
        if ok:
            cov.append(all(ok))
            badcov.append(all(bad))
            rows.append(dict(case_id=row['case_id'], joint_covered=all(ok), readable_fields=len(ok)))
    C = float(np.mean(cov))
    size = float(np.mean(sizes))
    sing = float(np.mean(single))
    bc = float(np.mean(badcov))
    assert bc < 0.85
    gate = bool(C >= 0.85 and size <= 2.0 and (sing >= 0.2) and (len(cov) >= 40))
    out = dict(primary_gate='PASS' if gate else 'FAIL', joint_patient_coverage=C, mean_field_set_size=size, singleton_fraction=sing, n_test_patients=len(cov), calibration_n=n, threshold_q=q, matched_control_threshold=control, method_gain='TIE', injected_bad_threshold=dict(q=-0.01, actual_joint_coverage=bc, gate='FAIL', refutes=True), rows=rows, external_referent=r['external_referent'], wall_seconds=time.perf_counter() - start, evaluation_scope=r['evaluation_scope'], meaning='Passing coverage with broad sets cannot establish an informative whole-bite diagnosis')
    write('raw/RESULTS_R4.json', out)
    (P / 'HANDOFF_R4.md').write_text(f"Joint bite uncertaintyR4: {out['primary_gate']}; coverage{C:.6f}, mean set size{size:.6f}, singleton{sing:.6f}, q{q:.6f}. Frozen criteria coverage>=.85,size<=2,singleton>=.2. Same-info joint conformalTIE. Clinical full-bite ability remainsOPEN. Next operation: independently measured per-tooth landmarks that improve class identifiability; avoid another conformal variant on unchanged probabilities. Same test patients seen in previous rounds, descriptive replay, no fresh validation.\\n")
    print('R4', out['primary_gate'], 'joint', C, 'size', size, 'singleton', sing, 'q', q)
if __name__ == '__main__':
    run()
