"""R4 changes the observation scope, preserving R1-R3 anterior-only endpoints."""
from dental_release.paths import expand as _release_expand
from common import *
from screening import partition, features, rates, threshold, roc_ci, noise_context
import sys, re, zipfile, time, pickle, resource
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline

def prereg():
    x = read(P / 'PREREG_R2.json')
    x.update(round='R4', frozen_utc=now(), findings=['any_open_bite'], changed_operation='Replace anterior-only vertical category by positive open-bite statements at any region, preserving anterior/lateral clause scope. Fit the binary any-region target using already frozen tooth-owned regional feature recipe. No R1-R3 label/gate is rewritten.', obstacle=_release_expand('@DENTAL_CASE_ID@ is anterior-negative yet lateral-positive: anterior category cannot represent the requested any-region open-bite finding.'), endpoint_definitions=dict(any_open_bite='Positive if anterior open category OR at least one explicit unnegated open-bite clause in any region. Negative if no positive open-bite clause and known non-open overbite category. Otherwise UNKNOWN. Near/possible/almost qualifiers retained as uncertainty context; no adjudicated clinical truth.'), test_exposure=_release_expand('R1-R3 scored and @DENTAL_CASE_ID@ source reviewed before R4 scope prereg. This is hypothesis-generating changed-scope replay, no independent validation.'), minimum_representation_extension='Region-indexed open-bite source polarity; local vertical/opposition features from each tooth, not one overbite summary.', strongest_equally_informed_control='Actually re-fit identical conventional HGB on same all-region labels/features; logistic same features/labels. Previous R2 anterior-only model rescored on new labels as scope ablation; this ablation is differently informed, not equally informed algorithm comparator.')
    freeze(P / 'PREREG_R4.json', x)
    freeze(P / 'DECOMPOSITION_R4.json', dict(idea='Any-region screening requires union of regional observations, not just anterior category.', mechanism='Local vertical relation -> regional report statement -> logical union -> binary score.', equation='Y_any=max(Y_anterior,Y_lateral_right,Y_lateral_left); one anterior category is insufficient.', operation=x['changed_operation'], representation='source clauses with region and polarity, linked to tooth-owned features; arch consumer aggregates.', assumptions='Known non-open vertical categories plus no positive open mention treated as report-negative; omission remains report-completeness closure, not anatomical absence.', leaves=[dict(leaf='Regional source statements', status='EXTERNALLY_MEASURED', stop='Clinician-authored report states region but no measured landmark distance; rater identity UNKNOWN.'), dict(leaf='Polarity/negation extraction', status='CONSTITUTIVE_CLOSURE', stop='Explicit scope rules; producer audit/injections. Blinded expert extraction review remains NOT_RUN.'), dict(leaf='Any-region union', status='DERIVED_UNDER_ASSUMPTIONS', stop='Boolean OR retains positive clauses; no physiological model needed.'), dict(leaf='Regional geometry', status='DERIVED_UNDER_ASSUMPTIONS', stop='Same X21 conditional meshes; physical pose, target FDI and named landmarks UNKNOWN.'), dict(leaf='Unreported posterior opening', status='UNKNOWN', stop='Adjudicated regional absence and repeated measurement required. Report absence does not prove physical absence.')]))
    state('R4_PREREG_FROZEN', _release_expand('Anterior scope failure documented @DENTAL_CASE_ID@'), 'Extract any-region source polarity; train and freeze before R4 scoring')

def label(text, prior):
    rr = []
    for sent in re.split('[.;!\\n]', text.lower()):
        for mm in re.finditer('open[ -]bite', sent):
            prefix = sent[max(0, mm.start() - 65):mm.start()]
            tail = sent[mm.end():]
            neg = bool(re.search('\\b(?:no|without|absence of|neither)\\s+(?:(?:a|the|anterior|lateral|posterior|right|left|bilateral|unilateral)\\s+){0,5}$', prefix) or re.match('\\s+(?:(?:is|was|are)\\s+)?(?:absent|not present|not observed|not seen)', tail))
            scope = 'lateral' if re.search('\\b(?:lateral|posterior)\\s*$', prefix) else 'anterior' if re.search('\\banterior\\s*$', prefix) else 'unspecified'
            rr.append(dict(positive=not neg, scope=scope, clause=sent.strip(), uncertain_qualifier=bool(re.search('\\b(?:almost|possible|seems|suspect|not completely|not.*certainty)\\b', sent))))
    value = 1 if prior['open_bite'] == 1 or any((r['positive'] for r in rr)) else 0 if prior['open_bite'] == 0 or (rr and all((not r['positive'] for r in rr))) else None
    return (value, rr)

def first(rs):
    return next((r['value'] for r in rs if r['value'] is not None), None)

def fit():
    st = time.perf_counter()
    m = partition()
    src = read(P / 'raw/LABELS.json')
    lab = {}
    zpath = read(X7 / 'raw/DATA_MANIFEST.json')['zip']
    pr = read(P / 'PREREG_R4.json')
    with zipfile.ZipFile(zpath) as z:
        for c in sum(m.values(), []):
            lab[c] = []
            for r in src[c]:
                b = z.read(r['member'])
                assert hashlib.sha256(b).hexdigest() == r['sha256']
                (v, clauses) = label(b.decode('utf8', errors='replace'), r['findings'])
                lab[c].append(dict(member=r['member'], sha256=r['sha256'], value=v, clauses=clauses, anterior_value=r['findings']['open_bite']))
    freeze(P / 'raw/LABELS_R4.json', lab)
    (feat, cols, loc) = features()
    ids = [c for c in m['train'] if first(lab[c]) is not None]
    y = np.array([first(lab[c]) for c in ids])
    X = np.array([[feat[c].get(k) for k in cols] for c in ids], float)
    qids = m['calibration'] + m['test']
    Xq = np.array([[feat[c].get(k) for k in cols] for c in qids], float)
    model = HistGradientBoostingClassifier(**pr['model']['HGB']).fit(X, y)
    p = model.predict_proba(Xq)[:, 1]
    control = HistGradientBoostingClassifier(**pr['model']['HGB']).fit(X, y).predict_proba(Xq)[:, 1]
    assert max(abs(p - control)) <= 1e-12
    logistic = make_pipeline(SimpleImputer(strategy='median', add_indicator=True, keep_empty_features=True), StandardScaler(), LogisticRegression(**pr['model']['logistic'])).fit(X, y).predict_proba(Xq)[:, 1]
    old = read(P / 'raw/SCORES_R2.json')['candidate']
    scores = dict(candidate=dict(zip(qids, p.tolist())), matched_control=dict(zip(qids, control.tolist())), logistic_control=dict(zip(qids, logistic.tolist())), anterior_scope_ablation={c: old[c]['open_bite'] for c in qids})
    freeze(P / 'raw/SCORES_R4.json', scores)
    mp = D / 'model_any_open_R4.pkl'
    mp.write_bytes(pickle.dumps(dict(model=model, columns=cols)))
    freeze(P / 'FROZEN_PREDICTIONS_R4.json', dict(frozen_utc=now(), round='R4', prediction_path=str(P / 'raw/SCORES_R4.json'), prediction_sha256=sha(P / 'raw/SCORES_R4.json'), labels_sha256=sha(P / 'raw/LABELS_R4.json'), prereg_sha256=sha(P / 'PREREG_R4.json'), code_sha256=sha(Path(__file__)), model_path=str(mp), model_sha256=sha(mp), regional_sources=loc, train_n=len(ids), train_pos=int(sum(y)), feature_count=len(cols), matched_probability_difference=float(max(abs(p - control))), wall_s=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, test_exposure='Previously viewed sources and outcomes; all-region source transformation after scope discovery is descriptive.'))

def freeze_threshold():
    m = partition()
    lab = read(P / 'raw/LABELS_R4.json')
    scores = read(P / 'raw/SCORES_R4.json')
    out = {}
    for (meth, ss) in scores.items():
        ids = [c for c in m['calibration'] if first(lab[c]) is not None]
        y = [first(lab[c]) for c in ids]
        s = [ss[c] for c in ids]
        out[meth] = {}
        for sp in [0.9, 0.95]:
            t = threshold(y, s, sp)
            out[meth][str(sp)] = dict(threshold=t, calibration_ids=ids, calibration=rates(y, np.asarray(s) >= t))
    freeze(P / 'raw/THRESHOLDS_R4.json', out)
    freeze(P / 'FROZEN_THRESHOLDS_R4.json', dict(frozen_utc=now(), threshold_sha256=sha(P / 'raw/THRESHOLDS_R4.json'), phase='calibration only before any-region test scoring'))

def evaluate():
    st = time.perf_counter()
    m = partition()
    lab = read(P / 'raw/LABELS_R4.json')
    sc = read(P / 'raw/SCORES_R4.json')
    th = read(P / 'raw/THRESHOLDS_R4.json')
    fr = read(P / 'FROZEN_PREDICTIONS_R4.json')
    assert sha(fr['prediction_path']) == fr['prediction_sha256']
    assert sha(P / 'raw/LABELS_R4.json') == fr['labels_sha256']
    ids = [c for c in m['test'] if first(lab[c]) is not None]
    y = np.array([first(lab[c]) for c in ids])
    pairs = []
    out = {}
    errors = []
    for c in ids:
        rr = [r for r in lab[c] if r['value'] is not None]
        if len(rr) >= 2:
            pairs.append(dict(case_id=c, report1=rr[0]['value'], report2=rr[1]['value'], members=[r['member'] for r in rr[:2]], hashes=[r['sha256'] for r in rr[:2]]))
    for (meth, ss) in sc.items():
        out[meth] = {}
        s = np.array([ss[c] for c in ids])
        for sp in [0.9, 0.95]:
            t = th[meth][str(sp)]['threshold']
            r = rates(y, s >= t)
            support = r['npos'] >= 20 and r['nneg'] >= 30
            specok = r['specificity'] >= sp - 1e-12
            sensok = r['sensitivity_wilson95'][0] >= 0.8
            gate = 'PASS' if support and specok and sensok else 'UNKNOWN_INSUFFICIENT_SUPPORT' if not support else 'FAIL'
            context = noise_context(pairs, ss, t, sp, 5400)
            if not specok and context['paired_comparison'] == 'REACHES_OBSERVED_REPORT_CONSISTENCY':
                context['paired_comparison'] = 'UNKNOWN_MODEL_TEST_SPECIFICITY_BELOW_TARGET'
            out[meth][str(sp)] = dict(target_specificity=sp, threshold=t, calibration=th[meth][str(sp)]['calibration'], fixed_threshold_test=r, empirical_test_roc=roc_ci(y, s, sp, 5400), primary_gate=gate, gates=dict(support=support, test_specificity=specok, sensitivity_lower_at_least_point8=sensok), report_noise_context=context)
            assert rates(y, np.ones(len(y), int))['specificity'] < sp
            if meth == 'candidate':
                for (c, yy, score) in zip(ids, y, s):
                    if int(score >= t) != yy:
                        errors.append(dict(case_id=c, finding='any_open_bite', target_specificity=sp, error='FN' if yy else 'FP', reference=int(yy), score=float(score), threshold=t, source_reports=lab[c]))
    changed = [dict(case_id=c, old_anterior=next((r['anterior_value'] for r in lab[c] if r['anterior_value'] is not None), None), new_any=first(lab[c]), source_reports=lab[c]) for c in ids if first(lab[c]) != next((r['anterior_value'] for r in lab[c] if r['anterior_value'] is not None), None)]
    rr = dict(round='R4', claim_type='information_link', external_referent=read(P / 'PREREG_R4.json')['external_referent'], finding='any_open_bite', methods=out, test_n=len(ids), test_rejected=len(m['test']) - len(ids), test_rejected_fraction=(len(m['test']) - len(ids)) / len(m['test']), scope_changed_test_cases=len(changed), anterior_reference_rounds_preserved=True, resolution='PER_ARCH', input_resolution='PER_TOOTH', timescale='SIMULTANEOUS', noise_floor_reached='UNKNOWN', inter_intra_noise_floor='UNKNOWN', evaluation_wall_s=time.perf_counter() - st, review_state='PENDING_INDEPENDENT_REVIEW', matched_control='Executed HGB parity max difference 0; supporting control only; no algorithm novelty.', injected_all_positive_rejected=True)
    write(P / 'raw/RESULTS_R4.json', rr)
    write(P / 'raw/ERRORS_R4.json', errors)
    write(P / 'raw/SCOPE_CHANGED_R4.json', changed)
    (P / 'HANDOFF_R4.md').write_text('R4 any-region open-bite screening; source scope changes and new error rows retained. Original anterior-only results unchanged. True inter/intra floor UNKNOWN, reused cohort descriptive only. Next: blinded region-specific assessment with named readers on new patients; repeat registered bite.\n')
    state('R4_DECIDED', 'Any-region scope scored; original anterior gates retained', 'Package per-finding results and region-specific acquisition contract')
    print(json.dumps({k: out['candidate'][k] for k in out['candidate']}, indent=2), flush=True)
if __name__ == '__main__':
    {'preregister': prereg, 'fit': fit, 'threshold': freeze_threshold, 'evaluate': evaluate}[sys.argv[1]]()
