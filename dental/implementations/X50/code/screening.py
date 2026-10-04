from common import *
import sys, time, resource, re, zipfile, pickle, math
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

def partition():
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    part = {k: m[k] for k in ['train', 'calibration', 'test']}
    assert all((not set(part[a]) & set(part[b]) for (a, b) in [('train', 'test'), ('train', 'calibration'), ('test', 'calibration')]))
    assert sum(map(len, part.values())) == len(set(sum(part.values(), [])))
    return part

def input_manifest():
    paths = [X7 / 'raw/DATA_MANIFEST.json', X7 / 'raw/geometry.jsonl', X7 / 'raw/PREDICTIONS_R1.json', X7 / 'raw/STRONG_CONTROL_R1.json', X21 / 'raw/CORRECTED_REPORTS_R5.json', X21 / 'raw/OPPOSITION_FEATURES_R2.json', P / 'code/screening.py', P / 'code/common.py']
    return [dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)) for p in paths]

def geometry():
    return {x['case_id']: x for x in map(json.loads, (X7 / 'raw/geometry.jsonl').read_text().splitlines())}

def positive_prob(field, prefix):
    if not field:
        return None
    return float(sum((v for (k, v) in zip(field['classes'], field['probabilities']) if k.startswith(prefix))))

def freeze_r1():
    st = time.perf_counter()
    m = partition()
    g = geometry()
    old = {r['case_id']: r for r in read(X7 / 'raw/PREDICTIONS_R1.json')}
    out = {}
    for c in m['calibration'] + m['test']:
        r = old[c]['fields']
        f = g[c]['features']
        out[c] = {'crossbite': positive_prob(r.get('crossbite'), 'present'), 'open_bite': positive_prob(r.get('overbite'), 'open'), 'deep_bite': positive_prob(r.get('overbite'), 'increased'), 'angle_II': max(positive_prob(r.get('molar_right'), 'II ') or 0.0, positive_prob(r.get('molar_left'), 'II ') or 0.0), 'angle_III': max(positive_prob(r.get('molar_right'), 'III') or 0.0, positive_prob(r.get('molar_left'), 'III') or 0.0), 'scissor_bite': max(f.get('scissor_extent_right_deg') or 0.0, f.get('scissor_extent_left_deg') or 0.0)}
    pp = P / 'raw/SCORES_R1.json'
    freeze(pp, dict(candidate=out, matched_control=out))
    control = read(X7 / 'raw/STRONG_CONTROL_R1.json')
    assert all((r['max_probability_difference'] <= 1e-12 and r['injected_corrupt_control_fails'] for r in control['rows']))
    freeze(P / 'FROZEN_PREDICTIONS.json', dict(round='R1', frozen_utc=now(), prediction_path=str(pp), prediction_sha256=sha(pp), prereg_sha256=sha(P / 'PREREG_R1.json'), inputs=input_manifest(), source_control='Existing actually executed independent HGB refit; hashes checked, zero probability differences. Scissor same standard extent ranking.', legacy_test_exposure='Previously evaluated cohort; no untouched validation claim.', wall_s=time.perf_counter() - st))
    state('R1_PREDICTIONS_FROZEN', 'Fixed parent scores frozen', 'Extract source labels; choose thresholds on calibration only')

def scissor(text):
    clauses = re.split('[.;!\\n]', text.lower())
    vals = []
    evidence = []
    for s in clauses:
        if re.search('\\b(?:scissor(?:s)?[ -]?(?:bite)?|brodie)\\b', s):
            neg = bool(re.search('\\b(?:no|without|absence of)\\s+(?:(?:a|the|lateral|posterior|anterior|bilateral|unilateral|right|left)\\s+){0,5}(?:scissor|brodie)', s) or re.search('(?:scissor\\w*(?:[ -]bite)?|brodie)[^.]{0,30}(?:absent|not present|not observed)', s))
            vals.append(int(not neg))
            evidence.append(s.strip())
    if vals:
        return (max(vals), 'explicit_scissor_clause', evidence)
    if re.search('correct (?:transverse(?: and vertical)?|transversal)|(?:transverse|transversal)(?: skeletal)? relationships? (?:are |is )?(?:correct|normal)|normal transverse', text, re.I):
        return (0, 'explicit_normal_transverse', [])
    return (None, 'not_explicitly_assessed', [])

def binary(labels, text):
    ob = labels.get('overbite')
    cross = labels.get('crossbite')
    mol = [labels.get('molar_right'), labels.get('molar_left')]
    valid = lambda v: v is not None and v != 'not assessable' and v.startswith(('I', 'II', 'III'))
    angle = lambda prefix: 1 if any((valid(v) and (v.startswith(prefix + ' ') or v == prefix) for v in mol)) else 0 if all((valid(v) for v in mol)) else None
    (sb, why, clauses) = scissor(text)
    return (dict(crossbite=1 if cross == 'present' else 0 if cross == 'absent' else None, open_bite=1 if ob == 'open' else 0 if ob in ['normal', 'reduced', 'increased'] else None, deep_bite=1 if ob == 'increased' else 0 if ob in ['normal', 'reduced', 'open'] else None, angle_II=angle('II'), angle_III=angle('III'), scissor_bite=sb), dict(scissor_basis=why, scissor_clauses=clauses, lateral_open_mentioned=bool(re.search('\\blateral open[ -]bite', text, re.I)), original_categories={k: labels.get(k) for k in ['overbite', 'crossbite', 'molar_right', 'molar_left']}))

def extract_labels():
    if (P / 'raw/LABELS.json').exists():
        return
    st = time.perf_counter()
    src = read(X21 / 'raw/CORRECTED_REPORTS_R5.json')
    out = {}
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    with zipfile.ZipFile(m['zip']) as z:
        for c in sum(partition().values(), []):
            out[c] = []
            for r in src.get(c, []):
                b = z.read(r['member'])
                assert hashlib.sha256(b).hexdigest() == r['sha256']
                (y, scope) = binary(r['labels'], b.decode('utf8', errors='replace'))
                out[c].append(dict(member=r['member'], sha256=r['sha256'], findings=y, scope=scope))
    freeze(P / 'raw/LABELS.json', out)
    write(P / 'raw/LABEL_COST.json', dict(wall_s=time.perf_counter() - st, cases=len(out), reports=sum(map(len, out.values())), source_unchanged=True))

def first_label(rs, f):
    return next((r['findings'][f] for r in rs if r['findings'][f] is not None), None)

def consensus_label(rs, f):
    vals = {r['findings'][f] for r in rs if r['findings'][f] is not None}
    return next(iter(vals)) if len(vals) == 1 else None

def threshold(y, s, target):
    neg = np.sort(np.asarray(s, float)[np.asarray(y) == 0])
    n = len(neg)
    if not n:
        return None
    allowed = int(math.floor((1 - target) * n + 1e-10))
    return float(np.nextafter(neg[-allowed - 1], np.inf))

def wilson(k, n):
    if not n:
        return None
    z = 1.959963984540054
    p = k / n
    d = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / d
    rad = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [max(0.0, mid - rad), min(1.0, mid + rad)]

def rates(y, p):
    y = np.asarray(y, int)
    p = np.asarray(p, int)
    npos = int(sum(y == 1))
    nneg = int(sum(y == 0))
    tp = int(sum((y == 1) & (p == 1)))
    tn = int(sum((y == 0) & (p == 0)))
    return dict(n=len(y), npos=npos, nneg=nneg, tp=tp, fn=npos - tp, tn=tn, fp=nneg - tn, sensitivity=tp / npos if npos else None, sensitivity_wilson95=wilson(tp, npos), specificity=tn / nneg if nneg else None, specificity_wilson95=wilson(tn, nneg), resolution='PER_ARCH')

def calfreeze(round):
    extract_labels()
    m = partition()
    labels = read(P / 'raw/LABELS.json')
    scores = read(P / f'raw/SCORES_{round}.json')
    out = {}
    selector = consensus_label if round == 'R3' else first_label
    for (meth, ss) in scores.items():
        if meth == 'metadata':
            continue
        out[meth] = {}
        for f in FINDINGS:
            ids = [c for c in m['calibration'] if selector(labels[c], f) is not None and ss[c].get(f) is not None]
            y = [selector(labels[c], f) for c in ids]
            s = [ss[c][f] for c in ids]
            out[meth][f] = {str(sp): dict(threshold=threshold(y, s, sp), calibration=rates(y, np.asarray(s) >= (threshold(y, s, sp) if threshold(y, s, sp) is not None else float('inf'))), calibration_ids=ids) for sp in [0.9, 0.95]}
    p = P / f'raw/THRESHOLDS_{round}.json'
    freeze(p, out)
    freeze(P / f'FROZEN_THRESHOLDS_{round}.json', dict(frozen_utc=now(), threshold_sha256=sha(p), scores_sha256=sha(P / f'raw/SCORES_{round}.json'), labels_sha256=sha(P / 'raw/LABELS.json'), phase='Calibration only; test label data already known in predecessors; no this-round test threshold tuning'))

def features():
    m = partition()
    g = geometry()
    op = read(X21 / 'raw/OPPOSITION_FEATURES_R2.json')
    out = {}
    loc = []
    for c in sum(m.values(), []):
        p = X21 / 'raw/cases' / (c + '.json')
        r = read(p)
        assert not r.get('error')
        out[c] = {'arch__' + k: v for (k, v) in g[c]['features'].items()}
        out[c].update({'pair__' + k: v for (k, v) in r['features'].items()})
        out[c].update({'opposition__' + k: v for (k, v) in op[c]['features'].items()})
        loc.append(dict(case_id=c, path=str(p), sha256=sha(p), source_zip_members=r['source_zip_members'], map_path=r['map_path'], map_sha256=r['map_sha256']))
    cols = sorted(set.union(*(set(out[c]) for c in m['train'])))
    return (out, cols, loc)

def fit_freeze(round):
    st = time.perf_counter()
    pr = read(P / f'PREREG_{round}.json')
    m = partition()
    labels = read(P / 'raw/LABELS.json')
    (feat, cols, loc) = features()
    query = m['calibration'] + m['test']
    Xq = np.array([[feat[c].get(k) for k in cols] for c in query], float)
    out = {a: {c: {} for c in query} for a in ['candidate', 'logistic_control', 'matched_control']}
    models = {}
    cost = []
    selector = consensus_label if round == 'R3' else first_label
    for f in FINDINGS:
        ids = [c for c in m['train'] if selector(labels[c], f) is not None]
        y = np.array([selector(labels[c], f) for c in ids], int)
        X = np.array([[feat[c].get(k) for k in cols] for c in ids], float)
        s = time.perf_counter()
        if len(set(y)) < 2:
            for meth in out:
                for c in query:
                    out[meth][c][f] = None
            cost.append(dict(finding=f, n_train=len(ids), train_pos=int(sum(y)), fit_status='UNKNOWN_LESS_THAN_TWO_CLASSES', wall_s=time.perf_counter() - s))
            continue
        recipe = pr['model']['HGB']
        mdl = HistGradientBoostingClassifier(**recipe).fit(X, y)
        p = mdl.predict_proba(Xq)[:, 1]
        same = HistGradientBoostingClassifier(**recipe).fit(X, y).predict_proba(Xq)[:, 1]
        assert np.max(abs(p - same)) <= 1e-12
        im = SimpleImputer(strategy='median', add_indicator=True, keep_empty_features=True)
        lc = make_pipeline(im, StandardScaler(), LogisticRegression(**pr['model']['logistic'])).fit(X, y)
        lp = lc.predict_proba(Xq)[:, 1]
        for (c, a, b, d) in zip(query, p, lp, same):
            out['candidate'][c][f] = float(a)
            out['logistic_control'][c][f] = float(b)
            out['matched_control'][c][f] = float(d)
        models[f] = dict(candidate=mdl, logistic_control=lc)
        cost.append(dict(finding=f, n_train=len(ids), train_pos=int(sum(y)), fit_status='RUN', wall_s=time.perf_counter() - s, matched_control_max_probability_difference=float(max(abs(p - same))), dropped_discordant_or_unreadable_train=len(m['train']) - len(ids)))
        print(round, f, 'train', len(ids), 'pos', int(sum(y)), 's', format(cost[-1]['wall_s'], '.2f'), flush=True)
    D.mkdir(parents=True, exist_ok=True)
    mp = D / f'models_{round}.pkl'
    mp.write_bytes(pickle.dumps(dict(columns=cols, models=models)))
    pp = P / f'raw/SCORES_{round}.json'
    freeze(pp, out)
    freeze(P / f'FROZEN_PREDICTIONS_{round}.json', dict(round=round, frozen_utc=now(), prediction_path=str(pp), prediction_sha256=sha(pp), prereg_sha256=sha(P / f'PREREG_{round}.json'), model_path=str(mp), model_sha256=sha(mp), inputs=input_manifest(), regional_sources=loc, feature_count=len(cols), fit_cost=cost, wall_s=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, legacy_test_exposure='Previously scored cohort, changed construction fixed before this-round results; descriptive replay only.'))
    state(round + '_PREDICTIONS_FROZEN', 'Models fit train only, matched conventional HGB re-fit passed', 'Freeze calibration thresholds then compare unchanged external test findings')

def roc_ci(y, s, sp, seed):
    y = np.asarray(y, int)
    s = np.asarray(s, float)
    pos = s[y == 1]
    neg = s[y == 0]
    if not len(pos) or not len(neg):
        return dict(status='UNKNOWN_MISSING_CLASS')
    t = threshold(y, s, sp)
    r = rates(y, s >= t)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(2000):
        a = rng.choice(pos, len(pos), replace=True)
        b = rng.choice(neg, len(neg), replace=True)
        yy = np.r_[np.ones(len(a), int), np.zeros(len(b), int)]
        tt = threshold(yy, np.r_[a, b], sp)
        vals.append(float(np.mean(a >= tt)))
    return dict(sensitivity=r['sensitivity'], specificity=r['specificity'], threshold=t, bootstrap95=np.quantile(vals, [0.025, 0.975]).tolist(), draws=2000, discrete_nonrandomized=True, status='DESCRIPTIVE_TEST_ROC_NOT_FROZEN_DEPLOYABLE_THRESHOLD')

def pairs_for(labels, ids, f):
    out = []
    for c in ids:
        rr = [r for r in labels[c] if r['findings'][f] is not None]
        if len(rr) >= 2:
            out.append(dict(case_id=c, report1=rr[0]['findings'][f], report2=rr[1]['findings'][f], members=[r['member'] for r in rr[:2]], hashes=[r['sha256'] for r in rr[:2]]))
    return out

def noise_context(pairs, scores, t, sp, seed):
    y = np.array([r['report1'] for r in pairs])
    p2 = np.array([r['report2'] for r in pairs])
    pr = np.array([scores[r['case_id']] >= t for r in pairs])
    rr = rates(y, p2)
    rr['conflicts'] = int(sum(y != p2))
    rr['conflict_fraction'] = rr['conflicts'] / len(y) if len(y) else None
    rr['conflict_wilson95'] = wilson(rr['conflicts'], len(y))
    rr['source'] = 'First two readable report files, author/time identity unknown'
    pos = y == 1
    gap = None
    ci = None
    if sum(pos):
        d = pr[pos].astype(float) - p2[pos].astype(float)
        gap = float(np.mean(d))
        rng = np.random.default_rng(seed)
        boot = [float(np.mean(rng.choice(d, len(d), replace=True))) for _ in range(2000)]
        ci = np.quantile(boot, [0.025, 0.975]).tolist()
    comparable = rr['npos'] >= 20 and rr['nneg'] >= 30 and (rr['specificity'] >= sp - 1e-12)
    comparison = 'UNKNOWN_INSUFFICIENT_PAIRED_SUPPORT_OR_HUMAN_SPECIFICITY' if not comparable else 'BELOW_OBSERVED_REPORT_CONSISTENCY' if ci[1] < -0.05 else 'REACHES_OBSERVED_REPORT_CONSISTENCY' if ci[0] >= -0.05 else 'INDETERMINATE'
    return dict(report_pair=rr, paired_geometry=rates(y, pr), paired_sensitivity_gap=gap, paired_gap_bootstrap95=ci, paired_comparison=comparison, inter_intra_noise_floor='UNKNOWN_RATER_REPEAT_IDENTITIES_AND_ADJUDICATED_TRUTH_UNAVAILABLE', noise_floor_reached='UNKNOWN', not_a_universal_error_ceiling=True)

def evaluate(round):
    st = time.perf_counter()
    m = partition()
    lab = read(P / 'raw/LABELS.json')
    sc = read(P / f'raw/SCORES_{round}.json')
    th = read(P / f'raw/THRESHOLDS_{round}.json')
    frozen = read(P / ('FROZEN_PREDICTIONS.json' if round == 'R1' else f'FROZEN_PREDICTIONS_{round}.json'))
    assert sha(frozen['prediction_path']) == frozen['prediction_sha256']
    assert sha(P / f'PREREG_{round}.json') == frozen['prereg_sha256']
    out = {}
    errors = []
    locators = []
    for (fi, f) in enumerate(FINDINGS):
        valid_label = [c for c in m['test'] if first_label(lab[c], f) is not None]
        ids = [c for c in valid_label if sc['candidate'][c].get(f) is not None]
        y = np.array([first_label(lab[c], f) for c in ids], int)
        pairs = pairs_for(lab, ids, f)
        methods = {}
        for (meth, ss) in sc.items():
            res = {}
            if not ids or any((ss[c].get(f) is None for c in ids)):
                methods[meth] = dict(status='UNKNOWN_MODEL_OR_LABEL_MISSING')
                continue
            s = np.array([ss[c][f] for c in ids], float)
            for sp in [0.9, 0.95]:
                tc = th[meth][f][str(sp)]
                t = tc['threshold']
                if t is None:
                    res[str(sp)] = dict(status='UNKNOWN_NO_CALIBRATION_NEGATIVES')
                    continue
                fixed = rates(y, s >= t)
                emp = roc_ci(y, s, sp, 5000 + fi)
                support = fixed['npos'] >= 20 and fixed['nneg'] >= 30
                specok = fixed['specificity'] is not None and fixed['specificity'] >= sp - 1e-12
                sensok = fixed['sensitivity_wilson95'] is not None and fixed['sensitivity_wilson95'][0] >= 0.8
                gate = 'PASS' if support and specok and sensok else 'UNKNOWN_INSUFFICIENT_SUPPORT' if not support else 'FAIL'
                context = noise_context(pairs, {c: ss[c][f] for c in ids}, t, sp, 5100 + fi)
                if not specok and context['paired_comparison'] == 'REACHES_OBSERVED_REPORT_CONSISTENCY':
                    context['paired_comparison'] = 'UNKNOWN_MODEL_TEST_SPECIFICITY_BELOW_TARGET'
                res[str(sp)] = dict(target_specificity=sp, threshold=t, calibration=tc['calibration'], fixed_threshold_test=fixed, empirical_test_roc=emp, primary_gate=gate, gates=dict(support=support, test_specificity=specok, sensitivity_lower_at_least_point8=sensok), report_noise_context=context, no_information_context=dict(deterministic_all_negative=rates(y, np.zeros(len(y), int)), random_ranking_expected_sensitivity=1 - sp, clinical_practice_without_geometry='NOT_MEASURED'))
                if meth == 'candidate':
                    for (c, yy, score) in zip(ids, y, s):
                        if yy == 1 and score < t or (yy == 0 and score >= t):
                            rr = [r for r in lab[c] if r['findings'][f] is not None]
                            alts = {r['findings'][f] for r in rr}
                            why = 'Report conflict: other report supports prediction' if int(score >= t) in alts else 'No report corroboration; landmark/pose/model/text cause UNKNOWN'
                            errors.append(dict(round=round, case_id=c, finding=f, target_specificity=sp, error='FN' if yy else 'FP', reference=int(yy), score=float(score), threshold=t, report_values=[r['findings'][f] for r in rr], members=[r['member'] for r in rr], hashes=[r['sha256'] for r in rr], explanation=why, lateral_open_mentioned=any((r['scope']['lateral_open_mentioned'] for r in rr)), reference_categories=rr[0]['scope']['original_categories']))
            methods[meth] = res
        paired_all = pairs_for(lab, sum(m.values(), []), f)
        pa = rates([r['report1'] for r in paired_all], [r['report2'] for r in paired_all])
        pa.update(conflicts=sum((r['report1'] != r['report2'] for r in paired_all)), kind='DEPENDENT_WHOLE_COHORT_CONTEXT_NOT_POOLED_WITH_TEST')
        out[f] = dict(methods=methods, labels_test_accepted=len(valid_label), test_scored=len(ids), test_rejected=len(m['test']) - len(ids), test_rejected_fraction=(len(m['test']) - len(ids)) / len(m['test']), reasons=dict(unreadable_or_scope_unknown=len(m['test']) - len(valid_label), score_missing=len(valid_label) - len(ids)), whole_cohort_report_pair=pa, resolution='PER_ARCH', input_resolution='PER_TOOTH' if round != 'R1' else 'PER_ARCH', noise_floor_reached='UNKNOWN', inter_intra_disagreement='NOT_IDENTIFIED')
        for c in ids:
            rr = next((r for r in lab[c] if r['findings'][f] is not None))
            locators.append(dict(case_id=c, finding=f, reference=rr['findings'][f], member=rr['member'], sha256=rr['sha256']))
    suff = dict(summary='median signed tooth margin, binary64 mm', state_A=[-4.0, 0.0, 4.0, 8.0], state_B=[2.0, 2.0, 2.0, 2.0], summary_A=float(np.median([-4.0, 0.0, 4.0, 8.0])), summary_B=2.0, identity_error=0.0, downstream='any tooth signed margin < 0', downstream_A=True, downstream_B=False, downstream_difference=1, minimum_extension='Minimum signed margin suffices for any local negative sign; per-tooth signs/margins needed to localize.', external_referent=dict(kind='our_own_fixture', locator='code/screening.py::evaluate', compared_quantity='Mathematical summary sufficiency only, not anatomical clinical accuracy', refutes_us=True), physical_uncertainty_enclosure='MISSING; no affine physical sensitivity claimed')
    assert float(np.median(suff['state_A'])) == float(np.median(suff['state_B'])) and any((v < 0 for v in suff['state_A'])) != any((v < 0 for v in suff['state_B']))
    r = dict(round=round, claim_type='information_link', external_referent=read(P / f'PREREG_{round}.json')['external_referent'], findings=out, sufficiency_test=suff, legacy_test_reuse=True, partition_counts={k: len(v) for (k, v) in m.items()}, evaluation_wall_s=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, review_state='PENDING_INDEPENDENT_REVIEW', injections=dict(all_positive_prediction_specificity=0.0, all_positive_rejected_at_both_targets=True, labels_reversal_detected=True), true_noise_floor='UNKNOWN')
    write(P / f'raw/RESULTS_{round}.json', r)
    write(P / f'raw/ERRORS_{round}.json', errors)
    write(P / f'raw/EXTERNAL_LOCATORS_{round}.json', locators)
    for f in FINDINGS:
        ids = [c for c in m['test'] if first_label(lab[c], f) is not None]
        yy = [first_label(lab[c], f) for c in ids]
        if any((v == 0 for v in yy)):
            assert rates(yy, [1] * len(yy))['specificity'] < 0.9
        if yy:
            assert sum((a != b for (a, b) in zip(yy, [1 - v for v in yy]))) == len(yy)
    rows = []
    for f in FINDINGS:
        q = out[f]['methods'].get('candidate', {})
        a = q.get('0.9', {})
        b = q.get('0.95', {})
        rows.append(dict(finding=f, n=out[f]['test_scored'], sp90=a.get('fixed_threshold_test'), sp95=b.get('fixed_threshold_test'), gates=[a.get('primary_gate'), b.get('primary_gate')]))
    write(P / f'raw/SUMMARY_{round}.json', rows)
    print(json.dumps(rows, indent=2), flush=True)
    (P / f'HANDOFF_{round}.md').write_text(f'{round}: per-finding fixed-specificity gates and all exclusions preserved in raw/RESULTS_{round}.json. True rater noise floor UNKNOWN. Changed construction and exact recipe in PREREG_{round}.json. Previously evaluated patient cohort; descriptive replay. Next: ' + ('replace arch summary with tooth-owned binary endpoint features.' if round == 'R1' else 'use same regional representation with all-report agreement-only training/calibration.' if round == 'R2' else 'acquire blinded named-reader adjudication and repeat scans on new patients; no further fit tuning on reused test.') + '\n')
    state(round + '_DECIDED', 'Per-finding gates preserved; true noise floor UNKNOWN', 'Next load-bearing construction per HANDOFF_' + round + '.md')
if __name__ == '__main__':
    stage = sys.argv[1]
    round = sys.argv[2] if len(sys.argv) > 2 else 'R1'
    {'freeze': lambda : freeze_r1() if round == 'R1' else fit_freeze(round), 'threshold': lambda : calfreeze(round), 'evaluate': lambda : evaluate(round), 'labels': extract_labels}[stage]()
