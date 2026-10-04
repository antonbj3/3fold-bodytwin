from common import *
from data import reports, load_geometry
from experiment import matrix, metrics, threshold, bootstrap
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import numpy as np, re, collections, sys, pickle

def labels(ids, g):
    ios = reports(ids, 'ios')
    photo = reports(ids)
    rows = []
    excluded = []
    for c in ids:
        rr = ios[c]
        pp = photo[c]
        missing = {}
        present = {}
        if rr:
            r = rr[0]
            text = r['source_text']
            sentences = re.split('(?<=[.!?])\\s+', text)
            for f in r['missing_teeth']:
                if f // 10 not in [1, 2, 3, 4]:
                    excluded.append(dict(case_id=c, fdi=f, reason='PRIMARY_ABSENCE_OUT_OF_SCOPE'))
                    continue
                clause = next((s for s in sentences if re.search('\\b' + str(f) + '\\b', s) and re.search('absence|missing|agenesis|absent|not present|lack', s, re.I)), None)
                if not clause:
                    excluded.append(dict(case_id=c, fdi=f, reason='NO_UNAMBIGUOUS_ABSENCE_CLAUSE'))
                    continue
                substitute = f + 40
                if re.search('decrown|root fragment|prosthe|replac|reshap', clause, re.I) or re.search('\\b' + str(substitute) + '\\b', text) or re.search('mixed dentition|primary dentition|deciduous dentition', text, re.I):
                    excluded.append(dict(case_id=c, fdi=f, reason='REPLACEMENT_REMNANT_OR_MIXED_DENTITION'))
                    continue
                missing[f] = dict(source_member=r['member'], source_sha256=r['sha256'], clause=clause)
            for s in sentences:
                if re.search('\\ball (?:permanent )?teeth are present\\b', s, re.I):
                    exceptions = {int(v) for v in re.findall('\\b([1-8][1-8])\\b', s)} if re.search('except', s, re.I) else set()
                    for q in range(1, 5):
                        for pos in range(1, 8):
                            f = q * 10 + pos
                            if f not in exceptions and f not in r['missing_teeth']:
                                present[f] = dict(source_member=r['member'], source_sha256=r['sha256'], clause=s)
        if pp:
            r = pp[0]
            for a in r['parsed']['assertions']:
                if a['polarity']:
                    for f in a['teeth']:
                        if rr and f in rr[0]['missing_teeth']:
                            excluded.append(dict(case_id=c, fdi=f, reason='PHOTO_PRESENT_IOS_MISSING_CONFLICT'))
                            missing.pop(f, None)
                            continue
                        present[f] = dict(source_member=r['member'], source_sha256=r['sha256'], clause=a['clause'])
        for (polarity, group) in [(1, missing), (0, present)]:
            for (f, s) in sorted(group.items()):
                if c not in g or f not in g[c]:
                    excluded.append(dict(case_id=c, fdi=f, reason='NO_GEOMETRY'))
                    continue
                if f in missing and f in present:
                    excluded.append(dict(case_id=c, fdi=f, reason='SOURCE_CONFLICT'))
                    continue
                rows.append(dict(case_id=c, fdi=f, y=polarity, **s, quantity='Surface tooth-site absent (visible substitute/remnant cases excluded)', resolution='PER_TOOTH', timescale='SIMULTANEOUS'))
    return (rows, excluded)

def fit():
    st = time.perf_counter()
    cpu = time.process_time()
    g = load_geometry()
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    (train, ex) = labels(m['train'], g)
    (cal, cex) = labels(m['calibration'], g)
    write('raw/R3_TRAIN_ROWS.json', train)
    write('raw/R3_CAL_ROWS.json', cal)
    write('raw/R3_TRAIN_CAL_EXCLUSIONS.json', ex + cex)
    cols = sorted({k for cc in g.values() for row in cc.values() for k in row})
    xx = matrix(train, g, cols)
    yy = np.array([r['y'] for r in train])
    calx = matrix(cal, g, cols)
    hgb = HistGradientBoostingClassifier(max_iter=150, max_depth=3, learning_rate=0.05, min_samples_leaf=20, l2_regularization=2, random_state=79)
    logit = make_pipeline(SimpleImputer(strategy='median', add_indicator=True), StandardScaler(), LogisticRegression(C=0.2, max_iter=2000, random_state=79))
    models = {'shape': hgb, 'logistic': logit}
    ths = {}
    for (name, model) in models.items():
        model.fit(xx, yy)
        ths[name] = float(threshold(cal, model.predict_proba(calx)[:, 1]))
    candidate = [dict(case_id=c, fdi=f) for c in m['test'] if c in g for f in sorted(g[c])]
    testx = matrix(candidate, g, cols)
    scores = {name: model.predict_proba(testx)[:, 1] for (name, model) in models.items()}
    pred = [dict(**r, scores={n: float(s[i]) for (n, s) in scores.items()}, predicted={**{n: int(s[i] >= ths[n]) for (n, s) in scores.items()}, 'any_vertex_rule': int(g[r['case_id']][r['fdi']]['label_count'] == 0), '100_vertex_rule': int(g[r['case_id']][r['fdi']]['label_count'] < 100)}) for (i, r) in enumerate(candidate)]
    write('raw/R3_PREDICTIONS.json', pred)
    modelp = DATA / 'absence_models.pkl'
    modelp.write_bytes(pickle.dumps(dict(models=models, columns=cols, thresholds=ths)))
    freeze('FROZEN_PREDICTIONS_R3.json', dict(frozen_utc=now(), predictions_path=str(ROOT / 'raw/R3_PREDICTIONS.json'), predictions_sha256=sha(ROOT / 'raw/R3_PREDICTIONS.json'), model_path=str(modelp), model_sha256=sha(modelp), thresholds=ths, prereg_sha256=sha(ROOT / 'PREREG_R3_TOOTH_ABSENCE_SPACING.json'), source_contract_sha256=sha(ROOT / 'code/absence.py'), parser_sha256=sha(ROOT / 'code/report_parser.py'), phase='Before numerical IOS-absence test comparison; older absence labels and exploratory text exposed, no untouched generalization'))
    write('raw/R3_FIT_COST.json', cost(st, cpu))
    print('R3 frozen', flush=True)

def evaluate():
    st = time.perf_counter()
    cpu = time.process_time()
    fr = read('FROZEN_PREDICTIONS_R3.json')
    assert sha(fr['predictions_path']) == fr['predictions_sha256']
    assert sha(ROOT / 'code/absence.py') == fr['source_contract_sha256']
    g = load_geometry()
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    (rows, ex) = labels(m['test'], g)
    pred = {(r['case_id'], r['fdi']): r for r in read(fr['predictions_path'])}
    for r in rows:
        r.update(pred[r['case_id'], r['fdi']])
    res = {}
    for name in ['shape', 'logistic', 'any_vertex_rule', '100_vertex_rule']:
        res[name] = metrics(rows, [r['predicted'][name] for r in rows])
        res[name]['bootstrap'] = bootstrap(rows, name)
    mm = res['shape']
    support = mm['positive_patients'] >= 10 and mm['negative_patients'] >= 20
    numeric = mm['balanced_accuracy'] >= 0.7 and mm['sensitivity'] >= 0.6 and (mm['specificity'] >= 0.8)
    result = dict(claim_type=['information_link', 'capability'], outcome='UNKNOWN_INSUFFICIENT_PATIENT_SUPPORT' if not support else 'PASS' if numeric else 'FAIL', metrics=res, numeric_gate=numeric, support_gate=support, external_referent=dict(kind='published_dataset', locator=str(ZIP) + '::*/reports_ios_en/*.txt + */reports_intraoral-photo_en/*.txt', compared_quantity='Named missing tooth sites versus explicit named visible treatment/all-teeth-present sites, replacement/remnant/mixed cases excluded', refutes_us=True), exclusions=ex, excluded_assertions=len(ex), exclusion_reasons=dict(collections.Counter((r['reason'] for r in ex))), included_assertions=len(rows), assertion_dropout_fraction=len(ex) / (len(ex) + len(rows)), evaluable_test_patients=len({r['case_id'] for r in rows}), patient_dropout_fraction=1 - len({r['case_id'] for r in rows}) / len(m['test']), cost=cost(st, cpu), limitations='Presence classes have different tooth/case spectrum; source-clause selection conservative. Existing X76 absence parser QA is producer review, not clinical adjudication.')
    write('raw/R3_TEST_ROWS.json', rows)
    write('raw/R3_RESULTS.json', result)
    write('raw/R3_ERROR_CASES.json', [r for r in rows if r['y'] != r['predicted']['shape']])
    (ROOT / 'HANDOFF_R3.md').write_text(f"R3 {result['outcome']}: absence sensitivity={mm['sensitivity']:.6f}, specificity={mm['specificity']:.6f}, BA={mm['balanced_accuracy']:.6f}. Named report source, not physical segmentation ground truth. Next: measure spacing positive-only agreement and refuse specificity where no explicit negatives exist; package independent material-label measurement port.\n")
    state('R3_ABSENCE_DECIDED', result['outcome'], 'Measure source support for spacing/wear, repeats, export figure and one-command demo')
    print(json.dumps({k: result[k] for k in ['outcome', 'metrics', 'excluded_assertions']}), flush=True)
if __name__ == '__main__':
    globals()[sys.argv[1]]()
