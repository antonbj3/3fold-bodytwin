"""Narrow negation-scope correction. Original X7 labels retained unchanged."""
from dental_release.paths import expand as _release_expand
import re, json, zipfile, hashlib, time, datetime, pickle
from pathlib import Path
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import balanced_accuracy_score, accuracy_score, roc_auc_score
from contact import P, D, sha, write
X7 = P.parent / _release_expand('X7')
MENTION = re.compile('cross[ -]?bites?', re.I)
MOD = '(?:(?:lateral|posterior|anterior|buccal|lingual|bilateral|unilateral|right|left|dental|a|the)\\s+){0,4}'
NEG = re.compile('(?:\\bno\\s+|\\bwithout\\s+|\\babsence\\s+of\\s+|\\bneither\\s+)' + MOD + 'cross[ -]?bites?', re.I)
AFTER = re.compile('^\\s+(?:(?:are|is|were|was)\\s+)?(?:not\\s+(?:observed|present|noted)|absent|not\\s+seen)', re.I)

def repair_text(text, original):
    decisions = []
    for sentence in re.split('[.;!\\n]', text):
        mm = list(MENTION.finditer(sentence))
        if not mm:
            continue
        negatives = list(NEG.finditer(sentence))
        exception = bool(re.search('\\b(?:except|however|but|excluding)\\b', sentence, re.I))
        for m in mm:
            neg = any((n.start() <= m.start() < n.end() for n in negatives)) or bool(AFTER.search(sentence[m.end():]))
            prefix = sentence[max(0, m.start() - 45):m.start()].lower()
            scopes = [s for s in ['lateral', 'posterior', 'anterior', 'right', 'left', 'bilateral', 'unilateral'] if re.search('\\b' + s + '\\b', prefix)]
            decisions.append(dict(mention=m.group(), explicitly_negated=neg, scope=scopes or ['unspecified'], exception_in_clause=exception, clause=sentence.strip()))
    allneg = bool(decisions) and all((x['explicitly_negated'] and (not x['exception_in_clause']) for x in decisions))
    corrected = 'absent' if allneg else original
    return (corrected, decisions, allneg)

def extract():
    st = time.perf_counter()
    original = json.load(open(X7 / 'raw/ALL_STRUCTURED_REPORTS.json'))
    manifest = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    out = {}
    changed = []
    no_lateral = []
    with zipfile.ZipFile(manifest['zip']) as z:
        for (c, reports) in original.items():
            out[c] = []
            for r in reports:
                blob = z.read(r['member'])
                assert hashlib.sha256(blob).hexdigest() == r['sha256']
                text = blob.decode('utf-8', errors='replace')
                lab = r['labels']
                old = lab.get('crossbite')
                (value, clauses, allneg) = repair_text(text, old)
                newlab = dict(lab)
                newlab['crossbite'] = value
                if allneg:
                    newlab['crossbite_teeth'] = []
                record = dict(member=r['member'], sha256=r['sha256'], labels=newlab, original_labels=lab, all_crossbite_mentions_negated=allneg, crossbite_scope_clauses=clauses, correction='all_mentions_explicitly_negated' if value != old else None)
                out[c].append(record)
                if value != old:
                    changed.append(dict(case_id=c, **record))
                if re.search('\\bno lateral cross[ -]?bites?', text, re.I):
                    no_lateral.append(dict(case_id=c, member=r['member'], sha256=r['sha256'], corrected=value, allneg=allneg, positive_injection_rejected=allneg and value == 'absent'))
    write(P / 'raw/CORRECTED_REPORTS_R5.json', out)
    write(P / 'raw/REPORT_CORRECTIONS_R5.json', changed)
    write(P / 'raw/NO_LATERAL_CONTROLS_R5.json', no_lateral)
    order = lambda r: hashlib.sha256(('X21-r5-manual:' + r['member']).encode()).hexdigest()
    negatives = sorted(changed, key=order)[:6]
    unchanged = sorted([dict(case_id=c, **r) for (c, rs) in out.items() for r in rs if not r['correction']], key=order)[:6]
    write(P / 'raw/REPORT_REVIEW_PANEL_R5.json', negatives + unchanged)
    write(P / 'raw/REPORT_REPAIR_COST_R5.json', dict(wall_s=time.perf_counter() - st, reports=sum(map(len, out.values())), changed_reports=len(changed), literal_no_lateral_reports=len(no_lateral)))
    print('reports', sum(map(len, out.values())), 'changes', len(changed), 'no lateral', len(no_lateral), flush=True)

def fit_freeze():
    st = time.perf_counter()
    m = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    part = m['inherited_partition']
    reports = json.load(open(P / 'raw/CORRECTED_REPORTS_R5.json'))
    rows = {c: json.load(open(P / 'raw/cases' / (c + '.json'))) for c in m['cases']}
    op = json.load(open(P / 'raw/OPPOSITION_FEATURES_R2.json'))
    cols = pickle.loads((D / 'models_R2.pkl').read_bytes())['columns']
    feat = {c: dict(rows[c]['features'], **op[c]['features']) for c in m['cases']}
    data = lambda ids: np.array([[feat[c].get(k) for k in cols] for c in ids], dtype=float)
    ids = [c for c in part['train'] if any((r['labels'].get('crossbite') is not None for r in reports.get(c, [])))]
    y = [next((r['labels']['crossbite'] for r in reports[c] if r['labels'].get('crossbite') is not None)) for c in ids]
    model = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_depth=3, min_samples_leaf=15, l2_regularization=1.0, random_state=21).fit(data(ids), y)
    qq = part['calibration'] + part['test']
    classes = model.classes_.tolist()
    prob = model.predict_proba(data(qq))
    pred = {c: dict(classes=classes, probabilities=pb.tolist(), point=classes[int(np.argmax(pb))]) for (c, pb) in zip(qq, prob)}
    scores = []
    for c in part['calibration']:
        values = {r['labels']['crossbite'] for r in reports.get(c, []) if r['labels'].get('crossbite') is not None}
        if values:
            scores.append(max((1 - pred[c]['probabilities'][classes.index(v)] if v in classes else 1.0 for v in values)))
    q = float(np.sort(scores)[min(len(scores), int(np.ceil((len(scores) + 1) * 0.9))) - 1])
    for c in pred:
        pred[c]['set'] = [v for (v, p) in zip(classes, pred[c]['probabilities']) if 1 - p <= q + 1e-12]
    write(P / 'raw/PREDICTIONS_R5.json', pred)
    (D / 'model_crossbite_R5.pkl').write_bytes(pickle.dumps(dict(model=model, columns=cols)))
    write(P / 'FROZEN_PREDICTIONS_R5.json', dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), claim_type='capability', prediction_sha256=sha(P / 'raw/PREDICTIONS_R5.json'), model_sha256=sha(D / 'model_crossbite_R5.pkl'), corrected_source_reports_sha256=sha(P / 'raw/CORRECTED_REPORTS_R5.json'), code_sha256=sha(Path(__file__)), prereg_sha256=sha(P / 'PREREG_R5.json'), before_corrected_test_comparison=True, legacy_test_reuse='Previously read/evaluated X7/X11/X21 source cohort, no untouched validation claim', scope='Reported crossbite mention with negation scope retained, not universal anatomy absence'))
    write(P / 'raw/FIT_COST_R5.json', dict(train_n=len(y), label_counts={k: y.count(k) for k in classes}, calibration_n=len(scores), q=q, fit_and_query_s=time.perf_counter() - st))
    print('fitR5', len(y), 'q', q, flush=True)

def discordance(reports, cases):
    n = 0
    k = 0
    rows = []
    for c in cases:
        vals = [r['labels']['crossbite'] for r in reports.get(c, []) if r['labels'].get('crossbite') is not None]
        if len(vals) > 1:
            n += 1
            k += int(len(set(vals)) > 1)
            rows.append(dict(case_id=c, values=vals, discordant=len(set(vals)) > 1))
    return dict(n_multi_report_cases=n, n_discordant_cases=k, discordance=k / n if n else None, rows=rows)

def evaluate():
    from report_capability import metrics, wilson
    st = time.perf_counter()
    frozen = json.load(open(P / 'FROZEN_PREDICTIONS_R5.json'))
    assert sha(P / 'raw/PREDICTIONS_R5.json') == frozen['prediction_sha256']
    assert sha(P / 'raw/CORRECTED_REPORTS_R5.json') == frozen['corrected_source_reports_sha256']
    pred = json.load(open(P / 'raw/PREDICTIONS_R5.json'))
    reports = json.load(open(P / 'raw/CORRECTED_REPORTS_R5.json'))
    old = json.load(open(X7 / 'raw/ALL_STRUCTURED_REPORTS.json'))
    m = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    test = m['inherited_partition']['test']
    r2 = json.load(open(P / 'raw/PREDICTIONS_R2.json'))
    x7 = {r['case_id']: r for r in json.load(open(X7 / 'raw/PREDICTIONS_R1.json'))}
    ids = [c for c in test if any((r['labels'].get('crossbite') is not None for r in reports.get(c, [])))]
    y = [next((r['labels']['crossbite'] for r in reports[c] if r['labels'].get('crossbite') is not None)) for c in ids]
    scores = {name: metrics(y, p) for (name, p) in [('corrected_label_model', [pred[c]['point'] for c in ids]), ('frozen_R2_model', [r2[c]['crossbite']['point'] for c in ids]), ('frozen_X7_model', [x7[c]['fields']['crossbite']['point'] for c in ids])]}
    covered = 0
    size = []
    rr = []
    negctrl = json.load(open(P / 'raw/NO_LATERAL_CONTROLS_R5.json'))
    for (c, truth) in zip(ids, y):
        values = [r['labels']['crossbite'] for r in reports[c] if r['labels'].get('crossbite') is not None]
        s = pred[c]['set']
        covered += int(set(values).issubset(set(s)))
        size.append(len(s))
        rr.append(dict(case_id=c, source_members=[r['member'] for r in reports[c]], source_hashes=[r['sha256'] for r in reports[c]], reference_first=truth, reference_all=values, prediction=pred[c], original_references=[r['labels'].get('crossbite') for r in old.get(c, [])], frozen_R2=r2[c]['crossbite']['point'], frozen_X7=x7[c]['fields']['crossbite']['point']))
    cov = covered / len(ids)
    ms = float(np.mean(size))
    all_rej = all((r['positive_injection_rejected'] for r in negctrl))
    gates = dict(balanced_accuracy=scores['corrected_label_model']['balanced_accuracy'] >= 0.65, set_coverage=cov >= 0.85, set_size=ms <= 2.5, explicit_negation_injection=all_rej)
    result = dict(round='R5', claim_type='capability', external_referent=json.load(open(P / 'PREREG_R5.json'))['external_referent'], scores=scores, all_report_set_coverage=cov, mean_set_size=ms, source_repairs=json.load(open(P / 'raw/REPORT_REPAIR_COST_R5.json')), test_report_disagreement_original=discordance(old, test), test_report_disagreement_corrected=discordance(reports, test), all_cohort_disagreement_original=discordance(old, m['cases']), all_cohort_disagreement_corrected=discordance(reports, m['cases']), gates=gates, primary_gate='PASS' if all(gates.values()) else 'FAIL', negation_controls=dict(n=len(negctrl), all_positive_injections_rejected=all_rej, remaining_ambiguous=[r for r in negctrl if not r['positive_injection_rejected']]), original_X7_X21_gates_preserved=True, source_scope='Explicit reported absence can be lateral/anterior; all scopes retained. No universal clinical absence conclusion.', physical_contact_validation='UNKNOWN', evaluation_wall_s=time.perf_counter() - st)
    write(P / 'raw/COMPARISON_ROWS_R5.json', rr)
    write(P / 'raw/RESULTS_R5.json', result)
    print(json.dumps({k: result[k] for k in ['round', 'scores', 'all_report_set_coverage', 'mean_set_size', 'source_repairs', 'gates', 'primary_gate', 'negation_controls']}, indent=2))
if __name__ == '__main__':
    import sys
    globals()[sys.argv[1]]()
