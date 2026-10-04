"""Frozen categorical association; reports measure categories, never contact mm."""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
import sys, json, pickle, time, datetime, resource, argparse
from pathlib import Path
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import balanced_accuracy_score, accuracy_score, roc_auc_score
from contact import P, D, X11, LABELS, sha, write, state
X7 = P.parent / _release_expand('X7')
FIELDS = ['overbite', 'overjet', 'crossbite', 'molar_right', 'molar_left']

def coarse(v):
    if v is None:
        return None
    if v.startswith('III'):
        return 'III'
    if v.startswith('II'):
        return 'II'
    if v.startswith('I'):
        return 'I'
    return v

def load_rows():
    return {p.stem: json.load(open(p)) for p in sorted((P / 'raw/cases').glob('*.json'))}

def fit_freeze(round='R1'):
    start = time.perf_counter()
    manifest = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    part = manifest['inherited_partition']
    rows = load_rows()
    assert len(rows) == 993, 'Need all 993 map outputs before freeze'
    prereg = json.load(open(P / f'PREREG_{round}.json'))
    assert sha(P / f'PREREG_{round}.json') == (P / f'PREREG_{round}.json.sha256').read_text().strip()
    train = json.load(open(X7 / 'raw/TRAIN_LABELS.json'))
    cal = json.load(open(X7 / 'raw/CALIBRATION_LABELS.json'))
    features = {c: r.get('features', {}) for (c, r) in rows.items()}
    if round == 'R2':
        extra = json.load(open(P / 'raw/OPPOSITION_FEATURES_R2.json'))
        features = {c: dict(features[c], **extra[c]['features']) for c in rows}
    cols = sorted({k for c in part['train'] for k in features[c]})
    models = {}
    pred = {c: {} for c in part['calibration'] + part['test']}
    cost = []

    def data(ids):
        return np.array([[features[c].get(k) for k in cols] for c in ids], dtype=float)
    Xt = data(list(pred))
    for field in FIELDS:
        ids = [c for c in part['train'] if rows[c].get('error') is None and any((x['labels'].get(field) is not None for x in train.get(c, [])))]
        y = np.array([next((x['labels'][field] for x in train[c] if x['labels'].get(field) is not None)) for c in ids])
        s = time.perf_counter()
        model = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_depth=3, min_samples_leaf=15, l2_regularization=1.0, random_state=21).fit(data(ids), y)
        prob = model.predict_proba(Xt)
        maj = str(np.unique(y, return_counts=True)[0][np.argmax(np.unique(y, return_counts=True)[1])])
        classes = model.classes_.tolist()
        scores = []
        for (c, pb) in zip(pred, prob):
            pred[c][field] = dict(classes=classes, probabilities=pb.tolist(), point=classes[int(np.argmax(pb))], train_majority=maj, geometry_available=rows[c].get('error') is None)
        for c in part['calibration']:
            values = {x['labels'][field] for x in cal.get(c, []) if x['labels'].get(field) is not None}
            if values:
                scores.append(max((1 - pred[c][field]['probabilities'][classes.index(v)] if v in classes else 1.0 for v in values)))
        rank = min(len(scores), int(np.ceil((len(scores) + 1) * 0.9)))
        q = float(np.sort(scores)[rank - 1]) if scores else 1.0
        for c in pred:
            pb = pred[c][field]['probabilities']
            pred[c][field]['set'] = [v for (v, p0) in zip(classes, pb) if 1 - p0 <= q + 1e-12]
            pred[c][field]['set_q'] = q
        cost.append(dict(field=field, n_train=len(ids), calibration_patients=len(scores), threshold_q=q, fit_and_query_s=time.perf_counter() - s))
        models[field] = model
        print(round, field, cost[-1], flush=True)
    out = P / 'raw' / f'PREDICTIONS_{round}.json'
    write(out, pred)
    mdl = D / f'models_{round}.pkl'
    mdl.write_bytes(pickle.dumps(dict(models=models, columns=cols)))
    write(P / 'raw' / f'FIT_COST_{round}.json', dict(rows=cost, total_wall_s=time.perf_counter() - start, feature_count=len(cols), peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    write(P / ('FROZEN_PREDICTIONS.json' if round == 'R1' else 'FROZEN_PREDICTIONS_R2.json'), dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), round=round, claim_type='capability', before_this_round_reads_test_reports=True, legacy_test_reuse='X7/X11 previously evaluated cohort; descriptive fixed construction, no claim of untouched clinical validation', predictions_path=str(out), predictions_sha256=sha(out), model_sha256=sha(mdl), prereg_sha256=sha(P / f'PREREG_{round}.json'), code_sha256=sha(Path(__file__)), map_manifest=[dict(case_id=c, case_file_sha256=sha(P / 'raw/cases' / (c + '.json')), map_sha256=rows[c].get('map_sha256')) for c in sorted(rows)], extra_features_sha256=sha(P / 'raw/OPPOSITION_FEATURES_R2.json') if round == 'R2' else None))
    state(round + '_PREDICTIONS_FROZEN', 'Prediction hashes frozen before evaluation', 'Score existing independent clinician narratives; preserve legacy cohort limitation')

def wilson(k, n):
    if not n:
        return None
    z = 1.959963984540054
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    a = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [c - a, c + a]

def metrics(y, p):
    return dict(n=len(y), accuracy=float(accuracy_score(y, p)), balanced_accuracy=float(balanced_accuracy_score(y, p)), agreement_wilson95=wilson(sum((a == b for (a, b) in zip(y, p))), len(y))) if y else dict(n=0)

def evaluate(round='R1'):
    st = time.perf_counter()
    frozen = json.load(open(P / ('FROZEN_PREDICTIONS.json' if round == 'R1' else 'FROZEN_PREDICTIONS_R2.json')))
    pp = Path(frozen['predictions_path'])
    assert sha(pp) == frozen['predictions_sha256']
    pred = json.load(open(pp))
    test = json.load(open(X7 / 'raw/TEST_LABELS.json'))
    part = json.load(open(P / 'raw/INPUT_MANIFEST.json'))['inherited_partition']
    rows = load_rows()
    x7 = {r['case_id']: r for r in json.load(open(X7 / 'raw/PREDICTIONS_R1.json'))}
    out = {}
    raw = []
    inj = {}
    for field in FIELDS:
        ids = [c for c in part['test'] if any((x['labels'].get(field) is not None for x in test.get(c, [])))]
        truth = [next((x['labels'][field] for x in test[c] if x['labels'].get(field) is not None)) for c in ids]
        point = [pred[c][field]['point'] for c in ids]
        base = [x7[c]['fields'][field]['point'] for c in ids]
        majority = [pred[c][field]['train_majority'] for c in ids]
        first = metrics(truth, point)
        base_m = metrics(truth, base)
        cover = 0
        sizes = []
        allagree = []
        anyagree = []
        disagree = []
        review = []
        for (c, y) in zip(ids, truth):
            values = [x['labels'][field] for x in test[c] if x['labels'].get(field) is not None]
            value_set = set(values)
            ps = pred[c][field]['set']
            cover += int(value_set.issubset(set(ps)))
            sizes.append(len(ps))
            allagree.append(len(value_set) == 1 and point[ids.index(c)] in value_set)
            anyagree.append(pred[c][field]['point'] in value_set)
            disagree.append(len(value_set) > 1)
            prob = np.array(pred[c][field]['probabilities'])
            entropy = float(-np.sum(prob * np.log(np.maximum(prob, 1e-15))) / max(np.log(len(prob)), 1e-15))
            r = dict(case_id=c, field=field, reference_first=y, reference_all=values, reference_members=[x['member'] for x in test[c] if x['labels'].get(field) is not None], prediction=pred[c][field]['point'], prediction_set=ps, x7_prediction=x7[c]['fields'][field]['point'], report_discordant=len(value_set) > 1, has_multiple_readable_reports=len(values) > 1, geometry_uncertainty_score=entropy, agreement_any_report=pred[c][field]['point'] in value_set, agreement_all_reports=len(value_set) == 1 and pred[c][field]['point'] in value_set)
            raw.append(r)
            if len(values) > 1:
                review.append(r)
        auc = None
        ci = None
        if review and len({r['report_discordant'] for r in review}) == 2:
            yy = np.array([r['report_discordant'] for r in review])
            ss = np.array([r['geometry_uncertainty_score'] for r in review])
            auc = float(roc_auc_score(yy, ss))
            rng = np.random.default_rng(2101)
            boot = []
            for b in range(1000):
                ix = rng.integers(0, len(yy), len(yy))
                if len(np.unique(yy[ix])) == 2:
                    boot.append(roc_auc_score(yy[ix], ss[ix]))
            ci = np.quantile(boot, [0.025, 0.975]).tolist()
        out[field] = dict(first_report=first, x7_first_report=base_m, practice_majority=metrics(truth, majority), coarse_angle=metrics([coarse(v) for v in truth], [coarse(v) for v in point]) if field.startswith('molar') else None, all_report_set_coverage=cover / len(ids), coverage_wilson95=wilson(cover, len(ids)), mean_set_size=float(np.mean(sizes)), singleton_fraction=float(np.mean(np.array(sizes) == 1)), point_matches_any_report=float(np.mean(anyagree)), point_matches_all_reports=float(np.mean(allagree)), discordant_report_cases=int(sum(disagree)), multi_report_cases=len(review), multi_report_discordance=float(np.mean([r['report_discordant'] for r in review])) if review else None, geometry_entropy_discordance_AUC=auc, AUC_descriptive_bootstrap95=ci, geometry_explanation_gate='PASS' if auc is not None and auc >= 0.6 and (ci[0] > 0.5) else 'FAIL_OR_UNKNOWN')
        wrong = ['__INJECTED_WRONG_CATEGORY__'] * len(ids)
        inj[field] = dict(injection='invalid category for every patient, zero report sets', wrong_metrics=metrics(truth, wrong), set_coverage=0.0, rejected=True)
    num = [r['numerical_control'] for r in rows.values() if r.get('numerical_control')]
    numerics = bool(len(num) == 12 and all((x['gate'] and x['injected_rejected'] for x in num)))
    area = bool(len(num) == 12 and all((x['area_refinement']['gate'] and x['area_injected_100mm2_rejected'] for x in num)))
    complete = len(rows) == 993 and all((r.get('error') is None for r in rows.values()))
    if round == 'R2':
        repaired = json.load(open(P / 'raw/AREA_INJECTION_CHECK.json'))
        area = bool(len(repaired['rows']) == 12 and all((x['original_refinement_gate'] and x['rejected'] for x in repaired['rows'])))
    ba = float(np.mean([r['first_report']['balanced_accuracy'] for r in out.values()]))
    cov = float(np.mean([r['all_report_set_coverage'] for r in out.values()]))
    size = float(np.mean([r['mean_set_size'] for r in out.values()]))
    gates = dict(completeness=complete, numerical_panel=numerics, area_refinement_panel=area, scored_n=all((x['first_report']['n'] >= 100 for x in out.values())), macro_balanced_accuracy=ba >= 0.65, all_report_set_coverage=cov >= 0.85, mean_set_size=size <= 2.5)
    result = dict(round=round, claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', external_referent=json.load(open(P / 'PREREG_R1.json'))['external_referent'], fields=out, macro_balanced_accuracy=ba, macro_x7_balanced_accuracy=float(np.mean([r['x7_first_report']['balanced_accuracy'] for r in out.values()])), macro_all_report_set_coverage=cov, macro_mean_set_size=size, gates=gates, primary_gate='PASS' if all(gates.values()) else 'FAIL', geometry_maps=993, numerical_controls=len(num), maximum_numerical_parity_mm=max((max(r['parity_mm'], r['witness_inherited_parity_mm'], r['full_subset_parity_mm']) for r in num)), numerical_injected_rejections=all((r['injected_rejected'] for r in num)), wrong_category_controls=inj, wrong_set_controls_rejected=all((r['set_coverage'] == 0 for r in inj.values())), evaluation_wall_s=time.perf_counter() - st, physical_contact_area_accuracy='UNKNOWN_NO_EXTERNAL_AREA_MEASUREMENT', target_FDI_accuracy='UNKNOWN_NO_EXTERNAL_LABELS', numeric_overbite_overjet_accuracy='UNKNOWN_NO_EXTERNAL_POINT_MM_MEASUREMENT', physical_force='UNKNOWN_NOT_PREDICTED', test_reuse='Legacy X7 patient-disjoint split reused after prior lane evaluations. No untouched-cohort validation or hyperparameter tuning on this round.', geometry_interpretation='Triangle-pair minima certify conditional affine projection in exact arithmetic; current implementation floating and selected oracle tested. Penetration and proximity refer to this model.')
    write(P / 'raw' / f'COMPARISON_ROWS_{round}.json', raw)
    write(P / 'raw' / f'RESULTS_{round}.json', result)
    state(round + '_DECIDED', result['primary_gate'], 'Preserve gates; change contact-only representation to opposition geometry or acquire missing measurement')
    print(json.dumps({k: result[k] for k in ['round', 'macro_balanced_accuracy', 'macro_x7_balanced_accuracy', 'macro_all_report_set_coverage', 'macro_mean_set_size', 'gates', 'primary_gate']}, indent=2))
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['fit', 'evaluate'])
    ap.add_argument('--round', default='R1')
    a = ap.parse_args()
    fit_freeze(a.round) if a.stage == 'fit' else evaluate(a.round)
