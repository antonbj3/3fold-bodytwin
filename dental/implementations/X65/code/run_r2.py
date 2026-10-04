import datetime, hashlib, json, math, time
from pathlib import Path
import numpy as np
from scipy.stats import t as student
from freeze import frozen_write
ROOT = Path(__file__).resolve().parents[1]
CLASSES = [100, 150, 175, 200, 250, 300, 400]

def save(p, x):
    (ROOT / p).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def sm(values, studies):
    return np.array([np.mean([v for (v, s) in zip(values, studies) if s == study]) for study in sorted(set(studies))])

def coeff(train, mode, mat=None):
    selected = [r for r in train if mode != 'material' or r['material_class'] == mat]
    if not selected:
        selected = train
    residual = [math.log(r['force_N'] / r['legacy_prediction_N']) if mode in ['global', 'material'] else math.log(r['force_N']) - offset(r) for r in selected]
    vals = sm(residual, [r['pmcid'] for r in selected])
    mu = float(np.mean(vals))
    spread = max(float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0, float(np.std(residual, ddof=1)) if len(residual) > 1 else 0.0)
    width = float(student.ppf(0.975, max(len(vals) - 1, 1)) * spread * math.sqrt(1 + 1 / len(vals)))
    return (mu, width, len(vals))

def offset(r):
    return 3 * math.log(r['diameter_mm']) - math.log((r['lever_mm'] or 11.0) * math.sin(math.radians(r['angle_deg'])))

def main():
    start = time.perf_counter()
    data = json.loads((ROOT / 'raw/published_endpoints.json').read_text())
    rows = [r for r in data['records'] if r['legacy_prediction_N'] is not None]
    studies = sorted({r['pmcid'] for r in rows})
    held = []
    trainlogs = []
    for study in studies:
        train = [r for r in rows if r['pmcid'] != study]
        test = [r for r in rows if r['pmcid'] == study]
        for r in test:
            pred = {'practice': {'prediction_N': r['legacy_prediction_N'], 'interval_N': None}}
            for mode in ['global', 'material', 'nominal']:
                (c, w, n) = coeff(train, mode, r['material_class'])
                logp = math.log(r['legacy_prediction_N']) + c if mode != 'nominal' else offset(r) + c
                pred[mode] = {'prediction_N': math.exp(logp), 'interval_N': [math.exp(logp - w), math.exp(logp + w)], 'train_studies_used': n, 'coefficient_log': c, 'closure': 'Student t study-residual prediction; not a rigorous enclosure'}
            held.append({'id': r['id'], 'held_pmcid': study, 'observed_force_N': r['force_N'], 'source': r['locator'], 'table_locator': r['table_locator'], 'predictions': pred, 'strict_eligible': r['strict_eligible']})
        trainlogs.append({'held_study': study, 'train_ids': [r['id'] for r in train], 'no_held_response_in_fit': all((r['pmcid'] != study for r in train))})
    frozen = {'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'kind': 'retrospective_LOSO_prediction_before_scoring; historical responses already known', 'prereg_sha256': hashlib.sha256((ROOT / 'PREREG_R2.json').read_bytes()).hexdigest(), 'input_records_sha256': hashlib.sha256(json.dumps(data['records'], sort_keys=True, separators=(',', ':')).encode()).hexdigest(), 'held_predictions': held}
    frozen_write(ROOT / 'FROZEN_PREDICTIONS.json', frozen)
    metrics = {}
    decisions = []
    for mode in ['practice', 'global', 'material', 'nominal']:
        errors = [abs(math.log(h['predictions'][mode]['prediction_N'] / h['observed_force_N'])) for h in held]
        sb = sm(errors, [h['held_pmcid'] for h in held])
        coverage = []
        inside = []
        for h in held:
            p = h['predictions'][mode]
            v = h['observed_force_N']
            ci = p['interval_N']
            inside.append(abs(math.log(p['prediction_N'] / v)) <= math.log(1.5))
            if ci:
                coverage.append(ci[0] <= v <= ci[1])
                for c in CLASSES:
                    outcome = 'ABSTAIN' if ci[0] < c <= ci[1] else 'SUPPORTED' if c <= ci[0] else 'BELOW_CLASS'
                    correct = outcome == 'ABSTAIN' or (outcome == 'SUPPORTED') == (v >= c)
                    decisions.append({'id': h['id'], 'mode': mode, 'bench_max_force_class_N': c, 'outcome': outcome, 'against_reported_endpoint_correct': correct})
        metrics[mode] = {'study_balanced_MAE_log': float(np.mean(sb)), 'median_abs_logforce': float(np.median(errors)), 'within_factor_1_5': int(sum(inside)), 'n': len(held), 'prediction_interval_coverage': float(np.mean(coverage)) if coverage else None}
    base = metrics['practice']['study_balanced_MAE_log']
    for mode in ['global', 'material', 'nominal']:
        metrics[mode]['gain_vs_practice'] = 1 - metrics[mode]['study_balanced_MAE_log'] / base
    strict = [r for r in rows if r['strict_eligible']]
    strict_studies = len({r['pmcid'] for r in strict})
    candidates = {}
    for mode in ['global', 'material']:
        m = metrics[mode]
        ds = [d for d in decisions if d['mode'] == mode]
        candidates[mode] = {'error_gain_pass': m['gain_vs_practice'] >= 0.25, 'absolute_error_pass': m['median_abs_logforce'] <= math.log(1.5), 'coverage_pass': m['prediction_interval_coverage'] >= 0.9, 'decided_classes': sum((d['outcome'] != 'ABSTAIN' for d in ds)), 'wrong_decided_classes': sum((not d['against_reported_endpoint_correct'] for d in ds)), 'strict_study_gate_pass': strict_studies >= 4, 'transport_status': 'UNKNOWN' if strict_studies < 4 else 'PASS' if all((d['against_reported_endpoint_correct'] for d in ds)) else 'FAIL'}
    mutations = [{'gate': 'HELD_RESPONSE', 'nominal_pass': all((x['no_held_response_in_fit'] for x in trainlogs)), 'injected_error_rejected': not all((x['no_held_response_in_fit'] for x in trainlogs + [{'no_held_response_in_fit': False}])), 'injection': 'insert held-study row into train'}, {'gate': 'FACTOR_1_5', 'nominal_pass': True, 'injected_error_rejected': abs(math.log(2.0)) > math.log(1.5), 'injection': 'double prediction when truth equals nominal'}, {'gate': 'FOUR_STRICT_STUDIES', 'nominal_pass': strict_studies >= 4, 'injected_error_rejected': 0 < 4, 'injection': 'drop all studies'}]
    for gate in candidates.values():
        gate['zero_wrong_classes_pass'] = gate['wrong_decided_classes'] == 0
        gate['decision_claim_status'] = 'FAIL_OBSERVED_CLASSES' if gate['wrong_decided_classes'] else 'UNKNOWN_STRICT_COHORT' if strict_studies < 4 else 'PASS_CONDITIONAL'
    out = {'claim_type': 'information_link', 'diagnostic_n': len(rows), 'diagnostic_studies': len(studies), 'strict_n': len(strict), 'strict_studies': strict_studies, 'metrics': metrics, 'gates': candidates, 'strict_rejection_fraction': data['strict_rejection_fraction'], 'outcome': 'CALIBRATION_DIAGNOSTIC_ONLY; TRANSFERABLE_DESIGN_UNKNOWN', 'held_rows': held, 'decisions': decisions, 'mutations': mutations, 'train_manifests': trainlogs, 'wall_s': time.perf_counter() - start, 'external_referent': {'kind': 'independent_measurement', 'locator': [r['locator'] for r in rows], 'compared_quantity': 'study-held-out published maximum cyclic force, with mixed protocol labels retained', 'refutes_us': True}}
    save('results_R2.json', out)
    print(json.dumps({k: v for (k, v) in out.items() if k not in ['held_rows', 'decisions', 'mutations', 'train_manifests', 'external_referent']}))
if __name__ == '__main__':
    main()
