"""Replay immutable local measurements, predictions, falsification and typed queries."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse, datetime as dt, hashlib, json, math, resource, time
from pathlib import Path
import numpy as np
import mpmath as mp
from fractions import Fraction as F
from scipy.optimize import brentq
import ports
from sources import load
ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'raw'
RAW.mkdir(exist_ok=True)

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()

def write(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def freeze_predictions(name, predictions, training):
    """Called before held-point evaluation. Historical source visibility is declared."""
    contract = json.loads((ROOT / ('PREREG_' + name + '.json')).read_text())
    payload = {'construction': name, 'prereg_sha256': sha(ROOT / ('PREREG_' + name + '.json')), 'training': training, 'predictions': predictions, 'data_visibility': contract['data_visibility']}
    h = hashlib.sha256(canonical(payload)).hexdigest()
    path = ROOT / ('FROZEN_PREDICTIONS_' + name + '.json')
    if path.exists():
        old = json.loads(path.read_text())
        if old['payload_sha256'] != h or old['payload'] != payload:
            raise ValueError('prediction drift ' + name)
    else:
        write(path.name, {'frozen_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'payload_sha256': h, 'payload': payload})
    return h

def sufficiency(data):
    isq = {g: [F(70) + F(str(p['mean_ISQ'])) - F(str(v[0]['mean_ISQ'])) for p in v] for (g, v) in data['ISQ'].items()}
    isq_error = abs(float(isq['OD'][0]) - float(isq['PISP'][0]))
    later = abs(float(isq['OD'][-1]) - float(isq['PISP'][-1]))
    wear = {g: [F(1) + F(str(p['source_mean_1e-4g'])) / 10 - F(str(v[0]['source_mean_1e-4g'])) / 10 for p in v] for (g, v) in data['WEAR'].items()}
    (g1, g2) = ('Filtek Z250', 'Admira fusion flow')
    wear_error = abs(float(wear[g1][0]) - float(wear[g2][0]))
    wlater = abs(float(wear[g1][-1]) - float(wear[g2][-1]))
    wstates = []
    for m in [5.0, 15.0]:
        with mp.workdps(80):
            gamma = mp.gamma(1 + mp.mpf(1) / int(m))
            scale = mp.mpf(500) / gamma
            mean = scale * gamma
            if float(mean) != 500.0:
                raise ValueError('mean witness not exactly500 at reported binary64 precision')
            q = scale * (-mp.log(mp.mpf('0.95'))) ** (1 / mp.mpf(int(m)))
            wstates.append({'m': m, 'sigma0_MPa_80digit': mp.nstr(scale, 80), 'mean_MPa': float(mean), 'mean_reconstruction_residual_80digit': mp.nstr(mean - mp.mpf(500), 8), 'Q05_MPa': float(q)})
    err = abs(wstates[0]['mean_MPa'] - wstates[1]['mean_MPa'])
    out = {'K31': {'summary': 'initial ISQ after explicitly constructed common-origin shift', 'identity_error_ISQ': isq_error, 'downstream_difference_day90_ISQ': later, 'states': {k: [float(a) for a in v] for (k, v) in isq.items()}, 'minimal_supported_extension': 'protocol-indexed temporal response knots; two actual equal-initial patients are not claimed', 'resolution': 'POPULATION', 'external_referent': 'PMC10092180 Table2 supplies changes, normalization is our witness'}, 'K34_LTD': {'summary': 'mean strength', 'identity_error_MPa': err, 'states': wstates, 'downstream_difference_Q05_MPa': abs(wstates[0]['Q05_MPa'] - wstates[1]['Q05_MPa']), 'minimal_extension_within_declared_Weibull_family': 'm plus mean/scale, and exact exposure key', 'resolution': 'PHENOMENOLOGICAL', 'external_referent': 'our_own_fixture for insufficiency only; source Table4 validates actual distributions'}, 'K34_WEAR': {'summary': 'first cumulative mass after common-origin shift', 'identity_error_mg': wear_error, 'downstream_difference_120k_mg': wlater, 'states': {k: [float(a) for a in wear[k]] for k in [g1, g2]}, 'minimal_supported_extension': 'protocol/cycle-indexed ordered cumulative mass history, not an early rate alone', 'resolution': 'POPULATION', 'external_referent': 'PMC10258403 Table3 supplies changes, normalization is our witness'}}
    assert isq_error == 0 and wear_error == 0 and (err == 0)
    write('raw/SUFFICIENCY_TESTS.json', out)
    return out

def exponential_fit(points):
    y0 = points[0]['mean_ISQ']
    t = np.array([p['day'] for p in points[1:]], float)
    d = y0 - np.array([p['mean_ISQ'] for p in points[1:]])
    tau = np.logspace(-1, 3, 2001)
    basis = -np.expm1(-t[None, :] / tau[:, None])
    A = np.maximum(0, basis @ d / np.sum(basis * basis, axis=1))
    sse = np.sum((basis * A[:, None] - d) ** 2, axis=1)
    j = int(np.argmin(sse))
    return (float(A[j]), float(tau[j]), y0)

def exponential_control(points):
    y0 = points[0]['mean_ISQ']
    best = None
    for tau in np.logspace(-1, 3, 2001):
        b = [-math.expm1(-p['day'] / float(tau)) for p in points[1:]]
        d = [y0 - p['mean_ISQ'] for p in points[1:]]
        A = max(0, sum((u * v for (u, v) in zip(b, d))) / sum((u * u for u in b)))
        sse = sum(((A * u - v) ** 2 for (u, v) in zip(b, d)))
        row = (sse, A, float(tau))
        if best is None or sse < best[0]:
            best = row
    return (best[1], best[2], y0)

def assay(name, points, pred, key, value, tolerance):
    rows = []
    for p in points:
        q = pred[str(p[key])]
        expected = p[value]
        e = abs(q - expected)
        tol = tolerance(p)
        rows.append({key: p[key], 'predicted': q, 'measured': expected, 'abs_error': e, 'tolerance': tol, 'pass': bool(e <= tol), 'resolution': 'POPULATION'})
    return rows

def run_isq(data, phase=False):
    name = 'K31_R2' if phase else 'K31_R1'
    reg = json.loads((ROOT / ('PREREG_' + name + '.json')).read_text())
    predictions = {}
    training = {}
    controls = []
    fits = {}
    query_details = {}
    for (g, points) in data['ISQ'].items():
        fitpts = [p for p in points if p['day'] in reg['training_days']]
        training[g] = fitpts
        if phase:
            ts = [p['day'] for p in fitpts]
            ys = [p['mean_ISQ'] for p in fitpts]
            details = {str(t): ports.interpolate_enclosed(ts, ys, t) for t in reg['held_days']}
            pred = {t: v['mean'] for (t, v) in details.items()}
            query_details[g] = details
            ctrl = {}
            for t in reg['held_days']:
                j = next((i for i in range(len(ts) - 1) if ts[i] <= t <= ts[i + 1]))
                w = (t - ts[j]) / (ts[j + 1] - ts[j])
                ctrl[str(t)] = ys[j] + (ys[j + 1] - ys[j]) * w
        else:
            (A, tau, y0) = exponential_fit(fitpts)
            (Ac, tc, yc) = exponential_control(fitpts)
            fits[g] = {'A_ISQ': A, 'tau_days': tau, 'y0_ISQ': y0, 'resolution': 'PHENOMENOLOGICAL', 'debt': 'same-implant RFA and spatial stiffness/BIC needed; not a material law'}
            pred = {str(t): y0 - A * -math.expm1(-t / tau) for t in reg['held_days']}
            ctrl = {str(t): yc - Ac * -math.expm1(-t / tc) for t in reg['held_days']}
        delta = max((abs(pred[k] - ctrl[k]) for k in pred))
        controls.append({'group': g, 'independent_control_max_error_ISQ': delta, 'pass': delta <= 1e-12})
        predictions[g] = pred
    frozen = freeze_predictions(name, predictions, training)
    validations = {}
    practice = {}
    for (g, points) in data['ISQ'].items():
        held = [p for p in points if p['day'] in reg['held_days']]
        validations[g] = assay(name, held, predictions[g], 'day', 'mean_ISQ', lambda p: 1.0)
        old = {str(t): points[0]['mean_ISQ'] for t in reg['held_days']}
        practice[g] = assay(name, held, old, 'day', 'mean_ISQ', lambda p: 1.0)
    passes = sum((row['pass'] for v in validations.values() for row in v))
    total = sum((len(v) for v in validations.values()))
    mutation = all((not abs(predictions[g][str(p['day'])] + 10 - p['mean_ISQ']) <= 1 for (g, v) in data['ISQ'].items() for p in v if p['day'] in reg['held_days']))
    refusal = []
    for qty in ['BIC', 'stiffness_MPa', 'individual_ISQ', 'patient_failure_probability']:
        try:
            ports.promote_isq(qty)
            refusal.append({'request': qty, 'refused': False})
        except ports.PortError as e:
            refusal.append({'request': qty, 'refused': True, 'reason': str(e)})
    external = dict(reg['external_referent'])
    external['refutes_us'] = passes < total
    return {'construction': name, 'chain': 'K31', 'claim_type': reg['claim_type'], 'resolution': 'POPULATION', 'timescale': 'HANDOVER', 'outcome': 'SCOPED_OBSERVATION_PASS' if passes == total else 'TEMPORAL_CLOSURE_REFUTED', 'external_referent': external, 'frozen_predictions_sha256': frozen, 'validation': validations, 'passes': passes, 'total': total, 'project_static_practice': practice, 'fits': fits, 'query_details': query_details, 'equally_informed_correctness_control': controls, 'injected_10ISQ_detected': mutation, 'promotion_refusals': refusal, 'empirical_scope': 'published group means only; no new subjects, no BIC, no calibrated physical heating-to-healing coupling', 'uncertainty': reg['uncertainty']}

def run_ltd(data):
    name = 'K34_LTD_R1'
    reg = json.loads((ROOT / ('PREREG_' + name + '.json')).read_text())
    p = reg['query_probability']
    pred = {}
    ctrl = []
    for (grade, states) in data['LTD'].items():
        pred[grade] = {s: ports.material_query(v, p) for (s, v) in states.items()}
        for (s, v) in states.items():
            q = pred[grade][s]['quantile_MPa']
            m = v['m']
            scale = v['sigma0_MPa']
            q2 = scale * math.exp(math.log(-math.log(1 - p)) / m)
            inv = brentq(lambda x: -math.expm1(-(x / scale) ** m) - p, 1e-08, scale * 2, xtol=1e-10)
            rel = max(abs(q - q2), abs(q - inv)) / q
            ctrl.append({'grade': grade, 'state': s, 'relative_error': rel, 'pass': rel <= 1e-12})
    frozen = freeze_predictions(name, pred, data['LTD'])
    rows = []
    for (grade, states) in data['LTD'].items():
        before = pred[grade]['unaged']['quantile_MPa']
        after = pred[grade]['aged']['quantile_MPa']
        rows.append({'grade': grade, 'unaged_Q05_MPa': before, 'aged_Q05_MPa': after, 'change_percent': 100 * (after / before - 1), 'resolution': 'POPULATION', 'mean_change_sign_not_sufficient': True, 'sampling_CI': 'UNKNOWN', 'empirical_effect_sign_admitted': False})
    refusals = []
    record = data['LTD']['T']['aged']
    for ctx in ['crown_FE_per_point', '5year_clinical_life', 'another_antagonist', 'different_sinter_lot']:
        try:
            ports.material_query(record, p, context=ctx)
            refusals.append({'request': ctx, 'refused': False})
        except ports.PortError as e:
            refusals.append({'request': ctx, 'refused': True, 'reason': str(e)})
    bad = dict(record)
    bad['sigma0_MPa'] *= 1.1
    mutation = abs(ports.material_query(bad, p)['quantile_MPa'] - pred['T']['aged']['quantile_MPa']) / pred['T']['aged']['quantile_MPa'] > 1e-12
    external = dict(reg['external_referent'])
    external['refutes_us'] = False
    return {'construction': name, 'chain': 'K34', 'claim_type': 'information_link', 'outcome': 'SCOPED_EXPOSURE_STATE_PASS', 'resolution': 'POPULATION', 'timescale': 'HANDOVER', 'external_referent': external, 'frozen_predictions_sha256': frozen, 'query_results': rows, 'material_states': pred, 'correctness_control': ctrl, 'source_transcription_pass': True, 'injected_scale_10percent_detected': mutation, 'promotion_refusals': refusals, 'protocol': reg['protocol'], 'uncertainty': reg['uncertainty'], 'clinical_aging_law': 'UNKNOWN', 'corrosion': 'NOT_BUILT', 'wear_geometry': 'NOT_BUILT', 'published_conclusion': 'Authors do not find statistically significant mean strength reduction after this exposure; fit-derived tail changes do not overturn that conclusion.'}

def affine_control(x, y, intercept):
    xa = np.array(x, float)
    ya = np.array(y, float)
    candidates = [(0.0, max(0.0, float(xa @ ya / (xa @ xa))))]
    if intercept:
        candidates.append((max(0.0, float(ya.mean())), 0.0))
        A = np.column_stack([np.ones(len(xa)), xa])
        (a, k) = np.linalg.lstsq(A, ya, rcond=None)[0]
        if a >= 0 and k >= 0:
            candidates.append((float(a), float(k)))
    return min(candidates, key=lambda ak: float(np.sum((ya - ak[0] - ak[1] * xa) ** 2)))

def run_wear(data, mode):
    name = 'K34_WEAR_R' + str(mode)
    reg = json.loads((ROOT / ('PREREG_' + name + '.json')).read_text())
    preds = {}
    train = {}
    controls = []
    fits = {}
    for (g, points) in data['WEAR'].items():
        fitpts = [p for p in points if p['cycles'] in reg['training_cycles']]
        train[g] = fitpts
        x = [p['cycles'] for p in fitpts]
        y = [str(F(str(p['source_mean_1e-4g'])) / 10) for p in fitpts]
        y = [F(v) for v in y]
        if mode == 3:
            preds[g] = {str(t): ports.wear_envelope(x, y, t) for t in reg['held_cycles']}
            delta = 0.0
            for t in reg['held_cycles']:
                j = next((i for i in range(len(x) - 1) if x[i] <= t <= x[i + 1]))
                (lo, hi) = (float(y[j]) - 0.005, float(y[j + 1]) + 0.005)
                got = preds[g][str(t)]['mass_interval_mg']
                delta = max(delta, abs(lo - got[0]), abs(hi - got[1]))
        else:
            ak = ports.fit_rate(x, y, intercept=mode == 2)
            (ac, kc) = affine_control(x, list(map(float, y)), intercept=mode == 2)
            preds[g] = {str(t): ports.rate_query(ak, t) for t in reg['held_cycles']}
            fits[g] = {'run_in_mg': float(ak[0]), 'rate_mg_per_cycle': float(ak[1]), 'exact_rationals': [str(v) for v in ak], 'resolution': 'PHENOMENOLOGICAL', 'debt': 'measure later mass/volume, sliding and hydration; extrapolation enclosure UNKNOWN'}
            delta = max((abs(preds[g][str(t)]['mean_mg'] - (ac + kc * t)) for t in reg['held_cycles']))
        controls.append({'material': g, 'max_arithmetic_difference_mg': delta, 'pass': delta <= 1e-10})
    frozen = freeze_predictions(name, preds, train)
    validation = {}
    mutations = []
    for (g, points) in data['WEAR'].items():
        rows = []
        for p in points:
            if p['cycles'] not in reg['held_cycles']:
                continue
            pred = preds[g][str(p['cycles'])]
            if mode == 3:
                (lo, hi) = pred['mass_interval_mg']
                passed = lo <= p['mean_mg'] <= hi
                row = {'cycles': p['cycles'], 'measured_mean_mg': p['mean_mg'], 'interval_mg': [lo, hi], 'width_mg': pred['width_mg'], 'pass': passed, 'resolution': 'POPULATION'}
                mutation = not lo <= p['mean_mg'] + 10 <= hi
            else:
                error = abs(pred['mean_mg'] - p['mean_mg'])
                tol = max(0.1, 1.96 * p['SEM_mg'])
                passed = error <= tol
                row = {'cycles': p['cycles'], 'predicted_mg': pred['mean_mg'], 'measured_mean_mg': p['mean_mg'], 'abs_error_mg': error, 'tolerance_mg': tol, 'pass': passed, 'resolution': 'POPULATION'}
                mutation = abs(pred['mean_mg'] - (p['mean_mg'] + 10)) > tol
            rows.append(row)
            mutations.append({'material': g, 'cycles': p['cycles'], 'injected_10mg_detected': mutation})
        validation[g] = rows
    refusals = []
    for g in data['WEAR']:
        try:
            ports.mass_to_volume(1)
            refusals.append({'material': g, 'request': 'volume_without_density', 'refused': False})
        except ports.PortError as e:
            refusals.append({'material': g, 'request': 'volume_without_density', 'refused': True, 'reason': str(e)})
        if mode == 3:
            pp = train[g]
            x = [p['cycles'] for p in pp]
            y = [p['mean_mg'] for p in pp]
            try:
                ports.wear_envelope(x, y, 150000)
                refusals.append({'material': g, 'request': 'cycle150000', 'refused': False})
            except ports.PortError as e:
                refusals.append({'material': g, 'request': 'cycle150000', 'refused': True, 'reason': str(e)})
    passes = sum((p['pass'] for v in validation.values() for p in v))
    total = sum((len(v) for v in validation.values()))
    external = dict(reg['external_referent'])
    external['refutes_us'] = passes < total
    return {'construction': name, 'chain': 'K34', 'claim_type': 'capability', 'resolution': 'POPULATION', 'timescale': 'HANDOVER', 'outcome': ('CONDITIONAL_MEASURED_STATE_PASS' if mode == 3 else 'SCOPED_RATE_PASS') if passes == total else 'RATE_TRANSFER_REFUTED', 'external_referent': external, 'validation': validation, 'passes': passes, 'total': total, 'predictions': preds, 'fits': fits, 'frozen_predictions_sha256': frozen, 'equally_informed_correctness_controls': controls, 'mutation_controls': mutations, 'promotion_refusals': refusals, 'uncertainty': reg['uncertainty'], 'full_antagonist_wear_chain': 'UNKNOWN: tip material and signed sliding unreported', 'full_clinical_geometry_chain': 'UNKNOWN: density, local spatial wear and matched anatomy absent'}

def update_state(stage, attempts):
    write('CURRENT_WORK_STATE.json', {'lane': 'X45-missing-chains-r2', 'updated_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'status': stage, 'latest_gate': {k: v['outcome'] for (k, v) in attempts.items()}, 'next_operation': {'initial': 'Freeze K31 temporal branches and K34 wear run-in successor', 'successor': 'Freeze measured monotone wear-state successor after failed rate transfer', 'envelope': 'Complete consumer query, figures, graph feedback and handoff', 'all': 'Independent review; same-specimen temporal stiffness/BIC and wear-volume measurements required'}[stage], 'limits': {'threads': 1, 'max_intermediates_bytes': 3000000000, 'gpu': False}, 'agents': 0})

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=['initial', 'successor', 'envelope', 'all'], default='all')
    a = ap.parse_args()
    start = time.perf_counter()
    cpu = time.process_time()
    data = load()
    suff = sufficiency(data)
    tasks = {'initial': [('K31_R1', lambda : run_isq(data)), ('K34_LTD_R1', lambda : run_ltd(data)), ('K34_WEAR_R1', lambda : run_wear(data, 1))], 'successor': [('K31_R2', lambda : run_isq(data, True)), ('K34_WEAR_R2', lambda : run_wear(data, 2))], 'envelope': [('K34_WEAR_R3', lambda : run_wear(data, 3))]}
    stages = list(tasks) if a.stage == 'all' else [a.stage]
    for stage in stages:
        for (name, fn) in tasks[stage]:
            before = time.perf_counter()
            result = fn()
            result['kernel_wall_seconds'] = time.perf_counter() - before
            write('raw/RESULT_' + name + '.json', result)
            print(name, result['outcome'], str(result.get('passes', '')) + '/' + str(result.get('total', '')))
    attempts = {}
    for p in sorted(RAW.glob('RESULT_*.json')):
        item = json.loads(p.read_text())
        attempts[item['construction']] = item
    elapsed = time.perf_counter() - start
    write('results.json', {'schema': 'X45-missing-chains-r2-v1', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'claim_type': ['information_link', 'capability'], 'external_referent': [v['external_referent'] for v in attempts.values()], 'constructed_chains': ['K31 observational time port', 'K34 exposure-conditioned coupon state and mass-wear history'], 'full_Northstar_chain_promotions': 0, 'attempts': attempts, 'sufficiency_tests': suff, 'source_manifest': data['manifest'], 'parser_controls': data['parser_controls'], 'cost': {'current_stage': a.stage, 'wall_seconds': elapsed, 'CPU_seconds': time.process_time() - cpu, 'max_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fit_discovery_validation_queries': 'All above timed from primary XML loading through output; no FE/GPU jobs', 'preparation_and_reasoning': 'COMMANDS + read manifest; model tokens and total discovery wall time UNKNOWN', 'inherited_discovery': 'UNKNOWN; X43/X44 and old cells read without rerun', 'physical_fallback': 'NOT_RUN, cost UNKNOWN'}, 'attrition': {'selected_chain_contracts': 7, 'new_chains_constructed_as_bounded_subports': 2, 'unbuilt_chain_contracts': 5, 'normalized_measurement_means': 44, 'source_strength_means_in_raw_Table3': 8, 'extracted_distribution_parameters': 16, 'source_historical_patient_screening': {'screened': 42, 'included': 27, 'excluded': 15, 'reason': 'source inclusion/exclusion criteria'}, 'source_failed_implants': {'initial': 54, 'in_function_after_followup': 53, 'failed_OD': 1}, 'note': 'Visit-level group denominator after failed implant is UNKNOWN; source SD not converted to mean CI.'}})
    write('ATTEMPTS.json', [{'construction': k, 'parent': json.loads((ROOT / ('PREREG_' + k + '.json')).read_text()).get('parent'), 'outcome': v['outcome'], 'evidence': 'raw/RESULT_' + k + '.json', 'next_operation': 'preserve failed closure, change state representation' if 'REFUTED' in v['outcome'] else 'independent review; physical promotion remains blocked'} for (k, v) in attempts.items()])
    update_state(a.stage, attempts)
    print('WALL', elapsed, 'RSS_MiB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
if __name__ == '__main__':
    main()
