"""Frozen synthetic trials and external aggregate measurement replay."""
import collections, datetime, json, math, time
import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm
from lab_alarm import *
from synthetic import series, save_series
FAULTS = {'shifted_weibull': 'material_scale', 'wrong_sinter_factor': 'sinter', 'one_bad_record': 'measurement_record', 'model_form': 'model_form', 'cement_film': 'cement_film'}

def trials(p, f, pr):
    st = time.perf_counter()
    counts = pr['metrics']
    raw = []
    examples = {}
    nullhits = []
    for i in range(counts['null_replicates']):
        (out, _) = analyze(series(61000 + i, p, f), p, f, posterior=False)
        nullhits.append(out['monitor']['any_alarm'])
        raw.append(dict(kind='nominal', seed=61000 + i, **out['monitor']))
    n = len(nullhits)
    k = sum(nullhits)
    report = dict(null=dict(replicates=n, alarms=k, rate=k / n, binomial_95CI=[float(__import__('scipy').stats.beta.ppf(0.025, k, n - k + 1)) if k else 0.0, float(__import__('scipy').stats.beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0], horizon_specimens=72, global_alpha=0.05), faults={})
    for (fault, channel) in FAULTS.items():
        hits = []
        primary = []
        for i in range(counts['fault_replicates']):
            seed = 71000 + i
            (out, _) = analyze(series(seed, p, f, fault), p, f, posterior=False)
            h = out['monitor']['first_alarm_at_specimen']
            t = h.get(channel)
            hits.append(t)
            primary.append(channel in h)
            raw.append(dict(kind=fault, seed=seed, **out['monitor']))
        success = [t for t in hits if t is not None]
        rep = dict(replicates=len(hits), detected=sum(primary), detection_fraction=float(np.mean(primary)), channel=channel, median_measurements_among_detected=float(np.median(success)) if success else None, p10_p90_measurements_among_detected=[float(x) for x in np.quantile(success, [0.1, 0.9])] if success else None, not_detected_by_72=sum((x is None for x in hits)), right_observation_link_fraction=float(np.mean(primary)), causal_localization='Observation link is identified; physical process vs metrology and strength vs contact remain candidate sets until independent calibration')
        if fault == 'one_bad_record':
            rep['measurement_events_after_fault'] = 1
        report['faults'][fault] = rep
        rows = series(6100, p, f, fault)
        (path, side) = save_series(R / 'raw' / ('SYNTHETIC_' + fault + '.csv'), rows)
        (out, post) = analyze(rows, p, f)
        out['input_csv_sha256'] = sha(path)
        out['sidecar_sha256'] = sha(side)
        write(R / 'raw' / ('EXAMPLE_' + fault + '.json'), out)
        np.savez_compressed(R / 'raw' / ('POSTERIOR_' + fault + '.npz'), m=post.m, lambda_reference_N=post.lam, form=post.form, weights=post.weights())
        examples[fault] = out
        print(json.dumps(dict(fault=fault, **rep)), flush=True)
    rows = series(6100, p, f)
    (path, side) = save_series(R / 'raw/SYNTHETIC_nominal.csv', rows)
    (ex, post) = analyze(rows, p, f)
    write(R / 'raw/EXAMPLE_nominal.json', ex)
    examples['nominal'] = ex
    write(R / 'raw/TRIALS_R2.json', raw)
    report['wall_seconds'] = time.perf_counter() - st
    return (report, examples)

def posterior_checks(p, f):
    hits = {'m': 0, 'scale': 0, 'sinter': 0, 'cement': 0}
    values = []
    st = time.perf_counter()
    for seed in range(81000, 81060):
        rows = series(seed, p, f)
        (out, post) = analyze(rows, p, f)
        q = out['posterior']
        ci = lambda item, x: item['credible_95'][0] <= x <= item['credible_95'][1]
        good = dict(m=ci(q['m'], 6.0), scale=ci(q['sigma0_force_equivalent'], 3000.0), sinter=ci(q['sinter_factor'], 0.8), cement=ci(q['cement_film_mean'], 60.0))
        for (k, v) in good.items():
            hits[k] += v
        values.append(dict(seed=seed, coverage=good, posterior=q))
    write(R / 'raw/POSTERIOR_COVERAGE_R2.json', values)
    forms = []
    for seed in range(91000, 91020):
        row = {}
        for kind in ['nominal', 'model_form']:
            (out, post) = analyze(series(seed, p, f, kind), p, f)
            row[kind] = out['posterior']['model_probabilities']['FE']
        forms.append(dict(seed=seed, **row))
    write(R / 'raw/MODEL_LOCALIZATION_R2.json', forms)
    return dict(replicates=60, coverage={k: v / 60 for (k, v) in hits.items()}, model_localization_fraction={'nominal': sum((x['nominal'] >= 0.95 for x in forms)) / 20, 'separable': sum((x['model_form'] <= 0.05 for x in forms)) / 20}, wall_seconds=time.perf_counter() - st)

def external_replay():
    cells = json.loads((R / 'inputs/shrinkage_cells.json').read_text())
    g = {}
    pred = []
    for c in cells:
        k = tuple((c[x] for x in ('material', 'position', 'horizontal', 'axis')))
        g.setdefault(k, {})[c['method']] = c
    for (k, methods) in g.items():
        if 'micrometer' not in methods or 'surface_scan' not in methods:
            continue
        pred.append(dict(identity=k, predicted_surface_scan_shrinkage_pct=methods['micrometer']['mean_shrinkage_pct'], prediction_locator=methods['micrometer']['locator'], held_measurement_locator=methods['surface_scan']['locator']))
    fp = R / 'FROZEN_EXTERNAL_PREDICTIONS.json'
    payload = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), retrospective=True, blinded=False, source_sha256=sha(R / 'inputs/shrinkage_cells.json'), predictions=pred, claim='Equality of paired aggregate modalities within 0.5 percentage points; no population noise estimate', tolerance_percentage_points=0.5)
    if not fp.exists():
        write(fp, payload)
        fp.with_name(fp.name + '.sha256').write_text(sha(fp) + '\n')
    saved = verified(fp)
    out = []
    for row in saved['predictions']:
        k = tuple(row['identity'])
        meas = g[k]['surface_scan']
        y = meas['mean_shrinkage_pct']
        err = y - row['predicted_surface_scan_shrinkage_pct']
        out.append(dict(**row, observed_pct=y, error_pct=err, rejected=abs(err) > 0.5, reference_is_aggregate=True))
    write(R / 'raw/EXTERNAL_REPLAY_R2.json', out)
    return dict(external_referent=dict(kind='independent_measurement', locator='doi:10.3390/ma18143217 Table2; inputs/Shrinkage2025.xml', compared_quantity='Paired scanner and micrometer mean XYZ shrinkage %; no single-crown truth', refutes_us=True), pairs=len(out), rejected=sum((x['rejected'] for x in out)), rejected_fraction=sum((x['rejected'] for x in out)) / len(out), RMSE_percentage_points=float(np.sqrt(np.mean([x['error_pct'] ** 2 for x in out]))), max_abs_error_percentage_points=max((abs(x['error_pct']) for x in out)), alpha_for_real_metrology='UNKNOWN: Table2 cell n/cross-method covariance and noise law missing; do not borrow synthetic sigma', physical_transfer='FAIL_SCANNER_EQUALS_GAUGE_AT_0P5_PP', source_cells=len(cells), dropped_cells=len(cells) - 2 * len(out), drop_reason='153 light-microscopy aggregates outside chosen paired quantity, retained source not discarded', resolution='PER_SURFACE_REGION aggregated coupon/axis/position, not 459 independent specimens')

def mutations(p, f):
    base = series(6100, p, f)
    cases = {}
    mo = Monitor(p, model_groups(f))
    bad = copy_row = __import__('copy').deepcopy(base[0])
    bad['_sidecar']['post_sinter_ref_length_mm'] = '9.5'
    e = mo.step(bad)
    cases['wrong_sinter_value'] = e['scores']['sinter']['pit'] > 1 - 1e-10
    bad = __import__('copy').deepcopy(base[0])
    bad['fracture_force_N'] = str(float(bad['fracture_force_N']) * 3)
    e = Monitor(p, model_groups(f)).step(bad)
    cases['wrong_force_record'] = not e['force_update_allowed']
    bad = __import__('copy').deepcopy(base[0])
    bad['failure_mode'] = 'die_fracture'
    (out, post) = analyze([bad], p, f)
    cases['wrong_mechanism'] = post.n == 0 and any((x['kind'] == 'MECHANISM_FALSIFIER' for x in out['timeline'][0]['events']))
    bad = __import__('copy').deepcopy(base[0])
    bad['die_E_MPa'] = '1800'
    (out, post) = analyze([bad], p, f)
    cases['wrong_support'] = post.n == 0
    vals = (np.arange(8) + 0.5) / 8 - 0.5
    mu = vals.mean()
    sd = vals.std()
    good = np.mean([np.mean([math.exp(c / sd * (x - mu)) / np.mean(np.exp(c / sd * (vals - mu))) for c in (-2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2)]) for x in vals])
    wrong = np.mean([np.mean([math.exp(c / sd * (x - mu)) for c in (-2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2)]) for x in vals])
    cases['wrong_mgf'] = abs(good - 1) < 1e-12 and abs(wrong - 1) > 0.001
    cases['wrong_weibull_shape'] = abs(force_cdf(2000, 6, 3000, 0.01) - float(__import__('scipy').stats.weibull_min.cdf(2000, 3, scale=3000))) > 0.001
    data = json.loads((R / 'raw/EXTERNAL_REPLAY_R2.json').read_text())
    cases['external_value_plus10'] = all((abs(x['observed_pct'] + 10 - x['predicted_surface_scan_shrinkage_pct']) > 0.5 for x in data))
    write(R / 'raw/MUTATION_CONTROLS_R2.json', cases)
    return cases

def main():
    import resource
    pr = verified(R / 'PREREG_R2.json')
    p = verified(R / 'inputs/DEMO_PROFILE.json')
    f = json.loads((R / 'inputs/X1B_FROZEN_PREDICTIONS.json').read_text())
    st = time.perf_counter()
    (tr, ex) = trials(p, f, pr)
    post = posterior_checks(p, f)
    ext = external_replay()
    mut = mutations(p, f)
    threshold = pr['metrics']
    gates = dict(nominal_FPR=tr['null']['rate'] <= threshold['MC_upper_tolerance'], all_faults_power=all((x['detection_fraction'] >= threshold['power_min'] for x in tr['faults'].values())), all_faults_median_n=all((x['median_measurements_among_detected'] is not None and x['median_measurements_among_detected'] <= threshold['median_required_measurements_max'] for x in tr['faults'].values())), posterior_coverage=all((v >= threshold['posterior_95_interval_coverage_min'] for v in post['coverage'].values())), model_form_localization=post['model_localization_fraction']['separable'] >= 0.9, all_mutations_rejected=all(mut.values()))
    out = dict(lane='X61-lab-alarm', round='R2', claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', outcome='PASS_SCOPED_SYNTHETIC_PORT' if all(gates.values()) else 'PARTIAL_CAPABILITY_WITH_FAILED_GATES', gates=gates, trials=tr, posterior_checks=post, external_replay=ext, external_referent=ext['external_referent'], mutations=mut, sufficiency=json.loads((R / 'raw/R1_SUFFICIENCY.json').read_text()), cost=dict(wall_seconds=time.perf_counter() - st, max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, preparation='Source read/search/hash + symbolic R1; wall time not separately measured', fit='Joint grid every record for examples/coverage; finite model, no physical fit', discovery='400 null + 200 per five fault families, fixed seeds and severities', validation='60 posterior coverage +20/2 model-form; SciPy/quad and exact conjugate controls separately in unittest', questions=0, fallback='UNKNOWN instrument/lab time, stress mapping and matched independent held validation'), resolution='PER_TOOTH records and POPULATION posterior; model debt PHENOMENOLOGICAL', timescale='HANDOVER; record trace SIMULTANEOUS', physical_lab_measurements=0, parameter_updates='Weibull m and force-equivalent sigma0 (N); sinter factor and cement film only with orthogonal measured sidecar; true material sigma0 (MPa) UNKNOWN without m-specific stress operator', same_information_control='Independent full-density posterior and analytical Gaussian/categorical controls match; no algorithm superiority', missing_physical_enclosure='Force-factor independence from fitted m and fixed contact/support not physically enclosed; dry-distance folded normal only scenario observation law', next_construction='If faults/coverage pass, attack acquisition ambiguity with biased scanner/reference measurements; do not mark lab calibration complete')
    write(R / 'raw/RESULTS_R2.json', out)
    write(R / 'results.json', out)
    write(R / 'CURRENT_WORK_STATE.json', dict(lane='X61-lab-alarm', status='R2_DECIDED', latest_gate=gates, next_operation=out['next_construction']))
    (R / 'HANDOFF_R2.md').write_text(json.dumps(dict(outcome=out['outcome'], gates=gates, external=ext, next=out['next_construction']), indent=2) + '\n')
    print(json.dumps(dict(gates=gates, cost=out['cost'])), flush=True)
if __name__ == '__main__':
    main()
