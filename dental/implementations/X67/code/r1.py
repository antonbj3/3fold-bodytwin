import time, numpy as np
from scipy.linalg import expm
from common import *
import legacy_drill_thermal as dt
SD = [2.8, 3.4, 3.8, 4.4]
OBS = {'none': [4.03, 2.46, 2.76, 2.45], '12': [1.19, 1.04, 1.05, 0.76], '30': [0.57, 0.43, 1.05, 0.37]}

def sindel(irr, **kw):
    p = dict(dt.PARAMS)
    p.update(kw.pop('params', {}))
    steps = [dict(D=d, d_prev=0 if i == 0 else SD[i - 1], depth=10.0, rpm=1000.0, feed=2.0, irrig=irr, Tcool=37.0, pause=0 if i == 0 else 30.0) for (i, d) in enumerate(SD)]
    r = dt.simulate(steps, t_cort_mm=1.5, T0=37.0, perfusion=False, h_mm=0.1, probes=[('tc', 3.0, 5.0)], cool_after_s=30.0, p=p, R_extra_mm=8.0, Z_extra_mm=6.0, **kw)
    return r

def per_drill(r, name):
    out = []
    t = r['t']
    y = r['probes'][name]
    for s in r['steps']:
        i0 = np.searchsorted(t, s['t_start_s'])
        i1 = np.searchsorted(t, s['t_end_drill_s'] + 30.0)
        out.append(float(np.nanmax(y[i0:i1]) - y[max(0, i0 - 1)]))
    return out

def abuhajar(**kw):
    p = dict(dt.PARAMS)
    p.update(kw.pop('params', {}))
    out = []
    for limiting in [True, False]:
        r = dt.simulate([dict(D=2.2, d_prev=0, depth=9, rpm=800, feed=2, irrig=True, Tcool=21)], t_cort_mm=1.5, T0=21, perfusion=False, h_mm=0.1, probes=[('c', 2.1, 1.5), ('s', 2.1, 7)], cool_after_s=30, p=p, wall_irrig_during=not limiting, R_extra_mm=8, Z_extra_mm=6, **kw)
        out.extend((float(np.nanmax(r['probes'][n]) - 21) for n in ['c', 's']))
    return out

def vec(s, a):
    return np.array(s['none'] + s['irr'] + s['irr'] + a)

def observed():
    return np.array(OBS['none'] + OBS['12'] + OBS['30'] + [0.82, 0.78, 0.3, 0.35])

def metrics(pred, y):
    lr = np.log(np.maximum(pred, 1e-12) / y)
    return dict(within_factor2_count=int(np.sum(np.abs(lr) <= np.log(2))), within_factor2_fraction=float(np.mean(np.abs(lr) <= np.log(2))), median_abs_log_error=float(np.median(abs(lr))), median_log_bias=float(np.median(lr)), mae_C=float(np.mean(abs(pred - y))))

def main():
    start = time.perf_counter()
    verify('PREREG_R1.json')
    old = __import__('json').loads((DENTAL / 'results/K3_drilling/anchors_thermal.json').read_text())
    y = observed()
    a = [old['abuhajar_central'][g][b] for g in ['limiting', 'nonlimiting'] for b in ['cortical', 'cancellous']]
    base = vec(old['sindel_central'], a)
    diagnostic = []

    def add(name, v, note):
        diagnostic.append(dict(mechanism=name, predicted_rise_C=v.tolist(), metrics=metrics(v, y), status='DIAGNOSIS_ONLY', note=note))
    add('inherited_h0.05', base, 'Actual recorded model; unknown feed was 2 mm/s')
    add('eta_physical_ceiling_1', 2 * base, 'Exact source-amplitude scaling for Tcool=T0 and fixed operator; eta=1 ceiling; numerical rounding enclosure absent')
    for (name, kwargs) in [('chi0', {'chi': 0.0}), ('chi1', {'chi': 1.0}), ('k0.53', {'params': {'k_cort': 0.53}}), ('h_irr25', {'params': {'h_irr': 25.0}})]:
        s = {k: per_drill(sindel(k == 'irr', **kwargs), 'tc') for k in ['none', 'irr']}
        v = vec(s, abuhajar(**kwargs))
        add(name, v, 'One changed mechanism; native 0.1 mm operator. No measured input repaired.')
        dump('raw/bisection_' + name + '.json', dict(sindel=s, abuhajar=v[-4:].tolist()))
        state('R1_BISECTION', name + ' recorded', 'Continue physical bisection then freeze LOSO')
    mult = float(np.exp(np.median(np.log(y / base))))
    add('unbounded_scalar_fit', mult * base, 'Posthoc descriptive best L1 log fit; not a validation or physical partition')
    C = np.array([1.0, 3.0])
    g = 0.3
    A = np.array([[-g / C[0], g / C[0]], [g / C[1], -g / C[1]]])
    EA = np.array([20.0, 0.0])
    EB = np.array([0.0, 20.0])
    ta = expm(A) * EA[None, :] / C[None, :]
    tb = expm(A) * EB[None, :] / C[None, :]
    wall_a = float(ta.sum(axis=1)[0])
    wall_b = float(tb.sum(axis=1)[0])
    suff = dict(kind='two_control_volume_illustration', external_referent='our_own_fixture', summary='total residual-bone energy J', identical_summary_A_J=float(EA.sum()), identical_summary_B_J=float(EB.sum()), identity_error_J=float(EA.sum() - EB.sum()), time_s=1.0, wall_rise_A_C=wall_a, wall_rise_B_C=wall_b, downstream_difference_C=wall_a - wall_b, minimum_extension='source location / local energy inventory; total energy alone is insufficient')
    assert suff['identity_error_J'] == 0 and suff['downstream_difference_C'] > 0.5
    extracted = load('raw/measurements.json')

    def get(s, g):
        return next((d['value_C'] for d in extracted if d['study'] == s and d['group'] == g))
    guide = dict(study='PMC10299697', locator='jcm-12-03944-t001 TC1 Mean; Groups A,C', compared_quantity='published table mean absolute temperature at 2 mm depth', same_nominal_flow_ml_min=35.0, A_entry_exit_C=get('PMC10299697', 'A_2'), C_handpiece_C=get('PMC10299697', 'C_2'), difference_C=get('PMC10299697', 'C_2') - get('PMC10299697', 'A_2'), resolution='POPULATION', support='PER_POINT', interpretation='Flow volume alone does not encode guide access; no identification of numeric h from this contrast')
    result = dict(claim_type='information_link', legacy_observed_rise_C=y.tolist(), legacy_predicted_rise_C=base.tolist(), scalar_fit_multiplier=mult, implied_eta_for_scalar_fit=0.5 * mult, diagnostic=diagnostic, sufficiency=suff, external_mechanism_measurement=guide, mechanism_owner='NOT_IDENTIFIED', reasons=['eta alone cannot absorb scalar bias without violating energy fraction', 'unknown feed, cooling access and sensor distance covary', 'one-at-a-time sweeps are not joint physical identification', 'Sindel is inherited extraction, primary fulltext absent'], strict_primary_physical_calibration_rows=load('raw/extraction_checks.json')['strict_complete_rows'], faults=dict(eta_greater_than_one_rejected=not 0 <= 0.5 * mult <= 1, primary_plus10_C_rejected=load('raw/extraction_checks.json')['injected_plus10_rejected']), runtime_s=time.perf_counter() - start)
    dump('attempts/R1_results.json', result)
    (ROOT / 'attempts/HANDOFF_R1.md').write_text(f'R1: fixed-source scalar multiplier {mult:.8g} would require eta={0.5 * mult:.8g}. Primary thermal source rows exist, but strict complete calibration rows=0. Native bisection retained; no unique physical owner identified. Total energy sufficiency identity error=0 J; wall difference={wall_a - wall_b:.8g} C in declared two-volume illustration. Next: freeze study-held-out source-location fit with explicit nominal-process debts.\n')
    state('R1_DECIDED', 'Physical ownership NOT_IDENTIFIED; scalar eta repair inadmissible', 'R2: leave-one-study-out bounded edge/wall source fit')
    print('R1', metrics(base, y), 'scalar eta', 0.5 * mult)
if __name__ == '__main__':
    main()
