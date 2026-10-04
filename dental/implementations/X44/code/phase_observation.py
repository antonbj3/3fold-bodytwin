"""K50 R2: empirical population observation response is not a constitutive material.
Work sign below is a conditional diagnostic. The article gives magnitudes but no
work-conjugate sign or equilibrated zero-energy initial state for this plotted mean.
"""
import numpy as np
from common import *
from calibration_port import SOURCE, frozen, alarm

def interpolate(knots, x):
    q = sorted(knots, key=lambda t: t['deflection_mm'])
    u = np.array([a['deflection_mm'] for a in q])
    f = np.array([a['F_mean_N'] for a in q])
    a = float(np.interp(x, u, f))
    return dict(F_mean_N=a, digitization_interval_N=[a - 0.3, a + 0.3], population_range_N=[float(np.interp(x, u, [v['F_range_lo_N'] for v in q])), float(np.interp(x, u, [v['F_range_hi_N'] for v in q]))])

def direct_affine(knots, x):
    q = sorted(knots, key=lambda t: t['deflection_mm'])
    for (a, b) in zip(q, q[1:]):
        if a['deflection_mm'] <= x <= b['deflection_mm']:
            w = (x - a['deflection_mm']) / (b['deflection_mm'] - a['deflection_mm'])
            return (1 - w) * a['F_mean_N'] + w * b['F_mean_N']
    raise ValueError('outside observed displacement support')

def observation_query(model, x, phase, quantity_unit='N'):
    if phase not in ['loading', 'unloading']:
        raise ValueError('phase required; force is not single-valued in displacement')
    if quantity_unit != 'N':
        raise ValueError('Table N/mm is ambiguous; figure force N selected with explicit provenance')
    if not 0 <= x <= 0.2:
        raise ValueError('outside observed deflection support')
    return interpolate(model[phase], x)

def run():
    pr = verify('PREREG_K50_R2.json')
    st = time.perf_counter()
    d = json.loads(SOURCE.read_text())
    rows = d['points']
    model = {phase: [rows[i] for i in pr['calibration_' + phase + '_k']] for phase in ['loading', 'unloading']}
    selected = pr['held_loading_k'] + pr['held_unloading_k']
    pred = [dict(k=i, phase=rows[i]['phase'], deflection_mm=rows[i]['deflection_mm'], prediction=observation_query(model, rows[i]['deflection_mm'], rows[i]['phase'])) for i in selected]
    fr = frozen('FROZEN_PREDICTIONS_R2.json', dict(round='R2', source_sha256=sha(SOURCE), prereg_sha256=sha(L / 'PREREG_K50_R2.json'), calibration_knots=model, predictions=pred, context='Historical descriptive split. Not blind, not new individual measurements.'))
    scored = []
    faults = []
    parity = []
    for pp in pred:
        r = rows[pp['k']]
        y = r['F_mean_N']
        yp = pp['prediction']['F_mean_N']
        ref = direct_affine(model[r['phase']], r['deflection_mm'])
        parity.append(abs(yp - ref))
        scored.append(dict(**r, prediction=pp['prediction'], error_N=abs(y - yp), tolerance_N=0.3 + 0.087 * y, pass_=not alarm(y, yp), resolution='POPULATION'))
        faults.append(dict(k=r['k'], kind='plus5N', rejected=bool(alarm(yp + 5, yp))))
    for (kind, args) in [('missing_phase', (0.1, None, 'N')), ('wrong_unit', (0.1, 'loading', 'N/mm')), ('out_of_support', (0.3, 'loading', 'N'))]:
        try:
            observation_query(model, *args)
            faults.append(dict(kind=kind, rejected=False))
        except ValueError as ex:
            faults.append(dict(kind=kind, rejected=True, reason=str(ex)))
    u = np.array([r['deflection_mm'] for r in rows])
    F = np.array([r['F_mean_N'] for r in rows])
    du = np.diff(u)
    coef = np.zeros(len(F))
    coef[:-1] += 0.5 * du
    coef[1:] += 0.5 * du
    work = float(coef @ F)
    error = float(0.3 * np.abs(coef).sum())
    ref_work = sum(((F[i] + F[i + 1]) * (u[i + 1] - u[i]) / 2 for i in range(len(F) - 1)))
    guard = dict(signed_work_if_plotted_magnitude_is_conjugate_Nmm=work, digitization_interval_Nmm=[work - error, work + error], work_quadrature_parity_Nmm=abs(work - ref_work), closed_displacement_cycle=bool(u[0] == u[-1]), passivity_if_equilibrated_and_sign_correct=bool(work + error >= 0), force_sign_known=False, initial_stored_energy_known=False, device_energy_known=False, material_promotion='REFUSED_MISSING_CONJUGATE_SIGN_AND_INITIAL_STATE', resolution='POPULATION', falsification_scope='If plotted force is a signed conjugate resistance and initial state is relaxed passive cycle, negative work would refute that interpretation. Without these inputs it is a protocol obstruction, not biological anti-passivity.')
    zero = rows[0]['F_mean_N']
    last = rows[-1]['F_mean_N']
    guard['residual_force_N'] = last - zero
    guard['residual_digitization_interval_N'] = [last - zero - 0.6, last - zero + 0.6]
    port = dict(id='healthy_incisor_group_0p5s_fig2b', anatomy='healthy_incisor_population', resolution='POPULATION', timescale='HANDOVER', source_locator=pr['external_referent']['locator'], source_sha256=sha(SOURCE), model=model, quantity_unit='N', input_unit='mm', phase_required=True, supported_deflection_mm=[0, 0.2], material_promotion=guard['material_promotion'], parameter_posterior_for_patient_E='UNKNOWN', phenomenological_debt=pr['phenomenological_debt'])
    write(L / 'OBSERVATION_PORT.json', port)
    gates = dict(held_observation=all((r['pass_'] for r in scored)), independent_affine=max(parity) <= 1e-10, injection_rejection=all((f['rejected'] for f in faults)), independent_work=abs(work - ref_work) < 1e-12, unsupported_material_promotion_refused=True)
    report = dict(chain='K50', round='R2', claim_type='capability', outcome='PASS_OBSERVATION_ONLY_MATERIAL_PROMOTION_REFUSED' if all(gates.values()) else 'FAIL', gates=gates, external_referent=pr['external_referent'], resolution='POPULATION', timescale='HANDOVER', frozen_predictions=fr, held_points=scored, held_pass_count=sum((r['pass_'] for r in scored)), held_count=len(scored), max_error_N=max((r['error_N'] for r in scored)), independent_affine_max_error_N=max(parity), protocol_guard=guard, fault_injections=faults, rejection=dict(source_points=len(rows), calibration=12, held=10, missing_points=0, material_promotions_considered=1, material_promotions_rejected=1, reason='Unmeasured force sign, initial state and matching geometry'), cost=cost(st), limitations=['Figure digitization ±0.3 N, no individual raw data; prior authors interpreted force, Table 2 labels N/mm: retained ambiguity', 'Two source branches are observed curves; not a new constitutive law or external patient validation', 'Initial energy and force sign missing; negative apparent cycle work cannot refute passive biology', 'All published points seen in legacy file before new prereg; split is descriptive reanalysis'], same_information_control='Piecewise affine interpolation matches; no algorithm superiority. Gain is preserving phase/protocol and refusing unsupported physical promotion.')
    write(L / 'raw/K50_R2.json', report)
    state('THREE_CHAINS_EXECUTED', gates, 'Complete 12-chain blockers and consumer demo with figure/graph feedback')
    print(json.dumps(dict(round='R2', gates={k: bool(v) for (k, v) in gates.items()}, held_pass=report['held_pass_count'], conditional_work=work, work_interval=guard['digitization_interval_Nmm'])))
if __name__ == '__main__':
    run()
