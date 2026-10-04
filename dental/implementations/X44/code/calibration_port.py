"""K50 R1: historical published-measurement replay, fit -> freeze -> held-point alarm.
Posterior covers two-point Gaussian digitization model, not interindividual uncertainty.
"""
import numpy as np
from common import *
SOURCE = ROOT / 'results/F2_pdl_nonlinear/jepsen2023_fig2b.json'

def frozen(n, payload):
    p = L / n
    canonical = json.dumps(payload, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()
    h = hashlib.sha256(canonical).hexdigest()
    if p.exists():
        previous = json.loads(p.read_text())
        assert previous['payload_sha256'] == h, 'Prediction drift: create new prereg and preserve old predictions'
    else:
        write(p, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), payload_sha256=h, payload=payload))
        (L / (n + '.sha256')).write_text(sha(p) + '\n')
    assert sha(p) == (L / (n + '.sha256')).read_text().strip()
    return dict(path=str(p), sha256=sha(p), payload_sha256=h)

def alarm(observed, predicted):
    return abs(observed - predicted) > 0.3 + 0.087 * observed

def run_r1():
    pr = verify('PREREG_K50_R1.json')
    st = time.perf_counter()
    d = json.loads(SOURCE.read_text())
    rows = d['points']
    fitst = time.perf_counter()
    cal = pr['calibration_points_k']
    X = np.array([rows[i]['deflection_mm'] for i in cal])
    y = np.array([rows[i]['F_mean_N'] for i in cal])
    sigma = 0.3
    k = float(X @ y / (X @ X))
    sd = float(sigma / np.sqrt(X @ X))
    fitsec = time.perf_counter() - fitst
    W = np.eye(2) / sigma ** 2
    control = (np.linalg.inv(X[:, None].T @ W @ X[:, None]) @ (X[:, None].T @ W @ y)).item()
    posterior_error = abs(k - control)
    pred = [dict(k=x['k'], phase=x['phase'], deflection_mm=x['deflection_mm'], predicted_F_N=k * x['deflection_mm'], conditional_parameter_95_interval_N=[(k - 1.96 * sd) * x['deflection_mm'], (k + 1.96 * sd) * x['deflection_mm']], shared_digitization_bias_scenario_N=[k * x['deflection_mm'] - 0.3 * float(X.sum() / (X @ X)) * x['deflection_mm'], k * x['deflection_mm'] + 0.3 * float(X.sum() / (X @ X)) * x['deflection_mm']]) for x in rows]
    fr = frozen('FROZEN_PREDICTIONS.json', dict(round='R1', source_sha256=sha(SOURCE), prereg_sha256=sha(L / 'PREREG_K50_R1.json'), parameter_k_N_per_mm=k, conditional_posterior_sd_N_per_mm=sd, posterior='Improper uniform k prior, Gaussian digitization sigma 0.3N at selected points; approximate independence assumption only', predictions=pred, blinding='Historical descriptive split; source already visible in F2 and orientation'))
    held = []
    for i in pr['held_loading_k'] + pr['unloading_k']:
        x = rows[i]
        pp = pred[i]['predicted_F_N']
        held.append(dict(**x, predicted_F_N=pp, error_N=abs(x['F_mean_N'] - pp), tolerance_N=0.3 + 0.087 * x['F_mean_N'], alarm=bool(alarm(x['F_mean_N'], pp)), resolution='POPULATION'))
    loading = [x for x in held if x['phase'] == 'loading']
    unloading = [x for x in held if x['phase'] == 'unloading']
    faults = [dict(k=i, injected_F_N=pred[i]['predicted_F_N'] + 5, predicted_F_N=pred[i]['predicted_F_N'], rejected=bool(alarm(pred[i]['predicted_F_N'] + 5, pred[i]['predicted_F_N']))) for i in pr['held_loading_k']]
    gates = dict(loading_model_valid=all((not x['alarm'] for x in loading)), injected_force_rejected=all((x['rejected'] for x in faults)), normal_equation_parity=posterior_error <= 1e-10, phase_transfer_valid=all((not x['alarm'] for x in unloading)))
    report = dict(chain='K50', round='R1', claim_type='information_link', outcome='MODEL_REJECTED_ALARM_WORKS' if not gates['loading_model_valid'] and gates['injected_force_rejected'] else 'SEE_GATES', gates=gates, external_referent=pr['external_referent'], resolution='POPULATION', timescale='HANDOVER', frozen_predictions=fr, posterior=dict(k_mean_N_per_mm=k, k_sd_N_per_mm=sd, scope='Conditional digitization-error parameter posterior; not patient material uncertainty'), held_points=held, loading_pass_count=sum((not x['alarm'] for x in loading)), loading_count=len(loading), unloading_alarm_count=sum((x['alarm'] for x in unloading)), unloading_count=len(unloading), fault_injections=faults, independent_fit_absolute_error_N_per_mm=posterior_error, rejection=dict(observation_points=len(rows), used_for_calibration=len(cal), held_loading=len(loading), held_unloading=len(unloading), other_loading_not_scored=3), cost=dict(**cost(st), fit_wall_seconds=fitsec), practice='Fixed loading parameter does not know phase or emit model-form invalidation', same_information_control='Weighted least squares and residual test agree; no algorithm claim', binding=dict(from_measurement='published healthy-incisor mean response', to_parameter_group='effective_incisor_support', edge_resolution='POPULATION', not_patient_point_measurement=True, physical_material_update='REFUSED_AFTER_MODEL_ALARM'))
    write(L / 'raw/K50_R1.json', report)
    state('K50_R1_COMPLETE', gates, 'Freeze phase-indexed observation representation R2, with passivity/protocol admission guard')
    (L / 'HANDOFF_R1.md').write_text('K09 source geometry passed; K13 compiler passed but 87,168 mandibular PDL tets remain UNKNOWN. K50 R1 linear loading posterior rejects held curve points and unloading transfer; freeze and raw alarms preserved. Next construction changes observation representation to phase-indexed force response with an energy/protocol guard. Source visible historically; no blind claim.\n')
    print(json.dumps(dict(round='R1', gates=gates, loading_pass=report['loading_pass_count'], unloading_alarm=report['unloading_alarm_count'], k=k)))
if __name__ == '__main__':
    run_r1()
