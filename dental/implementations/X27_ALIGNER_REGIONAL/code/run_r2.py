import json, hashlib, time, sys
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R.parents[1] / 'results/LANE_X19_ALIGNER_FORCE/code'))
from model import response

def put(n, x):
    (R / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')
COMP = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']

def main():
    start = time.perf_counter()
    pr = json.loads((R / 'PREREG_R2.json').read_text())
    assert hashlib.sha256((R / 'PREREG_R2.json').read_bytes()).hexdigest() == json.loads((R / 'PREREG_R2.sha256.json').read_text())['sha256']
    source = json.loads((R / 'inputs/CHO_SENSOR.json').read_text())
    pred = []
    controls = []
    for material in ['Zendura', 'Trioclear', 'Graphy']:
        for fdi in [12, 11]:
            calibration = [next((r for r in source if r['material'] == material and r['fdi'] == fdi and (r['activation_mm'] == 0.3) and (r['component'] == c))) for c in COMP]
            vector = np.array([r['mean'] for r in calibration])
            prediction = vector * 2
            tangent = np.linalg.lstsq(np.array([[0.3]]), vector[None, :], rcond=None)[0][0]
            controls.append(float(np.max(np.abs(tangent * 0.6 - prediction))))
            pred.append(dict(material=material, fdi=fdi, activation_mm=0.6, calibration_wrench=vector.tolist(), predicted_wrench=prediction.tolist(), calibration_locators=[r['locator'] for r in calibration], moment_origin=calibration[0]['moment_origin'], force_units='N', moment_units='Nmm', resolution_level='PER_TOOTH', time_scale='SIMULTANEOUS'))
    frozen = dict(frozen_utc=datetime.now(timezone.utc).isoformat(), prereg_sha256=hashlib.sha256((R / 'PREREG_R2.json').read_bytes()).hexdigest(), source_sha256=hashlib.sha256((R / 'inputs/CHO_SENSOR.json').read_bytes()).hexdigest(), predictions=pred)
    payload_hash = hashlib.sha256(json.dumps(pred, sort_keys=True).encode()).hexdigest()
    frozen['predictions_sha256'] = payload_hash
    if (R / 'FROZEN_PREDICTIONS_R2.json').exists():
        if json.loads((R / 'FROZEN_PREDICTIONS_R2.json').read_text())['predictions_sha256'] != payload_hash:
            raise RuntimeError('R2 drift')
    else:
        put('FROZEN_PREDICTIONS_R2.json', frozen)
    ge = pr['gates']
    scores = []
    injections = []
    for row in pred:
        ref = [next((r for r in source if r['material'] == row['material'] and r['fdi'] == row['fdi'] and (r['activation_mm'] == 0.6) and (r['component'] == c))) for c in COMP]
        truth = np.array([r['mean'] for r in ref])
        prediction = np.array(row['predicted_wrench'])
        sd = np.array([r['sd'] for r in ref])
        fe = float(np.linalg.norm(prediction[:3] - truth[:3]) / max(np.linalg.norm(truth[:3]), ge['force_norm_floor_N']))
        me = float(np.linalg.norm(prediction[3:] - truth[3:]) / max(np.linalg.norm(truth[3:]), ge['moment_norm_floor_Nmm']))
        scores.append(dict(**row, external_wrench=truth.tolist(), external_sd=sd.tolist(), force_vector_error=fe, moment_vector_error=me, force_gate=fe <= ge['force_vector_relative_error_max'], moment_gate=me <= ge['moment_vector_relative_error_max'], external_locators=[r['locator'] for r in ref], source_quantity='Group means of baseline-referenced sensor wrench; not specimen-matched contact'))
        corrupted = prediction.copy()
        corrupted[2] += 10
        corrupted[5] += 100
        cfe = float(np.linalg.norm(corrupted[:3] - truth[:3]) / max(np.linalg.norm(truth[:3]), ge['force_norm_floor_N']))
        cme = float(np.linalg.norm(corrupted[3:] - truth[3:]) / max(np.linalg.norm(truth[3:]), ge['moment_norm_floor_Nmm']))
        injections.append(dict(material=row['material'], fdi=row['fdi'], injected_force_error=cfe, injected_moment_error=cme, force_rejected=cfe > ge['force_vector_relative_error_max'], moment_rejected=cme > ge['moment_vector_relative_error_max']))
    arches = json.loads((R / 'inputs/ARCHES.json').read_text())
    practice = []
    for fdi in [12, 11]:
        values = []
        for a in arches:
            w = response(a, 12, -0.6, 2746.0, 0.75)['wrenches'].get(str(fdi), [0] * 6)
            F = np.array(w[:3])
            e = np.array(a['teeth'][str(fdi)]['buccal_unit'])
            z = np.array(a['frame']['occlusal_unit'])
            x = np.cross(e, z)
            x /= np.linalg.norm(x)
            values.append([float(F @ x), float(F @ e), float(F @ z)])
        vec = np.median(np.array(values), axis=0)
        for material in ['Zendura', 'Trioclear', 'Graphy']:
            truth = np.array(next((r['external_wrench'] for r in scores if r['material'] == material and r['fdi'] == fdi))[:3])
            practice.append(dict(material=material, fdi=fdi, force_N=vec.tolist(), vector_error=float(np.linalg.norm(vec - truth) / max(np.linalg.norm(truth), ge['force_norm_floor_N'])), scope='Nominal .75mm common-material ribbon proxy; translation used for source tipping because source root/geometry missing; no apex moment comparison'))
    gates = dict(force=all((r['force_gate'] for r in scores)), moment=all((r['moment_gate'] for r in scores)), numerical_control=max(controls) <= ge['numerical_prediction_parity_abs'], injected_errors_rejected=all((r['force_rejected'] and r['moment_rejected'] for r in injections)))
    result = dict(round='R2', claim_type='information_link', outcome='ONE_SEATED_CALIBRATION_SUFFICIENT' if gates['force'] and gates['moment'] else 'ONE_SEATED_CALIBRATION_INSUFFICIENT', external_referent=pr['external_referent'], rows=scores, practice_proxy=practice, gates=gates, injections=injections, force_groups_passing=sum((r['force_gate'] for r in scores)), moment_groups_passing=sum((r['moment_gate'] for r in scores)), total_groups=len(scores), acquisition_status='Unfit known group means, not independent held study; no new physical measurement', resolution_level='PER_TOOTH', time_scale='SIMULTANEOUS', wall_s=time.perf_counter() - start, negative_result=not (gates['force'] and gates['moment']), next_operation='Replace scalar/tangent force calibration by regional contact observability: measure shell housing pose plus six-axis wrench under multiple independent perturbations; design smallest full-rank probe set')
    put('round2/results.json', result)
    (R / 'round2/HANDOFF.md').write_text(f"# R2 handoff\n\nOne source-conditioned seated tangent calibration: {result['outcome']}. Force groups passing: {result['force_groups_passing']}/6; moment groups: {result['moment_groups_passing']}/6.\n\nSource: doi:10.4041/kjod25.003 T4/T5; .3 mm calibration and unfit .6 mm groups. Endpoint means were known; no blinded or independent-study validation. Different nominally identical shells can have changed active contacts.\n\nNext: {result['next_operation']}\n")
    put('CURRENT_WORK_STATE.json', dict(lane='X27-aligner-regional', state='R2_COMPLETE_R3_STARTING', latest_gate=result['outcome'], next_operation=result['next_operation']))
    print(json.dumps(dict(outcome=result['outcome'], force_pass=result['force_groups_passing'], moment_pass=result['moment_groups_passing'], rows=[dict(material=r['material'], fdi=r['fdi'], force_error=r['force_vector_error'], moment_error=r['moment_vector_error']) for r in scores]), indent=2))
if __name__ == '__main__':
    main()
