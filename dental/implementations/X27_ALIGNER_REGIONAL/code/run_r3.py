import itertools, json, hashlib, time, resource
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from contact_model import make_system, solve
from calibration_port import fit, matrix, predict
R = Path(__file__).resolve().parents[1]

def put(n, x):
    (R / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def candidates(e, z, H):
    out = []
    for amp in [0.1, 0.2, 0.4]:
        for (label, vec) in [('buccal+', e), ('buccal-', -e), ('occlusal+', z), ('occlusal-', -z)]:
            out.append(dict(id=f'{label}{amp}', tooth_pose=np.r_[amp * vec, np.zeros(3)].tolist()))
    axis = np.cross(e, z)
    axis /= np.linalg.norm(axis)
    for amp in [0.2, 0.4]:
        for sign in [-1, 1]:
            out.append(dict(id=f'tipping{sign:+d}_edge{amp}', tooth_pose=np.r_[np.zeros(3), sign * amp / H * axis].tolist()))
    return out

def observation(case, system, probe, prediction):
    patches = system['patch'][:5]
    return dict(specimen_id='OUR_OWN_FIXTURE_' + case, tooth_id=11, source_kind='our_own_fixture', baseline_seated=True, force_N=prediction['wrenches']['11'][:3], moment_Nmm=prediction['wrenches']['11'][3:], moment_origin='crown_centroid', tooth_pose=probe['tooth_pose'], shell_housing_pose=prediction['housing_state'][:6], gap_mm=[p['gap_mm'] for p in patches], region_signatures=[p['signature'] for p in patches], elapsed_s=0.0, temperature_C=37.0, probe_id=probe['id'], optical_housing_rigidity_residual_mm=0.0, resolution_level='PER_TOOTH', contact_resolution_level='PER_SURFACE_REGION', observation_note='Simulated exact states. No actual optical measurement, no physical validation.')

def design(obs, H, ge):
    scaled = np.array([1, 1, 1, 1 / H, 1 / H, 1 / H])
    X = [matrix(r)[0] * scaled[:, None] for r in obs]
    examined = {}
    best = None
    for cardinality in range(1, 6):
        tested = 0
        full = 0
        best_at = None
        for subset in itertools.combinations(range(len(obs)), cardinality):
            S = np.vstack([X[i] for i in subset])
            sv = np.linalg.svd(S, compute_uv=False)
            tested += 1
            rank = 0 if sv[0] == 0 else int(np.sum(sv > sv[0] * ge['rank_relative_tolerance']))
            if rank != ge['full_rank']:
                continue
            condition = float(sv[0] / sv[-1])
            full += 1
            if condition > ge['condition_number_max']:
                continue
            if best_at is None or condition < best_at['condition_number']:
                best_at = dict(indices=list(subset), condition_number=condition, singular_values=sv.tolist(), rank=rank)
        examined[str(cardinality)] = dict(subsets_examined=tested, full_rank_subsets=full)
        if best_at is not None:
            best = best_at
            break
    return dict(selected=best, enumeration=examined, search_scope='Minimum cardinality among the frozen 16 finite probes; not a universal physical lower bound')

def main():
    start = time.perf_counter()
    pr = json.loads((R / 'PREREG_R3.json').read_text())
    ge = pr['gates']
    assert hashlib.sha256((R / 'PREREG_R3.json').read_bytes()).hexdigest() == json.loads((R / 'PREREG_R3.sha256.json').read_text())['sha256']
    arches = json.loads((R / 'inputs/ARCHES.json').read_text())
    regions = json.loads((R / 'inputs/PARK_REGIONS.json').read_text())
    pred = []
    systems = {}
    for arch in arches:
        system = make_system(arch, 11, regions)
        systems[arch['case']] = system
        e = np.array(arch['teeth']['11']['buccal_unit'])
        z = system['z']
        e -= z * (e @ z)
        e /= np.linalg.norm(e)
        H = arch['teeth']['11']['crown_height_mm']
        probes = candidates(e, z, H)
        query_start = time.perf_counter()
        predictions = [solve(system, tooth_pose=np.array(p['tooth_pose'])) for p in probes]
        pred.append(dict(case=arch['case'], probes=probes, predictions=predictions, H_mm=H, fixture_stiffness_N_per_mm=[p['k_patch_N_per_mm'] for p in system['patch'][:5]]))
    frozen = dict(frozen_utc=datetime.now(timezone.utc).isoformat(), prereg_sha256=hashlib.sha256((R / 'PREREG_R3.json').read_bytes()).hexdigest(), source_kind='our_own_fixture', predictions=pred)
    h = hashlib.sha256(json.dumps(pred, sort_keys=True).encode()).hexdigest()
    frozen['predictions_sha256'] = h
    if (R / 'FROZEN_PREDICTIONS_R3.json').exists():
        if json.loads((R / 'FROZEN_PREDICTIONS_R3.json').read_text())['predictions_sha256'] != h:
            raise RuntimeError('R3 prediction drift')
    else:
        put('FROZEN_PREDICTIONS_R3.json', frozen)
    outputs = []
    allobs = []
    for row in pred:
        case = row['case']
        system = systems[case]
        H = row['H_mm']
        observations = [observation(case, system, p, s) for (p, s) in zip(row['probes'], row['predictions'])]
        allobs.extend(observations)
        selection = design(observations, H, ge)
        if selection['selected'] is None:
            outputs.append(dict(case=case, outcome='UNIDENTIFIABLE_ON_FROZEN_PROBES', design=selection))
            continue
        ids = selection['selected']['indices']
        chosen = [observations[i] for i in ids]
        identification = fit(chosen, H, tol=ge['scaled_model_residual_N'])
        k = np.array(identification['stiffness_N_per_mm'])
        reference = np.array(row['fixture_stiffness_N_per_mm'])
        kerror = float(np.max(np.abs(k - reference) / reference))
        holdouts = []
        for (i, obs) in enumerate(observations):
            if i in ids:
                continue
            p = predict(obs, k)
            truth = np.r_[obs['force_N'], obs['moment_Nmm']]
            holdouts.append(dict(probe_id=obs['probe_id'], wrench_max_abs=float(np.max(np.abs(p - truth))), source_kind='our_own_fixture', physical_validation=False))
        single = fit([chosen[0]], H, tol=ge['scaled_model_residual_N'])
        force_only = np.vstack([matrix(chosen[0])[0][:3]])
        source_rank = int(np.linalg.matrix_rank(force_only, tol=1e-09))
        missing = [dict(r) for r in chosen]
        missing[0].pop('shell_housing_pose')
        no_pose = fit(missing, H)
        corrupt = [dict(r) for r in chosen]
        corrupt[0]['force_N'] = list(corrupt[0]['force_N'])
        corrupt[0]['force_N'][0] += 10
        injected = fit(corrupt, H, tol=ge['scaled_model_residual_N'])
        origin = [dict(r) for r in chosen]
        origin[0]['moment_origin'] = 'root_apex'
        wrong_origin = fit(origin, H)
        gates = dict(full_rank=identification['rank'] == 5, conditioning=identification['condition_number'] <= ge['condition_number_max'], effective_stiffness_recovery=kerror <= ge['effective_stiffness_relative_error_max'], unfit_fixture_poses=max((h['wrench_max_abs'] for h in holdouts)) <= ge['unfit_probe_wrench_max_abs'], injected_force_rejected=injected['status'] == 'CONTACT_MODEL_REFUTED', missing_pose_refused=no_pose['status'] == 'NEEDS_MATCHED_MEASUREMENT', unknown_origin_refused=wrong_origin['status'] == 'NEEDS_EXACT_MOMENT_TRANSPORT')
        outputs.append(dict(case=case, outcome='CONDITIONAL_REGIONAL_CALIBRATION_DESIGNED' if all(gates.values()) else 'CALIBRATION_CONSTRUCTION_FAILED', minimum_perturbation_queries=len(ids), additional_seated_baseline_queries=1, total_physical_queries=len(ids) + 1, selected_probe_ids=[r['probe_id'] for r in chosen], design=selection, identification=identification, single_six_axis_probe=single, force_only_single_probe_rank=source_rank, stiffness_relative_error_max=kerror, unfit_fixture_poses=holdouts, gates=gates, injected_result=injected, missing_pose_result=no_pose, unknown_origin_result=wrong_origin, resolution_level='PER_SURFACE_REGION', time_scale='SIMULTANEOUS', physical_validation='UNKNOWN'))
    put('raw/R3_FIXTURE_OBSERVATIONS.json', allobs)
    result = dict(round='R3', claim_type='capability', outcome='REGIONAL_POSE_AND_WRENCH_ACQUISITION_PORT_EXECUTED_PHYSICAL_VALIDATION_UNKNOWN', external_referent=pr['external_referent'], empirical_referent=pr['empirical_referent'], fixture_referent=pr['fixture_referent'], rows=outputs, all_fixture_gates_pass=all((all(r.get('gates', {}).values()) and r.get('gates') for r in outputs)), wall_s=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, physical_measurements=0, independent_study_validation=False, costs=dict(forward_queries=48, physical_queries_executed=0, planned_queries_per_specimen=[r.get('total_physical_queries') for r in outputs], tokens=None), next_operation='Acquire matched finished-shell gap/thickness and 3D optical pose plus six-axis wrench, baseline plus selected independent tooth motions, then independent intermediate/target holdout on same specimen; unknown start branch needs shape scan')
    put('round3/results.json', result)
    template = {k: None for k in ['specimen_id', 'tooth_id', 'force_N', 'moment_Nmm', 'moment_origin', 'tooth_pose', 'shell_housing_pose', 'gap_mm', 'region_signatures', 'source_kind', 'baseline_seated', 'elapsed_s', 'temperature_C', 'optical_housing_rigidity_residual_mm']}
    template['source_kind'] = 'independent_measurement'
    put('inputs/MATCHED_LAB_MEASUREMENT_TEMPLATE.json', [template])
    (R / 'round3/HANDOFF.md').write_text('# R3 acquisition handoff\n\nHousing pose has become a measured input in a runnable regional inverse. Selected probe counts are model-conditional; unused pose predictions and recovery are OUR_OWN_FIXTURE, not empirical validation.\n\n' + '\n'.join((f"{r['case']}: {r.get('minimum_perturbation_queries')} independent perturbations + one seated baseline; {r.get('outcome')}." for r in outputs)) + '\n\n' + result['next_operation'] + '\n')
    put('CURRENT_WORK_STATE.json', dict(lane='X27-aligner-regional', state='R3_COMPLETE_PACKAGING', latest_gate=result['outcome'], next_operation='Verify inherited continuous-gap port, corruption checks and one-command replay; bind scoped coverage proposal'))
    print(json.dumps(dict(outcome=result['outcome'], rows=[dict(case=r['case'], minimum_probes=r.get('minimum_perturbation_queries'), selected=r.get('selected_probe_ids'), gates=r.get('gates')) for r in outputs]), indent=2))
if __name__ == '__main__':
    main()
