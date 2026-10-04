"""Assemble source-supported ports and explicit missing physical predictions."""
import collections, csv, datetime, hashlib, json, pathlib, time
import numpy as np
P = pathlib.Path(__file__).resolve().parents[1]
D = P.parent.parent

def put(n, x):
    (P / n).write_text(json.dumps(x, indent=2, allow_nan=False, ensure_ascii=False) + '\n')

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def run():
    rounds = {k: json.loads((P / f'results_{k}.json').read_text()) for k in ['R1', 'R2', 'R3']}
    rows = [json.loads(l) for l in (P / 'raw/PER_SITE_R2.jsonl').open()]
    old = {(r['case'], r['site']): r for r in map(json.loads, (P / 'raw/PER_SITE_R1.jsonl').read_text().splitlines())}
    suff = json.loads((P / 'SUMMARY_SUFFICIENCY.json').read_text())
    ports = []
    table = []
    rejections = collections.Counter()
    angle = []
    delta = []
    core_exclusions = []
    ray_level = []
    for r in rows:
        if r['status'] != 'MEASURED':
            continue
        q = r['rays'][1]
        core = r['trabecular_candidate_roi']
        prev = old[r['case'], r['site']]
        if q['envelope_chord_mm'] is not None:
            delta.append(q['envelope_chord_mm'] - prev['rays'][1]['envelope_chord_mm'])
            b1 = np.array(prev['provisional_buccal_axis_zyx'])
            b2 = np.array(r['provisional_buccal_axis_zyx'])
            angle.append(float(np.degrees(np.arccos(np.clip(abs(b1 @ b2), 0, 1)))))
        for qq in r['rays']:
            for (side, v) in qq['profiles'].items():
                rejections[v['rejection_reason'] or 'RETAINED_CONTRAST_ONLY'] += 1
                ray_level.append({'case': r['case'], 'site': r['site'], 'depth_mm': qq['depth_mm'], 'side_provisional': side, 'extent_mm': qq[side + '_extent_mm'], 'apparent_layer_mm_tau_04_05_06': v['apparent_layer_mm_tau_04_05_06'], 'rejection_reason': v['rejection_reason'], 'anatomical_cortex_mm': None, 'resolution': 'PER_POINT', 'raw_file': r['raw_file']})
        if core['n_voxels'] < 30:
            core_exclusions.append({'case': r['case'], 'site': r['site'], 'reason': 'LESS_THAN_30_INTERIOR_VOXELS', 'n_voxels': core['n_voxels']})
        import sys
        torque_dir = D / 'results/LANE_X68_INSERTION_TORQUE/code'
        if str(torque_dir) not in sys.path:
            sys.path.insert(0, str(torque_dir))
        import query as torque
        torque_result = torque.query({'same_regime_confirmed': False})
        if torque_result['status'] != 'UNKNOWN_REGIME_COMPATIBILITY':
            raise ValueError('X68 no longer rejects uncalibrated site transfer')
        thermal_dir = D / 'results/LANE_X67_DRILL_HEAT/code'
        if str(thermal_dir) not in sys.path:
            sys.path.insert(0, str(thermal_dir))
        import bench_calibrate
        try:
            bench_calibrate.check_metadata({'observer_model': 'CBCT_gray_profile'}, {'observer_model': 'CBCT_gray_profile'})
            raise RuntimeError('X67 admitted CBCT contrast as a thermal observation')
        except ValueError as exc:
            thermal_rejection = str(exc)
        port = {'case': r['case'], 'site': r['site'], 'resolution': 'PER_TOOTH', 'spatial_input_resolution': 'PER_POINT', 'time_relation': 'SIMULTANEOUS', 'source_pose_type': r['kind'], 'source_pose_preoperative_intent': 'UNKNOWN' if r['kind'] == 'J_postoperative' else 'existing algorithmic dentate plan', 'source_geometry': {'local_envelope_chords': r['rays'], 'axial_envelope_chord_mm': r['axial_envelope_chord_mm'], 'orientation': r['orientation']}, 'gray_input': core, 'spatial_raw': r['raw_file'], 'cortex_mm': {'crest': None, 'buccal': None, 'lingual': None}, 'crest_missing_reason': 'axis passes tooth or postoperative implant; no independently labeled inner cortex; side-ray bright layer cannot fill this port', 'density_g_cm3': None, 'density_calibration_error': 'UNBOUNDED', 'physical_error_enclosure': 'MISSING', 'X68': {'actual_operator_result': torque_result, 'input_status': 'UNKNOWN_LOCAL_MECHANICS_CALIBRATION', 'prediction_Ncm': None, 'explained_variance': None, 'final_bore_mm': None}, 'X67': {'actual_operator_admission': 'REJECTED', 'actual_operator_reason': thermal_rejection, 'input_status': 'UNKNOWN_INDEPENDENT_CAPACITY_AND_TRANSPORT', 'C1_C2_J_K': None, 'rho_cp_J_mm3_K': None, 'mechanical_work_J': None, 'temperature_C': None, 'explained_variance': None, 'drill_protocol': None}, 'debts': [{'quantity': 'apparent_bright_layer_is_cortex', 'resolution': 'PHENOMENOLOGICAL', 'replacement': 'registered inner/outer cortical boundary from independent microCT/histology'}, {'quantity': 'gray_to_density', 'resolution': 'PHENOMENOLOGICAL', 'replacement': 'same-acquisition independent low-density certified phantom with held-out insert/spatial positions'}, {'quantity': 'density_to_torque_heat', 'resolution': 'PHENOMENOLOGICAL', 'replacement': 'same-site material/mechanical response, actual drilling work, bore/thread geometry and regional sensor histories'}]}
        ports.append(port)
        gv = core['gray_quantiles_05_25_50_75_95']
        table.append({'case': r['case'], 'site': r['site'], 'kind': r['kind'], 'orientation_status': r['orientation']['status'], 'depth_mm': 4, 'envelope_chord_mm': q['envelope_chord_mm'], 'provisional_buccal_extent_mm': q['buccal_extent_mm'], 'provisional_lingual_extent_mm': q['lingual_extent_mm'], 'axial_envelope_chord_mm': r['axial_envelope_chord_mm'], 'core_n_voxels': core['n_voxels'], 'core_gray_q05': gv[0] if gv else None, 'core_gray_q50': gv[2] if gv else None, 'core_gray_q95': gv[4] if gv else None, 'anatomical_cortex_mm': 'UNKNOWN', 'density_g_cm3': 'UNKNOWN', 'torque_Ncm': 'UNKNOWN', 'temperature_C': 'UNKNOWN'})
    with (P / 'PER_SITE_MEASUREMENTS.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)
    (P / 'SITE_INPUT_PORTS.jsonl').write_text(''.join((json.dumps(v, allow_nan=False) + '\n' for v in ports)))
    (P / 'PER_POINT_LAYERS.jsonl').write_text(''.join((json.dumps(v, allow_nan=False) + '\n' for v in ray_level)))
    put('EXCLUSIONS.json', {'sites_retained': len(ports), 'profiles_total': len(ray_level), 'profile_reasons': dict(rejections), 'rejected_fraction': 1 - rejections['RETAINED_CONTRAST_ONLY'] / len(ray_level), 'interior_roi_rejections': core_exclusions, 'orientation_fallbacks': [{'case': r['case'], 'site': r['site'], 'reason': r.get('orientation')} for r in rows if r.get('orientation', {}).get('status') != 'LOCAL_ARCH_PROXY'], 'original_dataset_screen': {'candidate_cases': 24, 'original_counterparts': rounds['R3']['candidate_pairs'], 'rejected': rounds['R3']['rejections'], 'selected_pairs': 8, 'deferred_pairs': 13}, 'cortical_accuracy_rejections': {'candidate_profiles': 570, 'physical_validation_pairs': 0, 'reason': 'source has no independently measured inner boundary'}, 'response_admission': {'sites_considered': 95, 'sites_admitted': 0, 'reason': 'no matched density/mechanics/temperature truth'}})
    volume_bytes = sum((f['bytes'] for k in ['R1', 'R2'] for f in rounds[k]['data_files']))
    out = {'lane': 'X75-cbct-bone', 'claim_type': 'information_link', 'status': 'PARTIAL_OBSERVED_SITE_STATE_PHYSICAL_BONE_QUALITY_UNKNOWN', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'capability_delivered': 'runnable site geometry/gray profile extraction and spatially resolved proposed consumer ports; original-payload lineage comparison; enclosed mathematical summary test', 'external_referent': rounds['R2']['external_referent'], 'external_referents': [rounds['R2']['external_referent'], rounds['R3']['external_referent'], suff['heat']['external_referent']], 'rounds': rounds, 'sufficiency': suff, 'source_supported': {'sites': 95, 'case_files': 24, 'new_independent_density_refs': 0, 'interior_rois_observed': 94, 'arch_normal_sites': 94, 'gray_resolution': 'PER_POINT', 'pose_site_resolution': 'PER_TOOTH', 'population_proxy_resolution': 'POPULATION'}, 'descriptive_not_frozen_gate': {'median_R2_minus_R1_chord_mm': float(np.median(delta)), 'median_R1_R2_plane_angle_deg': float(np.median(angle))}, 'clinical_torque_temperature_protocol': 'UNKNOWN', 'variance_explained_torque': None, 'variance_explained_temperature': None, 'physical_cortex_density_accuracy': 'UNKNOWN', 'control_outcome': 'source ray/box and scalar gray replay agree within frozen tolerances; generic high-precision heat convolution agrees; no algorithm superiority claimed', 'uncertainty': {'source_numeric_geometry_mm': 'max3.20e-14 agreement with independent voxel-box enumeration; not anatomical boundary error', 'voxel_spacing_mm': 0.3, 'ray_sampling_step_mm': 0.15, 'subvoxel_sampling_caveat': 'interpolation does not improve original .3mm resolution', 'true_inner_cortex_density': 'UNBOUNDED', 'true_anatomical_plane': 'UNKNOWN', 'heat_model_numeric': 'outward Gaussian integral bounds in SUMMARY_SUFFICIENCY.json', 'physical_heat_remainder': 'MISSING'}, 'full_cost': {'numeric_wall_s_sum': sum((rounds[k].get('full_cost', rounds[k].get('cost', {}))['wall_s'] for k in rounds)) + suff['cost_wall_s'], 'preparation_discovery_wall_s': 'NOT_INSTRUMENTED; directed source/file/primary literature reads, no agents', 'fit': 'one fixed training mean per geometric representation; all thresholds frozen, no material fit', 'validation': '665 direct extraction controls per imaging round, all voxels in8original pairs, exact-summary and interval checks, injected corruption; no independent physical accuracy measurement', 'queries': 95, 'fallback_sites': 95, 'questions': 0, 'model_selection': 'adaptive R2 after R1, adaptive R3 after R2; no blind clinical validation', 'peak_rss_kib': max((rounds[k].get('full_cost', rounds[k].get('cost', {}))['peak_rss_kib'] for k in rounds)), 'data_intermediates_bytes': volume_bytes, 'max_per_lane_bytes': 3000000000, 'threads': 4, 'GPU': False, 'physical_measurements': 0}, 'limitations': ['case-disjoint IDs are not independently verified biological disjointness across source families', 'X8 poses are dentate planning candidates, not clinical implant measurements', 'J poses are postoperative and metal-contaminated, not preoperative intended cortex', 'tooth-filled envelope chords are geometry of an annotation closure, not solid bone width/height or clinically usable dimensions', 'no crestal inner boundary and only24 retained lateral contrast layers; no new cortex truth', 'original arrays8/8 identical and do not add independent calibration truth', 'source-to-physical moment/heat transfer missing', 'geometric reflection witness was exploratory after R2 freeze; only heat-histogram witness had predeclared setup/gates'], 'edges': [{'from': 'PER_POINT_LAYERS.jsonl', 'to': 'X68 proposed regional cortex/contact input', 'resolution': 'PER_POINT', 'timescale': 'SIMULTANEOUS', 'status': 'OBSERVED_CONTRAST_ONLY_MATERIAL_UNKNOWN'}, {'from': 'SITE_INPUT_PORTS.jsonl', 'to': 'X67 regional thermal capacity/material input', 'resolution': 'PER_POINT', 'timescale': 'SIMULTANEOUS', 'status': 'CAPACITY_UNIDENTIFIED'}, {'from': 'drilling temperature history', 'to': 'residual regional injury/osseointegration', 'resolution': 'PER_POINT', 'timescale': 'HANDOVER', 'status': 'UNKNOWN_SEPARATE_EMPIRICAL_LAW'}], 'large_array_manifest': [f for k in ['R1', 'R2'] for f in rounds[k]['data_files']], 'next_construction': 'independent registered ex-vivo cortex/low-density reference and source-to-mechanics/thermal test, not further fitting exposed CBCT gray'}
    put('results.json', out)
    put('CURRENT_WORK_STATE.json', {'lane': 'X75-cbct-bone', 'phase': 'THREE_ROUNDS_DECIDED_DELIVERY_VERIFY', 'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'latest_gate': '95 source geometry profiles; cortex/density/response UNKNOWN; summary tests refute scalars;8original pairs exact shared payload', 'next_operation': out['next_construction'], 'physical_measurement_needed': True, 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
    print('assembled', len(ports), 'ports; intermediates', volume_bytes, 'bytes')
if __name__ == '__main__':
    run()
