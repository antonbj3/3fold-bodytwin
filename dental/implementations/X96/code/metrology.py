from common import *
from full_geometry import cylinder_triangles, write_stl
pr = dict(round='R4', claim_type='information_link', capability='A lab can test frozen digital tool-envelope gaps in one registered case without refitting expected numbers', obstacle='Physical tool profile, registration and signed displacement not measured; point estimate from guide moments is insufficient', changed_operation='Freeze expected gap targets for first lexicographic case/FDI, export exact research tool cylinders and source voxel centers; typed metrology input accepts only matched units/frame/geometry', consumer='prototyping lab tool-profile and registered phantom measurement', selection='first lexicographic X8 case+FDI, all3 verified depth protocols, no outcome-based selection', metric=dict(absolute_gap_tolerance_mm=0.3, justification='predeclared bench pilot tolerance equal to one source voxel; no anatomical/clinical accuracy claim', exact_frame='TF2_INDEX_XYZ_MM', requires='independent measurement locator, same fixture/geometry hash and protocol'), decision='empty input => UNKNOWN; measured gap outside frozen +/-0.3mm => FAIL; no refit', strongest_equally_informed_control='Direct tolerance and schema validation; injected +1mm, wrong frame/unit/geometry each fail', falsifier='wrong data accepted, frozen expected value changes, or independent registered physical gap disagrees', external_referent=dict(kind='published_dataset', locator='PROTOCOLS.json and public TF2 voxel occupancy', compared_quantity='digital target for future independent phantom metrology; physical measurement not yet performed', refutes_us=True), full_cost=dict(preparation='existing measured digital source and protocols', fit='none', discovery='no new biological model', validation='three injected contract faults plus gap fault', queries='3 prospective lab targets', fallback='UNKNOWN no measurement'), resolution=dict(targets='PER_TOOTH', bench_tolerance='PHENOMENOLOGICAL', source='PER_POINT'), time_scale='SIMULTANEOUS')
freeze('PREREG_R4.json', pr)
rows = [json.loads(l) for l in (DATA / 'R2_GEOMETRY.jsonl').read_text().splitlines()]
pair = min(((r['case'], r['fdi']) for r in rows))
selected = [r for r in rows if (r['case'], r['fdi']) == pair]
protocols = load(ROOT / 'PROTOCOLS.json')['protocols']
targets = []
export = ROOT / 'export'
export.mkdir(exist_ok=True)
for p in protocols:
    r = next((r for r in selected if r['extra_mm'] == p['extra_mm']))
    f = export / (p['id'] + '.stl')
    tri = cylinder_triangles(np.array(r['entry_zyx_mm'])[::-1], np.array(r['axis_zyx'])[::-1], r['length_mm'] + r['extra_mm'], r['radius_mm'], 96)
    write_stl(f, tri)
    targets.append(dict(case=r['case'], fdi=r['fdi'], protocol=p['id'], extra_mm=r['extra_mm'], expected_gap_interval_mm=[r['union_lower_mm'], r['union_upper_mm']], unit='mm', frame='TF2_INDEX_XYZ_MM', geometry_sha256=sha(f), mesh_path=str(f), radial_polygon_sag_mm=r['radius_mm'] * (1 - np.cos(np.pi / 96)), polygon_sag_is_not_total_error=True, source_geometry_sha256=sha(DATA / 'R2_GEOMETRY.jsonl'), resolution='PER_TOOTH'))
import zipfile
first = selected[0]
if first['paired_revisions']:
    arr = []
    for src in first['point_sources']:
        with np.load(src['path']) as a:
            arr.extend((a[k].astype(float) * 0.3 for k in a.files if k.endswith('_voxels_zyx')))
    pts = np.unique(np.vstack(arr), axis=0)
else:
    with zipfile.ZipFile(tf2_io.ZIP) as z:
        (lab, sp, _) = tf2_io.read_mha_bytes(z.read(first['point_sources'][0]['member']))
        pts = np.argwhere(np.isin(lab, [3, 4])) * 0.3
np.savetxt(export / 'canal_voxel_centers_xyz_mm.csv', pts[:, ::-1], delimiter=',', header='x_mm,y_mm,z_mm', comments='', fmt='%.4f')
payload = dict(schema='X96-prospective-metrology-v1', targets=targets, absolute_gap_tolerance_mm=0.3, source_points_sha256=sha(export / 'canal_voxel_centers_xyz_mm.csv'), voxel_box_halfwidth_mm=0.15, unit='mm', frame='TF2_INDEX_XYZ_MM', claim='Expected digital full-radius research-envelope gap; actual bit shape and physical performance unvalidated', physical_measurements_performed=False)
f = ROOT / 'FROZEN_PREDICTIONS.json'
if f.exists():
    old = load(f)
    import copy
    original = copy.deepcopy(old)
    replayed = copy.deepcopy(payload)
    bindings = []
    for (prior, current) in zip(original['targets'], replayed['targets']):
        assert Path(prior['mesh_path']).name == Path(current['mesh_path']).name, 'Frozen mesh identity changed'
        assert sha(current['mesh_path']) == prior['geometry_sha256'], 'Frozen mesh bytes changed'
        bindings.append(dict(frozen_locator=prior['mesh_path'], replay_locator=current['mesh_path'], sha256=prior['geometry_sha256']))
        prior.pop('mesh_path')
        current.pop('mesh_path')
    assert all((original[k] == v for (k, v) in replayed.items())), 'Frozen scientific predictions changed'
    dump(ROOT / 'raw/REPLAY_PATH_BINDING.json', dict(bindings=bindings, frozen_predictions_sha256=sha(f), predictions_replaced=False))
else:
    payload['frozen_utc'] = now()
    dump(f, payload)
    f.with_suffix('.json.sha256').write_text(sha(f) + '\n')

def evaluate(m):
    if not m:
        return dict(status='UNKNOWN', reason='NO_INDEPENDENT_MEASUREMENT')
    verdict = []
    for r in m:
        t = next((t for t in targets if (r.get('case'), r.get('fdi'), r.get('protocol')) == (t['case'], t['fdi'], t['protocol'])), None)
        match = t is not None and all((r.get(k) == t[k] for k in ['unit', 'frame', 'geometry_sha256'])) and bool(r.get('independent_measurement_locator'))
        v = r.get('measured_gap_mm')
        valid = match and isinstance(v, (int, float)) and np.isfinite(v) and (v >= 0)
        passed = bool(valid and t['expected_gap_interval_mm'][0] - 0.3 <= v <= t['expected_gap_interval_mm'][1] + 0.3)
        verdict.append(dict(protocol=r.get('protocol'), valid_contract=bool(valid), passes_frozen_gap=passed))
    return dict(status='PASS' if all((r['passes_frozen_gap'] for r in verdict)) else 'FAIL', results=verdict)
t = targets[0]
fixture = {k: t[k] for k in ['case', 'fdi', 'protocol', 'unit', 'frame', 'geometry_sha256']}
fixture['independent_measurement_locator'] = 'OUR_OWN_INJECTED_FIXTURE_NOT_MEASUREMENT'
fixture['measured_gap_mm'] = float(np.mean(t['expected_gap_interval_mm']))
faults = []
for (key, val) in [('measured_gap_mm', fixture['measured_gap_mm'] + 1), ('unit', 'm'), ('frame', 'OTHER_FRAME'), ('geometry_sha256', '0' * 64)]:
    wrong = dict(fixture)
    wrong[key] = val
    faults.append(dict(injected_field=key, rejected=evaluate([wrong])['status'] == 'FAIL', external_referent_kind='our_own_fixture'))
if not (ROOT / 'raw/METROLOGY_INPUT.json').exists():
    dump(ROOT / 'raw/METROLOGY_INPUT.json', dict(measurements=[]))
result = evaluate(load(ROOT / 'raw/METROLOGY_INPUT.json')['measurements'])
dump(ROOT / 'raw/R4_METROLOGY.json', dict(faults=faults, measurement=result, schema_fixture_pass=evaluate([fixture])['status'] == 'PASS', actual_physical_reference='UNKNOWN', frozen_predictions_sha256=sha(f)))
e = np.zeros(3)
a = np.array([1.0, 0, 0])
pt = np.array([[5.0, 4.5, 0]])
from full_geometry import project_cylinder
gaps = [float(np.linalg.norm(pt - project_cylinder(pt, e + sign * np.array([0, 1.0, 0]), a, 10, 2))) for sign in [-1, 1]]
dump(ROOT / 'raw/SUFFICIENCY_R4_GUIDE.json', dict(summary=dict(entry_magnitude_mean_mm=[1.0, 1.0], apex_magnitude_mean_mm=[1.0, 1.0], magnitude_sd_mm=[0.0, 0.0], angle_deg=[0.0, 0.0]), identity_error_mm=0.0, downstream_gap_mm=gaps, downstream_difference_mm=abs(gaps[0] - gaps[1]), breach_fraction=[float(v < 2) for v in gaps], smallest_extension='signed achieved tool-to-canal clearance loss in registered frame', external_referent=dict(kind='our_own_fixture', locator='code/metrology.py', compared_quantity='guide moment sufficiency only', refutes_us=True)))
dump(ROOT / 'rounds/R4.json', dict(round='R4', claim_type='information_link', outcome='FROZEN_METROLOGY_PORT_EXECUTED_PHYSICAL_VALIDATION_UNKNOWN', gates={'fault_controls': all((r['rejected'] for r in faults)), 'physical_validation': result['status']}, external_referent=pr['external_referent'], targets=3, cost={'wall_time': 'included in run_all cost', 'manual_physical_measurement': 'NOT_PERFORMED'}))
(ROOT / 'HANDOFF_R4.md').write_text("R4 : frozen digital predictions and research STL is available in export/. Tom METROLOGY_INPUT.json gives UNKNOWN . Four injected errors fold. No physical measurements has been performed. The minimum next step is to measure actual tool profile and signed registered pose; outputs should be tested without readjustment to FROZEN_PREDICTIONS.json.\n")
state('R4_COMPLETE', 'METROLOGY_UNKNOWN', 'Finalize demo and pending-review feedback; next experiment requires independent tool metrology')
