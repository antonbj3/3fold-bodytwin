"""Verify frozen inputs and replay representative positive approach witnesses."""
import csv, sys
from common import *
from collision import *
from run_r1 import lib
from compare_cam import compare
from compare_stock import compare as compare_stock

def run():
    source_lock = read(ROOT / 'sources/PROXIMITY_SOURCE_LOCK.json')
    assert source_lock['same_bytes'] and sha(ROOT / 'sources/trimesh_proximity_4.11.3.py') == source_lock['published_source_sha256']
    assert sha(source_lock['installed_file']) == source_lock['installed_sha256']
    r1 = read(ROOT / 'RESULTS_R1.json')
    manifest = {r['id']: r for r in read(ROOT / 'INPUT_MANIFEST.json')['rows']}
    checked = 0
    for (file, digest0) in r1['raw_files'].items():
        assert sha(file) == digest0
        r = read(file)
        src = manifest[r['id']]
        assert sha(src['file']) == src['sha256']
        a = np.load(src['file'], allow_pickle=False)
        scene = Scene(a['vertices'], a['faces'])
        for v in r['points_and_tools']:
            for region in ['intaglio', 'exterior_and_rim']:
                positive = next((p for p in v['records'] if p['region'] == region and p['smallest_found_green_mm'] is not None), None)
                if positive is None:
                    continue
                t0 = next((t for t in positive['tools'] if t['status'] == 'FOUND'))
                tool = next((t for t in lib(v['library'], v['axes']) if t['id'] == t0['tool_id']))
                approach = float(np.linalg.norm(scene.bounds[1] - scene.bounds[0]) + 2 * tool['holder_diameter_mm'] + 2.0)
                assert pose_clearance(scene, np.asarray(t0['centre']), np.asarray(positive['normal']), np.asarray(t0['direction']), tool, approach)[0]
                checked += 1
    case = read(ROOT / 'LAB_CASES_20.json')['cases'][0]
    pred = read(case['prediction_file'])
    v = next((v for v in pred['points_and_tools'] if v['library'] == 'vhf' and v['axes'] == 5))
    index = next((i for (i, p) in enumerate(v['records']) if p['smallest_found_green_mm'] is not None))
    dest = ROOT / 'raw/INJECTED_CAM_OBSERVATION.csv'
    fields = ['case_id', 'library', 'axes', 'point_index', 'mesh_sha256', 'tool_library_sha256', 'observed_reachable', 'measured_assembly_and_frame', 'cam_locator']
    with open(dest, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow(dict(case_id=case['case_id'], library='vhf', axes=5, point_index=index, mesh_sha256=case['mesh_sha256'], tool_library_sha256=sha(ROOT / 'TOOL_LIBRARY.json'), observed_reachable='false', measured_assembly_and_frame='true', cam_locator='our_own_fixture: deliberately false reachability'))
    corrupted = compare(dest)
    assert corrupted['physical_validation'] == 'FAIL'
    missing = compare(ROOT / 'CAM_LOG_TEMPLATE.csv')
    assert missing['physical_validation'] == 'UNKNOWN' and sum(missing['counts'].values()) == 20
    r5 = read(DATA / 'R5' / (case['design_id'] + '.json'))
    index = next((i for (i, p) in enumerate(r5['records']) if p['local_component_stock_um'] is not None))
    fields = ['case_id', 'point_index', 'mesh_sha256', 'tool_library_sha256', 'frozen_predictions_sha256', 'observed_stock_um', 'measurement_bound_um', 'scene_conditioning_verified', 'measurement_locator', 'measurement_kind']
    injected_stock = ROOT / 'raw/INJECTED_STOCK_OBSERVATION.csv'
    with open(injected_stock, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow(dict(case_id=case['case_id'], point_index=index, mesh_sha256=case['mesh_sha256'], tool_library_sha256=sha(ROOT / 'TOOL_LIBRARY.json'), frozen_predictions_sha256=sha(ROOT / 'FROZEN_PREDICTIONS_R5.json'), observed_stock_um=r5['records'][index]['local_component_stock_um'] + 50, measurement_bound_um=1, scene_conditioning_verified='true', measurement_locator='our_own_fixture: injected50um', measurement_kind='our_own_fixture'))
    assert compare_stock(injected_stock)['scenario_transport'] == 'FAIL'
    empty_stock = compare_stock(ROOT / 'STOCK_LOG_TEMPLATE.csv')
    assert empty_stock['scenario_transport'] == 'UNKNOWN' and sum(empty_stock['counts'].values()) == 20
    freezes = []
    for path in ROOT.glob('FROZEN_PREDICTIONS*.json'):
        x = read(path)
        assert digest(x['payload']) == x['payload_sha256']
        assert sha(path) == path.with_suffix('.sha256').read_text().strip()
        freezes.append(str(path))
    for (name, result) in [('FROZEN_PREDICTIONS.json', 'RESULTS_R1.json'), ('FROZEN_PREDICTIONS_R2.json', 'RESULTS_R2.json'), ('FROZEN_PREDICTIONS_R3.json', 'RESULTS_R3.json'), ('FROZEN_PREDICTIONS_R4.json', 'RESULTS_R4.json'), ('FROZEN_PREDICTIONS_R4_V2.json', 'RESULTS_R4_V2.json'), ('FROZEN_PREDICTIONS_SURFACE_COVERAGE.json', 'RESULTS_SURFACE_COVERAGE.json')]:
        payload = read(ROOT / name)['payload']
        assert sha(ROOT / result) == payload.get('result_sha256', payload.get('results_sha256'))
        for (file, value) in payload.get('point_files', {}).items():
            assert sha(file) == value
    for (name, file) in [('FROZEN_PREDICTIONS_R2.json', 'run_r2.py'), ('FROZEN_PREDICTIONS_R3.json', 'run_r3.py'), ('FROZEN_PREDICTIONS_R4_V2.json', 'full_benchmark.py'), ('FROZEN_PREDICTIONS_SURFACE_COVERAGE.json', 'extend_surface_ports.py')]:
        assert sha(ROOT / 'code' / file) == read(ROOT / name)['payload']['code_sha256']
    payload = read(ROOT / 'FROZEN_PREDICTIONS_R5.json')['payload']
    assert sha(ROOT / 'RESULTS_R5.json') == payload['result_sha256']
    assert sha(ROOT / 'code/local_ray_branch.py') == payload['code_sha256']
    for (file, value) in payload['point_files'].items():
        assert sha(file) == value
    assert sha(ROOT / 'code/full_benchmark_v1_archived.py') == read(ROOT / 'FROZEN_PREDICTIONS_R4.json')['payload']['code_sha256']
    r4 = read(ROOT / 'RESULTS_R4_V2.json')
    counts = r4['counts']
    assert counts == {'ABSTAIN': 84, 'EXACT_SCENARIO_FULL_PLANAR_INTAGLIO': 488, 'UNKNOWN_SOURCE_SITE': 4}
    for row in r4['rows']:
        source = Path(row['source_design_file'])
        assert sha(source) == row['source_design_sha256']
        task = source.parents[1].parent / 'public' / source.name
        assert sha(task) == row['source_task_sha256']
    expected = read(ROOT / 'FROZEN_PREDICTIONS.json')['payload']['code'].get(str(ROOT / 'code/sufficiency.py'))
    reconcile = dict(path='code/sufficiency.py', original_sha256=expected, preserved_failed_version='code/sufficiency_v1_failed.py', preserved_sha256=sha(ROOT / 'code/sufficiency_v1_failed.py'), reason='Versioned repair of a failed fixture; collision/run_r1 numerical operators unchanged')
    if expected is None:
        reconcile['freeze_relation'] = 'Sufficiency fixture code added after R1 fields froze; retained failed version is independently hashed, not a changed field generator'
    else:
        assert expected == reconcile['preserved_sha256']
    for name in ['collision.py', 'run_r1.py', 'common.py']:
        p = ROOT / 'code' / name
        assert sha(p) == read(ROOT / 'FROZEN_PREDICTIONS.json')['payload']['code'][str(p)]
    out = dict(frozen_hashes_pass=True, replayed_positive_poses=checked, corrupted_CAM_observation_rejected=True, corrupted50um_stock_observation_rejected=True, unmeasured_20case_template_stays_UNKNOWN=True, unmeasured20case_stock_template_stays_UNKNOWN=True, freeze_files=freezes, code_reconciliation=reconcile, own_bytes=budget(), review_state='PENDING_INDEPENDENT_REVIEW')
    dump(ROOT / 'raw/DEMO_VERIFICATION.json', out)
    print(out)
if __name__ == '__main__':
    run()
