from dental_release.paths import expand as _release_expand
import csv
from common import *

def main():
    p = ROOT / 'PREREG_R2.json'
    if p.exists():
        assert sha(p) == p.with_suffix('.sha256').read_text().strip()
        return
    r1 = read(ROOT / 'RESULTS_R1.json')
    manifest = read(ROOT / 'INPUT_MANIFEST.json')
    selected = [r for r in manifest['rows'] if r['family'] == 'crown_loop'] + [r for r in manifest['rows'] if r['family'] == 'GenCAD_Foundation' and r['id'].endswith('_original')] + [r for r in manifest['rows'] if r['family'] == 'GenCAD_V2'][:7]
    cases = []
    for (i, r) in enumerate(selected):
        cases.append(dict(case_id=f'CAM{i + 1:02d}', design_id=r['id'], mesh_file=r['file'], mesh_sha256=r['sha256'], prediction_file=str(DATA / 'R1' / (r['id'] + '.json')), family=r['family']))
    assert len(cases) == 20
    dump(ROOT / 'LAB_CASES_20.json', dict(claim_type='capability', cases=cases, prospective=True, physical_CAM_jobs_run=0))
    with open(ROOT / 'CAM_LOG_TEMPLATE.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case_id', 'library', 'axes', 'point_index', 'mesh_sha256', 'tool_library_sha256', 'observed_reachable', 'measured_assembly_and_frame', 'cam_locator', 'remaining_stock_um', 'normal_gap_um', 'measurement_uncertainty_um'])
        for c in cases:
            w.writerow([c['case_id'], 'vhf', 5, 0, c['mesh_sha256'], sha(ROOT / 'TOOL_LIBRARY.json'), '', '', '', '', '', ''])
    payload = dict(round='R2', claim_type='capability', capability='Transfer the signed local stock/forced-overcut field to the nominal normal cement gap at the identical point, without hidden seating addition', obstacle='An access fraction or mean stock cannot identify gap obstruction because stock location and nominal gap co-vary', changed_operation='For 20 locked designs find safe normal stand-off poses where R1 found none; intersect finite union of cutting balls with each oriented normal ray; hand over signed material change at PER_POINT', consumer=['DENT-MFG-PROCESS-MODEL', _release_expand('X13')], resolution='PER_POINT', timescale='HANDOVER', metrics={'max_normal_standoff_mm': 0.4, 'offset_grid_mm': [0.05, 0.1, 0.2, 0.4], 'film_shift_for_reporting_um': 25.0, 'summary_identity_error_um_max': 0.0, 'summary_downstream_obstructed_fraction_min': 0.1, 'obstruction_scenario_gap_um': 40.0, 'sphere_ray_control_abs_um_max': 1e-06, 'external_gap_score_eligibility': 'Requires paired same-object nominal die/intaglio, installed CAM tool scene and local measured gap; none currently supplied'}, decision_criteria={'construction': 'Export 20 local field maps; exact same mean film states must differ by at least .1 in local obstruction fraction', 'physical_fit': 'UNKNOWN; no post hoc fit to published region means', 'safe_policy': 'No retained-solid intersection under the declared full-tool scenario; NOT_FOUND beyond .4 is UNKNOWN, not a measured binding height', 'forced_policy': 'Diagnostic union of balls at original targets, ignores access; explicit overcut scenario, not a CAM prediction'}, strongest_equal_information_control=dict(name='Scalar quadratic normal-ray sphere intersection and independent trimesh die ray intersections', same_information='same ball centres/radii and original die; no method-win claim'), practice='Same nominal spacer assigned everywhere regardless of spatial stock/overcut', falsifiers=['Injected +50um normal-ray root must fail scalar sphere control', 'Stock shift that removes all nominal normal gap must be marked potential interference', 'Any aggregate mean used to label full seating is rejected by the identical-mean pair', 'Published region gaps must remain ineligible without paired local geometry'], external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.4047/jap.2016.8.6.439 Table2', compared_quantity='regional manufactured internal gap in micrometres; transport currently UNKNOWN_UNPAIRED_GEOMETRY', refutes_us=True), full_cost={'prepare': 'Reuse frozen R1 meshes/witnesses and 20-case manifest', 'fit': 'none', 'discovery': 'Additional stand-off pose queries only where R1 had no zero-offset witness; ball union normal-ray intersections', 'validation': 'Scalar sphere root with injected error, exact mean-pair obstruction, frozen manifest and exported normal-ray die intersection', 'questions': 0, 'fallback': 'no normal die intersection or no safe cutter within400um -> UNKNOWN; preserve spatial output', 'physical': '20 CAM jobs + scans not executed; cost UNKNOWN'}, frozen_R1_sha256=sha(ROOT / 'FROZEN_PREDICTIONS.json'), lab_cases_sha256=sha(ROOT / 'LAB_CASES_20.json'))
    dump(p, payload)
    p.with_suffix('.sha256').write_text(sha(p) + '\n')
    dump(ROOT / 'DECOMPOSITION_R2.json', dict(leaves=[dict(leaf='sphere ray', status='DERIVED_UNDER_ASSUMPTIONS', equation='t_entry=n·(c-p)-sqrt(r²-|c-p|²+(n·(c-p))²)', stop_argument='Exact for the declared ball and unit normal; float rounding enclosure missing'), dict(leaf='normal gap', status='DERIVED_UNDER_ASSUMPTIONS', equation='first positive die intersection t0; stock entry s; h_new=t0-s', stop_argument='Exact on this fixed geometric ray. Does not solve seating/tilt or thin-film hydraulics.'), dict(leaf='cutting envelope', status='CONSTITUTIVE_CLOSURE', equation='finite set of final balls; swept cutting neck not reconstructed as as-built material', stop_argument='Sparse sphere cut union is a digital witness set; a real CAM stock field is needed to replace it.'), dict(leaf='physical fit', status='UNKNOWN', equation='registered measured intaglio minus die, same process and same object', stop_argument='The published regional measurements lack nominal CAD and local correspondence; must not be fitted as pixel truth.')]))
    state('R2_PREREG_FROZEN', sha(p), 'Construct safe local material-change port and test identical mean gap')
    print(sha(p))
if __name__ == '__main__':
    main()
