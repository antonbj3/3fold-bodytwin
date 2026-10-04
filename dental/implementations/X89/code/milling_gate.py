"""Read-only X49 extension: same triple plus local milling requirements.

This file is the new rule integration point; predecessor code stays unchanged.
"""
import argparse, importlib.util, copy
from codesign import *
from sliver_scene import SliverScene

def evaluate(crown_path, contract_path):
    start = time.perf_counter()
    c = read(contract_path)
    if c.get('units') != 'mm':
        raise ValueError('Explicit final crown mm required')
    if sha(crown_path) != c['crown_sha256']:
        raise ValueError('Crown/contract mismatch')
    p = c['geometry_file']
    if sha(p) != c['geometry_sha256']:
        raise ValueError('Labelled geometry drift')
    a = np.load(p, allow_pickle=False)
    scene = SliverScene(a['vertices'], a['faces'])
    import trimesh
    source = trimesh.load_mesh(crown_path, process=False)
    if source.triangles.shape != scene.tri.shape:
        raise ValueError('Geometry-to-export triangle mapping mismatch')
    parity = float(np.linalg.norm(source.triangles - scene.tri, axis=2).max())
    if parity > c['scenario']['roundoff_allowance_mm']:
        raise ValueError('Geometry-to-export coordinates mismatch')
    regions = {k: np.asarray(v, dtype=int) for (k, v) in c['region_faces'].items()}
    for ids in regions.values():
        if len(ids) == 0 or ids.min() < 0 or ids.max() >= len(scene.f):
            raise ValueError('Bad region facet labels')
    if sha(c['tool_card']) != c['tool_card_sha256']:
        raise ValueError('Tool card drift')
    m = c['scenario']
    tcard = read(c['tool_card'])['tools']
    rows = []
    for key in ('scale_green_to_final', 'modulus_N_mm2', 'force_N', 'film_error_tolerance_final_mm'):
        if not math.isfinite(m[key]) or m[key] <= 0:
            raise ValueError('Invalid ' + key)
    for key in ('holder_compliance_mm_N', 'other_geometry_error_final_mm', 'roundoff_allowance_mm'):
        if not math.isfinite(m[key]) or m[key] < 0:
            raise ValueError('Invalid ' + key)
    region_ids = {k: v[scene.valid[v]] for (k, v) in regions.items()}
    pts = samples(scene, region_ids, c['seed'], c['samples_per_region'])
    for point in pts:
        tools = [robust_search(scene, point['point'], point['normal'], t, m, c.get('installed_ports')) for t in tcard]
        vi = [t for t in tools if t['status'] == 'CONDITIONAL_LOCAL_WITNESS']
        rows.append(dict(**point, tools=tools, selected_tool=min(vi, key=lambda t: t['nominal_tip_diameter_mm'])['tool_id'] if vi else None))
    failed_mechanics = all((mechanical(t, m)['mechanical_status'] == 'SCENARIO_INFEASIBLE' for t in tcard)) and (not c.get('installed_ports'))
    all_witness = bool(rows) and all((x['selected_tool'] is not None for x in rows))
    required = ['same_tool_force_bound', 'measured_installed_response', 'measured_profile_envelope', 'fixture_CAM_swept_path', 'runout_dynamic_bound', 'lot_sinter_distortion_bound', 'continuous_surface_coverage']
    missing = [k for k in required if c.get('measurement_locators', {}).get(k) is None]
    return dict(status='FAIL' if failed_mechanics else 'UNKNOWN', scenario_verdict='NO_TOOL_MEETS_DECLARED_MECHANICS' if failed_mechanics else 'ALL_SAMPLED_POINTS_HAVE_CONDITIONAL_WITNESS' if all_witness else 'NO_LIBRARY_WITNESS_AT_SOME_POINTS', reason='Inverse local cutter/neck/holder requirements; physical admission unavailable', claim_type='capability', resolution='PER_POINT', time_scale='HANDOVER', crown_sha256=sha(crown_path), contract_sha256=sha(contract_path), scenario=m, geometry_to_export_max_mm=parity, geometry_to_export_resolution='PER_POINT', points=len(rows), found=sum((x['selected_tool'] is not None for x in rows)), records=rows, missing_measurements=missing, full_crown_status='UNKNOWN_SAMPLED_GEOMETRY_AND_NO_PAIRED_MEASUREMENT', minimum_cutting_diameter_global='UNKNOWN_CONTINUOUS_POSES_AND_CAM', alarm_no_library_witness=failed_mechanics or not all_witness, no_library_alarm_scope='Declared force/scenario and finite pose grid; not global machining impossibility', query_s=time.perf_counter() - start, rigorous_floating_enclosure=False)

def parent_check(paths, material, contract, units, sinter):
    sys.path.insert(0, str(DENTAL / 'results/LANE_X49_DESIGN_GATE/code'))
    spec = importlib.util.spec_from_file_location('x89_parent_gate', DENTAL / 'results/LANE_X49_DESIGN_GATE/code/design_gate.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.check(paths, material, contract, units, sinter)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('crown')
    ap.add_argument('--milling-contract', required=True)
    ap.add_argument('--prep')
    ap.add_argument('--antagonist')
    ap.add_argument('--contract', help='Existing X49 geometric triple contract')
    ap.add_argument('--material', default='3Y')
    ap.add_argument('--sinter-factor', type=float)
    ap.add_argument('--units', choices=['mm', 'um'])
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    try:
        milling = evaluate(args.crown, args.milling_contract)
        if args.prep and args.antagonist:
            out = parent_check(dict(prep=args.prep, antagonist=args.antagonist, crown=args.crown), args.material, args.contract, args.units, args.sinter_factor)
            old_verdict = out['verdict']
            out['rules']['milling_codesign'] = milling
            out['parent_verdict_preserved'] = old_verdict
            out['verdict'] = 'FAIL' if old_verdict == 'FAIL' or milling['status'] == 'FAIL' else 'UNKNOWN'
        else:
            out = dict(verdict=milling['status'], rules=dict(milling_codesign=milling), parent_status='NOT_RUN_NO_TRIPLE')
        dump(args.output, out)
        print(out['verdict'], milling['scenario_verdict'], milling['found'], '/', milling['points'])
        return {'PASS': 0, 'FAIL': 2, 'UNKNOWN': 3}[out['verdict']]
    except (ValueError, KeyError, TypeError, OSError, IndexError) as ex:
        dump(args.output, dict(verdict='UNKNOWN', input_rejected=True, reason=str(ex)))
        print(str(ex))
        return 4
if __name__ == '__main__':
    sys.exit(main())
