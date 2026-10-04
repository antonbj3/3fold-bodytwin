"""Falsifying controls on the actual external dataset and exported geometry."""
import json
from pathlib import Path
import numpy as np
from freeze import ROOT, DATA, write, sha
from planner import dist, load, reflection, extract_rail, dp_plan, exhaustive, planes, graft_mesh, graft_mesh_legacy, graft_levelset, unit
from constrained import cut_feasible
from implant_variant import envelope

def run():
    r1 = json.loads((ROOT / 'RESULTS_R1.json').read_text())
    r2 = json.loads((ROOT / 'RESULTS_R2.json').read_text())
    r3 = json.loads((ROOT / 'RESULTS_R3.json').read_text())
    cid = next((r['case'] for r in r2['cases'] if r['numerically_available']))
    raw = np.load(DATA / 'R2' / cid / 'evaluation_raw.npz')
    truth = raw['expert_roi'].astype(float)
    good = dist(truth, truth)
    wrong = dist(truth + np.array([20.0, 20.0, 20.0]), truth)
    external = dict(case=cid, source='published Post added-region samples', correct_reference_p95_mm=float(np.quantile(good, 0.95)), injected_translation_mm=[20.0, 20.0, 20.0], wrong_reference_p95_mm=float(np.quantile(wrong, 0.95)), correct_pass=bool(np.quantile(good, 0.95) <= 3), wrong_fails=bool(np.quantile(wrong, 0.95) > 3))
    plan = json.loads((DATA / 'R2' / cid / 'plan.json').read_text())
    detail = plan['grafts']['dp']
    nodes = np.array(detail['nodes_mm'])
    (mesh, check) = graft_mesh(nodes)
    nn = np.array(detail['cut_plane_normals'])
    ring = mesh[0][0][:48]
    p = nodes[0]
    norm = nn[0]
    residual = float(np.max(abs((ring - p) @ norm)))
    corrupt = float(np.max(abs((ring - (p + 20 * norm)) @ norm)))
    geometry = dict(correct_plane_residual_mm=residual, injected_plane_offset_mm=20.0, wrong_plane_residual_mm=corrupt, correct_pass=residual <= 1e-08, wrong_fails=corrupt > 1e-08)
    (legacy, _) = graft_mesh_legacy(nodes)
    wrong_edges = []
    for (i, (v, f)) in enumerate(legacy):
        axis = unit(nodes[i + 1] - nodes[i])
        edge = v[48:96] - v[:48]
        wrong_edges.append(float(np.linalg.norm(edge - (edge @ axis)[:, None] * axis, axis=1).max()))
    straight = dict(correct_max_transverse_edge_motion_mm=check['longitudinal_edge_transverse_motion_max_mm'], wrong_legacy_transverse_edge_motion_mm=max(wrong_edges), correct_pass=check['longitudinal_edge_transverse_motion_max_mm'] <= 1e-08, wrong_fails=max(wrong_edges) > 1e-08)
    (pre, a, _) = load(cid, 'Pre')
    (target, _) = reflection(pre, a)
    (curve, _) = extract_rail(target, pre, a)
    (_, candidate) = dp_plan(curve)
    control = exhaustive(curve, candidate['count'])
    disagreement = abs(candidate['max_chord_error_mm'] - control['max_chord_error_mm'])
    solver = dict(correct_disagreement_mm=disagreement, injected_objective_offset_mm=1.0, wrong_disagreement_mm=disagreement + 1.0, correct_pass=disagreement <= 1e-08, wrong_fails=disagreement + 1.0 > 1e-08)
    lengths = np.linalg.norm(np.diff(nodes, axis=0), axis=1)
    reported = np.array(detail['lengths_mm'])
    length_check = dict(correct_max_difference_mm=float(abs(lengths - reported).max()), wrong_10x_max_difference_mm=float(abs(lengths - 10 * reported).max()), correct_pass=bool(np.max(abs(lengths - reported)) <= 1e-08), wrong_fails=bool(np.max(abs(lengths - 10 * reported)) > 1e-08))
    simple = np.array([[-15.0, 0.0, 0.0], [15.0, 0.0, 0.0]])
    goodenv = envelope(simple, 8.0, np.array([0.0, 0.0, 1.0]), 5.0, prosthetic_space=0.0, length=10.0, implant_radius=2.0)
    badenv = envelope(simple, 3.0, np.array([0.0, 0.0, 1.0]), 5.0, prosthetic_space=0.0, length=10.0, implant_radius=2.0)
    implant = dict(kind='our_own_fixture_for_assay_only', correct_containment_upper_bound_mm=goodenv[0]['containment_upper_bound_mm'], injected_radius_mm=3.0, wrong_containment_upper_bound_mm=badenv[0]['containment_upper_bound_mm'], correct_pass=goodenv[0]['contained'], wrong_fails=not badenv[0]['contained'])
    folded = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    thickness = dict(correct_straight_accepted=cut_feasible(simple, 6.0)[0], wrong_folded_rejected=not cut_feasible(folded, 6.0)[0])
    aggregate = {}
    for round_id in ('R1', 'R2', 'R3', 'R4'):
        outcome = json.loads((ROOT / f'RESULTS_{round_id}.json').read_text())
        checks = {}
        for gate in outcome['gates']:
            if gate == 'all_cases_available':
                value = 117
                criterion = 'available==118'
                reject = value != 118
            elif gate in ('median_primary_p95_mm_max',):
                value = 999.0
                criterion = 'median<=3mm'
                reject = value > 3
            elif gate == 'solver_matches_exhaustive':
                value = 1.0
                criterion = 'disagreement<=1e-8mm'
                reject = value > 1e-08
            elif gate == 'all_miters_non_crossing':
                value = False
                criterion = 'all noncrossing'
                reject = not value
            elif gate == 'fraction_all_segment_envelopes_contained_min':
                value = 0.0
                criterion = 'fraction>=0.8'
                reject = value < 0.8
            elif gate == 'median_surface_penalty_vs_shape_only_mm_max':
                value = 999.0
                criterion = 'penalty<=2mm'
                reject = value > 2
            else:
                value = -999.0
                criterion = 'paired gain >= frozen positive minimum'
                reject = value < 0
            checks[gate] = dict(injected_value=value, criterion=criterion, injection_rejected=bool(reject))
        aggregate[round_id] = checks
    verdict = all((x['correct_pass'] and x['wrong_fails'] for x in (external, geometry, straight, solver, length_check, implant))) and all(thickness.values()) and all((c['injection_rejected'] for g in aggregate.values() for c in g.values()))
    field = json.loads((ROOT / 'FIELD_DEMO.json').read_text())
    field_invariant = field['pre_inside_voxels_lost'] == 0
    field_path = next((e['path'] for e in field['exports'] if Path(e['path']).name == 'fields.npz'))
    with np.load(field_path) as z:
        corrupted = z['union'].copy()
        inside = z['pre'] < 0
        point = tuple(np.argwhere(inside)[0])
        corrupted[point] = 1.0
        bad_count = int(np.sum(inside & (corrupted >= 0)))
    output = dict(passed=bool(verdict and field_invariant and (bad_count > 0)), external_surface_assay=external, cut_planes=geometry, straight_cylinder_side_generators=straight, independent_solver=solver, reported_lengths=length_check, conditional_implant_assay=implant, thickness=thickness, frozen_aggregate_gate_injections=aggregate, field_preservation=dict(correct_pre_voxels_lost=field['pre_inside_voxels_lost'], injected_one_voxel_lost=bad_count), interpretation='assay sensitivity and numerical contracts; does not validate anatomy, donor geometry, occlusion, mechanics or clinical outcomes', result_hashes={r: sha(ROOT / f'RESULTS_{r}.json') for r in ('R1', 'R2', 'R3', 'R4')})
    stl_checks = []
    dtype = np.dtype([('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attr', '<u2')])
    for rid in ('R3', 'R4'):
        for rr in json.loads((ROOT / f'FROZEN_PREDICTIONS_{rid}.json').read_text())['predictions']:
            pp = json.loads(Path(rr['path']).read_text())
            if pp['status'] != 'PREDICTED':
                continue
            detail = pp['grafts']['dp']
            nd = np.array(detail['nodes_mm'])
            for i in range(detail['count']):
                source = Path(rr['path']).parent / f'dp_segment_{i + 1}.stl'
                triangles = np.memmap(source, dtype=dtype, offset=84, mode='r')['vertices'].astype(float)
                side = triangles[::4]
                edges = side[:, 2] - side[:, 1]
                axis = unit(nd[i + 1] - nd[i])
                perp = edges - (edges @ axis)[:, None] * axis
                err = float(np.linalg.norm(perp, axis=1).max())
                stl_checks.append(dict(round=rid, case=pp['case'], segment=i + 1, path=str(source), sha256=sha(source), exported_transverse_edge_motion_mm=err, passes=err <= 0.0001))
    output['actual_STL_straightness'] = dict(tolerance_mm=0.0001, n=len(stl_checks), max_transverse_motion_mm=max((r['exported_transverse_edge_motion_mm'] for r in stl_checks)), all_pass=all((r['passes'] for r in stl_checks)), checks=stl_checks)
    output['passed'] = bool(output['passed'] and output['actual_STL_straightness']['all_pass'])
    write(ROOT / 'CONTROL_FALSIFICATION.json', output)
    if not output['passed']:
        raise AssertionError('A required control failed to reject its injected error')
    return output
if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
