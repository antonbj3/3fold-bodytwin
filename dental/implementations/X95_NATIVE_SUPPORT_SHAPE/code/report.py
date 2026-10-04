from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def run():
    r1 = read(R / 'raw/R1_NATIVE_PORT.json')
    r2 = read(R / 'raw/R2_SUPPORT_COVERAGE.json')
    r3 = read(R / 'raw/R3_COMPILED_SUPPORT.json')
    r4 = read(R / 'raw/R4_PORT_CERTIFICATE.json')
    ext = read(R / 'raw/EXTERNAL_FACIT.json')
    solver = read(R / 'raw/NATIVE_SOLVER.json')
    validation = read(R / 'raw/VALIDATION.json')
    p = np.load(DATA / 'NATIVE_PORT.npz')
    c = np.load(DATA / 'COMPILED_SUPPORT.npz')
    sl = p['slaves']
    dist = p['distances_mm']
    xyz = p['cV10'][sl]
    faces = p['slave_faces']
    mapping = {int(n): i for (i, n) in enumerate(sl)}
    R.joinpath('exports').mkdir(parents=True, exist_ok=True)
    R.joinpath('figures').mkdir(parents=True, exist_ok=True)
    sub = []
    for face in faces[:, 2:8]:
        (a, b, c0, d, e, f) = [mapping[int(n)] for n in face]
        sub.extend([[a, d, f], [d, b, e], [f, e, c0], [d, e, f]])
    lines = ['# vtk DataFile Version 3.0', 'X95 native support diagnostic; mm, no optimized specimen', 'ASCII', 'DATASET POLYDATA', f'POINTS {len(xyz)} double']
    lines.extend((' '.join((format(float(z), '.17g') for z in x)) for x in xyz))
    lines.append(f'POLYGONS {len(sub)} {4 * len(sub)}')
    lines.extend(('3 ' + ' '.join(map(str, t)) for t in sub))
    lines.append(f'POINT_DATA {len(xyz)}')
    for (name, vals) in [('projection_distance_mm', dist), ('within_006mm', dist <= 0.06), ('compiled_active', c['active_mask'])]:
        lines.extend([f'SCALARS {name} double 1', 'LOOKUP_TABLE default'])
        lines.extend((format(float(v), '.17g') for v in vals))
    lines.extend([f'CELL_DATA {len(sub)}', 'SCALARS source_crown_element_id_0based int 1', 'LOOKUP_TABLE default'])
    lines.extend((str(int(e)) for e in np.repeat(faces[:, 0], 4)))
    lines.extend(['SCALARS source_face_number int 1', 'LOOKUP_TABLE default'])
    lines.extend((str(int(e)) for e in np.repeat(faces[:, 1], 4)))
    (R / 'exports/native_support.vtk').write_text('\n'.join(lines) + '\n')
    import csv
    with (R / 'exports/native_support_nodes.csv').open('w') as f:
        w = csv.writer(f)
        w.writerow(['node_0based', 'x_mm', 'y_mm', 'z_mm', 'projection_mm', 'master_element_0based', 'master_face', 'compiled_active', 'within_006mm'])
        for (i, n) in enumerate(sl):
            w.writerow([int(n), *xyz[i], dist[i], *p['master_faces'][p['master_face'][i], :2], bool(c['active_mask'][i]), bool(dist[i] <= 0.06)])
    (fig, axs) = plt.subplots(2, 2, figsize=(11, 8), layout='constrained')
    ax = axs[0, 0]
    q = ax.scatter(xyz[:, 0], xyz[:, 2], c=dist, s=4, vmin=0, vmax=0.3, cmap='magma_r')
    fig.colorbar(q, ax=ax, label='Distance to native master [mm]')
    w = r1['worst_point']
    ax.scatter(*np.array(w['coordinate_mm'])[[0, 2]], marker='x', s=80, c='cyan')
    ax.plot([w['coordinate_mm'][0], w['closest_point_mm'][0]], [w['coordinate_mm'][2], w['closest_point_mm'][2]], c='cyan')
    ax.set(xlabel='Native frame x [mm]', ylabel='Native frame z [mm]', title='Declared support: local error retained')
    ax.set_aspect('equal')
    ax = axs[0, 1]
    d = np.sort(dist)
    ax.plot(d, np.arange(1, len(d) + 1) / len(d), label='All 5251 declared nodes')
    active = np.sort(dist[c['active_mask']])
    ax.plot(active, np.arange(1, len(active) + 1) / len(active), label='5114 solver-active nodes')
    ax.axvline(0.06, c='crimson', ls='--', label='Frozen 0.060 mm maximum')
    ax.set(xlabel='Projection distance [mm]', ylabel='Fraction of nodes', title=f'{int(np.sum(dist > 0.06))}/{len(dist)} = {100 * np.mean(dist > 0.06):.2f}% of native declared nodes exceed gate')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 0.55)
    ax = axs[1, 0]
    vals = r2['sufficiency']['strain_energy_N_mm']
    ax.bar(['Local load A', 'Local load B'], np.array(vals) * 1000000.0, color=['#336c97', '#ca804a'])
    ax.set(ylabel='Strain energy [micro N mm]', title=f"Exact same total load + support area\nenergy differs {r2['sufficiency']['relative_downstream_difference'] * 100:.2f}%")
    ax.text(0.02, 0.96, 'Conditional FE witness; not a physical measurement', transform=ax.transAxes, va='top', fontsize=8)
    ax = axs[1, 1]
    rr = ext['Chen_1mm_matched_source']['source_rows']
    ax.bar(['Chen 3Y, 1.0mm', 'Chen 5Y, 1.0mm'], [x['mean_N'] for x in rr], yerr=[x['sd_N'] for x in rr], color=['#336c97', '#ca804a'], capsize=4)
    ax.set(ylabel='Mean fracture force ± SD [N]', title='External Table 2: own premolar protocol')
    ax.text(0.02, 0.96, 'D1 absolute fracture force: UNKNOWN', transform=ax.transAxes, va='top', fontsize=9)
    fig.suptitle('X95: fixed native D1 support fails before shape optimization', fontsize=14)
    fig.savefig(R / 'figures/demo.png', dpi=180)
    fig.savefig(R / 'figures/demo.pdf')
    plt.close(fig)
    arrays = [{'path': str(f), 'sha256': sha(f), 'bytes': f.stat().st_size} for f in sorted(DATA.glob('*.npz'))]
    expr = 'distance(p,triangle union) >= min distance(p,triangle enclosing AABB)'
    outcome = {'claim_type': 'capability', 'outcome': 'FIXED_NATIVE_PORT_REJECTED_WITH_GEOMETRIC_CERTIFICATE; PHYSICAL_STRENGTH_UNKNOWN', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'case': 'D1 / L005_lo_M1_k6 thin', 'capability_answer': 'No defensible optimized/reference specimen pair can be selected under unchanged native support and the frozen .06mm maximum. Shape optimization was correctly blocked.', 'external_referent': {'kind': 'published_code', 'locator': 'https://www.dhondt.de/ccx_2.23.pdf section7.133 pp620-621; native gentiedmpc.f', 'compared_quantity': 'generation/non-generation of tied structural MPCs and native omission report, not empirical crown strength', 'refutes_us': True}, 'external_measurements': ext['external_referent'], 'closed_form_referent': {'kind': 'closed_form', 'locator': str(SRC / 'code/port_certificate.py'), 'compared_quantity': expr, 'refutes_us': True}, 'source_protocol_contrast': ext['Chen_1mm_matched_source'], 'rounds': {'R1': r1, 'R2': r2, 'R3': r3, 'R4': r4}, 'native_solver': {'saved_tensor_relative_L2_error': solver['relative_tensor_L2_error'], 'numerical_parity_pass': solver['native_numeric_parity_pass'], 'peak_RSS_MiB': solver['source_solver_peak_RSS_MiB'], 'native_full_MPC_owner_weight_parity': 'UNKNOWN_NOT_EXTRACTED', 'source_solver_scope': 'Unchanged native C3D10 source deck; actual original solve performed. Cached reuse if shared gate blocks is explicitly labeled; no empirical validation', 'native_replay_status': solver.get('replay_native_status', 'FRESH_NATIVE_SOLVE')}, 'absolute_fracture_force_N': None, 'absolute_force_status': 'UNKNOWN_UNMATCHED_D1_GEOMETRY_SUPPORT_LOAD_AND_MATERIAL_BATCH', 'physical_measurement_status': 'NOT_PERFORMED', 'shape_optimizer_runs': 0, 'candidate_specimen_status': 'NONE_PORT_GATE_BLOCKED', 'numerical_controls': {'same_Q_KKT_max_relative_response_error': max([x['displacement_relative_L2_error'] for x in r2['same_Q_full_FE_control']] + [x['same_Q_KKT_relative_displacement_error'] for x in r3['P1_diagnostics']]), 'algorithmic_advantage_claimed': False}, 'validation': validation, 'dropout': {'cohort': r1['dropout']['cases'], 'distance_qualified_support': r2['native_coverage'], 'compiled_source_support': r3['counts'], 'external_protocol_blocks': ext['original_X1C_bracket_recalculation']['dropout']}, 'resolution': {'local_projection_Q_and_loads': 'PER_POINT', 'support_patch_coverage': 'PER_SURFACE_REGION', 'elastic_energy': 'PER_TOOTH', 'source_fracture_group_means': 'POPULATION', 'volume_flaw_hazard': 'PHENOMENOLOGICAL', 'phenomenological_debt': 'Same-batch crown fracture and origin needed to replace borrowed Prott m; no hazard-based strength claim'}, 'edge': {'producer': 'native local surface/compiled support', 'consumer': 'researcher/material/cad shape feasibility', 'finest_common_resolution': 'PER_SURFACE_REGION', 'underlying_node_resolution': 'PER_POINT', 'timescale': 'SIMULTANEOUS'}, 'sufficiency': r2['sufficiency'], 'rigorous_enclosure': {'geometry': r4['rigorous_distance_enclosure_mm'], 'scope': 'one locked point and union of unchanged planar native master triangles', 'FE_shape_or_physical_enclosure': 'MISSING'}, 'full_cost': {'round_numerical_wall_seconds_sum': sum((r['cost']['wall_seconds'] for r in [r1, r2, r3, r4])) + ext['cost']['wall_seconds'] + solver['cost']['wall_seconds'], 'recorded_peak_RSS_MiB': max([r['cost']['peak_RSS_MiB'] for r in [r1, r2, r3, r4]] + [solver['source_solver_peak_RSS_MiB']]), 'FE_solves_P1': r2['FE_solve_calls'] + 4, 'native_C3D10_load_cases': 2, 'fit': 'None', 'upstream_acquisition': 'UNKNOWN', 'discovery_preparation_query_and_gate_wait': 'Elapsed lane timestamps and command logs retained; complete CPU/token attribution UNKNOWN', 'validation': 'Actual faults plus copied replay reported separately', 'physical_fit_instrument_manufacture_and_validation': 'UNKNOWN_NOT_PERFORMED', 'fallback': 'A conforming common-interface reconstruction from the unchanged preparation geometry and pointwise-certified support, then source-native MPC owner/weight export; physical seated support/force-displacement/fracture origins', 'thread_limit': 4, 'GPU_used': False, 'own_array_bytes': sum((x['bytes'] for x in arrays)), 'own_intermediate_limit_bytes': 3 * 1024 ** 3}, 'arrays': arrays, 'prereg_hashes': {p.name: sha(p) for p in sorted(SRC.glob('PREREG*.json'))}, 'next_construction': 'Change geometric interface: partition partially supported crown facets on die-SDF boundary and build a conforming shared interface. Certify .06mm and full required surface coverage before freezing any new minimax. Do not change material or force. Same-specimen support metrology is a separate prerequisite.'}
    outcome['dependency_lock_sha256'] = sha(SRC / 'DEPENDENCY_LOCK.json')
    if (R / 'raw/R4_DECIMAL_CHECK.json').exists():
        outcome['independent_geometry_enclosure_check'] = read(R / 'raw/R4_DECIMAL_CHECK.json')
    outcome['sufficiency']['minimum_extension_for_this_two_state_family'] = 'One bit specifying which of the two local load locations is used; both states retain identical force total and support area. For general loads retain the local load/Q relation; global minimal dimension UNKNOWN.'
    dump(R / 'results.json', outcome)
    if R == SRC and (not (SRC / 'FROZEN_PREDICTIONS.json').exists()):
        freeze(SRC / 'FROZEN_PREDICTIONS.json', {'frozen_utc': now(), 'prereg_hashes': outcome['prereg_hashes'], 'source_input_lock_sha256': sha(SRC / 'INPUT_LOCK.json'), 'case': outcome['case'], 'geometry_rejection': r4['rigorous_distance_enclosure_mm'], 'native_projection_max_mm': r1['interface_projection_mm']['max'], 'compiled_omitted_nodes': r3['counts']['compiled_omitted_nodes'], 'conditional_energies': r2['contrasts_P1_conditional'], 'arrays': arrays, 'candidate_force_N': None, 'candidate_pair': 'NONE_PORT_REJECTED', 'measurement_status': 'NO_PHYSICAL_MEASUREMENT_PERFORMED', 'future_measurements': 'registered seated die/cement/support map, load-deflection and origin; any reconstructed geometry is a new preregistration, not refit of these predictions'})
    print(json.dumps({'outcome': outcome['outcome'], 'result_sha256': sha(R / 'results.json'), 'validation_checks': validation['passed'], 'array_MiB': sum((x['bytes'] for x in arrays)) / 1024 ** 2}))
if __name__ == '__main__':
    run()
