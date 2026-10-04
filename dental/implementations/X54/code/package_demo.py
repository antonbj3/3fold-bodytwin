"""Reader-facing result, figure and region port. Every finding retains its scope."""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'review_cache/matplotlib'))
import csv, json, datetime, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from contact_model import P, D, BASE, sha, write, state
from run_r1 import area_shares

def main():
    r1 = json.loads((P / 'rounds/R1/results.json').read_text())
    r2 = json.loads((P / 'rounds/R2/results.json').read_text())
    r3 = json.loads((P / 'rounds/R3/results.json').read_text())
    r3b = json.loads((P / 'rounds/R3B/results.json').read_text())
    verify = json.loads((P / 'VERIFICATION.json').read_text())
    g = json.loads((P / 'raw/geometry_001_h02.json').read_text())
    query = json.loads((P / 'rounds/R2/case001_force100.json').read_text())
    a = json.loads((P / 'rounds/R2/case001_support_A.json').read_text())
    b = json.loads((P / 'rounds/R2/case001_support_B.json').read_text())
    relief = json.loads((P / 'rounds/R2/case001_relief30um.json').read_text())
    (fig, axs) = plt.subplots(2, 3, figsize=(16, 9), layout='constrained')
    colors = ['#164f86', '#dc7132', '#4d9279']
    types = ['I1', 'I2', 'C', 'P1', 'P2', 'M1', 'M2']
    xx = np.arange(7)
    for (ax, jaw, label) in [(axs[0, 0], 'upper', 'Ferrato 2017: upper arch'), (axs[0, 1], 'lower', 'Hattori 1996: lower arch')]:
        v = r2['external_comparison'][jaw]
        ax.bar(xx - 0.26, v['reference_pp'], 0.26, label='Published measurement', color='#777777')
        ax.bar(xx, v['rigid_area_pp'], 0.26, label='Rigid near-area rule', color=colors[1])
        ax.bar(xx + 0.26, v['candidate_pp'], 0.26, label='Conditional contact at 100 N', color=colors[0])
        ax.set_xticks(xx, types)
        ax.set_ylabel('Bilateral force / sensor share (%)')
        ax.set_title(label + '\n11-case aggregate; unmatched populations', fontsize=11)
        ax.grid(axis='y', alpha=0.2)
        ax.text(0.02, 0.97, f"MAE contact {v['MAE_pp']:.2f} pp / area {v['rigid_MAE_pp']:.2f} pp", transform=ax.transAxes, va='top', fontsize=9)
    axs[0, 0].legend(fontsize=8, loc='upper left', bbox_to_anchor=(0, 0.83))
    ax = axs[0, 2]
    keys = sorted(a['force_shares_pp']['upper'], key=int)
    xx = np.arange(len(keys))
    ax.plot(xx, [a['force_shares_pp']['upper'][t] for t in keys], 'o-', color=colors[1], label='Support A: 16 soft / 26 stiff')
    ax.plot(xx, [b['force_shares_pp']['upper'][t] for t in keys], 'o-', color=colors[0], label='Support B: 16 stiff / 26 soft')
    ax.set_xticks(xx, keys, rotation=45)
    ax.set_ylabel('Upper tooth normal-force share (%)')
    ax.set_title('Identical rigid gap and area: identity error 0\nCase 1; force-share difference 24.17 pp', fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)
    ax = axs[1, 0]
    rows = query['rows']
    xy = np.array([p['upper_point_mm'][:2] for p in rows])
    force = np.array([p['lambda_N'][0] for p in rows])
    plot = ax.scatter(xy[:, 0], xy[:, 1], c=force, s=12 + force * 14, cmap='viridis', edgecolors='white', linewidths=0.4)
    fig.colorbar(plot, ax=ax, label='Regional normal reaction (N)')
    for t in keys:
        ids = [i for (i, p) in enumerate(rows) if str(p['upper_fdi']) == t]
        if ids:
            pt = xy[ids].mean(0)
            ax.annotate(t, pt, xytext=(2, 3), textcoords='offset points', fontsize=8)
    ax.set_aspect('equal')
    ax.set_xlabel('Frame right (mm)')
    ax.set_ylabel('Frame anterior (mm)')
    ax.set_title('Case 1 registered-source regional contact\n100 N vertical resultant; hypothetical support', fontsize=11)
    ax = axs[1, 1]
    before = np.array([p['lambda_N'][0] for p in query['rows']])
    after = np.array([p['lambda_N'][0] for p in relief['rows']])
    order = np.argsort(before)[-8:]
    xx = np.arange(len(order))
    ax.bar(xx - 0.18, before[order], 0.36, color=colors[0], label='Original patch gaps')
    ax.bar(xx + 0.18, after[order], 0.36, color=colors[2], label='30 µm gap relief in peak region')
    ax.set_xticks(xx, [f"{rows[i]['upper_fdi']}:{rows[i]['region']}" for i in order], rotation=45)
    ax.set_ylabel('Regional normal reaction (N)')
    ax.set_title('Case 1 local crown-gap intervention\nFixed vertical force, area and normals', fontsize=11)
    ax.legend(fontsize=8)
    ax = axs[1, 2]
    x = np.arange(3)
    old = [max(row['max_share_difference_pp'].values()) for row in r2['refinement']]
    remaining = [row['remaining_area_difference_pp'] for row in r3['rows']]
    ax.bar(x - 0.18, old, 0.36, color=colors[1], label='Original raster 0.2 → 0.1 mm')
    ax.bar(x + 0.18, remaining, 0.36, color=colors[0], label='Continuous gap / remaining area')
    ax.axhline(2, color='darkred', ls='--', label='Frozen 2 pp gate')
    ax.set_xticks(x, ['Case 1', 'Case 2', 'Case 3'])
    ax.set_ylabel('Max tooth share difference (pp)')
    ax.set_title('Continuous affine gap removes sampled minimum\nOriginal named patches; full contact coverage unknown', fontsize=11)
    ax.legend(fontsize=8)
    fig.suptitle('Registered dental geometry → compliant Coulomb contact\nConditional simulation; individual load, preload and support are not measured', fontsize=15)
    fp = P / 'figures/deformable_contact.png'
    fig.savefig(fp, dpi=180)
    fig.savefig(P / 'figures/deformable_contact.pdf')
    plt.close(fig)
    gp = json.loads((P / 'raw/geometry_001_h02.json').read_text())
    R = np.asarray(gp['frame']['R'])
    center = np.asarray(gp['frame']['center'])
    export = P / 'regional_force_case001.csv'
    fields = ['case', 'upper_fdi', 'lower_fdi', 'region', 'upper_source_face', 'lower_source_face', 'x_mm', 'y_mm', 'z_mm', 'normal_force_N', 'tangential_force_N', 'vertical_force_N', 'mean_normal_traction_MPa', 'state', 'resolution', 'time_scale']
    points = []
    with export.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for p in query['rows']:
            xyz = np.asarray(p['upper_point_mm']) @ R.T + center
            points.append(xyz)
            writer.writerow(dict(case=1, upper_fdi=p['upper_fdi'], lower_fdi=p['lower_fdi'], region=p['region'], upper_source_face=p['upper_source_face'], lower_source_face=p['lower_source_face'], x_mm=xyz[0], y_mm=xyz[1], z_mm=xyz[2], normal_force_N=p['lambda_N'][0], tangential_force_N=np.linalg.norm(p['lambda_N'][1:]), vertical_force_N=p['vector_on_upper_N'][2], mean_normal_traction_MPa=p['mean_normal_traction_MPa'], state=p['state'], resolution='PER_SURFACE_REGION', time_scale='SIMULTANEOUS'))
    with (P / 'regional_force_case001.vtk').open('w') as f:
        f.write('# vtk DataFile Version 3.0\nX54 conditional source region force (mm,N)\nASCII\nDATASET POLYDATA\nPOINTS ' + str(len(points)) + ' float\n')
        for point in points:
            f.write(' '.join(map(str, point)) + '\n')
        f.write(f'VERTICES {len(points)} {2 * len(points)}\n')
        for i in range(len(points)):
            f.write(f'1 {i}\n')
        f.write(f'POINT_DATA {len(points)}\nSCALARS normal_reaction_N float 1\nLOOKUP_TABLE default\n')
        for p in query['rows']:
            f.write(str(p['lambda_N'][0]) + '\n')
    changes = {jaw: dict(changed=0, valid=0) for jaw in ('upper', 'lower')}
    for dec in r1['decisions']:
        s = json.loads((P / 'rounds/R1' / f"case{dec['case']:03d}_d0.05.json").read_text())
        if s['total_normal_N'] <= 1e-07:
            continue
        for jaw in changes:
            changes[jaw]['valid'] += 1
            changes[jaw]['changed'] += dec['changed'][jaw]
    disk = sum((f.stat().st_size for f in D.rglob('*') if f.is_file()))
    manifest = []
    for f in sorted(D.rglob('*')):
        if f.is_file() and f.suffix == '.npz':
            manifest.append(dict(path=str(f), bytes=f.stat().st_size, sha256=sha(f)))
    qrows = [s for row in r2['rows'] for s in list(row['force_targets'].values()) + list(row['scenarios'].values()) if 'numerics' in s]
    costs = dict(R1_wall_s=r1['full_cost']['this_run_wall_s'], R2_wall_s=r2['cost']['this_run_wall_s'], R3_wall_s=r3['cost']['wall_s'], R3B_wall_s=r3b['wall_s'], verification_s=verify['wall_s'], force_control_NCP_evaluations=r2['cost']['total_evaluations'], peak_RSS_MiB=max((s['cost']['peak_rss_MiB'] for s in qrows)), complete_preparation_fit_discovery_cost='UNKNOWN: inherited X11 fit and predecessor campaigns, human/provider costs not measured; not treated as free')
    result = dict(lane='X54-deformable-contact', claim_type=['information_link', 'capability'], status='PARTIAL_CONDITIONAL_CAPABILITY; INDIVIDUAL_PHYSICAL_FORCE_UNKNOWN', review_state='PENDING_INDEPENDENT_REVIEW', answer='Executable region-level compliant Coulomb contact, force-controlled crown-gap intervention, and exact-summary counterexample. No validated individual bite-force or crown-overcontact prediction.', external_referent=r2['external_referent'], external_comparison=r2['external_comparison'], rounds={r['round']: {'result_path': str(P / 'rounds' / r['round'] / 'results.json'), 'sha256': sha(P / 'rounds' / r['round'] / 'results.json')} for r in [r1, r2, r3]}, R3B=dict(path=str(P / 'rounds/R3B/results.json'), sha256=sha(P / 'rounds/R3B/results.json'), all_pass=r3b['all_pass']), frozen_predictions={n: sha(P / n) for n in ['FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R2.json', 'FROZEN_PREDICTIONS_R3.json']}, decision_change_valid_loaded_cases=changes, zero_force_case9_argmax='UNDEFINED; arbitrary first key in original R1 raw retained but excluded from decisions', sufficiency=dict(rigid_summary_identity_error_mm=0.0, maximum_downstream_share_difference_pp=max((max(x['downstream_max_share_difference_pp'].values()) for x in r2['sufficiency'])), counterexamples=r2['sufficiency'], smallest_extension='named tooth compliance and load/preload; regional force+area for regional traction; full pressure field for peak pressure'), intervention=dict(relief_um=30.0, successful=len([row for row in r2['rows'] if row.get('crown_edit', {}).get('gate')]), attempted=r2['central_accepted'], selected_region_force_reduction_N_range=[min((row['crown_edit']['reduction_N'] for row in r2['rows'] if 'crown_edit' in row)), max((row['crown_edit']['reduction_N'] for row in r2['rows'] if 'crown_edit' in row))], scope='conditional fixed-patch gap edit; no actual manufactured crown'), support_friction_ablation=dict(max_no_crown_effect_pp=max((row.get('scenario_max_share_difference_pp', {}).get('no_crown') or 0 for row in r2['rows'])), max_mu0_effect_pp=max((row.get('scenario_max_share_difference_pp', {}).get('mu0') or 0 for row in r2['rows'])), max_mu04_effect_pp=max((row.get('scenario_max_share_difference_pp', {}).get('mu04') or 0 for row in r2['rows']))), dropout=dict(force_roots=r2['force_root_counts'], attempted=12, central_rejected_case=9, unloaded_pose_compatible=r1['unloaded_pose_compatible'], normal_candidate_filter_mean_fraction=r1['dropout']['patch_normal_filter_mean_fraction'], omitted_source_faces_mean_fraction=float(np.mean([s['omitted_face_fraction'] for f in (P / 'raw').glob('geometry_*_h02.json') for s in json.loads(f.read_text())['sources']]))), resolution={'reaction_and_mean_traction': 'PER_SURFACE_REGION', 'source_witness': 'PER_POINT', 'force_share_argmax': 'PER_TOOTH', 'external_profiles': 'POPULATION of PER_TOOTH measurements', 'PDL_mu_nu_krotation': 'PHENOMENOLOGICAL'}, time_scale='SIMULTANEOUS', rigorous_enclosure='MISSING; NCP residuals and continuous affine minima are scoped float64 controls; no affine/linear sensitivity claim', phenomenological_debts=[dict(quantity='PDL6DOF support', replacement_measurement='per-tooth simultaneous force/wrench and6DOF motion at multiple independent directions and force/time levels'), dict(quantity='friction_mu', replacement_measurement='same surfaces and saliva condition tangential force/slip onset with normal force'), dict(quantity='preload and unloaded gap', replacement_measurement='paired registered scans at documented force with independent reference pose/preload'), dict(quantity='patch disk area/normal elasticity', replacement_measurement='regional pressure map and crown indentation or resolved FE comparison')], data=dict(directory=str(D), bytes=disk, limit_bytes=3000000000, manifest=manifest), figure=dict(path=str(fp), sha256=sha(fp)), region_port=dict(csv=str(export), csv_sha256=sha(export), vtk=str(P / 'regional_force_case001.vtk')), verification=dict(path=str(P / 'VERIFICATION.json'), count=verify['count'], all_pass=verify['all_pass']), cost=costs, physical_measurement_status='NOT_RUN', next_construction='Independent registered loaded/unloaded geometry + calibrated regional force observation and per-tooth support/preload; expand continuous contact candidate set before crown surface export. Frozen measurement comparator available.')
    r4 = json.loads((P / 'rounds/R4/results.json').read_text())
    result['rounds']['R4'] = dict(result_path=str(P / 'rounds/R4/results.json'), sha256=sha(P / 'rounds/R4/results.json'))
    result['frozen_predictions']['FROZEN_PREDICTIONS_R4.json'] = sha(P / 'FROZEN_PREDICTIONS_R4.json')
    result['contact_incidence_expansion'] = dict(new_regions=[r['expanded_patch_count'] - r['original_patch_count'] for r in r4['rows']], new_active_counts=[r['new_active_count'] for r in r4['rows']], force_share_differences_pp=[r['max_tooth_share_difference_pp'] for r in r4['rows']], all_whole_set_gates_pass=r4['all_whole_set_gates_pass'], scope=r4['rigorous_enclosure'])
    result['finite_scenario_ranges_pp'] = {}
    for row in r2['rows']:
        qs = [v for v in row['scenarios'].values() if 'force_control' in v]
        if 'force_control' in row['force_targets']['100']:
            qs.append(row['force_targets']['100'])
        if qs:
            result['finite_scenario_ranges_pp'][str(row['case'])] = {jaw: {t: [min((s['force_shares_pp'][jaw].get(t, 0) for s in qs)), max((s['force_shares_pp'][jaw].get(t, 0) for s in qs))] for t in qs[0]['force_shares_pp'][jaw]} for jaw in ('upper', 'lower')}
    result['uncertainty'] = 'Finite scenario ranges are descriptive, not rigorous parameter/population intervals. Individual physical uncertainty remains UNKNOWN.'
    costs['R4_wall_s'] = r4['cost']['wall_s']
    write(P / 'results.json', result)
    write(P / 'COSTS.json', costs)
    write(P / 'DATA_MANIFEST.json', dict(bytes=disk, files=manifest))
    write(P / 'ATTEMPTS.json', [dict(round=r['round'], parent_target='DENT-VAL-OCCLUSAL-REGIONAL-FORCES', claim_type=r['claim_type'], result_file=str(P / 'rounds' / r['round'] / 'results.json'), outcome=r.get('answer'), next_operation=n) for (r, n) in [(r1, 'force-controlled boundary and support sufficiency'), (r2, 'continuous regional source-triangle minima'), (r3, 'independent physical measurement/preload, contact-set expansion')]])
    state('ROUND_PACKAGE_COMPLETE', str(verify['count']) + ' numerical/interface checksPASS; physical/profile and original2pp refinement gates failed orUNKNOWN', 'Independent review, then loaded-pose regional-force/support acquisition; do not refit held-out force profiles', result_file='results.json', physical_measurement_status='NOT_RUN', review_state='PENDING_INDEPENDENT_REVIEW')
    print('Packaged', P / 'results.json', 'data', disk, 'bytes', flush=True)
if __name__ == '__main__':
    main()
