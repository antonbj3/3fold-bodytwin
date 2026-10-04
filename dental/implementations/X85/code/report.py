from common import *
import time, platform, csv

def run():
    import matplotlib.pyplot as plt
    r0 = read(ROOT / 'raw/R0_SUFFICIENCY.json')
    r1 = read(ROOT / 'raw/R1_UNIFORM_INVERSE.json')
    r2 = read(ROOT / 'raw/R2_ROWS.json')
    r3 = read(ROOT / 'raw/R3_ROWS.json')
    r4 = read(ROOT / 'raw/R4_ROWS.json')
    r5 = read(ROOT / 'raw/R5_ROWS.json')
    external = read(ROOT / 'raw/EXTERNAL_MEASUREMENT.json')
    unique = [r for r in r4 if r['participant'] == 'volume_1.0']
    roof = [r for r in unique if r.get('diagnostics', {}).get('manufacturing_status') == 'CONDITIONAL_ROOF_ONLY']
    patterns = [r for r in roof if r['contact_retention'] == 'PASS']
    positive = [r for r in patterns if r['contact_after']['reference']['area_mm2'] > 1e-12]
    representative = roof[0]
    inverse = representative['diagnostics']['inverse']
    force = np.array(inverse['at_selected_height']['force_N'])
    R = inverse['at_selected_height']['joint_l2_error_radius_N']
    tooth_sum = float(force[:2].sum())
    tooth_radius = float(np.sqrt(2) * R)
    controls = {f'R{i}': read(ROOT / f'raw/R{i}_CONTROLS.json') for i in [2, 4, 5]}
    byte_total = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    package_bytes = sum((p.stat().st_size for p in ROOT.rglob('*') if p.is_file()))
    result = dict(schema='regional-occlusal-interval-design-v1', lane='X85-occlusal-adjust', claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', status='CONDITIONAL_REGIONAL_DESIGN_AND_FULL_MESH_CONTACT_PROBE; PHYSICAL_FORCE_AND_SAVED_CHAIR_ADJUSTMENT_UNKNOWN', external_referent=read(ROOT / 'PREREG_R2.json')['external_referent'], additional_external_referents=[external['external_referent'], dict(kind='external_review', locator=str(ROOT.parent / 'LANE_XREVIEW_BATCH23/evidence/X82_INDEPENDENT.json'), compared_quantity='X82 uniform tooth-height conditional model intervals; not physical force', refutes_us=True), dict(kind='our_own_fixture', locator=str(ROOT / 'raw/probes'), compared_quantity='Regional simulated linear force vectors and interval calibration; not an external facit', refutes_us=False)], headline=dict(unique_roof_scenarios=16, source_case_clusters=4, roof_force_wall_antagonist_ideal_tool_passes=len(roof), also_external_contact_gate_passes=len(patterns), positive_native_contact_among_those=len(positive), actual_whole_crown_meshes=len(r5), whole_crown_wall_bound_passes=sum((r['continuous_wall_gate'] == 'PASS' for r in r5)), whole_crown_contact_masks_changed=sum((r['contact_change']['symdiff_mm2'] > 1e-09 for r in r5)), physical_crown_forces_validated=0, chairside_adjustment_saved_um=None), representative=dict(uid=representative['uid'], region_height_coefficients_um=np.array(inverse['coefficients_mm']) * 1000, maximum_surface_relief_um=representative['planned_relief_max_um'], area_mean_surface_relief_um=representative['planned_relief_mean_um'], unadjusted_region_force_N=[12, 8], adjusted_region_force_interval_N=inverse['at_selected_height']['force_interval_N'][:2], region_target_N=[7, 9], tooth_aggregate_interval_N=[tooth_sum - tooth_radius, tooth_sum + tooth_radius], joint_l2_error_radius_N=R, resolution='PER_SURFACE_REGION; geometry displacement PER_POINT; area average/aggregate PER_TOOTH', planned_relief_moved_to_CAD_um=representative['planned_relief_max_um'], model_remaining_adjustment_um=0.0, measured_chair_adjustment_saved_um=None, force_provenance='Own spring fixture; no matched loaded scan or physical crown measurement'), rounds={f'R{i}': read(ROOT / f'raw/R{i}_SUMMARY.json') for i in [2, 3, 4, 5]}, uniform_X82=dict(conditional_target_passes=r1['conditional_target_passes'], requested=r1['requested'], reuse='Existing inverse, not a new algorithm', region_transfer_passes=0), sufficiency=[r0, dict(summary='Same silicone count8 at tooth16', identity_error=external['contact_count_identity_error'], downstream_difference_pp=external['relative_force_difference_percentage_points'], smallest_extension='Same-session calibrated force observation, region-resolved for regional design', referent=external['external_referent'], resolution='PER_TOOTH')], controls=controls, control_outcome='Full-coordinate QP and unchanged GenCAD v4/v6 polygons agree at recorded tolerances; no method advantage claim', duplicate_accounting=dict(R2_rows=32, R2_available_surfaces=29, R4_rows=32, R4_unique_geometries=16, reason='Two predecessor labels collapse to the same new planar reference; count16 unique scenarios, not32 independent designs', R5_meshes=10, R5_case_clusters=len({r['case_key'] for r in r5})), attrition=dict(R2_missing_surfaces=3, R2_missing_fraction=3 / 32, R2_fixed_wall_inverse_rejections=25, R2_fixed_wall_rejection_fraction=25 / 29, R2_full_tool_certified=0, R3_failed=1, R4_unique_antagonist_or_geometry_rejections=16 - len(roof), R4_unique_external_contact_rejections_among_geometric_passes=len(roof) - len(patterns), R5_absolute_force_transfer_rejections=10, R5_absolute_force_transfer_rejection_fraction=1.0, external_absolute_reactions_admitted=0, reasons='Original3 missing X1B surfaces; active roof wall; conservative exterior tool proof; nominal source contact mismatch; no same-specimen force calibration'), resolution=dict(height='PER_POINT', force='PER_SURFACE_REGION', contact='PER_SURFACE_REGION', tooth_aggregate='PER_TOOTH', case_counts='POPULATION', support_and_targets='PHENOMENOLOGICAL'), edges=[dict(source='Measured IOS source/antagonist triangles', target='GenCAD projected gap and contact polygons', resolution='PER_SURFACE_REGION', timescale='SIMULTANEOUS', status='GEOMETRIC_ONLY'), dict(source='Same-geometry region height probes + calibrated reaction observations', target='Force interval inverse', resolution='PER_SURFACE_REGION', timescale='SIMULTANEOUS', status='SIMULATED; physical UNKNOWN'), dict(source='Regional design mesh, die and manufacturing plan', target='As-built lab geometry and heldout force test', resolution='PER_POINT', timescale='HANDOVER', status='PROSPECTIVE_NOT_RUN'), dict(source='Full crown force total', target='Regional patch force or point pressure', resolution='PER_SURFACE_REGION', timescale='SIMULTANEOUS', status='BLOCKED_NOT_IDENTIFIED')], phenomenological_debts=[dict(quantity='Linear reciprocal regional support', replacement_measurement='5 exact-operation force-vector probes plus independent heldout and declared load/rate/history; a physical curvature bound or restricted accepted corridor remains necessary'), dict(quantity='Exact heights and0.05N channel bound', replacement_measurement='Calibrated displacement, force and observation cross-talk bounds; real height-denominator enclosure required'), dict(quantity='Equal half-neighbor force target', replacement_measurement='Same-session spatial/native contact-force profile; not pooled population norm'), dict(quantity='Source gap vs fixture force binding', replacement_measurement='Loaded pose and spatial contact identity in the same specimen; present spring fixture is separate from IOS scan pose'), dict(quantity='Complete tool holder/fixture access', replacement_measurement='CAM sweep plus as-built scan of actual full crown')], limitations=['No absolute patient force or clinical recommendation', 'No measured chairside micrometres or adjustment time saved', 'R4 roofs are planar specimen shells, not complete anatomical crowns; baseline probe forces come from a separate declared spring fixture', 'Empty and nonempty source masks pass the same nominal contact threshold; six of nine distinct passing roofs have positive native contact', 'R5 full-mesh20um intervention is a geometric probe, not a force-selected design', 'All stated force-affine bounds are rigorous only under the explicitly declared linear full-contact law; physical linearity UNKNOWN', 'Formal floating-point, source geometry and physical nonlinear enclosures MISSING', 'R3 die editing invalidates physical calibration; numerical LP minimum checked, not a formal rounding certificate', 'Fullmesh topology self-intersection and full holder/fixture CAM remain UNKNOWN'], cost=dict(stage_seconds={f'R{i}': read(ROOT / f'raw/R{i}_SUMMARY.json')['seconds'] for i in [2, 3, 4, 5]}, R0_R1_seconds=r1['seconds'], peak_stage_rss_MiB=max((read(ROOT / f'raw/R{i}_SUMMARY.json').get('peak_rss_MiB', 0) for i in [2, 3, 4, 5])), fit='Simulated5 reaction-vector probes per distinct reference; physical fit0', validation='Full-force QP, source polygons, live injected faults and14 contract tests', preparation_upstream_cost='UNKNOWN', discovery_cost='Current stage times recorded; source-reading/token/historical costs UNKNOWN', questions='12 uniform inverse cases,32 inherited roof labels,29 available surfaces,16 unique new planes,10 actual fullmeshes', fallback='All failed/UNKNOWN cases retained', physical_acquisitions=0, GPU_seconds=0, threads_max=4), data_artifacts=dict(root=DATA, bytes=byte_total, limit_bytes=3000000000, files={str(p): dict(sha256=sha(p), bytes=p.stat().st_size) for p in sorted(DATA.glob('*.npz'))}), package_bytes=package_bytes, licenses=dict(Bits2Bites='CC BY-NC-SA, inherited local source binding; official author page verifies registered scans but not source archive licence text', Bite2Text='UNKNOWN in upstream archive/author page; private use only, no redistribution', STS_Tooth3D='CC BY4.0; X1B ancestor only, no new STS data exported', PMC10958416='CC BY4.0 as local XML licence; only factual selected numeric cells used', ToothFairy2='Unused', Teeth3DS='No mesh used; predecessor predicted FDI licence binding unresolved', mandible_defects='Unused'), source_direction_changes=[dict(source='Wiechens Table1', change='None to force law; confirmed that same count does not identify force and that T-Scan cannot silently provide absoluteN', own_calculation_before_web='raw/R0_SUFFICIENCY.json and raw/R1_UNIFORM_INVERSE.json')])
    if (ROOT / 'raw/TEST_RECEIPT.json').exists():
        result['test_receipt'] = read(ROOT / 'raw/TEST_RECEIPT.json')
    if (ROOT / 'raw/DEMO_RECEIPT.json').exists():
        result['demo_receipt'] = read(ROOT / 'raw/DEMO_RECEIPT.json')
    for i in [6, 7]:
        result['rounds'][f'R{i}'] = read(ROOT / f'raw/R{i}_SUMMARY.json')
        result['controls'][f'R{i}'] = read(ROOT / f'raw/R{i}_CONTROLS.json')
        result['cost']['stage_seconds'][f'R{i}'] = result['rounds'][f'R{i}']['seconds']
    result['bounded_height_recovery'] = dict(removal_only=read(ROOT / 'raw/R6_SUMMARY.json'), signed_generation=read(ROOT / 'raw/R7_SUMMARY.json'), representative=read(ROOT / 'raw/R7_REPRESENTATIVE.json'), height_error_budget_um=0.5, physical_metrology_budget='PHENOMENOLOGICAL, not measured CAD/CAM precision')
    result['headline'].update(bounded_height_removal_only_feasible=0, signed_generation_enclosed=120, signed_generation_worlds=120)
    result['edges'][0] = dict(source='Measured IOS source/antagonist triangles', target='GenCAD signed gap at each projected point', resolution='PER_POINT', timescale='SIMULTANEOUS', status='GEOMETRIC_ONLY')
    result['edges'].insert(1, dict(source='Projected point-gap field', target='Pointwise contact indicator; consumer assembles polygons/regions', resolution='PER_POINT', timescale='SIMULTANEOUS', status='GEOMETRIC_ONLY; regional aggregation belongs to consumer'))
    result['cost']['validation'] = 'Full-force QP, source polygons, live injected faults and17 contract tests'
    dump(ROOT / 'results.json', result)
    (fig, ax) = plt.subplots(2, 2, figsize=(11, 8))
    b = np.array([12.0, 8.0])
    iv = np.array(inverse['at_selected_height']['force_interval_N'])[:2]
    cent = iv.mean(1)
    ax[0, 0].bar(np.arange(2) - 0.16, b, 0.3, label='Unadjusted simulated')
    ax[0, 0].bar(np.arange(2) + 0.16, cent, 0.3, yerr=(iv[:, 1] - iv[:, 0]) / 2, label='Designed interval')
    ax[0, 0].axhspan(7, 9, color='green', alpha=0.15)
    ax[0, 0].set_xticks([0, 1], ['Region A', 'Region B'])
    ax[0, 0].set_ylabel('N · PER_SURFACE_REGION')
    ax[0, 0].set_title('Conditional spring fixture; no physical force facit')
    ax[0, 0].legend(fontsize=8)
    with np.load(representative['artifact_path'], allow_pickle=False) as a:
        xy = a['xy']
        f = a['faces']
        change = (a['unadjusted_z'] - a['adjusted_z']) * 1000
    im = ax[0, 1].tripcolor(xy[:, 0], xy[:, 1], f, change, shading='gouraud')
    fig.colorbar(im, ax=ax[0, 1], label='Removal µm · PER_POINT')
    ax[0, 1].set_aspect('equal')
    ax[0, 1].set_xlabel('x mm')
    ax[0, 1].set_ylabel('y mm')
    ax[0, 1].set_title('Regional surface operation')
    vals = [r['contact_change']['symdiff_mm2'] for r in r5]
    ax[1, 0].bar(np.arange(len(vals)) + 1, vals)
    ax[1, 0].set_xlabel('Full crown geometry probe 1–10')
    ax[1, 0].set_ylabel('Contact-mask change mm² · PER_SURFACE_REGION')
    ax[1, 0].set_title('Actual mesh edit; same source/antagonist')
    ax[1, 1].bar(['Original\nsculpted', 'Plane\ngeometry', 'Plane +\nsource contact', 'Physical\nforce'], [0, len(roof), len(patterns), 0])
    ax[1, 1].set_ylabel('Distinct scenarios · POPULATION')
    ax[1, 1].set_title('Separate gates; no clinical crown acceptance')
    fig.tight_layout()
    (ROOT / 'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT / 'figures/demo.png', dpi=160)
    plt.close(fig)
    with (ROOT / 'raw/DESIGN_TABLE.csv').open('w', newline='') as ff:
        writer = csv.writer(ff)
        writer.writerow(['uid', 'geometry_force_tool', 'native_contact', 'positive_native_area_mm2', 'mean_relief_um', 'max_relief_um', 'measured_chair_saved_um'])
        for r in unique:
            writer.writerow([r['uid'], r.get('diagnostics', {}).get('manufacturing_status'), r.get('contact_retention'), r.get('contact_after', {}).get('reference', {}).get('area_mm2'), r.get('planned_relief_mean_um'), r.get('planned_relief_max_um'), 'UNKNOWN'])
    state('DEMO_PACKAGED', 'Conditional regional design; fullmesh contact probe; measured chair saving UNKNOWN', 'Independent review, then same-specimen regional height/force acquisition under LAB_PROTOCOL.md')
    return result
if __name__ == '__main__':
    run()
