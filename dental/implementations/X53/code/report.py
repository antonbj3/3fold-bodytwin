from dental_release.paths import expand as _release_expand
import csv, collections, resource
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def report():
    r1 = read(ROOT / 'RESULTS_R1.json')
    r2 = read(ROOT / 'RESULTS_R2.json')
    r3 = read(ROOT / 'RESULTS_R3.json')
    verify = read(ROOT / 'raw/DEMO_VERIFICATION.json')
    source = read(ROOT / 'raw/SOURCE_CHECKS.json')
    r5 = read(ROOT / 'RESULTS_R5.json')
    r4 = read(ROOT / 'RESULTS_R4_V2.json')
    coverage = read(ROOT / 'RESULTS_SURFACE_COVERAGE.json')
    stats = []
    bigger = 0
    ball_lost = 0
    N = 0
    gain = 0
    for family in ['crown_loop', 'GenCAD_Foundation', 'GenCAD_V2']:
        rows = [r for r in r1['rows'] if r['family'] == family]
        vals = {}
        for axes in [4, 5]:
            v = [next((v for v in r['variants'] if v['library'] == 'vhf' and v['axes'] == axes))['by_region']['intaglio'] for r in rows]
            vals[str(axes)] = float(np.mean([r['found_fraction'] for r in v]))
            vals['ball'] = float(np.mean([r['ball_only_fraction'] for r in v]))
        stats.append(dict(family=family, designs=len(rows), intaglio_points=len(rows) * 64, **vals, resolution='PER_SURFACE_REGION', population_mean_scope='Unweighted mean of design-level sample fractions, not whole-population prevalence'))
    per_point = []
    for file in r1['raw_files']:
        a = read(file)
        v4 = next((v for v in a['points_and_tools'] if v['library'] == 'vhf' and v['axes'] == 4))
        v5 = next((v for v in a['points_and_tools'] if v['library'] == 'vhf' and v['axes'] == 5))
        for (i, p) in enumerate(v5['records']):
            per_point.append(dict(design_id=a['id'], point_index=i, region=p['region'], x_mm=p['point'][0], y_mm=p['point'][1], z_mm=p['point'][2], smallest_found_green_diameter_mm=p['smallest_found_green_mm'], largest_found_green_diameter_mm=p['largest_found_green_mm'], status=p['status'], resolution='PER_POINT'))
            if p['region'] == 'intaglio':
                N += 1
                t = p['tools']
                bigger += t[0]['status'] != 'FOUND' and any((x['status'] == 'FOUND' for x in t[1:]))
                ball_lost += any((x['ball_only'] for x in t)) and p['smallest_found_green_mm'] is None
                gain += p['smallest_found_green_mm'] is not None and v4['records'][i]['smallest_found_green_mm'] is None
    with open(ROOT / 'PER_POINT_TOOL_MAP.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(per_point[0]))
        w.writeheader()
        w.writerows(per_point)
    with open(ROOT / 'PER_POINT_CEMENT_PORT.csv', 'w', newline='') as f:
        fields = ['design_id', 'point_index', 'x_mm', 'y_mm', 'z_mm', 'normal_x', 'normal_y', 'normal_z', 'normal_gap_nominal_um', 'stock_safe_ball_union_um', 'forced_signed_stock_um', 'normal_gap_after_safe_union_um', 'normal_gap_after_forced_union_um', 'stock_status', 'fit_status', 'resolution', 'timescale']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for path in r2['raw_files']:
            a = read(path)
            for (i, p) in enumerate(a['records']):
                row = {k: p.get(k) for k in fields if k in p}
                row.update(design_id=a['id'], point_index=i, x_mm=p['point'][0], y_mm=p['point'][1], z_mm=p['point'][2], normal_x=p['normal'][0], normal_y=p['normal'][1], normal_z=p['normal'][2], resolution='PER_POINT', timescale='HANDOVER')
                w.writerow(row)
    export_rows = [dict(r, access_hoeffding_95_halfwidth=r['hoeffding_95_halfwidth'], stock_fraction_confidence_bound='MISSING_DEPENDENT_SPARSE_UNION') for r in coverage['rows']]
    local_summaries = {(r['id'], r['region']): r for r in r5['rows']}
    for r in export_rows:
        r.pop('hoeffding_95_halfwidth')
        for k in ['stock_gt25um_fraction', 'safe_union_overcut_gt25um_fraction', 'forced_overcut_gt25um_fraction', 'unknown_stock_fraction']:
            r['legacy_global_first_' + k] = r.pop(k)
        r.update(local_summaries[r['id'], r['region']])
    fields = sorted(set().union(*(r.keys() for r in export_rows)))
    with open(ROOT / 'PER_DESIGN_SURFACE_FRACTIONS.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(export_rows)
    fields = ['legacy_global_first_stock_um', 'design_id', 'point_index', 'region', 'x_mm', 'y_mm', 'z_mm', 'normal_x', 'normal_y', 'normal_z', 'normal_gap_nominal_um', 'stock_safe_ball_union_um', 'forced_signed_stock_um', 'normal_gap_after_safe_union_um', 'normal_gap_after_forced_union_um', 'stock_status', 'fit_status', 'resolution', 'timescale']
    with open(ROOT / 'PER_POINT_STOCK_PORT.csv', 'w', newline='') as allf, open(ROOT / 'PER_POINT_CEMENT_PORT.csv', 'w', newline='') as innerf:
        wa = csv.DictWriter(allf, fieldnames=fields)
        wi = csv.DictWriter(innerf, fieldnames=fields)
        wa.writeheader()
        wi.writeheader()
        for path in r5['raw_files']:
            a = read(path)
            for (i, p) in enumerate(a['records']):
                row = {k: p.get(k) for k in fields if k in p}
                row.update(design_id=a['id'], point_index=i, x_mm=p['point'][0], y_mm=p['point'][1], z_mm=p['point'][2], normal_x=p['normal'][0], normal_y=p['normal'][1], normal_z=p['normal'][2], stock_safe_ball_union_um=p['local_component_stock_um'], forced_signed_stock_um=p['local_component_forced_stock_um'], normal_gap_after_safe_union_um=p['local_component_normal_gap_um'], normal_gap_after_forced_union_um=p['local_component_forced_normal_gap_um'], stock_status=p['local_component_status'], resolution='PER_POINT', timescale='HANDOVER')
                wa.writerow(row)
                if p['region'] == 'intaglio':
                    wi.writerow(row)
    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(221)
    ys = np.arange(len(stats))
    width = 0.23
    for (key, lab, col, shift) in [('ball', 'Ideal ball only', '#aaaaaa', -width), ('4', 'Whole tool, 4 axes', '#2878b5', 0), ('5', 'Whole tool, 5 axes', '#e88729', width)]:
        ax.barh(ys + shift, [100 * s[key] for s in stats], height=width, label=lab, color=col)
    ax.set_yticks(ys, [s['family'] for s in stats])
    ax.set_xlim(0, 100)
    ax.set_xlabel('Found intaglio sample fraction (%)')
    ax.set_title('Declared centre-neck scenario, scale 0.8')
    ax.legend(fontsize=9)
    ax.text(0.01, -0.32, '64 area-uniform points/region/design; per-design 95% sampling half-width ±17 pp.\nNo physical machine or full-surface certificate.', transform=ax.transAxes, fontsize=9)
    ax2 = fig.add_subplot(222, projection='3d')
    a = read(DATA / 'R1/crownloop_D1.json')
    v = next((v for v in a['points_and_tools'] if v['library'] == 'vhf' and v['axes'] == 5))
    pp = [p for p in v['records'] if p['region'] == 'intaglio']
    P = np.array([p['point'] for p in pp])
    d = np.array([p['smallest_found_green_mm'] or np.nan for p in pp])
    known = np.isfinite(d)
    sc = ax2.scatter(*P[known].T, c=d[known], vmin=0.6, vmax=2, cmap='viridis', s=25)
    ax2.scatter(*P[~known].T, c='black', marker='x', label='No grid witness', s=30)
    ax2.set_title('D1: minimum found nominal cutter')
    ax2.set_xlabel('x (mm)')
    ax2.set_ylabel('y (mm)')
    ax2.set_zlabel('z (mm)')
    fig.colorbar(sc, ax=ax2, pad=0.12, shrink=0.65, label='Tool diameter, green mm')
    ax2.legend(fontsize=8)
    ax3 = fig.add_subplot(223)
    lengths = sorted((r['upper_green_mm'] for r in r3['rows'] if r['upper_green_mm'] is not None))
    ax3.plot(np.arange(1, len(lengths) + 1), lengths, color='#2878b5', label='Found-pose upper bracket')
    ax3.axhline(3, color='#d84b4b', ls='--', label='Published nominal small-tool L3')
    ax3.set_xlabel('Larger-only points sorted by required reach')
    ax3.set_ylabel('Required neck reach, green mm')
    ax3.set_title('Preserve the intaglio: invert the tool requirement')
    ax3.legend(fontsize=9)
    ax3.text(0.01, -0.27, '111 digital points; bracket width ≤0.01 mm.\nNeck taper/holder/fixture and cutting rigidity remain UNKNOWN.', transform=ax3.transAxes, fontsize=9)
    ax4 = fig.add_subplot(224)
    A = np.array([60.0, 60.0])
    B = np.array([20.0, 100.0])
    x = np.arange(2)
    ax4.bar(x - 0.18, A, 0.36, label='Arrangement A', color='#2878b5')
    ax4.bar(x + 0.18, B, 0.36, label='Arrangement B', color='#e88729')
    ax4.axhline(40, color='gray', ls=':', label='40 µm scenario threshold')
    ax4.set_xticks(x, ['Point class 1', 'Point class 2'])
    ax4.set_ylabel('Local normal gap (µm)')
    ax4.set_title('Identical mean 60 µm; identity error 0.0')
    ax4.legend(fontsize=9)
    ax4.text(0.01, -0.27, 'Constructed fields: restricted fraction 0% versus 50%.\nIllustrates summary loss; no empirical cement seating prediction.', transform=ax4.transAxes, fontsize=9)
    fig.suptitle('X53: pointwise tool access → unchanged-fit tool requirement', fontsize=15)
    fig.subplots_adjust(left=0.14, right=0.96, hspace=0.62, wspace=0.36, top=0.91, bottom=0.11)
    fig.savefig(ROOT / 'figures/MILLING_TOOLS.png', dpi=170)
    fig.savefig(ROOT / 'figures/MILLING_TOOLS.pdf')
    plt.close(fig)
    (fig2, ax) = plt.subplots(figsize=(9, 3.6))
    ax.broken_barh([(-0.375, 0.125), (0.125, 0.125)], (1.1, 0.5), facecolors=['#aaaaaa', '#2878b5'])
    ax.broken_barh([(-0.375, 0.625)], (0.1, 0.5), facecolors='#e88729')
    ax.axvline(0, color='black', ls='--', label='Target t=0')
    ax.axvline(-0.375, color='gray', ls=':', label='Same global first root -0.375 mm')
    ax.scatter([0.125, -0.375], [1.35, 0.35], marker='x', s=100, color='red', label='Selected local boundary')
    ax.set_yticks([0.35, 1.35], ['Connected B', 'Disconnected A'])
    ax.set_xlabel('Oriented normal-ray coordinate (mm)')
    ax.set_title('Identical global root; local boundary differs by 0.500 mm')
    ax.set_xlim(-0.45, 0.32)
    ax.set_ylim(-0.1, 1.85)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=8)
    fig2.text(0.17, 0.02, 'Dyadic sphere intervals checked by exact rational oracle; not manufactured stock.', fontsize=9)
    fig2.subplots_adjust(left=0.19, bottom=0.33, top=0.85)
    fig2.savefig(ROOT / 'figures/LOCAL_RAY_COMPONENT.png', dpi=170)
    fig2.savefig(ROOT / 'figures/LOCAL_RAY_COMPONENT.pdf')
    plt.close(fig2)
    limits = ['No independently measured installed CAM/tool/holder/fixture or physical milling truth; all physical gates UNKNOWN', 'Minimum found on a finite pose grid is an upper bound on the true minimum feasible library diameter; NOT_FOUND is not proof of impossibility', 'Tool neck/taper/holder and .8 green scale have labelled unknown/closure fields. Ceramill35deg is a scenario, not an OEM number. VHF actual cutter set is .6/1/2mm, not2.5mm', '60 crown surface fractions are Monte Carlo estimates; R4 certifies only the planar intaglio of488 roof designs by exact rational support, excluding exterior/axial/fixture', 'Crown-loop source SDF pitch120um and intaglio classifier can affect local results; no guaranteed geometric error enclosure', 'Sparse final-ball union is not full CAM swept stock; fixed-normal gap omits subsequent seating/tilt/cement/sinter errors', 'Rigorous floating-point enclosure of the collision/ray operation is MISSING; dense independent geometric controls passed stated1e-6mm tolerance', 'Stock fraction confidence bounds MISSING: sparse ball unions are built from the same sample; R1 fixed-oracle Hoeffding bound cannot transfer', 'Source reference-frame gate failed: frozen centre neck=L3 omits tip-to-centre conversion; source-amended tool card is provided for new scenes', 'No cutting force, long-neck stiffness, wear or approved material/tool compatibility established']
    output = dict(claim_type='capability', status='DIGITAL_POINT_QUERY_AND_INVERSE_TOOL_SPECIFICATION_DELIVERED', review_state='PENDING_INDEPENDENT_REVIEW', external_referent=dict(kind='published_code', locator='https://github.com/mikedh/trimesh/blob/4.11.3/trimesh/proximity.py', compared_quantity='continuous segment-to-retained-triangle minimum distance checked against independent dense proximity plus Lipschitz enclosure; not physical milling access', refutes_us=True), external_referents=[dict(kind='published_dataset', locator='https://doi.org/10.5281/zenodo.10597292', compared_quantity='Crown-loop source anatomy and derived meshes; no CAM measurement', refutes_us=True), dict(kind='published_dataset', locator='https://ditto.ing.unimore.it/bits2bites/', compared_quantity='Source paired arch meshes in GenCAD V2; geometry only', refutes_us=True), dict(kind='independent_measurement', locator='https://doi.org/10.4047/jap.2016.8.6.439 Table2', compared_quantity='Manufactured regional gap in um', refutes_us=True, comparison_status='UNKNOWN_UNPAIRED_DESIGN_AND_MACHINE', eligible_rows=0, excluded_rows=3)], requested_designs=61, accepted_designs=60, dropout=dict(count=1, fraction=1 / 61, reason='V2 missing construction', excluded_unflared_replicas_outside_denominator=31), per_point_queries=60 * 128, per_intaglio_queries=N, larger_only_count=bigger, ball_only_without_wholetool_witness_count=ball_lost, five_axis_extra_witness_count=gain, found_fraction_summary=stats, minimum_neck=dict(count=r3['bracketed'], range_green_mm=r3['required_neck_green_mm_range'], median_green_mm=r3['required_neck_green_mm_median'], bracket_mm_max=r3['max_bracket_mm'], resolution='PER_POINT', condition='Legacy centre-neck scenario; nominal L3 must add D1/2 and new source-converted fields are not recomputed. Finite-pose grid bracket only: lower has no grid witness, upper has a saved witness; no continuous-pose minimum or viable cutting-load prediction'), cement_handover=dict(cases=20, points=1280, new_standoff_witnesses=sum((r['stand_off_new_witnesses'] for r in r2['rows'])), unknown_stock_points=int(sum((r['unknown_stock_fraction'] * r['n'] for r in r2['rows']))), forced_overcut_gt25um_points=int(sum((r['forced_overcut_gt25um_fraction'] * r['n'] for r in r2['rows']))), resolution='PER_POINT', timescale='HANDOVER', physical_validation=False), sufficiency_tests=[read(ROOT / 'raw/SUFFICIENCY_R1.json'), r2['controls']], failed_sufficiency_attempt=read(ROOT / 'raw/SUFFICIENCY_R1_ATTEMPT1.json'), control_outcomes=dict(unit_tests=8, unit_test_failures=0, parent_replay=read(ROOT / 'raw/PARENT_REPLAY.json')['all_match'], original_parent_verdicts=24, source_parameter_checks=5, demo_verification=verify, inverse_upper_witness_replays=r3['all_upper_witnesses_replayed'], inverse_corruption_rejected=r3['all_short_neck_injections_rejected']), full_cost=dict(preparation_search_and_reading='UNMEASURED; source code/mesh prep not hidden in speedup', fit_seconds=0, discovery_wall_seconds=dict(R1=r1['wall_seconds_this_invocation'], R2=r2['seconds_this_invocation'], R3=r3['seconds']), parent_replay_seconds=read(ROOT / 'raw/PARENT_REPLAY.json')['seconds'], R1_peak_RSS_MiB=r1['peak_rss_MiB'], R2_peak_RSS_MiB=r2['peak_rss_MiB'], physical_CAM_jobs=0, physical_scan_cost='UNKNOWN_NOT_RUN', questions=0, GPU_seconds=0, fallback='Conditional witness / grid NOT_FOUND / physical UNKNOWN', own_bytes=budget()), phenomenological_debts=[dict(parameter='installed neck radial profile + gauge + holder', replacement_measurement='tool microscopy/profile scan and installed gauge/holder scan'), dict(parameter='A/B mounting and fixture', replacement_measurement='CAM transforms plus fixture/stock scan'), dict(parameter='green shrink .8', replacement_measurement='actual blank batch factor and paired green/sintered part')], limits=limits, artifacts=dict(figure='figures/MILLING_TOOLS.png', point_tools='PER_POINT_TOOL_MAP.csv', cement_port='PER_POINT_CEMENT_PORT.csv', lab_cases='LAB_CASES_20.json', lab_protocol='LAB_PROTOCOL_20.md'), raw_data=[dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size) for p in DATA.rglob('*') if p.is_file()], frozen_predictions=[dict(path=str(p), sha256=sha(p)) for p in ROOT.glob('FROZEN_PREDICTIONS*.json')], graph=dict(target='DENT-MFG-ASBUILT-DEVIATION', dispatch_status='FAILED_EXISTING_HISTORY_RESULT_HASH_MISMATCH', source_graph_modified=False))
    output['full_GenCAD_benchmark'] = dict(requested_tasks=576, counts=r4['counts'], exact_scope='All points on serialized planar inner triangles only; no axial/exterior/fixture certificate', fraction_without_design_certificate=88 / 576, reasons_without_certificate={'ABSTAIN_SOURCE_CONSTRUCTION': 84, 'UNKNOWN_SOURCE_SITE': 4}, physical_gate='UNKNOWN', resolution='PER_POINT', rational_support_enclosure=True, all_numerical_controls_pass=r4['all_numerical_controls_pass'], all_injections_rejected=r4['all_injections_rejected'], evidence='RESULTS_R4_V2.json', original_wrong_source_labels_preserved='raw/RESULTS_R4_V1_SOURCE_ELIGIBILITY_FAILED.json')
    output['surface_reporting'] = dict(designs=60, regions=120, points=7680, intaglio_film_points=3840, region_resolution='PER_SURFACE_REGION', point_resolution='PER_POINT', same20_R2_intaglios_unchanged=True, physical_gate='UNKNOWN', stock_interpretation='Separate sparse final-ball union per region; threshold25um and unknown fraction retained; forced unions are diagnostics, not CAM compensation', unknown_stock_points=int(round(sum((r['unknown_stock_fraction'] * r['n'] for r in coverage['rows'])))), fractions_file='PER_DESIGN_SURFACE_FRACTIONS.csv', raw_result='RESULTS_SURFACE_COVERAGE.json')
    output['local_ray_branch'] = dict(points=7680, changed_safe_points=r5['changed_safe_points'], intaglio_changes=0, exterior_changes=1, old_safe_overcut_count=1, new_safe_overcut_count=0, scope='Local connected component, still sparse ball union', controls=r5['controls'], physical_gate='UNKNOWN', evidence='RESULTS_R5.json')
    output['sufficiency_tests'].append(r5['controls'])
    output['external_referents'].append(r5['external_referent'])
    output['full_cost']['discovery_wall_seconds']['R5'] = r5['seconds']
    output['nominal_tool_reference_frame'] = read(ROOT / 'raw/TIP_REFERENCE_CHECK.json')
    output['published_control_source'] = read(ROOT / 'sources/PROXIMITY_SOURCE_LOCK.json')
    output['surface_reporting']['uncertainty'] = read(ROOT / 'REPORTING_INTERPRETATION.json')
    for debt in output['phenomenological_debts']:
        debt['resolution'] = 'PHENOMENOLOGICAL'
    output['external_referents'].append(r4['external_referent'])
    output['full_cost']['discovery_wall_seconds']['R4'] = r4['seconds']
    output['full_cost']['discovery_wall_seconds']['requested_surface_coverage'] = coverage['seconds']
    output['full_cost']['validation_wall_seconds'] = 'UNMEASURED total; parent replay separately recorded'
    output['full_cost']['per_on_demand_query_seconds'] = 'UNMEASURED; finite tool/pose enumeration without fit'
    output['full_cost']['R4_peak_memory'] = 'UNMEASURED; no heavy run, small serialized roof tasks'
    output['artifacts'].update(stock_port='PER_POINT_STOCK_PORT.csv', surface_fractions='PER_DESIGN_SURFACE_FRACTIONS.csv', all_original_tasks='RESULTS_R4_V2.json', local_branch_figure='figures/LOCAL_RAY_COMPONENT.png')
    output['preserved_failure_versions'] = ['raw/SUFFICIENCY_R1_ATTEMPT1.json', 'raw/VERIFY_DEMO_V1_FAILED.log', 'raw/R4_RUN_V1_FAILED.log', 'raw/RESULTS_R4_V1_SOURCE_ELIGIBILITY_FAILED.json', 'SOURCE_CORRECTIONS.json', 'RESULTS_SURFACE_COVERAGE.json: one falsely local exterior overcut corrected by R5']
    dump(ROOT / 'results.json', output)
    table = '| Design family | Designs | Ideal ball | Whole tool, 4 axes | Whole tool, 5 axes |\n|---|---:|---:|---:|---:|\n' + '\n'.join((f"| {s['family']} | {s['designs']} | {100 * s['ball']:.1f}% | {100 * s['4']:.1f}% | {100 * s['5']:.1f}% |" for s in stats))
    readme = _release_expand(f"# Query which tool reaches, then preserve the local fit\n\nRun `./run_all.sh` from this directory. It verifies the archived input/prediction hashes, runs the controls, replays positive approaches, and writes the table, point maps and figure. First construction takes about seven minutes on these local inputs. Subsequent runs preserve frozen results and verify them. Python/library and data paths are in `RUNTIME_LOCK.json` and `INPUT_MANIFEST.json`; no network or GPU is required by the demo.\n\nThe change is a point × tool × machine-pose relation that retains the neck, shank and holder. A smallest cutter diameter by itself cannot answer access. The inverse query now specifies a neck-length bracket for the small cutter at zero stock offset, rather than compensating away the intaglio.\n\n{table}\n\nFractions are intaglio Monte Carlo queries, PER_SURFACE_REGION, in a declared centre-neck vhf tool/holder scenario at green-to-design scale0.8. Each design has64 random area-uniform intaglio points; its95% Hoeffding half-width is17.0 percentage points. Rows aggregate the same queries; they are not clinical population statistics.\n\nOf{N} intaglio queries, {bigger} have a larger-tool witness without a0.6mm witness. {ball_lost} fit an ideal ball but have no whole-tool pose found. The fifth axis adds{gain} point witnesses to the fourth-axis set. The{r3['bracketed']} inverse small-tool requirements have found-pose upper values{r3['required_neck_green_mm_range'][0]:.3f}–{r3['required_neck_green_mm_range'][1]:.3f}mm and bracket width≤0.01mm (PER_POINT, green scale); the published nominal small-tool neck L3 is3mm. This is geometric reach, without measured long-neck stiffness or approved tool compatibility.\n\nAll60 designs have per-region sampled remaining-stock/overcut fractions in `PER_DESIGN_SURFACE_FRACTIONS.csv`, with unknowns retained in the denominator. The ±17 percentage-point bound applies only to the fixed zero-offset access query; the sampled ball-union stock fractions have no established area confidence bound. `PER_POINT_STOCK_PORT.csv` has7680 points;3840 intaglio points carry signed local stock into a **fixed normal-ray** gap (`PER_POINT_CEMENT_PORT.csv`). Two constructed states have exactly the same mean film60µm, identity error0.0, yet the local fraction below the40µm scenario threshold changes0%→50%. The smallest sufficient extension for this pair is the spatially registered nominal film plus signed stock. No cement seating height or hydraulics is predicted from the mean.\n\nThe separate full GenCAD V2 task audit covers576 original tasks:488 designs have an **exact rational supporting-halfspace certificate for their entire planar intaglio**,84 source constructions abstain and4 source sites remain UNKNOWN. The certificate uses a vertical approach and proves every modeled tool capsule stays below the retained-mesh supporting plane. It does not certify the roof exterior, clinical axial walls or a mounted blank. Same-information numerical checks pass at four inner vertices/design and an injected1mm green gauge fails every exact certificate. `RESULTS_R4_V2.json` retains all task/output hashes and inner-vertex tool maps. The earlier R4 version merged the four missing-site labels into abstention; its source-eligibility failure and frozen version are preserved.\n\n![Pointwise tool access and inverse neck requirement](figures/MILLING_TOOLS.png)\n\nA fifth construction fixes the **ray branch** representation. The first-root reduction once labelled a collision-free exterior pose as overcut: it selected a disconnected cutting ball beyond a thin wall. R5 preserves the connected interval component at the target. It changes1/7680 safe-stock answers (one exterior point), removes that false overcut and leaves every intaglio safe-stock answer unchanged. A pair has identical first root−0.375mm with identity error0.0, yet local boundaries+0.125/−0.375mm differ by0.500mm. An exact rational sphere-interval oracle derived from the [sphere equation](https://mathworld.wolfram.com/Sphere.html) agrees, and injected50µm errors fail. The active CSVs use R5 local components; frozen R2/global-root fields remain available as legacy evidence. Physical stock remains UNKNOWN.\n\n![Why a global first root loses the local surface](figures/LOCAL_RAY_COMPONENT.png)\n\nWhat holds against an external referent: independent [trimesh proximity code](https://github.com/mikedh/trimesh/blob/4.11.3/trimesh/proximity.py) brackets the continuous geometry distances using a dense-sample Lipschitz enclosure. Five nominal input checks match the original [VHF catalogue, p19](https://henryschein.com.au/Documents/Product%20Docs/VHF/Brochure-Tools.pdf), [VHF machine specs](https://www.vhf.com/en-us/products/dental/dental-milling-machines/k5/) and [Ha & Cho methods](https://doi.org/10.4047/jap.2016.8.6.439). All24 original GenCAD ideal-ball verdicts replay. Eight controls and injected clearance, holder, impossible cavity, short-neck and CAM errors are rejected where specified. These validate the declared digital operation and source parameters.\n\nWhat does not hold: **physical milling and fit remain UNKNOWN**. The independent measured regional gaps in [Ha & Cho Table2](https://doi.org/10.4047/jap.2016.8.6.439) lack paired CAD, installed assembly and local correspondence:0/3 cells are eligible as physical facit, all three exclusions retained. No CAM job or scan was run. Neck taper, gauge, holder, fixtures, mounting and scale remain measured-data debts. Ceramill±35° is an explicit scenario; only its five axes and0.6/1.0/2.5mm cutters are sourced. The vhf diameter set here is0.6/1.0/2.0mm. A source-frame check fails the frozen neck interpretation: p19 L3 ends at the tool tip, whereas frozen centre-to-neck lengths were set equal to L3. `NOMINAL_TIP_REFERENCE_CARD.json` supplies L3−D1/2 for new queries; old point fields remain the declared centre-neck scenario. Do not treat them as source-converted nominal-tool or measured-machine predictions. The drawing also restricts short diamond tools to specific machine models, including K4 edition rather than original K4. A missing grid witness cannot certify global impossibility. Surface coverage is sampled; collision/ray arithmetic has no rigorous floating-point enclosure. Crown-loop SDF pitch120µm is resolution, not an error guarantee. Forced-overcut fields and sparse final-ball unions are diagnostics, not executable CAM stock simulation.\n\nUse `LAB_PROTOCOL_20.md`, frozen hashes and `@DENTAL_PYTHON@ code/compare_cam.py CAM_LOG_TEMPLATE.csv` for the twenty-case lab handover. The empty template returns UNKNOWN; an injected wrong observation at a positive witness returns FAIL. Query a new local mesh with `@DENTAL_PYTHON@ code/query_mesh.py mesh.stl --point x y z --normal nx ny nz --axes 5 --library vhf`; for signed stock use `STOCK_LOG_TEMPLATE.csv` and `@DENTAL_PYTHON@ code/compare_stock.py STOCK_LOG_TEMPLATE.csv`; it retains discrepancy intervals and rejects an injected50µm observation. Optional `--fixture` and`--tool-card` retain measured geometry inputs. Axes are mounted about local x/y; verify that transform on the actual machine.\n\nData licences: STS-Tooth3D [Zenodo](https://doi.org/10.5281/zenodo.10597292), CC BY4.0; GenCAD V2 reuses local [Bits2Bites](https://ditto.ing.unimore.it/bits2bites/) and Bite2Text-derived meshes under the parent's licence lock (Bits2Bites CC BY-NC-SA; Bite2Text terms remain as recorded by the parent). Teeth3DS-derived Foundation geometry uses [the source terms](https://osf.io/xctdy/), with CC BY-SA reported in the parent. These are attributed source-derived research designs, not released person records. ToothFairy2 and mandibular-defect data were not used. Article/OEM source documents retain their individual rights and are archived privately. No new dataset or code publication is authorized. `SOURCE_CORRECTIONS.json` withdraws an unverified guessed V2 DOI while preserving the frozen original metadata.\n")
    for (a, b) in [('Of' + str(N), 'Of ' + str(N)), ('Of3840', 'Of 3840')]:
        readme = readme.replace(a, b)
    import re
    chunks = re.split('(\\[[^\\]]*\\]\\([^)]*\\)|`[^`]*`)', readme)
    for i in range(0, len(chunks), 2):
        chunks[i] = re.sub('(?<=[a-z])(?=\\d)', ' ', chunks[i])
        chunks[i] = re.sub('(?<=\\d)(?=(?:mm|µm|um|random|intaglio|point|cells|original|inverse|physical))', ' ', chunks[i])
    readme = ''.join(chunks).replace('Bits 2Bites', 'Bits2Bites').replace('Bite 2Text', 'Bite2Text').replace('Teeth 3DS', 'Teeth3DS').replace('STS-Tooth 3D', 'STS-Tooth3D').replace('ToothFairy 2', 'ToothFairy2')
    readme = re.sub('(?<=[;:,])(?=\\d)', ' ', readme)
    readme = readme.replace('optional`', 'optional `').replace('and`', 'and `').replace('CC BY4.0', 'CC BY 4.0')
    (ROOT / 'README_DEMO.md').write_text(readme)
    (ROOT / 'RESULTS.md').write_text(f"# Pointwise access and inverse cutter requirement\n\nThe executable capability now provides local compatible-tool witnesses on60 existing designs, and a zero-stock neck-length requirement on111 points where the nominal small cutter had no found pose. It is checked against published nominal tool/machine sources and independent triangle-proximity code. The576-task GenCAD audit additionally certifies the complete planar intaglio of488 roof models under exact rational halfspace assumptions, retaining84 abstentions and4 missing-source UNKNOWNs. Physical milling/fit remains UNKNOWN because no paired CAM scene or scan was available.\n\n{table}\n\nEach region has64 sample points per design. PER_POINT results, witnesses, rejections and signed cement handovers are retained in the raw data and CSVs; no aggregate certifies an entire surface. See README_DEMO.md for the figure and validity limits.\n\nR1: {bigger}/{N} larger-only intaglio witnesses, {ball_lost}/{N} ball-only false positives relative to the declared full-tool grid, {gain} extra five-axis witnesses. R2:199 new stand-off witnesses,50/1280 unknown sparse-ball stock points,4/1280 forced-ball overcut points above25µm. R3:{r3['bracketed']} neck brackets, range{r3['required_neck_green_mm_range'][0]:.3f}–{r3['required_neck_green_mm_range'][1]:.3f}mm, maximum width{r3['max_bracket_mm']:.6f}mm.\n\nRequested coverage: all60 meshes have both sampled regions,7680 signed stock points and3840 intaglio film handovers. The original20 R2 intaglio maps remain unchanged. The17pp R1 access bound is not a stock-fraction bound because the sparse union depends on the sampled points. Fractions above25µm and UNKNOWN rates are in PER_DESIGN_SURFACE_FRACTIONS.csv; they describe sparse ball unions, not measured CAM stock. R4:488 exact planar-intaglio certificates;1952 numerical pose checks pass and488 injected early-holder cases fail. The unmodeled exterior/axial/fixture remains UNKNOWN.\n\nR5 fixes a failed representation gate: one exterior point had a negative global root from a disconnected far-side cut sphere. Local connected-component selection removes this false-overcut label; all3840 intaglio safe-stock answers are unchanged. Same-first-root dyadic pair: identity error0.0mm, local-boundary difference0.5mm, exact rational oracle error0.0mm; injected50µm rejected. Generated CSVs now use local components and preserve legacy root fields.\n\nSource-frame negative result: the frozen neck card used published tip-referenced L3 as a centre-referenced reach. Old fields are preserved as scenarios. NOMINAL_TIP_REFERENCE_CARD.json converts L3−D1/2; three exact identity checks pass and injected unconverted fields fail. Source numeric equality alone did not validate the coordinate interpretation.\n\nControls: eight tests pass;24/24 unchanged parent verdicts replay; five nominal source checks pass; all111 inverse upper poses replay and injected3mm replacements fail. Equal-information control is verification, not a method-superiority claim. Three published fit cells excluded for missing paired local geometry/machine; dropout100%. Design dropout1/61,1.64%, missing V2 construction. Thirty-one unflared source replicas were outside the denominator.\n\nFailed sufficiency v1: identical .6mm tip and planar curvature,2mm versus6mm pockets; downstream difference0.0, gateFAIL. Capsule-rounded shank caused both to collide. The retained v2 state pair1mm/6mm has identity error0.0 and downstream access-fraction difference1.0. Film pair: identical60µm mean to machine precision, downstream restricted-fraction difference0.5. These fixtures are mathematical instruments and not independent physical facit.\n\nDiscovery wall time:{r1['wall_seconds_this_invocation']:.2f}s R1, {r2['seconds_this_invocation']:.2f}s R2, {r3['seconds']:.2f}s R3, {r4['seconds']:.2f}s R4, {coverage['seconds']:.2f}s requested surface coverage, {r5['seconds']:.2f}s R5. Preparation/search cost unmeasured; no speedup claim. R1 peak RSS{r1['peak_rss_MiB']:.2f}MiB; GPU0. All actual artefact hashes and costs are in results.json. Pending independent review.\n")
    result_text = (ROOT / 'RESULTS.md').read_text()
    result_text = re.sub('(?<=[a-z])(?=\\d)', ' ', result_text)
    result_text = re.sub('(?<=\\d)(?=(?:mm|µm|MiB|points|intaglio))', ' ', result_text)
    result_text = re.sub('(?<=[;:,])(?=\\d)', ' ', result_text)
    (ROOT / 'RESULTS.md').write_text(result_text)
    packet = read(ROOT / 'GRAPH_PACKET.json')
    feedback = dict(packet.get('result_binding_template', {}))
    feedback.update(target_id='DENT-MFG-ASBUILT-DEVIATION', result_file='results/LANE_X53_MILLING_TOOLS/results.json', sha256=sha(ROOT / 'results.json'), review_state='PENDING_INDEPENDENT_REVIEW', measured_quantity='PER_POINT declared full-tool access and inverse neck upper bracket; signed fixed-ray stock/film handover', units='mm, um, count, fraction', uncertainty='Finite orientation set and64-point region sampling; 95% Hoeffding half-width0.170; physical geometry and rigorous arithmetic enclosure UNKNOWN', population_regime='60 source-derived crown meshes;3840 local film points;111 inverse points;576 GenCAD tasks including488 exact planar-intaglio certificates; conditional neck/holder/.8 scale', preregistered_gate='PREREG_R1/R2/R3/R4_V2/R5 and PREREG_COVERAGE_R2 immutable hashes; preserved R4v1 source-eligibility gate failed; scoped digital controls PASS; physical gates UNKNOWN', baseline='Original ideal-ball GenCAD checker replay24/24; independent published trimesh proximity; no algorithm-superiority contest', outcome='DIGITAL_CAPABILITY_DELIVERED_PHYSICAL_UNKNOWN', negative_result=True, claim_type='capability', resolution='PER_POINT', timescale='HANDOVER', negative_findings=['112 ball-only points lack whole-tool grid witness', '3/3 independent physical fit cells ineligible', 'Original sufficiency fixture gate failed', 'R4v1 missed four missing-source labels; unchanged metrics in R4V2', 'R2 global ray root mislabelled one exterior point; R5 local components correct it', 'Nominal L3 coordinate interpretation failed; legacy fields are scenarios, converted card provided separately', 'dispatch blocked by existing history Result hash mismatch'], proposed_edges=[dict(from_artifact='PER_POINT_TOOL_MAP.csv', to_target='DENT-MFG-ASBUILT-DEVIATION', resolution='PER_POINT', timescale='SIMULTANEOUS'), dict(from_artifact='PER_POINT_CEMENT_PORT.csv', to_consumer=_release_expand('X13'), resolution='PER_POINT', timescale='HANDOVER', status='conditional geometric input only; hydraulic closure UNKNOWN')])
    dump(ROOT / 'GRAPH_FEEDBACK.json', feedback)
    dump(ROOT / 'ATTEMPTS.json', [dict(round='R1', changed_operation='Whole-tool finite-pose query', outcome='Digital capability; original sufficiency pair failed, revised state pair rejects radius summary', next_operation='R2 local stock/film handover', evidence='RESULTS_R1.json'), dict(round='R2', changed_operation='Signed ball-union normal-ray handover', outcome='20 local maps, identical mean insufficient; physical fit UNKNOWN', next_operation='R3 invert required neck rather than change cavity', evidence='RESULTS_R2.json'), dict(round='R3', changed_operation='Inverse nested shank-envelope neck bracket', outcome='111 found-pose upper brackets and111 short-neck injections rejected', next_operation='R4 exact supporting-halfspace certificate on all576 original tasks', evidence='RESULTS_R3.json'), dict(round='R4', changed_operation='Exact rational halfspace support on original576 task outputs', outcome='488 complete planar-intaglio certificates,84 abstentions,4 source UNKNOWN; initial source-eligibility label failure preserved', next_operation='R5 local connected ray components', evidence='RESULTS_R4_V2.json'), dict(round='R5', changed_operation='Local connected normal-ray cut component', outcome='One false exterior overcut removed; identical global-root pair fails sufficiency by.5mm; mathematical oracle passes', next_operation='R6 measured scene and new prediction freeze before20 actual CAM jobs', evidence='RESULTS_R5.json')])
    state('FIVE_DIGITAL_CONSTRUCTIONS_DELIVERED', 'Scoped digital scenario gates PASS; nominal tip/centre reference gate FAILED for frozen library; physical geometry/CAM/fit UNKNOWN', 'R6: measure actual radial neck profile/gauge/holder/fixture and A/B mount; freeze a new scene before running20 CAM jobs')
    print('reported', output['status'], bigger, N, r3['required_neck_green_mm_range'])
if __name__ == '__main__':
    report()
