from dental_release.paths import expand as _release_expand
from cadlib import *
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def report():
    r1 = read(ROOT / 'RESULTS_R1.json')
    r2 = read(ROOT / 'RESULTS_R2.json')
    r3 = read(ROOT / 'RESULTS_R3.json')
    r5 = read(ROOT / 'raw/R5_SCORED.json')
    g = [r for r in r5 if r['status'] == 'GENERATED']
    r4 = read(ROOT / 'raw/R4B_SCORED.json')
    g4 = [r for r in r4 if r['status'] == 'GENERATED']
    ctrl = read(ROOT / 'raw/CONTROLS.json')
    suff = read(ROOT / 'raw/SUFFICIENCY.json')
    material = {m: dict(requested=sum((r['material'] == m for r in r5)), generated=sum((r['material'] == m for r in g)), wall_pass=sum((r['material'] == m and r['wall']['status'] == 'PASS' for r in g)), contact_pass=sum((r['material'] == m and all(r['V6_gates'].values()) for r in g)), technical_cap_pass=sum((r['material'] == m and r['digital_cap_status'] == 'PASS' for r in g)), shape_pass=sum((r['material'] == m and r['shape']['gate'] == 'PASS' for r in g))) for m in sorted({r['material'] for r in r5})}
    lookup = {r['uid']: r for r in g4}
    regress = [r['uid'] for r in g if all(lookup[r['uid']]['V6_gates'].values()) and (not all(r['V6_gates'].values()))]
    gains = [r['uid'] for r in g if lookup[r['uid']]['wall']['status'] != 'PASS' and r['wall']['status'] == 'PASS']
    shapevals = [r['shape']['p95_mm'] for r in g]
    e2 = [r['margin_error']['max_um'] for r in read(ROOT / 'raw/R2_SCORED.json') if r['cohort'] == 'X52' and 'margin_error' in r]
    ex = dict(read(ROOT / 'PREREG_R5.json')['external_referent'])
    ex['refutes_us'] = True
    costs = {}
    for name in ['R1', 'R2', 'R3', 'R4', 'R4B', 'R5']:
        p = ROOT / ('FROZEN_PREDICTIONS_' + name + '.json')
        if p.exists():
            a = read(p)
            costs[name] = {k: a[k] for k in ['seconds', 'peak_rss_MiB'] if k in a}
    result = dict(claim_type='capability', status='PARTIAL_CAPABILITY_PENDING_INDEPENDENT_REVIEW', capability='Executable automatic digital CAD proposal and falsifying geometry chain on source-linked conical preparations; clinical margin and physical release unresolved', external_referent=ex, external_referents=[dict(r1['external_referent'], refutes_us=True), suff['external_referent'], dict(kind='independent_measurement', locator=str(BASE / 'LANE_X13_CEMENT_GAP/FACIT.csv'), compared_quantity='Predecessor physical gap calibration findings constrain interpretation; not remeasured here', refutes_us=False)], summary=dict(R1_margin=r1['summary'], R2_margin=r2['summary'], R2_X52_max_error_median_um=float(np.median(e2)), R2_X52_error_range_um=[float(min(e2)), float(max(e2))], requested_material_cases=len(r5), generated=len(g), unique_generated_teeth=len({r['key'] for r in g}), case_clusters=len({r['case_key'] for r in g}), wall_pass=sum((r['wall']['status'] == 'PASS' for r in g)), nominal_min_max_film_pass=sum((r['min_film']['status'] == 'PASS' and r['max_film']['status'] == 'PASS' for r in g)), contact_pass=sum((all(r['V6_gates'].values()) for r in g)), technical_cap_pass=sum((r['digital_cap_status'] == 'PASS' for r in g)), technical_and_shape_pass=sum((r['digital_cap_status'] == 'PASS' and r['shape']['gate'] == 'PASS' for r in g)), shape_pass=sum((r['shape']['gate'] == 'PASS' for r in g)), shape_p95_median_mm=float(np.median(shapevals)), shape_p95_range_mm=[float(min(shapevals)), float(max(shapevals))], full_no_manual_qualified_count=0, full_no_manual_physical_rate='UNKNOWN_NOT_MEASURED', milling_point_witnesses=sum((r['milling']['found'] for r in g)), milling_sampled_points=sum((r['milling']['points'] for r in g)), export_3mf_count=len(g), by_material=material), changed_geometry=dict(wall_new_pass=gains, contact_regressions=regress, interpretation='Actual neighbour and ring deformation; wall improvement does not establish anatomical or physical validity'), dropout=dict(requested=len(r5), generated=len(g), unavailable=len(r5) - len(g), unavailable_fraction=(len(r5) - len(g)) / len(r5), reasons=dict(Counter((r['status'] for r in r5 if r['status'] != 'GENERATED'))), R2_reference_missing=2, commercial_adjustment_sources=dict(candidate_records=2, admitted=0, rejected_fraction=1.0, reasons=['Acceptance is not adjustment-free fraction', 'No retrieved primary adjustment-fraction table'])), sufficiency=dict(identity_error=suff['identity_error'], downstream=suff['downstream'], minimum_extension=suff['minimum_extension']), controls=dict(count=ctrl['count'], all_pass=ctrl['all_pass'], file='raw/CONTROLS.json'), cost=dict(numerical_rounds=costs, source_preparation=read(ROOT / 'INPUT_LOCK_R3.json')['seconds'], fit='Included in generation runtime; not separately instrumented', discovery_reasoning='UNKNOWN_NOT_INSTRUMENTED', questions=0, physical_validation='NOT_PERFORMED', fallback='R1 and R2 semantic failures preserved; R3 STL support failure; R4 LP scheduler failure; R4B runtime repair and R5 new boundary operation', data_bytes=budget(), lane_bytes=sum((p.stat().st_size for p in ROOT.rglob('*') if p.is_file())), GPU=0, resource_note='Initial runs set numerical env<=4; uninstrumented initial HiGHS native pool count UNKNOWN. Final runs initialise a one-thread HiGHS and pin 4 CPUs affinity.'), resolutions=dict(curve_and_collision='PER_POINT', contact_and_role='PER_SURFACE_REGION', case='PER_TOOTH', aggregate_counts='POPULATION descriptive only', material_gap_force_tool='PHENOMENOLOGICAL until matched measurements'), phenomenological_debts=[dict(quantity=q, replacement_measurement=m) for (q, m) in [('CAD80 µm gap to physical film', 'same-part die/crown scan plus dry/wet spatial film and process/material/batch'), ('height to regional force', 'X82/X85 same-specimen response with registered geometry/frame/basis, force and height metrology'), ('X89 mounted tool response', 'same tool load-deflection/force bound, holder and full CAM scene'), ('scan geometry', 'traceable scanner uncertainty and independent marked preparation finish line')]], limits=['R1/R2 margin references are simulated, not annotated clinical preparations', 'R4B/R5 require source-linked exact conical CAD; unsupported generic preparations refuse', 'Native reference morphology includes actual IOS target surfaces; bilateral shape thresholds are inherited engineering queries', 'All58 shape gates fail; technical cap passes exclude morphology and do not imply complete crown success', 'Minimum wall excludes annular cervical closing rim; ring thickening alters cervical contour', 'Self-intersections and whole-arch insertion sweep remain UNKNOWN', 'GenCADV6 contact is a represented PL field; source/roundoff enclosures missing', 'Local milling samples plus1N scenario are not whole-CAM or LS2 process qualification', 'No clinical, tenfold speed, no-adjustment rate or manufactured-part claim'], next_construction='Replace synthetic cervical-to-occlusal exterior with full mirrored donor axial surface and a source-bound finish-line collar; couple continuous wall and registered-neighbour exclusion. Acquire independent prepared-scan margin annotation and same-part gap/force measurements before release.', reproduction=dict(command='./run_all.sh', full_command='./run_all.sh --recompute', case_command='./run_all.sh --case examples/source_linked_case.json', scope='Default fresh2material reconstruction, selected margin replay, all58 geometry scores and26 controls; --recompute regenerates all58 and all130available margins.'))
    if (ROOT / 'RESULTS_R6.json').exists():
        result['full_donor_pilot'] = read(ROOT / 'RESULTS_R6.json')
    dump(ROOT / 'results.json', result)
    fig = plt.figure(figsize=(14, 4.6))
    ax = fig.add_subplot(131)
    labels = ['R1 feature', 'R2 cut']
    ax.bar(np.arange(2) - 0.18, [3, 99], 0.36, label='Closed proposals', color='#5596b5')
    ax.bar(np.arange(2) + 0.18, [0, 0], 0.36, label='Correct <=25 um', color='#de9a48')
    ax.set_xticks(range(2), labels)
    ax.set_ylim(0, 108)
    ax.set_title('X52: closure does not locate margin')
    ax.set_ylabel('Synthetic cases / 99')
    ax.legend(fontsize=8)
    a = np.load(g[0]['mesh_path'])
    ax = fig.add_subplot(132, projection='3d')
    cols = np.array(['#9fb8cd', '#e19775', '#7da387'])
    tri = a['vertices'][a['faces']]
    poly = Poly3DCollection(tri, facecolors=cols[a['roles']], linewidth=0.1, edgecolors='#34485c', alpha=0.9)
    ax.add_collection3d(poly)
    v = a['vertices']
    ax.set_xlim(v[:, 0].min(), v[:, 0].max())
    ax.set_ylim(v[:, 1].min(), v[:, 1].max())
    ax.set_zlim(v[:, 2].min(), v[:, 2].max())
    ax.set_box_aspect(np.ptp(v, axis=0))
    ax.set_title('Generated research shell\nExterior / intaglio / rim')
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')
    ax.set_zlabel('mm')
    ax.view_init(22, 125)
    ax = fig.add_subplot(133)
    counts = [len(g), result['summary']['wall_pass'], result['summary']['contact_pass'], result['summary']['technical_cap_pass'], result['summary']['shape_pass'], 0]
    labs = ['Generated', 'Wall', 'V6 contact', 'Technical cap', 'Native shape', 'Full qualified']
    ax.barh(labs[::-1], counts[::-1], color=['#bd6b66', '#bd6b66', '#de9a48', '#de9a48', '#5596b5', '#5596b5'])
    ax.set_xlim(0, 64)
    ax.set_xlabel('Material-case requests / 64')
    ax.set_title('Separate gates; not a success funnel')
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/demo.png', dpi=170)
    fig.savefig(ROOT / 'figures/demo.pdf')
    plt.close(fig)
    table = '\n'.join((f"| {m} | {x['generated']}/{x['requested']} | {x['wall_pass']} | {x['contact_pass']} | {x['technical_cap_pass']} | {x['shape_pass']} |" for (m, x) in material.items()))
    text = _release_expand(f"# Automatic CAD proposals with explicit limits\n\nThe tool generates source-linked research crowns from an unlabelled preparation die, a mirrored same-subject opposite tooth and registered surrounding geometry. It detects a closed margin, derives an insertion witness, constructs a nominal offset cavity, adjusts the exterior and exports millimetre 3MF with face regions. **It does not deliver a qualified crown without manual/experimental review: 0 complete passes are established.**\n\nRun `./run_all.sh`. Full regeneration: `./run_all.sh --recompute`. One editable example manifest: `./run_all.sh --case examples/source_linked_case.json`. The case command writes a research part under the lane data disk; its physical release is UNKNOWN.\n\n| Material scenario | Generated / requested | Wall PASS | V6 contact PASS | Technical cap PASS | Original-shape PASS |\n|---|---:|---:|---:|---:|---:|\n{table}\n\nCounts are PER_TOOTH/material scenario, aggregated descriptively over {len({r['key'] for r in g})} teeth in {len({r['case_key'] for r in g})} case records, not independent clinical trials. Technical cap means margin/topology/wall/nominal film/static antagonist/V6 predicates; it excludes shape and physical qualification. All {len(g)} original-shape comparisons fail against independent original IOS triangles (median sampled bidirectional p95 = {np.median(shapevals):.3f} mm). The full outcome is therefore negative. Wall repairs add {len(gains)} passes and lose {len(regress)} previous contact passes; the unchanged thresholds expose this tradeoff.\n\nOn the known synthetic preparations, sharp-cycle detection recovers 31/31 available cone rims (1 missing / 32 requests). On X52 it yields 3 closed but incorrect curves / 99. The global minimum-cut construction yields 99 closed curves, 97 reference comparisons fail 25 µm and 2 references are unavailable. Median sampled maximum error among the 97 is {np.median(e2):.1f} µm. These are fixture comparisons: there is no independent clinical finish-line annotation. A closure proof does not establish margin semantics.\n\n[Figure](figures/demo.png). Exact per-case values, triangle witnesses, contact regions, hashes and rejection reasons are in `results.json`, `raw/R5_SCORED.json` and the data paths it names. All earlier failed rounds remain. Nominal 40–120 µm film constraints are proved on declared convex caps; 80 µm is a declared CAD scenario for both product scenarios, not a material-specific measured setting. Product wall values retain the V2 local manufacturer-source binding.\n\nWhat does not hold: actual cement film and force intervals are UNKNOWN without X13/X82/X85 calibration; {result['summary']['milling_point_witnesses']}/{result['summary']['milling_sampled_points']} local tool witnesses are sampled 1 N scenarios, not CAM qualification. Source uncertainty, self-intersections, full-arch insertion motion and clinical contour remain unresolved. Cervical rim is excluded from minimum wall. Initial STL rounding breaks convex support at fractions of a micrometre; source-CAD recovery is allowed only after facet bijection and a 0.01 µm serialization bound. It is unavailable for arbitrary clinical scans. These limitations are not hidden behind the {result['summary']['technical_cap_pass']} technical-cap passes.\n\nExact-summary counterexample: equal endpoints, height and mean radius with identity error 0 produce feasible versus impossible insertion; the second shape has a 0.5 mm waist rebound. The minimal distinguishing extra scalar for this pair is the positive radius rebound; arbitrary 3D cases need the oriented normal field and incidence. The unchanged exact Farkas verifier and a conventional same-information LP both run. 26 controls include deliberately wrong units, region IDs, open curves, Farkas multipliers, wall thickness, film limit, contact value and tool force; all reject.\n\nCommercial comparison: [3Shape reports more than 94% design acceptance](https://www.3shape.com/en/services/automate), with review and optional editing. This is not an adjustment-free fraction. No comparable commercial no-manual/no-chairside-adjustment percentage was established; see `LITERATURE.md`. No 10× speed claim is made.\n\nData/licences: Bits2Bites CC BY-NC-SA inherited local notice; Bite2Text exact archive licence UNKNOWN, local private research only; STS-derived X52 geometry CC BY 4.0. Predicted FDI labels are inherited. Teeth3DS, ToothFairy2 and mandibular-defect meshes are not newly used here. No personal records are read and nothing is published. Native geometry stays under `@DENTAL_WORK_ROOT@/PROOF_LANE_CAD_10X/`; parent sources are read-only. Region colours are display/region labels within one intended physical material.\n\nDependencies are pinned in `DEPENDENCIES.json`; the installed Newton Python environment is required. No GPU; final runs pin 4 CPUs and initialise a one-thread LP. Physical lab specification is `LAB_PROTOCOL.md`. Status PENDING_INDEPENDENT_REVIEW.\n")
    if (ROOT / 'RESULTS_R6.json').exists():
        text += '\nFull-donor pilot R6: three predetermined cases replace the roof with the complete reflected donor mesh. All three refuse at the source-topology gate (multiple/branching boundaries or component loss), before any shape score. These are retained failures, not three inaccurate scored crowns. The next construction must recover a source-supported tooth boundary from the original arch before stitching. See RESULTS_R6.json and HANDOFF_R6.md.\n'
    (ROOT / 'README_DEMO.md').write_text(text)
    (ROOT / 'RESULTS.md').write_text(text)
    return result
if __name__ == '__main__':
    report()
