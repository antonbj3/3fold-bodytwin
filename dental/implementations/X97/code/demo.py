import sys, os
from pathlib import Path
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R))
import json, hashlib, time, resource
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from occlusion_module import analyze, load_demo_input, query_height
from occlusion_module.api import clean
from occlusion_module._vendor import contact, quality, reachability

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, x):
    p = R / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(clean(x), indent=2, allow_nan=False) + '\n')

def run():
    phase = os.environ.get('X97_DEMO_ROUND', 'R2')
    latest = phase == 'R2'
    start = time.perf_counter()
    cohort = json.loads((R / 'inputs/cohort.json').read_text())
    old = json.loads((R / 'inputs/X95_ROUND1.json').read_text())['rows']
    old2 = json.loads((R / 'inputs/X95_ROUND2.json').read_text())['rows']
    old3 = json.loads((R / 'inputs/X95_ROUND3.json').read_text())['rows']
    rows = []
    parities = []
    for (row, truth, tr2, tr3) in zip(cohort, old, old2, old3):
        (b, c) = load_demo_input(R, row)
        out = analyze(b, c, {'local_boundary': latest})
        target = tr3 if latest else truth
        adj = out['proposed_adjustment']
        m = out['contact_map']['metrics']
        after = adj['contact']
        pos = out['contacts_requiring_buildup']
        field = np.array(out['contact_map']['gap_mm'], dtype=float)
        xy = np.array(b['xy_mm'])
        ff = np.array(b['faces'])
        ind = quality.contact_map(xy, ff, field, b['reference_gap_mm'])
        diffs = dict(original_mm2=abs(m['symdiff_mm2'] - truth['original']['symdiff_mm2']), after_mm2=abs(after['symdiff_mm2'] - target['outputs']['regional']['contact']['symdiff_mm2']), missing_patch_mm2=abs(pos['within_removal_cap']['unavoidable_missing_contact_mm2'] - tr2['nominal']['unavoidable_missing_contact_mm2']), max_relief_mm=abs(adj['gates']['max_relief_mm'] - target['outputs']['regional']['wall_and_removal']['max_relief_mm']), independent_polygon_mm2=abs(m['symdiff_mm2'] - ind['contact_symdiff_mm2']), protected_identity_mm=adj['gates']['protected_identity_error_mm'])
        good = (not latest or adj['gates']['local_vertical_interval_order_pass']) and all((v <= 1e-07 for (k, v) in diffs.items() if k.endswith('mm2'))) and (diffs['max_relief_mm'] <= 1e-08) and (diffs['protected_identity_mm'] == 0)
        path = R / 'exports' / phase / (row['key'] + '.npz')
        path.parent.mkdir(exist_ok=True, parents=True)
        np.savez_compressed(path, vertices_mm=np.array(adj['edited_vertices_mm']), faces=np.array(adj['faces']), height_change_mm=np.array(adj['height_change_per_vertex_mm']), xy_mm=xy, query_faces=ff, gap_mm=field, reference_gap_mm=b['reference_gap_mm'])
        minimal = {k: v for (k, v) in out.items() if k not in ['contact_map', 'proposed_adjustment']}
        minimal.update(contact_map={'metrics': m, 'uncertainty': out['contact_map']['uncertainty'], 'spatial_artifact': str(path.relative_to(R))}, proposed_adjustment={'status': adj['status'], 'contact': after, 'gates': adj['gates'], 'optimization': adj['optimization'], 'regional_faces': len(adj['regional_height_adjustments']), 'height_field_artifact': str(path.relative_to(R))})
        dump('raw/' + phase + '/' + row['key'] + '.json', minimal)
        rows.append(dict(key=row['key'], case_id=row['case_key'], source_fdi=row['source_fdi'], dataset=row['dataset'], before_symdiff_mm2=m['symdiff_mm2'], after_symdiff_mm2=after['symdiff_mm2'], unrecoverable_without_cap_mm2=pos['without_depth_cap']['unavoidable_missing_contact_mm2'], max_relief_mm=adj['gates']['max_relief_mm'], geometry_status=adj['gates'].get('removal_status_R3', adj['gates']['removal_status']), physical_force='UNKNOWN', parity=diffs, parity_pass=good, query_triangles=out['coverage']['query_triangles'], retained_triangles=out['coverage']['retained_triangles'], reference_supported_triangles=out['coverage']['reference_supported_triangles'], reference_common_area_fraction=out['coverage']['reference_common_area_fraction'], artifact=str(path.relative_to(R)), artifact_sha256=sha(path), artifact_bytes=path.stat().st_size))
        print(row['key'], 'parity', good, 'unrecoverable', pos['without_depth_cap']['unavoidable_missing_contact_mm2'], flush=True)
        dump('raw/DEMO_PROGRESS.json', {'rows': rows, 'completed': len(rows), 'requested': len(cohort)})
    w = reachability.sufficiency(contact)
    dump('raw/SUFFICIENCY.json', w)
    s = json.loads((R / 'inputs/F5367_R4_state.json').read_text())
    binding = dict(case_id=s['case'], response_geometry_sha256=s['geometry_sha256'], pose_id='calibration0', calibration_pose_id='calibration0', matched_load_rate_history=True)
    h = query_height(binding, s, dict(operation='uniform_tooth_height', fdi=16, height_change_mm=-0.02))
    i = s['predicted_fdi'].index(16)
    dump('raw/HEIGHT_QUERY.json', h)
    frozen = json.loads((R / 'FROZEN_PREDICTIONS.json').read_text())
    result = dict(schema='X97-results-v1', round=phase, claim_type='capability', outcome='COMPOSED_RESEARCH_API; PHYSICAL_FORCE_UNKNOWN', review_state='PENDING_INDEPENDENT_REVIEW', external_referent=dict(kind='published_dataset', locator='https://ditto.ing.unimore.it/bits2bites/', compared_quantity='Retrospective spatial source-gap-band contact, common-support polygon area, not loaded contact pressure', refutes_us=False), cohort_requested=10, cohort_retained=len(rows), rejected_records=0, rejected_fraction=0, region_dropout=dict(requested_triangles=sum((x['query_triangles'] for x in rows)), retained_triangles=sum((x['retained_triangles'] for x in rows)), reason='Missing common crown/antagonist/source field support'), rows=rows, r1_preserved=json.loads((R / 'raw/R1_RESULTS.json').read_text())['gate_pass'], gate_pass=all((x['parity_pass'] for x in rows)), crowns_requiring_positive_height=sum((x['unrecoverable_without_cap_mm2'] > 1e-07 for x in rows)), height_query=dict(case=s['case'], fdi=16, edit_mm=-0.02, force_interval_N=h['force_interval_N'][i], evidence_kind='SIMULATED_PROSPECTIVE_ACQUISITION', physical_validity='UNKNOWN', resolution='PER_TOOTH'), sufficiency=w, controls=dict(direct_upstream='Archived original operator calls match; no algorithm advantage claimed', independent_polygon_max_mm2=max((x['parity']['independent_polygon_mm2'] for x in rows)), part_test_log='raw/TEST_R1.stderr'), limitations=['R1 global vertical separator incompleteness preserved; latest local constraints change four proposals and remain conditional on original solid validity.', 'Source native scanner/pose error and all floating enclosures MISSING.', 'No physiological or manufacturing validation; operational 100 micrometre band is not a pressure map.', 'No force can be inferred from contact area alone.'], cost=dict(wall_seconds=time.perf_counter() - start, max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, source_fit='None in X97; source constitutive closures retained', source_preparation='X95 23MB derivative reused; prior source extraction/discovery cost UNKNOWN', threads=1, GPU=False, physical_measurements=0), prereg_sha256=sha(R / ('PREREG_' + phase + '.json')), source_manifest_sha256=sha(R / 'SOURCE_MANIFEST.json'), frozen_predictions_sha256=sha(R / ('FROZEN_PREDICTIONS_R2.json' if latest else 'FROZEN_PREDICTIONS.json')), measurement_protocol='MEASUREMENT_PROTOCOL.md')
    result['region_dropout']['grid_rejected_fraction'] = 1.0 - result['region_dropout']['retained_triangles'] / result['region_dropout']['requested_triangles']
    result['region_dropout']['reference_supported_triangles'] = sum((x['reference_supported_triangles'] for x in rows))
    result['region_dropout']['reference_rejected_triangle_fraction'] = 1.0 - result['region_dropout']['retained_triangles'] / result['region_dropout']['reference_supported_triangles']
    result['region_dropout']['source_area_coverage_range'] = [min((x['reference_common_area_fraction'] for x in rows)), max((x['reference_common_area_fraction'] for x in rows))]
    result['external_negative_referent'] = dict(kind='independent_measurement', locator='doi:10.1111/jopr.13838 Table 1; inputs/X94_EXTERNAL_VALIDATION.json/source_quantile_table', compared_quantity='178 participant MIP force-share population quantiles; patient/crown force identification not supported', refutes_us=True)
    result['inherited_negative_results'] = {'X95_robust_geometric_advantage': 'FAIL preserved; no physical or method-superiority upgrade', 'X54_physical_force_transfer': 'FAIL/UNKNOWN preserved', 'X94_geometry_force_identification': 'FAIL preserved'}
    dump('results.json', result)
    dump('raw/' + phase + '_RESULTS.json', result)
    (fig, axes) = plt.subplots(1, 2, figsize=(11, 4))
    x = np.arange(len(rows))
    axes[0].bar(x - 0.2, [r['before_symdiff_mm2'] for r in rows], 0.4, label='Original crown')
    axes[0].bar(x + 0.2, [r['after_symdiff_mm2'] for r in rows], 0.4, label='Regional removal')
    axes[0].set_ylabel('Source contact symmetric difference (mm²)')
    axes[0].set_xlabel('Crown index (10 crowns, 4 source cases)')
    axes[0].legend()
    axes[1].bar(x, [r['unrecoverable_without_cap_mm2'] for r in rows])
    axes[1].set_ylabel('Source patch requiring positive height (mm²)')
    axes[1].set_xlabel('Crown index')
    fig.suptitle('Occlusion API: retrospective geometry only; physical force UNKNOWN')
    fig.tight_layout()
    fig.savefig(R / 'figures/demo.png', dpi=160)
    plt.close(fig)
    assert result['gate_pass'], 'Frozen parity gate failed; preserve raw failures'
    assert w['summary_identity_error'] == 0 and w['bit_identical'] and (w['downstream_difference_mm2'] > 0)
    dump('CURRENT_WORK_STATE.json', dict(lane='X97-occlusion-module', phase=phase + '_DECIDED', latest_gate={'parity': result['gate_pass'], 'retained': len(rows), 'requires_positive_height': result['crowns_requiring_positive_height']}, next_operation='Package and isolated replay; next physical construction is same-specimen calibrated regional height probes and held-out measurement', review_state='PENDING_INDEPENDENT_REVIEW'))
    print(json.dumps({'gate_pass': result['gate_pass'], 'retained': len(rows), 'requires_positive_height': result['crowns_requiring_positive_height'], 'seconds': result['cost']['wall_seconds']}))
if __name__ == '__main__':
    run()
