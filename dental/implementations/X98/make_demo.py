"""One-command replay; source snapshot read-only and portable."""
import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[k] = '1'
import argparse
import csv
import datetime
import hashlib
import json
import resource
import subprocess
import sys
import time
from pathlib import Path
from implant_safety import SafetyModule
from implant_safety.module import sha, load
import validation
ROOT = Path(__file__).resolve().parent

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(exist_ok=True, parents=True)
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def state(stage, gate, next_operation):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X98-implant-safety-module', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), stage=stage, latest_gate=gate, next_operation=next_operation))

def run(output):
    start = time.perf_counter()
    output.mkdir(exist_ok=True, parents=True)
    module = SafetyModule(ROOT)
    pr = load(ROOT / 'PREREG_R1.json')
    rows = []
    table = []
    for site in module.sites.values():
        for system in pr['selection']['systems']:
            for guide in pr['selection']['guides']:
                q = dict(site_id=site['site_id'], cbct_pose=module.default_pose(site['site_id']), guide_type=guide, implant_system=system, length_mm=site['pose']['length_mm'], confidence=0.95)
                r = module.query(**q)
                rows.append(r)
                g = r['terms']['guide_combination']
                geo = r['digital_geometry']
                c = r['combined_conditional']
                table.append(dict(site_id=site['site_id'], guide=guide, system=system, length_mm=q['length_mm'], nominal_gap_mm=r['reference_rule']['nominal_interval_mm'][0], annotation_loss_upper_mm=r['terms']['annotation_revision_loss']['interval_mm'][1], extra_depth_mm=r['terms']['drill_extra_depth']['maximum_or_example_mm'], tool_loss_upper_mm=r['terms']['directed_tool_gap_loss']['interval_mm'][1], swept_gap_mm=geo['tool_envelope']['interval_mm'][0], guide_entry_mm=g['entry_mm'], guide_apex_mm=g['apex_mm'], guide_angle_deg=g['angle_deg'], guide_rotation_mm=g['rotation_mm'], guide_combined_mm=g['combined_displacement_mm'], conditional_gap_lower_mm=c['interval_mm'][0], conditional_gap_upper_mm=c['interval_mm'][1], conditional_margin_lower_mm=c['margin_to_2mm_interval_mm'][0], rule_class=r['reference_rule']['class_vs_2mm'], tool_class=geo['tool_class_vs_2mm'], conditional_class=c['class_vs_2mm'], physical_status=r['physical_safety']['status'], resolution='PER_TOOTH/PHENOMENOLOGICAL'))
    state('R1_NUMERICS_COMPLETE', 'DIGITAL_QUERIES_COMPUTED', 'Execute independent controls and fault injections')
    (output / 'QUERIES.jsonl').write_text(''.join((json.dumps(x, ensure_ascii=False, allow_nan=False) + '\n' for x in rows)))
    with (output / 'TABLE.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)
    replay = []
    for s in module.sites.values():
        r = next((r for r in rows if r['site_id'] == s['site_id']))
        for (name, src) in [('nominal', s['published_nominal_interval']), ('revision_union', s['published_union_interval'])]:
            b = r['digital_geometry'][name]['interval_mm']
            err = max(abs(b[0] - src[0]), abs(b[1] - src[1]))
            assert err <= 1e-06
            replay.append(dict(site_id=s['site_id'], quantity=name, replay_max_error_mm=err))
    ctr = validation.controls(module, rows)
    suf = validation.sufficiency()
    checks = validation.contract_tests(module, rows[0]['query'] | {'site_id': rows[0]['site_id']})
    dump(output / 'CONTROLS.json', ctr)
    dump(output / 'SUFFICIENCY.json', suf)
    dump(output / 'CONTRACT_TESTS.json', checks)
    dump(output / 'SOURCE_REPLAY.json', replay)
    predictions = dict(schema='implant-safety-frozen-v1', claim_type='capability', scope='Digital envelope and conditional guide scenario only; no anatomical or injury prediction', source_manifest_sha256=sha(ROOT / 'SOURCE_MANIFEST.json'), prereg_sha256=sha(ROOT / 'PREREG_R1.json'), rows=table)
    frozen = ROOT / 'FROZEN_PREDICTIONS_R2.json'
    if frozen.exists():
        assert sha(frozen) == frozen.with_suffix('.json.sha256').read_text().strip()
        old = load(frozen)
        assert {k: v for (k, v) in old.items() if k != 'frozen_utc'} == predictions
    else:
        predictions['frozen_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        dump(frozen, predictions)
        frozen.with_suffix('.json.sha256').write_text(sha(frozen) + '\n')
    demo = next((r for r in rows if r['site_id'] == 'ToothFairy2P_083/FDI34' and r['query']['guide_type'] == 'fully_guided' and (r['query']['implant_system'] == 'ProofLane_EV_Guided_2017')))
    dump(ROOT / 'examples/query.json', dict(site_id=demo['site_id'], **demo['query']))
    dump(output / 'EXAMPLE_RESULT.json', demo)
    summary = []
    for system in pr['selection']['systems']:
        ss = [r for r in rows if r['query']['implant_system'] == system and r['query']['guide_type'] == 'fully_guided']
        n = sum((r['reference_rule']['class_vs_2mm'] == 'AT_OR_ABOVE_REFERENCE' for r in ss))
        tool = sum((r['reference_rule']['class_vs_2mm'] == 'AT_OR_ABOVE_REFERENCE' and r['digital_geometry']['tool_class_vs_2mm'] == 'BELOW_REFERENCE' for r in ss))
        rev = sum((r['reference_rule']['class_vs_2mm'] == 'AT_OR_ABOVE_REFERENCE' and r['digital_geometry']['revision_union']['interval_mm'][1] < 2 for r in ss))
        summary.append(dict(system=system, sites=len(ss), nominal_ge2mm=n, changed_by_revision=rev, changed_after_tool_envelope=tool, resolution='POPULATION', scope='22 virtual sites, four image groups; not a clinical cohort'))
    subprocess.run([sys.executable, '-s', str(ROOT / 'plot_demo.py'), str(output)], check=True)
    results = dict(schema='dental-result-v1', lane='X98-implant-safety-module', claim_type='capability', outcome='SOURCE_RESOLVED_QUERY_AVAILABLE; PHYSICAL_SAFETY_UNKNOWN', review_state='PENDING_INDEPENDENT_REVIEW', sites=len(module.sites), image_groups=len({s['image_group'] for s in module.sites.values()}), queries=len(rows), summary=summary, controls=dict(geometry_queries=len(ctr['geometry']), guide_controls=len(ctr['guide']), geometry_max_error_mm=ctr['geometry_max_error_mm'], source_replay_max_error_mm=max((x['replay_max_error_mm'] for x in replay)), actual_distance_injection_rejections=len(ctr['geometry']), contract_checks=len(checks), outcome='PASS'), sufficiency=suf, external_referent=dict(kind='published_dataset', locator='https://doi.org/10.1016/j.media.2026.104095; SOURCE_MANIFEST.json; data/sites.json external_source_bindings', compared_quantity='Minimum generic cylinder-to-published canal voxel occupancy separation at the frozen CBCT pose', refutes_us=False), external_protocol_referent=dict(kind='external_review', locator='evidence/UPSTREAM_X96_REVIEW.json C01/C04/C05 and primary manual locators in data/protocols.json', compared_quantity='Extra depth and datum arithmetic; source moments; numeric verifier behavior', refutes_us=True), rejected=dict(all_X8_sites=2581, same_pose_dense_bone_sites=22, excluded_from_demo=2559, excluded_fraction=2559 / 2581, reason='No combined dense revision + same-pose X75 bone source; no data imputed', bone_profile_rows=95, X8_pose_rows=80, other_pose_rows_not_joined=15, protocol_sources_attempted_upstream=5, protocol_sources_retained=3, protocol_sources_rejected=2, restricted_T3_excluded_from_main_table='datum example supports actual12.6mm only; evaluated separately in contract tests', physical_certificates=0, physical_UNKNOWN_queries=len(rows)), resolution=dict(raw_geometry='PER_POINT', site_query='PER_TOOTH', guide_budget='PHENOMENOLOGICAL', counts='POPULATION', time_scale='SIMULTANEOUS'), cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, threads=1, gpu=False, preparation_reasoning_seconds='UNKNOWN; not zero', fit=0, human_queries=0, physical_fallback_queries=len(rows), source_snapshot_bytes=sum((p.stat().st_size for p in (ROOT / 'data').iterdir())), inherited_geometry_cost='X73 total known1618.73s; X96 coldR1 153.722s; X75 allround numeric361s; these are shared prior costs, not independent new measurements'), limitations=['Digital masks are releases, not independently delineated true nerve wall', 'Generic 4mm cylinders on virtual dentate axes; no commercial tool/implant shape validation', 'Guide moments and implant-to-drill transport are explicit constitutive closures', 'No rigorous IEEE float enclosure or affine sensitivity claim', 'No clinical recommendation, injury risk or patient-level coverage'], artifacts=[dict(path=str(p.relative_to(ROOT)), sha256=sha(p), bytes=p.stat().st_size) for p in sorted(output.iterdir()) if p.is_file() and p.name in ['QUERIES.jsonl', 'TABLE.csv', 'CONTROLS.json', 'SUFFICIENCY.json', 'CONTRACT_TESTS.json', 'SOURCE_REPLAY.json', 'EXAMPLE_RESULT.json', 'FIGURE.png', 'FIGURE.pdf']])
    dump(output / 'RESULTS_R2_BASE.json', results)
    if output == ROOT / 'raw':
        dump(ROOT / 'results.json', results)
    state('R2_BASE_COMPLETE', results['outcome'], 'Complete registered-source, portable release and final-verifier mutation tests')
    print(json.dumps(dict(outcome=results['outcome'], sites=results['sites'], queries=len(rows), summary=summary, controls=results['controls'])))
    return results
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output-dir', type=Path, default=ROOT / 'raw')
    a = p.parse_args()
    run(a.output_dir.resolve())
