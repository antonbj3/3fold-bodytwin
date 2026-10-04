from dental_release.paths import expand as _release_expand
from common_local import *
import time, resource, subprocess
start = time.perf_counter()
check_frozen()
xf = pathlib.Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk/R5_FULL_GEOMETRY.jsonl'))
paths = [xf, X73 / 'raw/LINEAGE_SITES.json', X75 / 'raw/FROZEN_POSES.json', X75 / 'raw/SOURCE_MANIFEST_R2.json', X75 / 'SITE_INPUT_PORTS.jsonl', X75 / 'raw/PER_SITE_R2.jsonl', X87 / 'inputs/guide_profiles.json', X87 / 'science.py']
parents = [_release_expand('X73'), _release_expand('X87'), _release_expand('X8'), _release_expand('X75'), _release_expand('X68'), _release_expand('X67'), _release_expand('X65'), _release_expand('X56'), _release_expand('X80')]
paths += [D / l / 'results.json' for l in parents]
paths += [D / 'LANE_X68_INSERTION_TORQUE/code/query.py', D / 'LANE_X68_INSERTION_TORQUE/code/calibration.py', D / 'LANE_X68_INSERTION_TORQUE/code/enclosures.py', D / 'LANE_X67_DRILL_HEAT/code/bench_calibrate.py', D / 'LANE_X65_ISO14801_FATIGUE/code/endpoint_consumer.py', D / 'LANE_X56_MICROMOTION/code/query_motion.py']
register_inputs(paths)
x8 = lines(xf)
dense = load(X73 / 'raw/LINEAGE_SITES.json')
poses = load(X75 / 'raw/FROZEN_POSES.json')
bone = lines(X75 / 'SITE_INPUT_PORTS.jsonl')
bm = {r['case']: r for r in load(X75 / 'raw/SOURCE_MANIFEST_R2.json')}
key = lambda s: (s['case'], int(s['fdi']))
xb = {key(s): s for s in x8}
db = {key(s): s for s in dense}
bp = {key(s): s for s in poses if s['kind'].startswith('X8')}
bb = {(s['case'], int(s['site'].rsplit('_', 1)[1])): s for s in bone if s['source_pose_type'].startswith('X8')}
assert len(xb) == len(x8) == 2581 and len(db) == len(dense)
geo = module('x87_science', X87 / 'science.py')
profiles = load(X87 / 'inputs/guide_profiles.json')['profiles']
tq = module('x68_query', D / 'LANE_X68_INSERTION_TORQUE/code/query.py')
thermal = module('x67_bench', D / 'LANE_X67_DRILL_HEAT/code/bench_calibrate.py')
fatigue = module('x65_endpoint', D / 'LANE_X65_ISO14801_FATIGUE/code/endpoint_consumer.py')
ports = {}
records = []
controls = []
plan = []
source_replay_errors = []
for (k, s) in xb.items():
    d = db.get(k)
    b = bb.get(k)
    if d:
        source_replay_errors.append(abs(d['distances']['tf2']['lower_mm'] - s['lower_mm']))
        assert d['pose'] == {q: s[q] for q in d['pose']}
        assert next((q['member_sha256'] for q in d['source_bindings'] if q['dataset'] == 'tf2')) == s['source_member_sha256']
    if b:
        p = bp[k]
        assert p['source_label_member_sha256'] == s['source_member_sha256'] == bm[k[0]]['label_member_sha256']
        for (name1, name2) in [('entry', 'entry_zyx_mm'), ('axis', 'axis_zyx')]:
            assert np.array_equal(p[name1], s[name2])
        assert p['length_mm'] == s['length_mm'] and p['radius_mm'] == s['radius_mm']
    image_group = d['image_group'] if d else s['case']
    record = dict(case=s['case'], fdi=s['fdi'], image_group=image_group, pose={q: s[q] for q in ['entry_zyx_mm', 'axis_zyx', 'length_mm', 'radius_mm']}, source_label_sha256=s['source_member_sha256'], bone_image_member_sha256=bm[k[0]]['image_member_sha256'] if b else None, nominal_gap_mm=[s['lower_mm'], s['upper_mm']], dense_release_union_gap_mm=[d['union']['lower_mm'], d['union']['upper_mm']] if d else None, segment=d['segment'] if d else 'UNKNOWN_NO_DENSE_PAIR', gray=b['gray_input'] if b else None, regional_geometry=b['source_geometry'] if b else None, resolution='PER_TOOTH', physical_recommended_length_mm=None, physical_recommended_diameter_mm=None, connection=None, material=None, drill_protocol=None, insertion_torque_Ncm=None, temperature_C=None, local_micromotion_um=None, fatigue_survival=None, nerve_injury_probability=None, evidence_status='DIGITAL_GEOMETRY_MEASURED; PHYSICAL_PLAN_UNKNOWN')
    if b:
        tr = tq.query({'same_regime_confirmed': False})
        fr = fatigue.endpoint_query({'case': s['case'], 'fdi': s['fdi']})
        try:
            thermal.check_metadata({'observer_model': 'uncalibrated_cbct_gray'}, {'observer_model': 'uncalibrated_cbct_gray'})
            raise AssertionError('Gray accepted as temperature')
        except ValueError as exc:
            hr = {'status': 'REJECTED', 'reason': str(exc)}
        record['operator_calls'] = {'X68': tr, 'X67': hr, 'X65': fr}
    plan.append(record)
    if not d:
        continue
    for p in profiles:
        for target in [0.9, 0.95]:
            g = geo.guide(p, 1 - target)
            controls.append(geo.controls(p, 1 - target, g))
            margin = g['guide_budget_mm']
            nominal = s['lower_mm'] >= 2
            union = d['union']
            passed = union['lower_mm'] >= margin
            failed = union['upper_mm'] < margin
            records.append(dict(case=s['case'], fdi=s['fdi'], image_group=image_group, segment=d['segment'], guide=p['id'], target=target, nominal_gap_mm=s['lower_mm'], union_gap_lower_mm=union['lower_mm'], union_gap_upper_mm=union['upper_mm'], conditional_budget_mm=margin, remaining_gap_lower_mm=union['lower_mm'] - margin, fixed_2mm_accept=nominal, conditional_accept=passed, conditional_reject=failed, conditional_unresolved=not (passed or failed), fixed_accept_not_certified=nominal and (not passed), fixed_reject_conditional_accept=not nominal and passed, release_flip=nominal and union['upper_mm'] < 2, physical_insufficient='UNKNOWN', physical_excess='UNKNOWN', resolution='PER_TOOTH', time_scale='SIMULTANEOUS', evidence='PHENOMENOLOGICAL exact-population-moment assumption; observed release union only'))
nested = sum((any((key(x) == key(p) for x in x8)) for p in poses if p['kind'].startswith('X8')))
assert nested == len(bb)
planned = 2.0
toward = np.array([-0.5, -1.5])
away = -toward
summaries = [[planned, float(np.mean(abs(q))), float(np.var(abs(q)))] for q in [toward, away]]
assert np.array_equal(*summaries)
witness = dict(summary=summaries, identity_error=0.0, summary_bitwise_identical=np.array_equal(*summaries), actual_clearance_mm=[(planned + q).tolist() for q in [toward, away]], below_2mm_fractions=[float(np.mean(planned + q < 2)) for q in [toward, away]], downstream_difference=1.0, minimum_extension='Signed clearance loss for the fixed cylinder/canal query; full pose and canal geometry for new queries', resolution='PHENOMENOLOGICAL', external_referent={'kind': 'our_own_fixture', 'locator': 'code/run_r1.py', 'compared_quantity': 'Exact counterexample to scalar-summary sufficiency, not anatomy', 'refutes_us': True})
summary = []
for p in profiles:
    for target in [0.9, 0.95]:
        for seg in sorted({r['segment'] for r in records}):
            selected = [r for r in records if r['guide'] == p['id'] and r['target'] == target and (r['segment'] == seg)]
            grouped = {}
            for r in selected:
                grouped.setdefault((r['image_group'], r['fdi']), []).append(r)
            rr = [dict(image_group=k[0], **{field: any((s[field] for s in rows)) for field in ['fixed_accept_not_certified', 'fixed_reject_conditional_accept', 'release_flip']}) for (k, rows) in grouped.items()]
            summary.append(dict(guide=p['id'], target=target, segment=seg, resolution='POPULATION', **{f: boot(rr, f) for f in ['fixed_accept_not_certified', 'fixed_reject_conditional_accept', 'release_flip']}))
for p in profiles:
    for target in [0.9, 0.95]:
        if not any((r['segment'] == 'posterior_endpoint_10mm_proxy' and r['guide'] == p['id'] and (r['target'] == target) for r in summary)):
            summary.append(dict(guide=p['id'], target=target, segment='posterior_endpoint_10mm_proxy', resolution='POPULATION', fixed_accept_not_certified=boot([], 'x'), fixed_reject_conditional_accept=boot([], 'x'), release_flip=boot([], 'x')))
(P / 'raw/SITE_PLANS_R1.jsonl').write_text(''.join((json.dumps(s, ensure_ascii=False, allow_nan=False) + '\n' for s in plan)))
csvout(P / 'raw/PER_SITE_GUIDE_R1.csv', records)
dump(P / 'raw/POPULATION_R1.json', summary)
dump(P / 'raw/SUFFICIENCY_R1.json', witness)
intersection = [r for r in dense if key(r) in bb]
dump(P / 'raw/DENSE_BONE_INTERSECTION.json', [dict(case=r['case'], fdi=r['fdi'], image_group=r['image_group']) for r in intersection])
res = dict(claim_type=['information_link', 'capability'], status='DIGITAL_SITE_COMPOSITION_COMPLETE_PHYSICAL_PLAN_UNIDENTIFIED', counts=dict(x8_sites=len(x8), x8_cases=len({r['case'] for r in x8}), dense_poses=len(dense), dense_unique_image_fdi=len({(r['image_group'], r['fdi']) for r in dense}), bone_same_pose_sites=len(bb), dense_and_bone_same_pose_sites=len(intersection), dense_and_bone_cases=len({r['case'] for r in intersection}), full_physical_plans=0), dropout=dict(dense_missing=len(x8) - len(dense), bone_missing=len(x8) - len(bb), bone_other_pose_excluded=len(bone) - len(bb)), identity_control=dict(nested_join_count=nested, hash_and_pose_mismatches=0, source_distance_replay_max_error_mm=max(source_replay_errors), source_replay_gate_1e_10='PASS' if max(source_replay_errors) <= 1e-10 else 'FAIL_PRESERVED'), controls=dict(guide_controls=len(controls), max_error_mm=max((c['error_mm'] for c in controls)), injected_budget_plus1_rejected=all((c['injected_plus1mm_rejected'] for c in controls))), physical_too_little_or_too_much='UNKNOWN; a conservative bound failing cannot establish harm', external_referent=load(P / 'PREREG_R1.json')['external_referent'], cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, upstream_reuse='Inherited computations not rerun; source cost retained in parent results.json', fit=0, questions=0))
dump(P / 'rounds/R1.json', res)
state('R1_COMPLETE', res['status'], 'Freeze R2 dimension inverse on exact-pose dense/bone intersection')
(P / 'HANDOFF_R1.md').write_text('# R1: gemensamma platsidentiteter\n\n' + json.dumps(res['counts'], indent=2) + "\n\nNo full physical plan is identified. The next design changes the dimensions and remeasures the digital geometry ; the fixed pose distance must not be reused for second dimensions.\n")
print(json.dumps(res, indent=2))
