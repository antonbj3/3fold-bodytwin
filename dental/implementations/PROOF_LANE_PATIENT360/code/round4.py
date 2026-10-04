from dental_release.paths import expand as _release_expand
import sys, time, copy
from collections import Counter
import numpy as np
from shared import *
run = Path(sys.argv[1])
out = run / 'R4'
out.mkdir(exist_ok=True)
start = time.perf_counter()
x15 = module('measure_r1', RESULTS / 'LANE_X15_BONE_DEHISCENCE/measure_r1.py')
x15r2 = module('x15_exact_control', RESULTS / 'LANE_X15_BONE_DEHISCENCE/measure_r2.py')
pr = load(R / 'PREREG_R4.json')
tol = pr['gates']['exact_boundary_original_X15_R2_tolerance_mm']

def compare(a, b):
    if len(a) != len(b) or any((u[2] != v[2] for (u, v) in zip(a, b))):
        return dict(accept=False, reason='LABEL_SEQUENCE_MISMATCH', maximum_error_mm=None)
    err = max((abs(u[j] - v[j]) for (u, v) in zip(a, b) for j in [0, 1]), default=0)
    return dict(accept=err <= tol, reason='ENDPOINT_COMPARISON', maximum_error_mm=err)
summary = []
for case in [_release_expand('@DENTAL_CASE_A@'), _release_expand('@DENTAL_CASE_B@')]:
    q = np.load(DATA / run.name / case / 'labels.npz')
    lab = q['labels']
    sp = q['spacing']
    src = load(run / case / 'X15_BONE.json')
    frame = load(run / case / 'FRAME.json')
    direction = np.array(frame['direction'])
    offset = np.array(frame['offset_xyz_mm'])
    rows = []
    sample_errors = []
    statuses = Counter()
    errors = []
    mutations = []
    sections = Counter((r['status'] for r in src['rows']))
    freeze(out / (case + '_FROZEN_PREDICTIONS.json'), dict(prereg_R4_sha256=sha(R / 'PREREG_R4.json'), patient_id=case, source_sections=artifact(run / case / 'X15_BONE.json'), label_file=artifact(DATA / run.name / case / 'labels.npz'), frame=artifact(run / case / 'FRAME.json'), prediction='All exact grid-face label intervals agree with independent cube-slab control within 1e-8mm; sampled gate remains FAIL; physical anatomy UNKNOWN', physical_measurement='NOT_RUN'))
    for r in src['rows']:
        if r['status'] != 'MEASURED_SECTION':
            continue
        origin = np.array(r['origin_xyz_mm'])
        b = np.array(r['buccal_direction_xyz'])
        ex = x15.exact_ray(lab, origin, b, sp)
        slab = x15r2.slab_control(lab, origin, b, sp)
        v = compare(ex, slab)
        stored = compare(ex, r['exact_runs'])
        errors.append(v)
        errors.append(stored)
        wrong = copy.deepcopy(slab)
        wrong[0][1] += 1.0
        wronglabel = copy.deepcopy(slab)
        wronglabel[0][2] += 100
        mutations.append(dict(endpoint_plus_1mm_rejected=not compare(ex, wrong)['accept'], label_change_rejected=not compare(ex, wronglabel)['accept']))
        item = dict(patient_id=case, fdi=r['tooth'], jaw=r['jaw'], height_from_apex_mm=r['height_from_apex_mm'], resolution='PER_SURFACE_REGION', time_scale='SIMULTANEOUS', origin_world_xyz_mm=direction @ origin + offset, direction_world_xyz=direction @ b, frame_sha256=sha(run / case / 'FRAME.json'), exact_runs=ex, slab_runs=slab, full_interval_control=v, stored_exact_replay=stored)
        for (side, name) in [(1, 'buccal'), (-1, 'lingual')]:
            m = x15.side_measure(ex, r['tooth'], r['jaw'], side)
            s = x15.side_measure(slab, r['tooth'], r['jaw'], side)
            item[name] = dict(primary_exact=m, slab_control=s, original_sampled=r.get(name, {}))
            statuses[m['status']] += 1
            if r.get(name, {}).get('control_error_mm') is not None:
                sample_errors.append(r[name]['control_error_mm'])
        rows.append(item)
    passed = all((e['accept'] for e in errors))
    faults = all((m['endpoint_plus_1mm_rejected'] and m['label_change_rejected'] for m in mutations))
    rec = dict(patient_id=case, claim_type='capability', resolution='PER_SURFACE_REGION', sections_attempted=len(src['rows']), section_status_counts=sections, sections_evaluated=len(rows), section_rejected=len(src['rows']) - len(rows), side_status_counts=statuses, sides_evaluated=2 * len(rows), sides_measured=statuses['MEASURED_ANNOTATION_SURROGATE'], sides_rejected=2 * len(rows) - statuses['MEASURED_ANNOTATION_SURROGATE'], sampled_max_error_mm=max(sample_errors), sampled_tolerance_mm=0.051, sampled_gate='PASS' if max(sample_errors) <= 0.051 else 'FAIL', exact_max_error_mm=max((e['maximum_error_mm'] or 0 for e in errors)), exact_gate='PASS' if passed else 'FAIL', fault_injections_rejected=faults, primary_port='EXACT_ANNOTATION_GEOMETRY' if passed and faults else 'ABSTAIN', physical='UNKNOWN_ANATOMY_AND_MOVEMENT', external_referent=load(run / case / 'CBCT_RESULTS.json')['external_referent'])
    dump(out / (case + '_EXACT_BONE.json'), dict(summary=rec, rows=rows, faults=mutations))
    summary.append(rec)
    print(case, 'sampled', rec['sampled_gate'], rec['sampled_max_error_mm'], 'exact', rec['exact_gate'], rec['exact_max_error_mm'], 'sides', rec['sides_measured'], flush=True)
    del lab, q
dump(out / 'RESULTS_R4.json', dict(claim_type='capability', patients=summary, prereg=artifact(R / 'PREREG_R4.json'), runtime_s=time.perf_counter() - start, external_referent=pr['external_referent']))
