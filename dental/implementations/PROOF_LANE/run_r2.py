import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[k] = '4'
import copy, json, time
import numpy as np
import trimesh
from gencad_bench.io import ROOT, dump, sha, digest, state, now
from gencad_bench.geometry import cap, planes_from_cap
from gencad_bench.generators.rounded import generate, scale_profile
from gencad_bench.generators.baselines import parametric
from gencad_bench.checks.milling import check as mill
from gencad_bench.evaluate import evaluate

def volume(profile):
    (v, f) = cap(profile)
    t = v[f]
    return float(np.einsum('ij,ij->i', t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() / 6)

def gap_diagnostic(inner, prep):
    (v, _) = cap(inner)
    (pv, pf) = cap(prep)
    m = trimesh.Trimesh(pv, pf, process=False)
    d = trimesh.proximity.closest_point(m, v)[1]
    return dict(max_vertex_gap_mm=float(d.max()), mean_vertex_gap_mm=float(d.mean()), scope='sampled intaglio-to-preparation distance; maximum is a lower bound on global maximum', upper_film_requirement='NOT SPECIFIED; this is a limitation, not an approval')

def main():
    tic = time.perf_counter()
    if sha(ROOT / 'PREREG_R2.json') != (ROOT / 'PREREG_R2.sha256').read_text().strip():
        raise ValueError('R2 prereg drift')
    taskset = json.loads((ROOT / 'data/TASKSET_FROZEN.json').read_text())['payload']
    if sha(ROOT / 'data/tasks.json') != taskset['tasks_sha256']:
        raise ValueError('Task drift')
    tasks = json.loads((ROOT / 'data/tasks.json').read_text())
    private = taskset['references']
    cache = {}
    controls = {}
    rows = []
    trace = []
    state('R2_RUNNING', 'R1 failures and UNKNOWN preserved', 'Compare rounded cavity against uniformly expanded cavity')
    for t in tasks:
        gid = t['geometry_id']
        r = t['requirements']['mill_radius_mm']
        tol = t['requirements']['mill_allowance_mm']
        if gid not in controls:
            for scale in np.round(np.arange(1.0, 2.01, 0.1), 1):
                p = scale_profile(t['public']['intaglio_template'], float(scale))
                (v, _) = cap(p)
                st = time.perf_counter()
                out = mill(planes_from_cap(p), v.tolist(), r, tol)
                trace.append(dict(geometry_id=gid, scale=float(scale), milling=out, wall_s=time.perf_counter() - st))
                if out['status'] == 'PASS':
                    controls[gid] = (p, float(scale))
                    break
            if gid not in controls:
                controls[gid] = (None, None)
        ref = private[gid]
        if sha(ref['path']) != ref['sha256']:
            raise ValueError('Reference drift')
        with np.load(ref['path']) as f:
            reference = {k: f[k] for k in ['vertices', 'faces']}
        for method in ['rounded_access_IFU', 'uniform_scale_control']:
            start = time.perf_counter()
            if method == 'rounded_access_IFU':
                d = generate(t)
            else:
                (p, scale) = controls[gid]
                if p is None:
                    rows.append(dict(task_id=t['task_id'], geometry_id=gid, method=method, status='UNKNOWN_NO_SCALE_CERTIFICATE'))
                    continue
                tt = copy.deepcopy(t)
                tt['public']['intaglio_template'] = p
                d = parametric(tt)
                d['generator'] = method
            res = evaluate(t, d, reference, cache)
            rows.append(dict(task_id=t['task_id'], geometry_id=gid, material=t['material']['product'], method=method, intaglio_volume_mm3=volume(d['intaglio']), baseline_intaglio_volume_mm3=volume(t['public']['intaglio_template']), added_intaglio_volume_mm3=volume(d['intaglio']) - volume(t['public']['intaglio_template']), gap_diagnostic=gap_diagnostic(d['intaglio'], t['preparation']['profile']), scale=controls[gid][1] if method == 'uniform_scale_control' else None, wall_s=time.perf_counter() - start, design=d, **res))
        dump(ROOT / 'rounds/R2/partial.json', dict(rows=rows, control_trace=trace))
        if len(rows) % 6 == 0:
            print('R2 designs', len(rows), flush=True)
    old = json.loads((ROOT / 'rounds/R1/results.json').read_text())
    base = {r['task_id']: r for r in old['rows'] if r['generator'] == 'parametric_IFU'}
    a = [r for r in rows if r['method'] == 'rounded_access_IFU']
    b = [r for r in rows if r['method'] == 'uniform_scale_control']
    gained = sum((base[r['task_id']]['checks']['milling']['status'] == 'FAIL' and r['checks']['milling']['status'] == 'PASS' for r in a))
    regressed = sum((base[r['task_id']]['checks']['milling']['status'] == 'PASS' and r['checks']['milling']['status'] != 'PASS' for r in a))
    failures = sum((any((r['checks'][k]['status'] != 'PASS' for k in ['cement', 'wall'])) for r in a))
    out = dict(round='R2', run_utc=now(), prereg_sha256=sha(ROOT / 'PREREG_R2.json'), rows=rows, control_trace=trace, gates=dict(new_milling_passes=gained, previous_pass_regressions=regressed, cement_wall_nonpasses=failures, construction_pass=gained >= 1 and regressed == 0 and (failures == 0)), complete_level1_passes=sum((r['level1']['status'] == 'PASS' for r in a)), wall_s=time.perf_counter() - tic, outcome='Digital geometry repair; clinical fit and full task eligibility remain UNKNOWN', clinical_fit='Minimum film permits overlarge cement space; no upper-film/retention claim', external_referent=json.loads((ROOT / 'PREREG_R2.json').read_text())['external_referent'])
    dump(ROOT / 'rounds/R2/results.json', out)
    state('R2_COMPLETE', str(out['gates']), 'Finalize reviewable package and bind pending graph feedback')
    print(json.dumps(out['gates']), flush=True)
if __name__ == '__main__':
    main()
