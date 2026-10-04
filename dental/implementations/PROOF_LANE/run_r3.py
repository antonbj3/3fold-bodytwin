import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[k] = '4'
import json, time
from gencad_bench.io import ROOT, dump, sha, state, now
from gencad_bench.geometry import cap, planes_from_cap
from gencad_bench.checks import cement, milling, cement_max
from gencad_bench.generators.rounded import round_intaglio
from gencad_bench.generators.rolling_ball import profile

def main():
    tic = time.perf_counter()
    if sha(ROOT / 'PREREG_R3.json') != (ROOT / 'PREREG_R3.sha256').read_text().strip():
        raise ValueError('R3 drift')
    if sha(ROOT / 'PREREG_R3_V2.json') != (ROOT / 'PREREG_R3_V2.sha256').read_text().strip():
        raise ValueError('R3 V2 drift')
    tasks = json.loads((ROOT / 'data/tasks.json').read_text())
    rows = []
    cache = {}
    for task in tasks[::3]:
        (P, _) = cap(task['preparation']['profile'])
        for (method, p) in [('original', task['public']['intaglio_template']), ('R2_rounded', round_intaglio(task['public']['intaglio_template'])), ('classical_rolling_ball', profile(task['preparation']['profile']))]:
            started = time.perf_counter()
            mx = cement_max.check(task['preparation']['profile'], p, 0.12)
            try:
                pl = planes_from_cap(p)
                (v, _) = cap(p)
                convex = True
                low = cement.check(P.tolist(), pl, 0.03)
                mill = milling.check(pl, v.tolist(), 0.5, 0.05)
            except ValueError as e:
                convex = False
                low = mill = dict(status='UNKNOWN', reason=str(e))
            rows.append(dict(geometry_id=task['geometry_id'], method=method, maximum_film=mx, minimum_film=low, milling=mill, exact_convex_cap=convex, profile=p, wall_s=time.perf_counter() - started))
        dump(ROOT / 'rounds/R3/partial.json', rows)
        print('R3', task['geometry_id'], flush=True)
    original = [r for r in rows if r['method'] == 'original']
    changed = [r for r in rows if r['method'] == 'R2_rounded']
    gates = dict(original_maximum_film_PASS=sum((r['maximum_film']['status'] == 'PASS' for r in original)), changed_maximum_film_FAIL=sum((r['maximum_film']['status'] == 'FAIL' for r in changed)))
    out = dict(round='R3_V2', prereg_v2_sha256=sha(ROOT / 'PREREG_R3_V2.json'), run_utc=now(), rows=rows, gates=gates, representation_gate_pass=all((r['exact_convex_cap'] for r in rows if r['method'] == 'classical_rolling_ball')), gate_pass=gates['original_maximum_film_PASS'] >= 1 and gates['changed_maximum_film_FAIL'] >= 1, n_anatomies=8, task_replicas_per_anatomy=3, new_scenario_maximum_mm=0.12, original_benchmark_verdicts_changed=False, clinical_threshold='UNKNOWN; 0.12 mm is an explicitly modeled adversarial test requirement', wall_s=time.perf_counter() - tic, external_referent=json.loads((ROOT / 'PREREG_R3.json').read_text())['external_referent'])
    dump(ROOT / 'rounds/R3/results.json', out)
    state('R3_COMPLETE', str(gates), 'Finalize package; next need paired clinical preparation and independently measured fit')
    print(json.dumps(gates), flush=True)
if __name__ == '__main__':
    main()
