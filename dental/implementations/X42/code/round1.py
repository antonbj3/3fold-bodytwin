from common import *
import resource, collections
METHODS = {'population': population, 'parametric': parametric, 'x1b_morphology': crown_loop, 'x18_antagonist': antagonist, 'constraint_optimizer': standard_lp, 'x42_contact': contact_r1}

def run():
    t0 = time.perf_counter()
    rows = []
    cost = collections.Counter()
    fails = []
    fcounts = collections.Counter()
    for (ii, key) in enumerate(split()['dev_round1']):
        refs = dev_reference(key)
        ss = {}
        for t0r in tasks(key):
            t = t0r
            if t['status'] == 'READY':
                if t['geometry_file'] not in ss:
                    ss[t['geometry_file']] = scene(t)
                t = dict(ss[t['geometry_file']])
                t.update(t0r)
                t['preparation_z'] = np.full(len(t['xy']), t['preparation_height_mm'])
                t['obstacle_b'] = ss[t['geometry_file']]['obstacle_b']
                fcounts[feasibility(t)['status']] += 1
            else:
                fcounts['UNKNOWN_SITE'] += 1
            for (name, fn) in METHODS.items():
                row = dict(case_key=key, family=t['family'], level=t['level'], method=name)
                if t['status'] != 'READY':
                    row.update(verdict='UNKNOWN_SITE')
                else:
                    tic = time.perf_counter()
                    try:
                        d = fn(t)
                    except Exception as e:
                        cost[name] += time.perf_counter() - tic
                        fails.append(dict(task_id=t['task_id'], method=name, exception=repr(e)))
                        row.update(verdict='NOT_RUN', reason=repr(e))
                        rows.append(row)
                        continue
                    cost[name] += time.perf_counter() - tic
                    ch = all_checks(t, d)
                    row.update(verdict=ch['validity'], checks=ch['checks'])
                    if d['status'] == 'DESIGN':
                        row.update(metric(t, np.asarray(d['outer_vertices'])[:, 2], refs[t['family']]))
                rows.append(row)
        dump(ROOT / 'raw/R1_ROWS.json', rows)
        state('R1_DEV_RUNNING', {'completed_cases': ii + 1}, 'Next preregistered dev case; hidden not opened')
        print('R1 dev', ii + 1, 'seconds', round(time.perf_counter() - t0), flush=True)
    out = dict(round='R1', claim_type='algorithm', summary=summarize(rows), feasibility_counts=dict(fcounts), cost_seconds=dict(cost), wall_seconds=time.perf_counter() - t0, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, failed_executions=fails, hidden_queries=0, prereg_sha256=sha(ROOT / 'PREREG_R1.json'), external_referent=json.loads((ROOT / 'PREREG_R1.json').read_text())['external_referent'])
    dump(ROOT / 'rounds/R1.json', out)
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    run()
