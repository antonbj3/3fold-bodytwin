from common import *
import resource, traceback, argparse

def main():
    import integrity, remesh_r7, physics, verify, report
    args = argparse.ArgumentParser()
    args.add_argument('--rebuild', action='store_true')
    q = args.parse_args()
    started = time.perf_counter()
    cpu = time.process_time()
    check = integrity.lock()
    state('FINAL_DEMO_RUNNING', 'Input and frozen protocol hashes pass', 'Geometry, conditional physics, source/fault checks, leaderboard and figure')
    cache = P / 'raw/R7_GEOMETRY.json'
    fresh = {}
    if q.rebuild or not cache.exists():
        rows = remesh_r7.run()
        mode = 'ALL_R7_REGENERATED'
    else:
        rows = read(cache)['rows']
        mode = 'VERIFIED_CACHE_PLUS_FRESH_CASE'
        for r in rows:
            if r['status'] == 'ADAPTED':
                if sha(r['mesh_path']) != r['mesh_sha256'] or sha(r['field_path']) != r['field_sha256']:
                    raise ValueError('Cached geometry drift ' + r['uid'])
        first = next((r for r in rows if r.get('geometry_pass')))
        rec = next((r for r in cohort() if r['uid'] == first['uid']))
        new = remesh_r7.construct(rec)
        keys = ['mesh_sha256', 'field_sha256', 'gates', 'margin_error_mm', 'stl_sha256']
        fresh = dict(uid=first['uid'], compared_keys=keys, pass_gate=all((new[k] == first[k] for k in keys)))
        if not fresh['pass_gate']:
            dump(P / 'raw/FRESH_REPLAY_FAILURE.json', dict(before=first, after=new))
            raise ValueError('Fresh representative geometry does not reproduce')
    dump(P / 'raw/FRESH_GEOMETRY_REPLAY.json', dict(mode=mode, **fresh))
    state('GEOMETRY_DECIDED', sum((r.get('geometry_pass', False) for r in rows)), 'Evaluate every available crown through inherited physical chains')
    rows = physics.run(rows)
    tests = verify.run(rows)
    state('NUMERICAL_CONTROLS_DECIDED', tests['all_pass'], 'Write quality vectors, conditional ranking, exports and external-reference report')
    sources = integrity.verify()
    art = integrity.manifest()
    account = dict(mode=mode, wall_seconds_before_reporting=time.perf_counter() - started, cpu_seconds_before_reporting=time.process_time() - cpu, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, threads=4, gpu=False, prior_full_geometry_seconds={f'R{k}': read(P / f'raw/R{k}_GEOMETRY.json')['seconds'] for k in range(1, 8)}, source_integrity=sources)
    result = report.run(rows, tests, account, art)
    account.update(wall_seconds_total=time.perf_counter() - started, cpu_seconds_total=time.process_time() - cpu, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    result['resources'] = account
    dump(P / 'results.json', result)
    dump(P / 'raw/RUN_ACCOUNTING.json', account)
    integrity.verify()
    state('DEMO_COMPLETE', 'Numerical checks pass; calibrated physical force/film UNKNOWN', 'Independent review; matched specimen metrology and fracture measurement')
    print(json.dumps(dict(complete=result['complete_panel'], physical=result['physical'], ranking=result['ranking'], faults=tests['count'], resources=account), indent=2))
if __name__ == '__main__':
    try:
        main()
    except Exception:
        (P / 'raw').mkdir(exist_ok=True)
        (P / 'raw/LAST_FAILURE.txt').write_text(traceback.format_exc())
        raise
