from score import *

def run():
    st = time.perf_counter()
    inputs = {r['key']: r for r in read(ROOT / 'FROZEN_INPUTS_B.json')['records']}
    rows = []
    for rec in read(ROOT / 'FROZEN_PREDICTIONS_C.json')['rows']:
        if rec['method'] != 'exact_margin':
            continue
        r = dict(key=rec['key'], family=rec['family'])
        target = npz(inputs[rec['key']]['private_path'])['target']
        m = npz(rec['mesh_path'])
        inner = m['vertices'][m['faces'][m['roles'] == 1]]
        ww = wall(target, inner)
        r.update(original_exterior_wall=ww, original_sampled_failure=ww['sampled_min_mm'] < 0.5, original_geometric_wall_pass=ww['continuous_lower_mm'] >= 0.5, necessary_same_point_max_displacement_mm=max(0.0, 0.5 - ww['sampled_min_mm']), bound_scope='A necessary maximum pointwise change for preserving all original surface correspondences; not a p95 lower bound')
        rows.append(r)
        dump(ROOT / 'raw/PREPARATION_ORACLE_CHECKPOINT.json', rows)
        print(r['key'], r['original_sampled_failure'], flush=True)
    result = dict(rows=rows, requested=18, source_sampled_failure=sum((r['original_sampled_failure'] for r in rows)), source_geometric_pass=sum((r['original_geometric_wall_pass'] for r in rows)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, external_referent=REFERENT)
    dump(ROOT / 'raw/PREPARATION_ORACLE.json', result)
    return result
if __name__ == '__main__':
    run()
