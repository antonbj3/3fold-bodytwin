from exclusion import *

def validate_one(w, prep, point):
    p = tuple((Q(x) for x in w['preparation_point']))
    s = tuple((Q(x) for x in w['source_point']))
    t = tuple((rational_point(v) for v in prep))
    (on, _) = point_triangle(p, t)
    ds = sum(((x - y) ** 2 for (x, y) in zip(s, p)))
    return on == 0 and s == rational_point(point) and (ds == Q(w['distance_squared_mm2'])) and (ds < Q(w['threshold_squared_mm2']))

def run():
    st = time.perf_counter()
    out = read(R / 'RESULTS_D.json')
    source = {r['key']: (r, rr) for (r, rr) in inputs()}
    n = 0
    fail = []
    wrong_rejected = False
    for row in out['rows']:
        (rec, rr) = source[row['key']]
        m = npz(rec['mesh_path'])
        prep = m['vertices'][m['faces'][m['roles'] == 1]]
        points = scorer.sample(npz(rr['private_path'])['target'], 8192)
        w = read(row['exact_witness_path'])['witnesses']
        assert sha(row['exact_witness_path']) == row['exact_witness_sha256']
        assert len(w) == row['exact_exclusion_count']
        for witness in w:
            tri = prep[witness['preparation_triangle_index']]
            point = points[witness['native_probe_index']]
            okay = validate_one(witness, tri, point)
            n += 1
            if not okay:
                fail.append((row['key'], witness['native_probe_index']))
            if not wrong_rejected:
                bad = dict(witness, distance_squared_mm2=str(Q(witness['distance_squared_mm2']) + 1))
                wrong_rejected = not validate_one(bad, tri, point)
        (_, dist, _) = fast_nearest(prep, points)
        bad_count = int(np.sum(np.maximum(0, 0.1 - 0.06 - dist) > row['shape_threshold_mm']))
        assert bad_count == 0
    assert n and (not fail) and wrong_rejected
    from controls import tile
    reference = tile(0, 0.625)
    a = tile(0, 0.75)
    b = tile(0, 0.5)
    inner = tile(0, 0.125)
    ma = scorer.metrics(a, reference, 128)
    mb = scorer.metrics(b, reference, 128)
    err = max((abs(ma[k] - mb[k]) for k in ['p95_mm', 'rms_mm', 'sampled_max_mm']))
    wa = distance(inner, sample(a, 128)).min()
    wb = distance(inner, sample(b, 128)).min()
    assert err == 0 and ma['p95_mm'] == 0.125 and (wa == 0.625) and (wb == 0.375)
    report = dict(verified_rational_witnesses=n, failures=fail, all_pass=True, injected_wrong_squared_distance_rejected=wrong_rejected, injected_wall0p1_removes_all_exclusions=True, shape_sufficiency=dict(identity_error_mm=err, shape_summaries=[ma, mb], wall_min_mm=[float(wa), float(wb)], difference_mm=float(wa - wb), minimum_extension='For this pair, the minimum wall distance. A shape quantile does not determine full-crown feasibility.', resolution='PER_POINT to PER_TOOTH', source='our_own_fixture; exact Euclidean parallel-plane check'), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    save(R / 'raw/D_VERIFICATION.json', report)
    print(n, 'exact witnesses verified; injected distance rejected; shape identity0, wall difference', wa - wb)
    return report
if __name__ == '__main__':
    run()
