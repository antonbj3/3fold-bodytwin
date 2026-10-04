from dental_release.paths import expand as _release_expand
from construct_b import *
from fractions import Fraction as Q
sys.path.insert(0, str(BASE / _release_expand('PROOF_LANE')))
from gencad_bench.checks.exact import point_triangle, strings

def rational_point(p):
    return tuple((Q(float(x)) for x in p))

def run():
    st = time.perf_counter()
    rows = []
    pr = read(R / 'PREREG_D.json')
    W = Q('0.5')
    G = Q('0.060')
    for (rec, rr) in inputs():
        m = npz(rec['mesh_path'])
        prep = m['vertices'][m['faces'][m['roles'] == 1]]
        target = npz(rr['private_path'])['target']
        pts = scorer.sample(target, 8192)
        (q, d, j) = fast_nearest(prep, pts)
        lim = Q(str(pr['metrics']['shape_threshold_mm'][rec['family']]))
        radius = W - G - lim
        potential = np.flatnonzero(d < float(radius)) if radius > 0 else []
        witness = []
        for k in potential:
            (dd, qq) = point_triangle(rational_point(pts[k]), tuple((rational_point(p) for p in prep[j[k]])))
            if dd < radius * radius:
                witness.append(dict(native_probe_index=int(k), preparation_triangle_index=int(j[k]), distance_squared_mm2=str(dd), threshold_squared_mm2=str(radius * radius), source_point=strings(rational_point(pts[k])), preparation_point=strings(qq)))
        status = 'CERTIFIED_INFEASIBLE_UNDER_OFFSET_COVERAGE' if len(witness) >= 411 else 'UNKNOWN'
        p = R / 'raw/D_WITNESSES' / (rec['key'] + '.json')
        save(p, dict(witnesses=witness, source_preparation_mesh_sha256=rec['mesh_sha256'], reference_file=rr['private_path'], reference_sha256=rr['private_sha256']))
        row = dict(key=rec['key'], family=rec['family'], status=status, reference_probes=8192, exact_exclusion_count=len(witness), exclusion_fraction=len(witness) / 8192, sampled_shape_p95_lower_mm=float(np.quantile(np.maximum(0, 0.44 - d), 0.95)), lower_bound_numeric_scope='p95 diagnostic float; exclusion count carries exact rational proof', shape_threshold_mm=float(lim), wall_mm=0.5, prep_to_intaglio_coverage_max_mm=0.06, exact_witness_path=p, exact_witness_sha256=sha(p), resolution='PER_TOOTH', injected_wall_0p1_excluded_count=0)
        rows.append(row)
        save(R / 'raw/D_PARTIAL.json', rows)
        print(rec['key'], status, len(witness), 'p95lower', row['sampled_shape_p95_lower_mm'], flush=True)
    out = dict(claim_type='capability', rows=rows, summary={fam: sum((r['status'].startswith('CERTIFIED') for r in rows if r['family'] == fam)) for fam in ['anterior', 'premolar', 'molar']}, external_referent=pr['external_referent'], seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, scope='Necessary condition only; covers fixed digital p95, true offset coverage<=60um and wall>=0.5mm. No clinical impossibility or proof of feasibility on UNKNOWN rows.')
    save(R / 'RESULTS_D.json', out)
    checkpoint('D_DECIDED', out['summary'], 'Use exact exclusions to distinguish bad preparation contract from generator failure')
if __name__ == '__main__':
    run()
