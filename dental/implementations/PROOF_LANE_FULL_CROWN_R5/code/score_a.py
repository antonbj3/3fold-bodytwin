from dental_release.paths import expand as _release_expand
from construct_a import *
from fractions import Fraction as Q

def exact_failure(t1, t2):
    sys.path.insert(0, str(BASE / _release_expand('PROOF_LANE')))
    from gencad_bench.checks.exact import triangle_distance, vec, strings
    (d, a, b) = triangle_distance(tuple((tuple((Q(float(x)) for x in p)) for p in t1)), tuple((tuple((Q(float(x)) for x in p)) for p in t2)))
    return dict(fails=d < Q(1, 4), distance_squared_mm2=str(d), outer_point=strings(a), inner_point=strings(b), rigorous_stored_coordinate_arithmetic=True, scope='Exact rational triangle pair; not scan uncertainty')

def evaluate(tag='A'):
    st = time.perf_counter()
    fr = read(R / f'FROZEN_PREDICTIONS_{tag}.json')
    idx = {r['key']: (rec, r) for (rec, r) in inputs()}
    rows = []
    pr = read(BASE / 'PROOF_LANE_FULL_CROWN_R4/PREREG_B.json')['metrics']
    for rec in fr['rows']:
        t0 = time.perf_counter()
        rr = idx[rec['key']][1]
        p = npz(rr['public_path'])
        target = npz(rr['private_path'])['target']
        m = npz(rec['mesh_path'])
        (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
        ext = v[f[roles == 0]]
        inner = v[f[roles == 1]]
        shape = scorer.metrics(ext, target)
        wall = scorer.wall(ext, inner)
        mm = trimesh.Trimesh(v, f, process=False)
        om = trimesh.Trimesh(v, f[roles == 0], process=False)
        curve = om.vertices[loops(om)[0]]
        margin = curve_error(curve, p['margin_curve'])
        (q, d, j) = closest(inner, ext.mean(1))
        i = int(d.argmin())
        ex = exact_failure(ext[i], inner[j[i]])
        gap = distance(m['prep_triangles'], sample(inner, 4096))
        (xy, ff) = scorer.grid(target)
        ceil = scorer.height(p['antagonist'], xy, True)
        pg = ceil - scorer.height(ext, xy)
        rg = ceil - scorer.height(target, xy)
        cm = scorer.contact.compare(xy, ff, pg, rg)
        cg = scorer.contact.gates(cm, pr['v6'])
        gates = dict(shape=shape['p95_mm'] <= pr['shape_threshold_mm'][rec['family']], margin=margin <= 0.025, closed=bool(mm.is_watertight and mm.is_winding_consistent and (mm.volume > 0)), wall=wall['continuous_lower_mm'] >= 0.5, offset=rec['offset']['max_plane_offset_residual_mm'] <= 0.01)
        row = dict(rec, shape=shape, wall=wall, exact_wall_witness=ex, gap_unsigned_sampled_mm=dict(min=float(gap.min()), median=float(np.median(gap)), max=float(gap.max())), margin_mm=margin, gates=gates, technical_joint=all((gates[k] for k in ['shape', 'margin', 'closed', 'wall'])), complete_with_offset=all(gates.values()), v6_contact=cm, v6_gates=cg, v6_pose_status='UNVERIFIED_RAW_POSE_NOT_FUNCTIONAL_FACIT', score_seconds=time.perf_counter() - t0, resolution='PER_TOOTH', status='SCORED')
        rows.append(row)
        save(R / f'raw/{tag}_SCORE.json', rows)
        print(rec['key'], rec['method'], gates, 'p95', shape['p95_mm'], 'wall', wall['sampled_min_mm'], flush=True)
        checkpoint(tag + '_SCORING', len(rows), 'Complete frozen controls and candidates, preserve all failed gates')
    out = dict(claim_type='capability', rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, summary={method: {fam: {k: sum((r['gates'][k] for r in rows if r['method'] == method and r['family'] == fam)) for k in ['shape', 'margin', 'closed', 'wall', 'offset']} for fam in ['anterior', 'premolar', 'molar']} for method in sorted({r['method'] for r in rows})})
    save(R / f'RESULTS_{tag}.json', out)
    checkpoint(tag + '_DECIDED', out['summary'], 'Change offset representation; keep same preparation and exterior targets')
if __name__ == '__main__':
    evaluate(sys.argv[1] if len(sys.argv) > 1 else 'A')
