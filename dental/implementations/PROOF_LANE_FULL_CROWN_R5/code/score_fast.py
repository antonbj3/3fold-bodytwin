from construct_b import *
import geometry as geometry_r4
from score_a import exact_failure
slow_closest = geometry_r4.closest
geometry_r4.closest = fast_nearest
scorer.closest = fast_nearest

def evaluate(tag):
    st = time.perf_counter()
    fr = read(R / f'FROZEN_PREDICTIONS_{tag}.json')
    idx = {r['key']: (rec, r) for (rec, r) in inputs()}
    pr = read(BASE / 'PROOF_LANE_FULL_CROWN_R4/PREREG_B.json')['metrics']
    rows = []
    cache = {}
    for rec in fr['rows']:
        t0 = time.perf_counter()
        if 'mesh_path' not in rec:
            rows.append(dict(rec, technical_joint=False, complete_with_offset=False))
            save(R / f'raw/{tag}_SCORE.json', rows)
            continue
        rr = idx[rec['key']][1]
        p = npz(rr['public_path'])
        target = npz(rr['private_path'])['target']
        m = npz(rec['mesh_path'])
        (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
        ext = v[f[roles == 0]]
        inner = v[f[roles == 1]]
        probes = sample(ext, 128)
        (_, dc, _) = fast_nearest(inner, probes)
        (_, ds, _) = slow_closest(inner, probes)
        parity = float(np.max(abs(dc - ds)))
        shape = scorer.metrics(ext, target)
        wall = scorer.wall(ext, inner)
        mm = trimesh.Trimesh(v, f, process=False)
        om = trimesh.Trimesh(v, f[roles == 0], process=False)
        curve = om.vertices[loops(om)[0]]
        margin = curve_error(curve, p['margin_curve'])
        (q, d, j) = fast_nearest(inner, ext.mean(1))
        i = int(d.argmin())
        ex = exact_failure(ext[i], inner[j[i]])
        ex.update(outer_triangle=ext[i], inner_triangle=inner[j[i]])
        gp = sample(inner, 4096)
        (q, gap, _) = fast_nearest(m['prep_triangles'], gp)
        (pv, pf) = compact(m['vertices'], m['faces'][roles == 1])
        original = npz(idx[rec['key']][0]['mesh_path'])
        (pv, pf) = compact(original['vertices'], original['faces'][original['roles'] == 1])
        prepmesh = trimesh.Trimesh(pv, pf, process=False)
        prepcur = pv[loops(prepmesh)[0]]
        desired = boundary_gap(q, prepcur)
        gres = abs(gap - desired)
        (xy, ff) = scorer.grid(target)
        ceil = scorer.height(p['antagonist'], xy, True)
        pg = ceil - scorer.height(ext, xy)
        rg = ceil - scorer.height(target, xy)
        cm = scorer.contact.compare(xy, ff, pg, rg)
        cg = scorer.contact.gates(cm, pr['v6'])
        v4 = scorer.contact.Q.contact_map(xy, ff, pg, rg)
        gates = dict(shape=shape['p95_mm'] <= pr['shape_threshold_mm'][rec['family']], margin=margin <= 0.025, closed=bool(mm.is_watertight and mm.is_winding_consistent and (mm.volume > 0)), wall=wall['continuous_lower_mm'] >= 0.5, offset=bool(gres.max() <= 0.01 and rec['offset']['max_plane_offset_residual_mm'] <= 0.01))
        row = dict(rec, shape=shape, wall=wall, exact_wall_witness=ex, gap_unsigned_sampled_mm=dict(min=float(gap.min()), median=float(np.median(gap)), max=float(gap.max()), max_field_residual_mm=float(gres.max()), scope='4096 area-centroid probes; not a continuous certificate or measured seated film'), margin_mm=margin, gates=gates, technical_joint=all((gates[k] for k in ['shape', 'margin', 'closed', 'wall'])), complete_with_offset=all(gates.values()), v4_contact=v4, v6_contact=cm, v6_gates=cg, v6_pose_status='UNVERIFIED_RAW_POSE_NOT_FUNCTIONAL_FACIT', score_seconds=time.perf_counter() - t0, resolution='PER_TOOTH', status='SCORED', distance_backend_parity_max_mm=parity)
        rows.append(row)
        save(R / f'raw/{tag}_SCORE.json', rows)
        print(tag, rec['key'], rec['method'], gates, 'p95', shape['p95_mm'], 'wall', wall['sampled_min_mm'], 'gaperr', gres.max(), flush=True)
        checkpoint(tag + '_SCORING', len(rows), 'Retain every result and execute X49 post-export control')
    out = dict(claim_type='capability', rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, summary={method: {fam: dict(requested=6, scored=sum((r.get('status') == 'SCORED' for r in rows if r['method'] == method and r['family'] == fam)), technical_joint=sum((r.get('technical_joint', False) for r in rows if r['method'] == method and r['family'] == fam)), complete_with_offset=sum((r.get('complete_with_offset', False) for r in rows if r['method'] == method and r['family'] == fam)), gates={k: sum((r.get('gates', {}).get(k, False) for r in rows if r['method'] == method and r['family'] == fam)) for k in ['shape', 'margin', 'closed', 'wall', 'offset']}) for fam in ['anterior', 'premolar', 'molar']} for method in sorted({r['method'] for r in rows})})
    save(R / f'RESULTS_{tag}.json', out)
    checkpoint(tag + '_DECIDED', out['summary'], 'Preserve outcomes and export only qualifying candidates')
if __name__ == '__main__':
    evaluate(sys.argv[1])
