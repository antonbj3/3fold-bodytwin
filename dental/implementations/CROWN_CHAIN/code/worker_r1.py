from common_chain import *
import time
import numpy as np
import trimesh
from scipy.optimize import linprog

def main(index):
    start = time.perf_counter()
    pr = read(ROOT / 'PREREG_R1.json')
    inp = read(ROOT / 'FROZEN_INPUTS.json')
    r = inp['cohort'][index]
    ir = next((x for x in inp['inputs'] if x['key'] == r['key']))
    for (path, digest) in [(r['mesh_path'], r['mesh_sha256']), (ir['public_path'], ir['public_sha256']), (ir['private_path'], ir['private_sha256'])]:
        if sha(path) != digest:
            raise ValueError('Source data drift ' + path)
    with np.load(r['mesh_path'], allow_pickle=False) as a:
        (v, f, roles) = (a['vertices'], a['faces'], a['roles'])
    with np.load(ir['public_path'], allow_pickle=False) as a:
        p = dict(a)
    with np.load(ir['private_path'], allow_pickle=False) as a:
        target = a['target']
    (co, geom, score) = scoped_modules(BASE / 'PROOF_LANE_FULL_CROWN_R4/code', ['common', 'geometry', 'score'])
    ext = v[f[roles == 0]]
    inner = v[f[roles == 1]]
    shape = geom.metrics(ext, target)
    wall = score.wall(ext, inner)
    outer = trimesh.Trimesh(v, f[roles == 0], process=False)
    outer.remove_unreferenced_vertices()
    curve = outer.vertices[geom.loops(outer)[0]]
    margin = geom.curve_error(curve, p['margin_curve'])
    crown = trimesh.Trimesh(v, f, process=False)
    centers = ext.mean(1)
    (cp, dd, ids) = geom.closest(inner, centers)
    i = int(np.argmin(dd))
    wall_candidate = dict(exterior_face=int(np.flatnonzero(roles == 0)[i]), intaglio_face=int(np.flatnonzero(roles == 1)[ids[i]]), exterior_point_mm=centers[i], inner_point_mm=cp[i], distance_mm=float(dd[i]))
    cone = module('chain_cone', BASE / 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py')
    normals = -crown.face_normals[roles == 1]
    normal_valid = np.linalg.norm(normals, axis=1) > 0
    con = cone.classify(normals[normal_valid])
    verified = cone.recheck(con) if con['certificates'] else con['status'] == 'YES_STRICT'
    lp = linprog([0, 0, 0, -1], A_ub=np.column_stack([-normals[normal_valid], np.ones(normal_valid.sum())]), b_ub=np.zeros(normal_valid.sum()), bounds=[(-1, 1)] * 3 + [(None, None)], method='highs')
    (xco, collision, xrun) = scoped_modules(BASE / 'LANE_X53_MILLING_TOOLS/code', ['common', 'collision', 'run_r1'])
    valid = np.linalg.norm(np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]]), axis=1) / 2 >= 1e-14
    scene = collision.Scene(v, f[valid])
    newroles = roles[valid]
    pts = xrun.sample(scene, np.flatnonzero(newroles == 1), seed=5300 + index, n=pr['metrics']['tool_sample_count_per_region'])
    tools = xrun.lib('vhf', 5)
    toolrows = []
    for pt in pts:
        rr = []
        for t in tools:
            a = collision.search(scene, pt['point'], pt['normal'], t, 5, offsets=(0.0,))
            rr.append(dict(tool_id=t['id'], **a))
        toolrows.append(dict(**pt, tools=rr, has_witness=any((q['status'] == 'FOUND' for q in rr))))
    (x85co, x85geom, x85) = scoped_modules(BASE / 'LANE_X85_OCCLUSAL_ADJUST/code', ['common', 'geometry', 'round5'])
    fq = x85.full_crown_force_query(dict(mesh_sha256=r['mesh_sha256'], frame=r['case_key']))
    try:
        x85.full_crown_force_query(dict(mesh_sha256=r['mesh_sha256'], frame=r['case_key']), dict(geometry_sha256='other-specimen', frame='other-frame'))
        wrong_binding = False
    except ValueError:
        wrong_binding = True
    sys.path.insert(0, str(BASE / 'LANE_X49_DESIGN_GATE/code'))
    gate = module('chain_x49', BASE / 'LANE_X49_DESIGN_GATE/code/design_gate.py')
    wall_status = gate.classify_interval([wall['continuous_lower_mm'], wall['sampled_min_mm']], 0.5)
    rules = dict(shape=dict(status='PASS' if shape['p95_mm'] <= pr['metrics']['shape_p95_mm'][r['family']] else 'FAIL'), margin=dict(status='PASS' if margin <= 0.025 else 'FAIL'), wall=dict(status=wall_status), closed=dict(status='PASS' if crown.is_watertight and crown.is_winding_consistent and (crown.volume > 0) else 'FAIL'), insertion=dict(status='FAIL' if con['status'] == 'NO_NONZERO_DIRECTION' else 'UNKNOWN'), milling=dict(status='UNKNOWN'), force=dict(status=fq['status']), fracture=dict(status='UNKNOWN'), cement=dict(status='UNKNOWN'), lab=dict(status='UNKNOWN'))
    agg = gate.finish(dict(rules=rules), time.perf_counter(), time.process_time())
    out = dict(key=r['key'], case_key=r['case_key'], fdi=ir['fdi'], family=r['family'], mesh_path=r['mesh_path'], mesh_sha256=r['mesh_sha256'], unit='mm', frame=dict(kind='R4 local tooth frame', case=r['case_key'], origin_world=p['origin_world'], bite='UNKNOWN'), shape=shape, natural_floor=r['natural_floor'], shape_family_gate=rules['shape']['status'], within_own_natural_floor=shape['p95_mm'] <= r['natural_floor']['p95_mm'], margin_max_mm=margin, wall=wall, wall_candidate=wall_candidate, cone=con, cone_verified=verified, lp_control=dict(success=lp.success, strict_margin=lp.x[3] if lp.success else None), tool_points=toolrows, tool_dropped_degenerate_facets=int((~valid).sum()), tool_scope='Legacy declared scaled envelope, finite5axis poses, sampled points; no actual OEM assembly/blank/CAM certificate', force=fq, wrong_force_binding_rejected=wrong_binding, design_gate=agg, replay_delta_mm=dict(shape=abs(shape['p95_mm'] - r['in_situ']['p95_mm']), wall=abs(wall['sampled_min_mm'] - r['wall']['sampled_min_mm'])), cost=dict(wall_s=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), float_enclosure='MISSING; exact Farkas certificates apply to stated finite normal model only')
    if max(out['replay_delta_mm'].values()) > 1e-07:
        raise ValueError('Frozen replay tolerance failed')
    dump(ROOT / 'raw/R1' / f"{r['key']}.json", out)
    print(r['key'], agg['verdict'], con['status'], sum((q['has_witness'] for q in toolrows)), 'tool points', flush=True)
if __name__ == '__main__':
    main(int(sys.argv[1]))
