from common_chain import *
import numpy as np
import trimesh, time

def data_path(p):
    prefix = str(BASE / 'PROOF_LANE_GENCAD_V4/payload')
    return Path(str(p).replace(prefix, '/tmp/chain_v4_payload', 1)) if str(p).startswith(prefix) else Path(p)

def run():
    verify_freeze(ROOT / 'PREREG_R3.json')
    pr = read(ROOT / 'PREREG_R3.json')
    t0 = time.perf_counter()
    for (p, h) in pr['files'].items():
        if sha(data_path(p)) != h:
            raise ValueError('R3 source drift ' + p)
    (co, geom, score) = scoped_modules(BASE / 'PROOF_LANE_FULL_CROWN_R4/code', ['common', 'geometry', 'score'])
    cone = module('r3_cone', BASE / 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py')
    (xco, collision, xrun) = scoped_modules(BASE / 'LANE_X53_MILLING_TOOLS/code', ['common', 'collision', 'run_r1'])
    (x85co, x85geom, x85) = scoped_modules(BASE / 'LANE_X85_OCCLUSAL_ADJUST/code', ['common', 'geometry', 'round5'])
    rows = []
    for (ix, r) in enumerate(pr['cohort']):
        st = time.perf_counter()
        with np.load(r['artifact_path'], allow_pickle=False) as a:
            (v, f, roles) = (a['vertices'], a['faces'], a['face_roles'])
        with np.load(data_path(BASE / 'PROOF_LANE_GENCAD_V4/payload/whole_private' / r['key'] / 'reference.npz'), allow_pickle=False) as a:
            source = dict(a)
        with np.load(data_path(BASE / 'PROOF_LANE_GENCAD_V4/payload/whole_inputs' / r['key'] / 'preparation.npz'), allow_pickle=False) as a:
            prep = dict(a)
        ext = v[f[roles == 0]]
        inner = v[f[roles == 1]]
        target = source['source_triangles']
        target = target[target.mean(1)[:, 2] >= float(prep['margin_z'])]
        shape = geom.metrics(ext, target)
        w = score.wall(ext, inner)
        cm = trimesh.Trimesh(v, f, process=False)
        norms = -cm.face_normals[roles == 1]
        validnorm = np.linalg.norm(norms, axis=1) > 0
        ins = cone.classify(norms[validnorm])
        ins_verified = cone.recheck(ins) if ins['certificates'] else ins['status'] == 'YES_STRICT'
        valid = cm.area_faces >= 1e-14
        scene = collision.Scene(v, f[valid])
        rr = roles[valid]
        pts = xrun.sample(scene, np.flatnonzero(rr == 1), seed=6300 + ix, n=8)
        tools = xrun.lib('vhf', 5)
        trows = []
        for pt in pts:
            tool = [dict(tool_id=t['id'], **collision.search(scene, pt['point'], pt['normal'], t, 5, offsets=(0.0,))) for t in tools]
            trows.append(dict(**pt, tools=tool, has_witness=any((t['status'] == 'FOUND' for t in tool))))
        (xy, faces) = geom.grid(target)
        ceil = score.height(source['antagonist_triangles'], xy, True)
        reference_gap = ceil - score.height(target, xy)
        gap = ceil - score.height(ext, xy)
        contact = score.contact.compare(xy, faces, gap, reference_gap)
        cg = score.contact.gates(contact, read(BASE / 'PROOF_LANE_FULL_CROWN_R4/PREREG_B.json')['metrics']['v6'])
        fq = x85.full_crown_force_query(dict(mesh_sha256=r['artifact_sha256'], frame=r['key']))
        known_lower = r['all_exterior_wall_lower_mm']
        probe_lower = known_lower - 0.025
        row = dict(key=r['key'], case_key=r['case_key'], dataset=r['dataset'], fdi=r['source_fdi'], mesh_path=r['artifact_path'], mesh_sha256=r['artifact_sha256'], frame=r['key'], unit='mm', shape=shape, shape_legacy_gate='PASS' if shape['p95_mm'] <= 0.35 else 'FAIL', contralateral_generation='NOT_THIS_AUXILIARY_ARM; inherited R2 generated exterior', own_natural_floor='UNKNOWN_NOT_MEASURED_FOR_THIS_ARM', wall_recomputed=w, wall_inherited_lower_mm=known_lower, wall_probe_lower_mm=probe_lower, probe_wall_gate='PASS' if probe_lower >= 0.5 else 'UNKNOWN', probe_scope='All exterior25um bounded-motion over same intaglio; analytic perturbation with prior nominal bound; full rigorous arithmetic/source enclosure MISSING', cone=ins, cone_verified=ins_verified, tool_points=trows, tool_dropped_degenerate_facets=int((~valid).sum()), contact=contact, contact_gates=cg, bite_pose='SOURCE_REGISTERED; independent loaded pose error UNKNOWN', force=fq, field_force_classification='UNKNOWN_NO_COMPLETE_SAME_MESH_GAP_GRAPH_AND_SUPPORT_OBSERVATION', field_force_interval_N=None, physical_full_chain='UNKNOWN', design_chain='FAIL' if shape['p95_mm'] > 0.35 or ins['status'] == 'NO_NONZERO_DIRECTION' else 'UNKNOWN', cost=dict(wall_s=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        rows.append(row)
        dump(ROOT / 'raw/R3_ROWS.json', rows)
        print(r['key'], row['design_chain'], row['probe_wall_gate'], ins['status'], flush=True)
    result = dict(claim_type='capability', round='R3', rows=rows, attrition=pr['attrition'], seconds=time.perf_counter() - t0, external_referent=pr['external_referent'])
    dump(ROOT / 'raw/R3_ALL.json', result)
    return result
if __name__ == '__main__':
    run()
