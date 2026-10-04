"""Runnable case adapter: existing mm scene + preparation + opposite tooth -> research3MF.

No source graph or predecessor file is written. Physical force/gap/CAM remain unknown.
"""
import thread_guard
from cadlib import *
from readonly import parents, scope, force_refusal
from pipeline_r4 import canonical_prep, offset_cap
import pipeline, three_mf, trimesh
P = parents()
pipeline.offset_cap = offset_cap
with scope(BASE / 'PROOF_LANE_FULL_CROWN_R4/code', ['common', 'geometry', 'generate', 'generate_c', 'generate_d']):
    AD = module('generate_d', BASE / 'PROOF_LANE_FULL_CROWN_R4/code/generate_d.py')

def generate_case(record, profile):
    t = P['tasks'].load_task(record['task_path'])
    if t.get('units') != 'mm':
        raise ValueError('Explicit mm public scene required')
    if record.get('units', 'mm') != 'mm':
        raise ValueError('Preparation input must declare mm')
    (prep, ser) = canonical_prep(record)
    data = np.load(record['public_geometry'])
    result = pipeline.create(prep, t, profile, P, data['donor'])
    result['serialization'] = ser
    if result['status'] != 'GENERATED':
        return result
    ext = trimesh.Trimesh(result['ov'], result['of'], process=False)
    nn = data['neighbors']
    mid = nn.mean(1)[:, 0]
    (e, info) = AD.adapt(ext, dict(mesial=nn[mid < 0], distal=nn[mid >= 0]))
    base = e.vertices.copy()
    z0 = float(result['iv'][:, 2].min())
    ob = np.flatnonzero(abs(base[:, 2] - z0) < 1e-07)
    attempts = []
    chosen = None
    for alpha in np.arange(9) * 0.25:
        v = base.copy()
        v[ob, :2] += alpha * v[ob, :2] / np.linalg.norm(v[ob, :2], axis=1)[:, None]
        w = P['cap'].wall_check(v, e.faces, result['iv'], result['inf'], profile['wall_mm'])
        attempts.append(dict(extra_flare_mm=float(alpha), wall=w))
        if w['status'] == 'PASS':
            chosen = (v, float(alpha))
            break
    if chosen is None:
        ov = base
        alpha = 0.0
        status = 'NO_WALL_FEASIBLE_IN_FROZEN_GRID'
    else:
        (ov, alpha) = chosen
        status = 'WALL_FEASIBLE_DISCRETE_RING'
    (m, roles) = pipeline.stitch(ov, e.faces, result['iv'], result['inf'], np.arange(24), ob)
    result.update(mesh=m, roles=roles, ov=ov, of=e.faces, neighbour_adaptation=info, ring_repair=dict(status=status, extra_flare_mm=alpha, attempts=attempts), force_port=force_refusal(record['key']))
    return result

def cli():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    manifest = read(args.manifest)
    if manifest.get('units') != 'mm':
        raise ValueError('Manifest units must be mm')
    profile = next((x for x in read(ROOT / 'MATERIAL_PROFILES.json')['profiles'] if x['id'] == manifest['material']))
    r = generate_case(manifest['record'], profile)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    if r['status'] == 'GENERATED':
        m = r['mesh']
        three_mf.write(out / 'crown.3mf', m.vertices, m.faces, r['roles'], dict(material=profile['id'], physical_release='UNKNOWN', source=manifest['record']['key']))
        np.savez_compressed(out / 'crown.npz', vertices=m.vertices, faces=m.faces, roles=r['roles'])
        m.export(out / 'crown.stl')
    report = {k: v for (k, v) in r.items() if k not in ['mesh', 'roles', 'ov', 'of', 'iv', 'inf', 'outer_roof', 'prior', 'undercut_dot', 'undercut_faces']}
    report['physical_release'] = 'UNKNOWN'
    dump(out / 'CAD_REPORT.json', report)
    print(r['status'], out)
if __name__ == '__main__':
    cli()
