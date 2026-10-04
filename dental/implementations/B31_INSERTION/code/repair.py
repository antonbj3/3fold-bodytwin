import pathlib, json, datetime, hashlib, os, time
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[k] = '4'
import numpy as np, trimesh
from sweep import *
R = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def freeze():
    j = json.load(open(R / 'rounds/R1/COHORT.json'))
    row = next((r for r in j['rows'] if r['full_path_status'] == 'COLLISION'))
    p = R / 'PREREG_R2.json'
    p.write_text(json.dumps(dict(round='R2', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), claim_type='capability', capability='Determine whether any inside-only repair can remove a source-bound path collision at frozen exterior/margin', obstacle='Collision witness lies on a repair-invariant triangle', changed_operation='Invariant-surface impossibility certificate plus explicit0.2mm intaglio-only apex expansion and path rerun', key=row['key'], prediction='Same exact witness survives; no feasible inside-only repair for fixed straight path', repair={'apex_delta_mm': [0, 0, 0.2], 'protected_margin_distance_mm': 0.3, 'fixed_roles': [0, 2]}, gate={'margin_max_mm': 0.01, 'exact_invariant_triangle_error': 0, 'same_witness_validation': True, 'export_roundtrip_margin_max_mm': 0.01}, strongest_equally_informed_control='Direct barycentric sweep recheck on full changed mesh; independent SciPy layout', falsifier='Changing exterior witness vertex by1mm must invalidate saved witness; moving an exterior/margin vertex must fail protected-repair contract', full_cost={'fit': 0, 'preparation': 'source mesh loading', 'discovery': 'no search, invariant proof', 'validation': 'full path+independent LP', 'query': 'two neighbour paths', 'fallback': 'R3axis search', 'physical': 'UNKNOWN'}, resolution='PER_POINT', timescale='SIMULTANEOUS'), indent=2) + '\n')
    p.with_suffix('.sha256').write_text(sha(p) + '\n')
    (R / 'FROZEN_PREDICTIONS_R2.json').write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(p), candidate='inside-only intaglio apex+0.2mm', predicted_path='COLLISION same exact invariant witness; minimal feasible inside-only repair does not exist in declared class'), indent=2) + '\n')
    print('R2FROZEN', row['key'])

def run(outdir):
    tick = time.perf_counter()
    pre = json.load(open(R / 'PREREG_R2.json'))
    co = json.load(open(R / 'FROZEN_COHORT.json'))
    rec = next((x for x in co['rows'] if x['key'] == pre['key']))
    old = np.load(rec['mesh_path'])
    p = np.load(rec['public_path'])
    v = old['vertices'].copy()
    f = old['faces']
    roles = old['roles']
    fixed = np.unique(f[np.isin(roles, [0, 2])])
    editable = np.setdiff1d(np.unique(f[roles == 1]), fixed)
    assert len(editable) == 1
    apex = int(editable[0])
    margin = p['margin_curve']
    d = np.linalg.norm(margin - v[apex], axis=1).min()
    assert d > pre['repair']['protected_margin_distance_mm']
    v[apex] += np.array(pre['repair']['apex_delta_mm'])
    assert np.array_equal(v[fixed], old['vertices'][fixed])
    outdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(outdir / 'rejected_inside_candidate.npz', vertices=v, faces=f, roles=roles)
    mesh = trimesh.Trimesh(v, f, process=False)
    mesh.export(outdir / 'rejected_inside_candidate.stl')
    mesh.export(outdir / 'rejected_inside_candidate.3mf')
    paths = {side: check_surface(v[f], Obstacle(p[side]), np.array([0, 0, 40.0]), face_ids=np.arange(len(f)), roles=roles, budget_s=180) for side in ['mesial', 'distal']}
    oldrow = next((x for x in json.load(open(R / 'rounds/R1/COHORT.json'))['rows'] if x['key'] == pre['key']))
    invariants = []
    for (side, q) in oldrow['neighbours'].items():
        w = q['path'].get('witness')
        if w:
            A = v[f[w['source_face_id']]]
            B = p[side][w['obstacle_face_id']]
            err = float(np.max(np.abs(A - np.array(w['source_triangle_mm']))))
            valid = validate_witness(A, B, np.array(w['translation_mm']), w)
            assert err == 0 and valid
            corrupt = A.copy()
            corrupt[0, 0] += 1.0
            invariants.append(dict(side=side, triangle_identity_error=err, same_exact_witness_survives=valid, corrupted_triangle_rejected=not validate_witness(corrupt, B, np.array(w['translation_mm']), w)))
    assert all((x['corrupted_triangle_rejected'] for x in invariants))
    from scipy.spatial import cKDTree
    reread = trimesh.load(outdir / 'rejected_inside_candidate.stl', process=False)
    marginerr = float(cKDTree(reread.vertices).query(old['vertices'][fixed])[0].max())
    assert marginerr <= 0.01
    three = trimesh.load(outdir / 'rejected_inside_candidate.3mf', process=False)
    three = three.to_geometry() if isinstance(three, trimesh.Scene) else three
    err3 = float(cKDTree(three.vertices).query(old['vertices'][fixed])[0].max())
    assert err3 <= 0.01
    r = dict(claim_type='capability', key=pre['key'], status='NO_FEASIBLE_INSIDE_ONLY_REPAIR_FOR_FIXED_PATH', inside_apex_change_mm=0.2, protected_margin_band_mm=0.3, native_margin_error_mm=0.0, native_exterior_identity_error_mm=0.0, STL_protected_vertex_roundtrip_max_mm=marginerr, three_mf_roundtrip_max_mm=err3, invariant_witnesses=invariants, paths=paths, minimal_feasible_inside_change='does not exist for unchanged witness triangle, regardless of amount; this is a negative certificate, not global crown-design impossibility', export_status='REJECTED_REPAIR_RESEARCH_ONLY', exports={x.name: sha(x) for x in outdir.glob('rejected_inside_candidate.*')}, seconds=time.perf_counter() - tick, unknown='wall validity/selfintersection, physical preparation/force and fabrication qualification unchanged UNKNOWN/previous failed gates')
    (outdir / 'REPAIR.json').write_text(json.dumps(r, indent=2) + '\n')
    print('R2inside0.2mm exported; protected margin0error; same exact collisions survive')
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--freeze', action='store_true')
    ap.add_argument('--out', default=str(R / 'rounds/R2'))
    a = ap.parse_args()
    freeze() if a.freeze else run(pathlib.Path(a.out))
