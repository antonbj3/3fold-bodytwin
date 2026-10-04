"""Supplemental engineering positive/negative control, not a new crown score."""
import sys, time
from pathlib import Path
if __name__ == '__main__':
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
    from util import *
    from integrity import Bundle
    import trimesh
    b = Bundle()
    pr = b.json('PREREG_R6.json')
    threshold = pr['metrics']['exterior_p95_distance_max_mm']
    freeze(ROOT / 'FROZEN_NATIVE_DISTANCE_CONTROL.json', dict(claim_type='capability', purpose='Verify native-distance gate admits actual source vertices and rejects displaced vertices; retrospective engineering control, no changed crown decision', benchmark_sha256=sha(ROOT / 'BENCHMARK_LOCK.json'), script_sha256=sha(Path(__file__)), cases=pr['cases'], tolerance_mm=threshold, shift_mm=[0.0, 0.0, 2.0], samples_per_case=128))
    rows = []
    start = time.perf_counter()
    for key in pr['cases']:
        a = b.npz('payload/preparations_r6/' + key + '/fields.npz')
        mesh = trimesh.Trimesh(a['source_vertices'], a['source_faces'], process=False)
        indices = np.linspace(0, len(mesh.vertices) - 1, min(128, len(mesh.vertices)), dtype=int)
        pts = mesh.vertices[indices]
        native = trimesh.proximity.closest_point(mesh, pts)[1]
        shift = trimesh.proximity.closest_point(mesh, pts + [0, 0, 2.0])[1]
        exact = trimesh.proximity.closest_point_naive(mesh, pts[:16])[1]
        zero = float(np.quantile(native, 0.95))
        bad = float(np.quantile(shift, 0.95))
        rows.append(dict(case_key=key, native_p95_mm=zero, shifted_p95_mm=bad, positive=zero <= threshold, shift_rejected=bad > threshold, independent_distance_max_error_mm=float(max(abs(native[:16] - exact))), resolution='PER_SURFACE_REGION', point_resolution='PER_POINT'))
    b.verify()
    out = dict(claim_type='capability', all_positive_and_negative=all((r['positive'] and r['shift_rejected'] and (r['independent_distance_max_error_mm'] <= 1e-08) for r in rows)), rows=rows, seconds=time.perf_counter() - start, external_referent=pr['external_referent'], scope='Native source identity and a displaced surface query only. Does not rescue any failed R4/R6 crown.')
    dump(ROOT / 'raw/NATIVE_DISTANCE_CONTROL.json', out)
    print(out['all_positive_and_negative'])
    if not out['all_positive_and_negative']:
        raise RuntimeError('source-distance positive/negative control failed')
