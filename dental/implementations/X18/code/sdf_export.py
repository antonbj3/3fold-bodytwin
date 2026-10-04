"""Euclidean mesh SDF samples for an exported roof shell; interpolation UNKNOWN."""
import time
from scipy.optimize import minimize
from geometry import *
from distance import closest

def run():
    pfile = H / 'raw/PREDICTIONS_R2.json'
    if not pfile.exists():
        return
    t = time.perf_counter()
    preds = json.loads(pfile.read_text())
    p = next((p for p in preds if p['status'] == 'DESIGNED' and p['fdi'] % 10 == 6))
    z = np.load(p['file'])
    (vv, ff) = shell(z['xy'], z['z_informed'], z['faces'], 1.2)
    tri = vv[ff]
    step = 0.5
    lo = vv.min(0) - 0.5
    hi = vv.max(0) + 0.5
    axes = [np.arange(lo[k], hi[k] + step / 2, step) for k in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1)
    points = grid.reshape(-1, 3)
    distance = np.full(len(points), np.inf)
    owners = np.full(len(points), -1, int)
    for first in range(0, len(points), 32):
        P = np.repeat(points[first:first + 32], len(tri), axis=0)
        T = np.tile(tri, (min(32, len(points) - first), 1, 1))
        (dd, _) = closest(P, T)
        dd = dd.reshape(-1, len(tri))
        ix = np.argmin(dd, axis=1)
        distance[first:first + len(ix)] = dd[np.arange(len(ix)), ix]
        owners[first:first + len(ix)] = ix
    (roof, _) = query_height(vv[ff[:len(z['faces'])]], points[:, :2], False)
    inside = np.isfinite(roof) & (points[:, 2] <= roof) & (points[:, 2] >= roof - 1.2)
    signed = distance * np.where(inside, -1, 1)
    path = D / 'demo_roof_sdf.npz'
    np.savez_compressed(path, sdf_mm=signed.reshape(grid.shape[:3]).astype(np.float32), origin_mm=lo, spacing_mm=np.full(3, step), shape=grid.shape[:3], nearest_triangle=owners.reshape(grid.shape[:3]).astype(np.int32))
    k = len(points) // 2
    T = tri[owners[k]]
    P = points[k]
    a = T[0]
    E = (T[1:] - a).T
    opt = minimize(lambda w: np.sum((a + E @ w - P) ** 2), np.array([1 / 3, 1 / 3]), jac=lambda w: 2 * E.T @ (a + E @ w - P), constraints=[dict(type='ineq', fun=lambda w: w, jac=lambda w: np.eye(2)), dict(type='ineq', fun=lambda w: 1 - w.sum(), jac=lambda w: -np.ones(2))], method='SLSQP', options=dict(ftol=1e-13, maxiter=500))
    err = abs(np.sqrt(opt.fun) - distance[k])
    fault = abs(np.sqrt(opt.fun) - (distance[k] + 1)) > 1e-07
    dump(H / 'raw/SDF_EXPORT.json', dict(case=p['case'], fdi=p['fdi'], path=str(path), sha256=sha(path), bytes=path.stat().st_size, spacing_mm=step, shape=grid.shape[:3], wall_seconds=time.perf_counter() - t, external_referent=dict(kind='closed_form', locator='https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf', compared_quantity='Euclidean point-to-triangle minimum distance at exported mesh-SDF grid points', refutes_us=True), independent_triangle_oracle_error_mm=err, injected1mm_distance_rejected=fault, status='PASS' if err < 1e-07 and fault else 'FAIL', scope='Euclidean distances evaluated at grid samples; sign from closed vertical-offset roof domain. Trilinear SDF interpolation/topology error unmeasured. Shape modification used exact roof fields and triangle constraints, not this coarse SDF. Not a whole preparation/crown field.'))
    assert err < 1e-07 and fault
if __name__ == '__main__':
    run()
