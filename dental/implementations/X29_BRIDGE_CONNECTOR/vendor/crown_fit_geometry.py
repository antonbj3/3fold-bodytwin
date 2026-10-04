'crown_fit_geometry.py — K2 geometry : dental frame, precise redistance (A143 -safe), synthetic preparation, crown design.\n\nPart of the K2 chain (DENT-MFG-PROCESS-MODEL / DENT-MFG-ASBUILT-DEVIATION):\n  designed crown -> milled green body (crown_milling_access ) -> sintering ( crown_sinter_shrinkage )\n  -> manufactured crown -> cement gap / fit ( crown_cement_fit ).\n\nUnits: length in mm internally; columns/deviations are reported in µm of calling cells.\nRepresentation: signed distance field (SDF , mm , negative inside) on a regular grid in the frame of the tooth\n(z = tooth length axis, positive towards occlusal ).\n\nA143 (COMPUTE_CELL_INVENTORY_cs_engines.md , C37 /§ queries : offset of a Union/Boolean combination or of a\nvariable offset, is NOT a distance field ; to offset it again gives distance error. Followed by every boolean\nOperation, variable offset, erosion and dilation of `redistance()` , which recalculates the exact distance to\nthe surface of the zero level (zero throughputs on grating edges -> kd-tree -> point/plane-distance ). Custom code, no import of\nField engine offset ( field_csg_fused_v1 ). Verified against analytical union by two spheres in\nresults/K2_crown_fit/run_a143_check.py.\n\nData source: OpenMandible base model (git e 1 f 8 cef, DOI 10.1016/j.dental.2021.01.009 ), Tooth_L6 = FDI 36\n( enamel + dentin + pulp ). Data license UNKNOWN -> is processed local , derived geometry is not spread.\n'
from dental_release.paths import expand as _release_expand
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree
import os
WORKERS = int(os.environ.get('THREADS', '2'))
OM_STL = _release_expand('@DENTAL_DATA_ROOT@/biomech/OpenMandible/00 OpenMandible Base model/materials/stl')
OM_SHA256 = {'Tooth_L6_Enamel.stl': 'af12b89f568ce9743d187bd6762f8d43e755566facd63db53dcd7a9a69642ff8', 'Tooth_L6_Dentin.stl': '7d0b326f167f9021cbfb8e1110f91e517b02c9fc6d60bf4f938a7c6b77dfd028', 'Tooth_L6_Pulp.stl': 'eb982db7f1478fc55dd7ef495a0597ca64e704c5ac90daf191f2eba18fd90890'}

class Grid:
    """Regular grid : point (i,j,k) is in origin + h*(i,j,k). Field has shape shape (x,y,z)."""

    def __init__(self, lo, hi, h):
        self.h = float(h)
        self.origin = np.asarray(lo, float)
        self.shape = tuple((int(np.ceil((hi[d] - lo[d]) / h - 1e-06)) + 1 for d in range(len(lo))))

    def axes(self):
        return [self.origin[d] + self.h * np.arange(self.shape[d]) for d in range(len(self.shape))]

    def mesh(self):
        return np.meshgrid(*self.axes(), indexing='ij')

    def to_index(self, pts):
        return ((np.asarray(pts) - self.origin) / self.h).T

    def sample(self, field, pts, order=1, cval=None):
        """Trilinear interpolation of fields in arbitrary points (N ,dim)."""
        if cval is None:
            cval = float(np.max(field))
        return ndimage.map_coordinates(field, self.to_index(pts), order=order, mode='constant', cval=cval)

def zero_crossings(phi, grid):
    """Zero throughputs on lattice edges (linear interpolation) + normals (gradient, interpolated)."""
    pts = []
    ax = grid.axes()
    dim = phi.ndim
    for d in range(dim):
        a = [slice(None)] * dim
        b = [slice(None)] * dim
        a[d] = slice(0, -1)
        b[d] = slice(1, None)
        (pa, pb) = (phi[tuple(a)], phi[tuple(b)])
        m = (pa > 0) != (pb > 0)
        idx = np.nonzero(m)
        t = pa[idx] / (pa[idx] - pb[idx])
        coords = [ax[e][idx[e]] for e in range(dim)]
        coords[d] = coords[d] + t * grid.h
        pts.append(np.stack(coords, 1))
    P = np.concatenate(pts, 0)
    g = np.gradient(phi.astype(np.float32), grid.h)
    N = np.stack([grid.sample(gi, P, order=1, cval=0.0) for gi in g], 1)
    N /= np.linalg.norm(N, axis=1, keepdims=True) + 1e-12
    return (P, N)

def distance_to_surface(X, P, N, tree=None, workers=WORKERS):
    """Unsigned distance from X to the surface that is sampled by points P with normals N .
    Near the surface (within 2 point distance) point-plan-distance , otherwise euclidic point distance is used."""
    if tree is None:
        tree = cKDTree(P)
    (d, i) = tree.query(X, k=1, workers=workers)
    dp = np.abs(np.einsum('ij,ij->i', X - P[i], N[i]))
    return np.where(d < 2.5 * _spacing(P, tree), dp, d)
_SP_CACHE = {}

def _spacing(P, tree):
    key = (len(P), float(P[0, 0]))
    if key not in _SP_CACHE:
        s = P[::max(1, len(P) // 2000)]
        (dd, _) = tree.query(s, k=2)
        _SP_CACHE[key] = float(np.median(dd[:, 1]))
    return _SP_CACHE[key]

def point_triangle_distance(X, A, B, C):
    """Precisely distance point-triangle (Ericson, Real-Time Collision Detection § 5.1.5 ), vectorised."""
    (ab, ac, ap) = (B - A, C - A, X - A)
    d1 = np.einsum('ij,ij->i', ab, ap)
    d2 = np.einsum('ij,ij->i', ac, ap)
    bp = X - B
    d3 = np.einsum('ij,ij->i', ab, bp)
    d4 = np.einsum('ij,ij->i', ac, bp)
    cp = X - C
    d5 = np.einsum('ij,ij->i', ab, cp)
    d6 = np.einsum('ij,ij->i', ac, cp)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    den = va + vb + vc
    den = np.where(np.abs(den) < 1e-30, 1e-30, den)
    v = vb / den
    w = vc / den
    Q = A + ab * v[:, None] + ac * w[:, None]
    t_ab = np.clip(d1 / np.where(d1 - d3 == 0, 1e-30, d1 - d3), 0, 1)
    t_ac = np.clip(d2 / np.where(d2 - d6 == 0, 1e-30, d2 - d6), 0, 1)
    t_bc = np.clip((d4 - d3) / np.where(d4 - d3 + (d5 - d6) == 0, 1e-30, d4 - d3 + (d5 - d6)), 0, 1)
    m = (d1 <= 0) & (d2 <= 0)
    Q[m] = A[m]
    m2 = (d3 >= 0) & (d4 <= d3) & ~m
    Q[m2] = B[m2]
    m |= m2
    m2 = (d6 >= 0) & (d5 <= d6) & ~m
    Q[m2] = C[m2]
    m |= m2
    m2 = (vc <= 0) & (d1 >= 0) & (d3 <= 0) & ~m
    Q[m2] = (A + ab * t_ab[:, None])[m2]
    m |= m2
    m2 = (vb <= 0) & (d2 >= 0) & (d6 <= 0) & ~m
    Q[m2] = (A + ac * t_ac[:, None])[m2]
    m |= m2
    m2 = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0) & ~m
    Q[m2] = (B + (C - B) * t_bc[:, None])[m2]
    return np.linalg.norm(X - Q, axis=1)

def surface_mesh(phi, grid):
    """"Marching cubes" triangulation in grating coordinates -> (V, F)."""
    from skimage.measure import marching_cubes
    (V, F, _, _) = marching_cubes(phi, 0.0, spacing=(grid.h,) * 3, allow_degenerate=False)
    return (V + grid.origin, F)

def mesh_distance(X, V, F, k=8, tree=None, workers=WORKERS, chunk=1000000):
    """Exactly distance from X to Triangle Meshe (V,F): min over the k triangles whose centre of gravity is closest."""
    Cc = V[F].mean(1)
    if tree is None:
        tree = cKDTree(Cc)
    out = np.empty(len(X))
    for s0 in range(0, len(X), chunk):
        x = X[s0:s0 + chunk]
        (_, ids) = tree.query(x, k=k, workers=workers)
        best = np.full(len(x), np.inf)
        for j in range(k):
            f = F[ids[:, j]]
            best = np.minimum(best, point_triangle_distance(x, V[f[:, 0]], V[f[:, 1]], V[f[:, 2]]))
        out[s0:s0 + chunk] = best
    return out

def signed_edt(phi, grid):
    """Pitch-corrected signed EDT from coating (phi< = 0 ): d = d_out - h/ 2 outside, -( d_in - h/ 2 ) inside.
    Same construction as the field engine faltkarna_v1_mesh_to_sdf d_out/d_in ; recommended before erosion of
    Union under the A143 notification (3fold-motion-engine/_private/romi_collab/lanes/DENTAL_ANSWERS_20260922.md § 2 ,
    build/A143/RESULTS.md p. 3,5,9). Errors ≤ ~h/ 2 ; provides the right topology."""
    occ = phi <= 0
    d_in = ndimage.distance_transform_edt(occ, sampling=grid.h)
    d_out = ndimage.distance_transform_edt(~occ, sampling=grid.h)
    return np.where(occ, -(d_in - 0.5 * grid.h), d_out - 0.5 * grid.h).astype(np.float32)

def redistance_lowmem(phi, grid, surface=None, workers=None, k=12, levels=(0.0,), shell=None):
    "Signed distance to the zero level of phi (or to the given triangle surface=(V,F)), character from phi.\n    Step 1 (topology, entire grid): pitch-corrected signed EDT (signed_edt, A143 recommendation).\n    Step 2 (accuracy): in thin shells |d_edt - L| < shell around each level L in `levels` (0 = surface area, and levels that\n    the next operation thresholds on, e.g. -r before erosion) is replaced by the exact point–triangle-distance to\n    the triangulation of the zero level (marching cubes; failure against smooth surface O(h^ 2 *scrow)) or to the given mesh.\n    Candidate triangles: the nearest (weight point) to the voxel's nearest seed voxel (sign change), via the EDT index.\n    Reference check: results/K2_crown_fit/run_a143_check.py."
    workers = workers or WORKERS
    if shell is None:
        shell = 3 * grid.h + 0.02
    (V, F) = surface_mesh(phi, grid) if surface is None else surface
    out = signed_edt(phi, grid)
    pos = phi > 0
    seed = np.zeros(phi.shape, bool)
    for d in range(phi.ndim):
        a = [slice(None)] * phi.ndim
        b = [slice(None)] * phi.ndim
        a[d] = slice(0, -1)
        b[d] = slice(1, None)
        diff = pos[tuple(a)] != pos[tuple(b)]
        seed[tuple(a)] |= diff
        seed[tuple(b)] |= diff
    si = np.flatnonzero(seed)
    ssub = np.unravel_index(si, phi.shape)
    Xs = np.stack([grid.origin[d] + grid.h * ssub[d] for d in range(phi.ndim)], 1)
    Cc = V[F].mean(1)
    (_, cand) = cKDTree(Cc).query(Xs, k=k, workers=workers)
    rank = np.full(phi.size, -1, np.int64)
    rank[si] = np.arange(len(si))
    del Xs, ssub
    m = np.zeros(phi.shape, bool)
    for L in levels:
        m |= np.abs(out - L) < shell
    sel = np.flatnonzero(m.ravel())
    del m
    (_, inds) = ndimage.distance_transform_edt(~seed, return_indices=True)
    del seed
    CH = 1000000
    for s0 in range(0, len(sel), CH):
        ii = sel[s0:s0 + CH]
        sub = np.unravel_index(ii, phi.shape)
        X = np.stack([grid.origin[d] + grid.h * sub[d] for d in range(phi.ndim)], 1)
        sf = np.ravel_multi_index(tuple((inds[d][sub] for d in range(phi.ndim))), phi.shape)
        ids = cand[rank[sf]]
        best = np.full(len(ii), np.inf)
        for j in range(ids.shape[1]):
            f = F[ids[:, j]]
            best = np.minimum(best, point_triangle_distance(X, V[f[:, 0]], V[f[:, 1]], V[f[:, 2]]))
        out.ravel()[ii] = np.where(pos.ravel()[ii], best, -best)
    return out

def erode(phi, grid, r, nxt=(0.0,), **kw):
    """Phil must be exactly close to the level -r (see levels at the previous redistance)."""
    return redistance_lowmem(phi + r, grid, levels=tuple(nxt), **kw)

def dilate(phi, grid, r, nxt=(0.0,), **kw):
    return redistance_lowmem(phi - r, grid, levels=tuple(nxt), **kw)

def opening(phi, grid, r, nxt=(0.0,), **kw):
    """Phi exactly close -r required. Erosion -> exactly close to +r -> dilatation."""
    return dilate(erode(phi, grid, r, nxt=(0.0, r), **kw), grid, r, nxt=nxt, **kw)

def closing(phi, grid, r, nxt=(0.0,), **kw):
    """phi exactly close to +r required."""
    return erode(dilate(phi, grid, r, nxt=(0.0, -r), **kw), grid, r, nxt=nxt, **kw)

def load_tooth(tooth='L6', check_sha=True):
    import hashlib
    import os
    import trimesh
    ms = {}
    for part in ('Enamel', 'Dentin', 'Pulp'):
        fn = f'Tooth_{tooth}_{part}.stl'
        p = os.path.join(OM_STL, fn)
        if check_sha and fn in OM_SHA256:
            h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
            assert h == OM_SHA256[fn], f'sha 256 differs for {fn}'
        ms[part] = trimesh.load(p)
    return ms

def tooth_frame(ms):
    'Frame : z = main axis (PCA of enamel + dentin ), positive to the centre of gravity of the enamel; origin at the PCA mean.'
    V = np.vstack([ms['Dentin'].vertices, ms['Enamel'].vertices])
    c = V.mean(0)
    (w, v) = np.linalg.eigh((V - c).T @ (V - c))
    ez = v[:, -1]
    if (ms['Enamel'].vertices.mean(0) - ms['Dentin'].vertices.mean(0)) @ ez < 0:
        ez = -ez
    ex = v[:, 0] - v[:, 0] @ ez * ez
    ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    R = np.stack([ex, ey, ez], 0)
    Ve = (ms['Enamel'].vertices - c) @ R.T
    c2 = c + R.T @ np.array([Ve[:, 0].mean(), Ve[:, 1].mean(), 0.0])
    return (R, c2)

def union_boundary_mesh(ms, R, c, occ, grid):
    """Triangles on the Union's edge ( enamel  ∪  dentin  ∪  pulp ) in dental frame: a triangle is retained if the point  1,5  h
    externally along its normal distance from the centre of gravity lies outside the Union's surface (internal interfaces removed)."""
    (Vs, Fs, off) = ([], [], 0)
    for part in ('Enamel', 'Dentin', 'Pulp'):
        m = ms[part]
        V = (m.vertices - c) @ R.T
        Nf = m.face_normals @ R.T
        Cc = V[m.faces].mean(1)
        probe = Cc + 1.5 * grid.h * Nf
        ijk = np.round(grid.to_index(probe)).astype(int)
        ok = np.all((ijk >= 0) & (ijk < np.array(grid.shape)[:, None]), 0)
        inside = np.zeros(len(Cc), bool)
        inside[ok] = occ[ijk[0, ok], ijk[1, ok], ijk[2, ok]]
        Vs.append(V)
        Fs.append(m.faces[~inside] + off)
        off += len(V)
    return (np.concatenate(Vs), np.concatenate(Fs))

def occupancy_from_meshes(ms, R, c, grid):
    """the coating (bool) of the Union by slicing along z and the even-odd rule per disc;"""
    from matplotlib.path import Path
    (X, Y) = np.meshgrid(*grid.axes()[:2], indexing='ij')
    XY = np.stack([X.ravel(), Y.ravel()], 1)
    occ = np.zeros(grid.shape, bool)
    zs = grid.axes()[2]
    for part in ('Enamel', 'Dentin', 'Pulp'):
        m = ms[part].copy()
        m.apply_translation(-c)
        T = np.eye(4)
        T[:3, :3] = R
        m.apply_transform(T)
        for (k, z) in enumerate(zs):
            if z < m.bounds[0, 2] or z > m.bounds[1, 2]:
                continue
            s = m.section(plane_origin=[0, 0, z], plane_normal=[0, 0, 1])
            if s is None:
                continue
            inside = np.zeros(len(XY), np.int8)
            for poly in s.discrete:
                pth = Path(poly[:, :2])
                inside ^= pth.contains_points(XY).astype(np.int8)
            occ[:, :, k] |= inside.reshape(X.shape).astype(bool)
    return occ

def tooth_sdf(ms, R, c, grid, levels=(0.0,)):
    """SDF for tooth : characters from coating (sliced mesh), exact point–triangle-distance to original meshens
    Union stripe in shell around `levels` (e.g. 0 and - occlusal reduction)."""
    occ = occupancy_from_meshes(ms, R, c, grid)
    phi0 = np.where(occ, -1.0, 1.0).astype(np.float32)
    VF = union_boundary_mesh(ms, R, c, occ, grid)
    return (redistance_lowmem(phi0, grid, surface=VF, levels=levels), VF)

def slice_sdf2d(phi3, grid, z):
    """2D - SDF of the cross section of the tooth at height z (refused in 2D )."""
    zs = grid.axes()[2]
    k = int(round((z - zs[0]) / grid.h))
    s = phi3[:, :, k].astype(np.float32)
    g2 = Grid(grid.origin[:2], grid.origin[:2] + grid.h * (np.array(grid.shape[:2]) - 1), grid.h)
    (P, N) = zero_crossings(s, g2)
    tree = cKDTree(P)
    X = np.stack([a.ravel() for a in g2.mesh()], 1)
    d = distance_to_surface(X, P, N, tree).reshape(s.shape)
    return (np.where(s > 0, d, -d).astype(np.float32), float(zs[k]))

def make_preparation(phi_T, grid, z_m, occl_red=2.0, axial_red=1.0, half_angle_deg=3.0, line_angle_radius=0.5, chamfer_radius=None, final_levels=(0.0, 0.05, 0.1)):
    """Synthetic whole chandelier repair (documented, non-clinical scan):
    - margin: plane at z_m  (≈ 1   mm  above enamel lowest point), chamfer depth =  axial_red .
    - axial walls: the outline of the cross section at z_m displaced towards axial_red , convergence half_angle_deg per wall .
    - occlusal : anatomical reduction occl_red (anatomy displaced downwards along the axis, following the cups).
    - line angles are rounded (opening, radius line_angle_radius); chamfer-hole heel (closure, radius chamfer_radius).
    Returns SDF for prepared tooth (including root under z_m) and the margin curve points ."""
    if chamfer_radius is None:
        chamfer_radius = 0.9 * axial_red
    (X, Y, Z) = grid.mesh()
    (phi2, zk) = slice_sdf2d(phi_T, grid, z_m)
    t = np.tan(np.radians(half_angle_deg))
    lat = phi2[:, :, None] + axial_red + np.maximum(Z - z_m, 0.0) * t
    n = int(round(occl_red / grid.h))
    occ = np.empty_like(phi_T)
    occ[:, :, :-n] = phi_T[:, :, n:]
    occ[:, :, -n:] = np.abs(phi_T[:, :, -1:]) + 1.0
    above = np.maximum(lat, occ)
    slab = Z - z_m
    P1 = np.minimum(above, slab).astype(np.float32)
    P1 = redistance_lowmem(P1, grid, levels=(0.0, -line_angle_radius))
    if line_angle_radius > 0:
        P1 = opening(P1, grid, line_angle_radius)
    prep = np.where(Z > z_m + 1e-09, P1, phi_T).astype(np.float32)
    prep = redistance_lowmem(prep, grid, levels=(0.0, chamfer_radius))
    if chamfer_radius > 0:
        prep = closing(prep, grid, chamfer_radius, nxt=final_levels)
    (Pz, Nz) = zero_crossings(phi2, Grid(grid.origin[:2], grid.origin[:2] + grid.h * (np.array(grid.shape[:2]) - 1), grid.h))
    margin = np.column_stack([Pz, np.full(len(Pz), zk)])
    return (prep, margin)

def spacer_field(prep, grid, margin_pts, s_marg, s_int, s_occ=None, b0=0.0, b1=1.0, nz_occ=0.7):
    """Set cement gap s (x) [mm ]: s_marg within b 0 from the margin curve, linear ramp to s_int at b 1 ,
    s_occ where the preparation normal points occlusal (n_z > nz_occ ), softly mixed."""
    if s_occ is None:
        s_occ = s_int
    (X, Y, Z) = grid.mesh()
    tree = cKDTree(margin_pts)
    sel = np.abs(prep) < 0.6
    d = np.full(prep.shape, 10.0, np.float32)
    Q = np.stack([X[sel], Y[sel], Z[sel]], 1)
    d[sel] = tree.query(Q, workers=WORKERS)[0]
    w = np.clip((d - b0) / max(b1 - b0, 1e-06), 0, 1)
    gz = np.gradient(prep, grid.h, axis=2)
    wo = np.clip((gz - (nz_occ - 0.1)) / 0.2, 0, 1)
    s_body = (1 - wo) * s_int + wo * s_occ
    return ((1 - w) * s_marg + w * s_body).astype(np.float32)

def design_cavity(prep, grid, s, levels=(0.0,)):
    """Avsedd intaglio (kavitet) = {prep - s(x) <= 0}; variabel offset -> redistans (A143).
    Levels: levels at which the thresholds milling cell (e.g. - r_eff , - r_c )."""
    return redistance_lowmem((prep - s).astype(np.float32), grid, levels=levels)
