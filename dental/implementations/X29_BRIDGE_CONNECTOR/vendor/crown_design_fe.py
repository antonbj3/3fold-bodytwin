"crown_design_fe.py — DENT-DESIGN-CROWN: FE (CalculiX C3D10) for crown tied to die + contact load + Weibull.\n\nKedja: krona-/die-SDF (crown_design_geometry) -> STL (isotropisk ommeshning) -> gmsh tet4 (crown_design_gmsh.py,\n/usr/bin/python3) -> C3D10 genom rakkantade mittnoder (samma metod som cells/physics/fe_reference_core.py) ->\n*TIE krona(intaglio + marginring, slav-noder) mot die-yta (master) = bundet cementskikt (ANTAGANDE: cementet har\ndie-module, thickness s is included in die) -> fixed RV on die-bottom (z_m - 3 mm ) -> load case as * STEP :\n  the sphere R is set along the load direction d in the the crown's outer surface; contact points ( distance < 20 µm) is clustered;\n  Candidate contacts within 150 µm ; force fractions from the Hertz model (steel sphere displaced δ at F_ref = 1500 N ,\n  f_i ∝ (δ n_i · d − g_i )_+^ 1,5 ), contact radius a_i = sqrt(R · penetration) (floor a_min 0,25 mm ; fixed a to rubber disc/film);\n  pressure p ∝ sqrt(1 -r² /a² ) per surface facet, consistent node loads for 6 -node triangle (corner 0 , centre nodes A/ 3 );\n  friction-free: force along the normal sphere, standardized so that the component along d = 1 N (lateral residue is reported);\n  \"Training\" (directed load): force along d.\nPost-treatment: nodal tensions (.frd) -> ytfacetternas σ1 (tangential main voltage, facet mid) and\nelementens σ1 (centroid) -> Weibull-integraler I_S = Σ A (σ1+)^m, I_V = Σ V (σ1+)^m per lastfall (vid 1 N skalas\nlinear: σ ∝ F). The contact zone (within 3 a from the loading centre) is shown separately (Hertz consnips are not modelled).\nEnheter: mm, N, MPa.\n"
import os
import re
import sys
import json
import time
import subprocess
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import nnls
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'solvers'))
sys.path.insert(0, os.path.join(HERE, '..', 'manufacturing'))
import ccx_highres
import crown_fit_geometry as G
GMSH_PY = '/usr/bin/python3'
GMSH_SCRIPT = os.path.join(HERE, 'crown_design_gmsh.py')
STEEL = (200000.0, 0.3)
ZONE = 3.0

def tet_mesh(stl, out_npz, size_max, size_min=0.0):
    r = subprocess.run([GMSH_PY, GMSH_SCRIPT, stl, out_npz, str(size_max), str(size_min)], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS='4'))
    if r.returncode != 0:
        raise RuntimeError('gmsh: ' + r.stderr[-2000:])
    d = np.load(out_npz)
    (V, T) = (d['V'], d['T'])
    vol = np.einsum('ij,ij->i', np.cross(V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 0]]), V[T[:, 3]] - V[T[:, 0]])
    T = T.copy()
    neg = vol < 0
    T[neg] = T[neg][:, [0, 2, 1, 3]]
    return (V, T)

def to_tet10(V, T):
    "Centre nodes on straight edges. Abaqus/CalculiX C3D10-ordning: 5(1-2) 6(2-3) 7(3-1) 8(1-4) 9(2-4) 10(3-4)."
    pairs = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    E = np.concatenate([np.sort(T[:, list(p)], 1) for p in pairs])
    (Eu, inv) = np.unique(E, axis=0, return_inverse=True)
    mid = 0.5 * (V[Eu[:, 0]] + V[Eu[:, 1]])
    V10 = np.vstack([V, mid])
    M = inv.reshape(6, -1).T + len(V)
    return (V10, np.hstack([T, M]))
FACES10 = {1: [0, 1, 2, 4, 5, 6], 2: [0, 3, 1, 7, 8, 4], 3: [1, 3, 2, 8, 9, 5], 4: [2, 3, 0, 9, 7, 6]}

def boundary_faces10(T10):
    "Randfacetter (C3D10): (elements, face no. 1.. 4, 6 nodes [3 vertices, 3 mid])."
    F = []
    for (f, loc) in FACES10.items():
        F.append(np.column_stack([np.arange(len(T10)), np.full(len(T10), f), T10[:, loc]]))
    F = np.concatenate(F)
    key = np.sort(F[:, 2:5], 1)
    (_, idx, cnt) = np.unique(key, axis=0, return_index=True, return_counts=True)
    return F[idx[cnt == 1]]

def face_geom(V, Fb):
    (A, B, C) = (V[Fb[:, 2]], V[Fb[:, 3]], V[Fb[:, 4]])
    n = np.cross(B - A, C - A)
    area = 0.5 * np.linalg.norm(n, axis=1)
    n /= 2 * area[:, None] + 1e-300
    return ((A + B + C) / 3.0, n, area)

class Model:
    "Crown + die as two C3D10 - network tied with * TIE ."

    def __init__(self, crown_VT, die_VT, sdf_crown, sdf_die, grid, z_m, z_bot, tie_tol=None):
        (Vc, Tc) = to_tet10(*crown_VT)
        (Vd, Td) = to_tet10(*die_VT)
        (self.Vc, self.Tc, self.Vd, self.Td) = (Vc, Tc, Vd, Td)
        self.nc = len(Vc)
        (self.grid, self.z_m, self.z_bot) = (grid, z_m, z_bot)
        tol = tie_tol or 1.2 * grid.h
        Fc = boundary_faces10(Tc)
        (cc, ncv, ac) = face_geom(Vc, Fc)
        dd = grid.sample(sdf_die, cc, cval=10.0)
        tied = np.abs(dd) < tol
        (self.Fc_all, self.Fc_tied, self.Fc_free) = (Fc, Fc[tied], Fc[~tied])
        self.tied_mask = tied
        (self.cc_free, self.n_free, self.a_free) = (cc[~tied], ncv[~tied], ac[~tied])
        Fd = boundary_faces10(Td)
        (cd, nd, ad) = face_geom(Vd, Fd)
        dc = grid.sample(sdf_crown, cd, cval=10.0)
        self.Fd_master = Fd[(np.abs(dc) < tol) & (cd[:, 2] > z_m - 0.3)]
        self.fix_nodes = np.nonzero(Vd[:, 2] < z_bot + 0.001 + 0.02)[0]
        self.face_class = np.where(tied, 'intaglio', 'outer')
        (Vs, _) = G.surface_mesh(sdf_crown, grid)
        self._outer_pts = Vs[np.abs(grid.sample(sdf_die, Vs, cval=10.0)) > tol]

    def outer_points(self):
        "Points on the the crown's free outer surface of the SDF:s zero (net independent) level for spherical setting and contacts."
        return self._outer_pts

def seat_sphere(P, R, d, xy0, perp=None):
    "Sphere with centre c = c0 + s d is carried along d until first contact with the point cloud P.\n    c0 = startpunkt (xy0 in the plane perpendicular to d, long 'above'). Returns the center and p."
    d = d / np.linalg.norm(d)
    c0 = np.asarray(xy0, float) - 30.0 * d
    q = P - c0
    along = q @ d
    perp2 = np.einsum('ij,ij->i', q, q) - along ** 2
    ok = perp2 < R * R
    if not ok.any():
        return (c0, -1000000000.0)
    s_hit = along[ok] - np.sqrt(R * R - perp2[ok])
    s = float(s_hit.min())
    return (c0 + s * d, s)

def fossa_seat(P, R, d, center_guess, step0=0.3, step_min=0.005):
    "The sphere 'rolls' to local deepest position along the d starting in center_guess (mountain climbing in the plane\n    perpendicular to d, step halved): stable fossa position (typically three cup contacts). No global search,\n    So the sphere can't slip off the crown."
    d = d / np.linalg.norm(d)
    u = np.cross(d, [1.0, 0, 0])
    if np.linalg.norm(u) < 1e-06:
        u = np.cross(d, [0, 1.0, 0])
    u /= np.linalg.norm(u)
    v = np.cross(d, u)
    x = np.asarray(center_guess, float)
    (c, s) = seat_sphere(P, R, d, x)
    step = step0
    while step > step_min:
        moved = False
        for (du, dv) in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)):
            xn = x + step * (du * u + dv * v)
            (cn, sn) = seat_sphere(P, R, d, xn)
            if sn > s + 1e-07:
                (x, c, s, moved) = (xn, cn, sn, True)
                break
        if not moved:
            step /= 2
    return c

def contacts(P, C, R, tol=0.02, link=0.8, return_gaps=False):
    "Contact Candidates = local minima of gap g(p) = |p-C|-R among surface points with g < tol: a point is retained if\n    It has the lowest gap within the radius 'link' mm (non-maximum pressure). Returns points (and gap [mm ])."
    dist = np.linalg.norm(P - C, axis=1) - R
    m = dist < tol
    (Q, gq) = (P[m], dist[m])
    if len(Q) == 0:
        return (np.zeros((0, 3)), np.zeros(0)) if return_gaps else np.zeros((0, 3))
    o = np.argsort(gq)
    (Q, gq) = (Q[o], gq[o])
    tree = cKDTree(Q)
    taken = np.zeros(len(Q), bool)
    (pts, gaps) = ([], [])
    for i in range(len(Q)):
        if taken[i]:
            continue
        pts.append(Q[i])
        gaps.append(max(0.0, gq[i]))
        taken[tree.query_ball_point(Q[i], link)] = True
    (pts, gaps) = (np.array(pts), np.array(gaps))
    return (pts, gaps) if return_gaps else pts

def hertz_share(gaps, normals, d, F_ref, R, E_crown, nu_crown, E_ind=None, nu_ind=None):
    "Followive contact distribution: rigid sphere displaced δ along d; contact in with initial gap g_i carries\n    f_i = k (δ_n ,i - g_i )_+^ 1,5 (Hertz spherical plane, k = 4 / 3 E* sqrt(R)), δ_n ,i = δ ( n_i · d).\n    Σ f_i (n_i·d) = F_ref. Provides continuous dependence on the gap (a third kusp 30 or 35 µm removes almost\n    The same proportion) instead of a hard tolerance that alternates with the surface discretisation."
    E_ind = E_ind or STEEL[0]
    nu_ind = nu_ind or STEEL[1]
    Es = 1.0 / ((1 - nu_crown ** 2) / E_crown + (1 - nu_ind ** 2) / E_ind)
    k = 4.0 / 3.0 * Es * np.sqrt(R)
    c = np.clip(normals @ d, 0.001, 1.0)

    def total(delta):
        pen = np.clip(delta * c - gaps, 0, None)
        return (float(np.sum(k * pen ** 1.5 * c)), k * pen ** 1.5)
    (lo, hi) = (0.0, 1.0)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if total(mid)[0] < F_ref:
            lo = mid
        else:
            hi = mid
    (tot, f) = total(hi)
    return (f / tot, hi)

def hertz_a(f, R, E_crown, nu_crown, E_ind=STEEL[0], nu_ind=STEEL[1]):
    Es = 1.0 / ((1 - nu_crown ** 2) / E_crown + (1 - nu_ind ** 2) / E_ind)
    return (3.0 * f * R / (4.0 * Es)) ** (1.0 / 3.0)

def fossa_point(P, center_guess, radius=2.5, depth=3.0):
    "Central Foss = lowest point of the crown's occlusive surface area: points within the radius mm (xy) of the occlusive\n    centroid and maximum'depth' mm below the crown's peak (excludes axial walls and margin)."
    r = np.linalg.norm(P[:, :2] - np.asarray(center_guess)[:2], axis=1)
    Q = P[(r < radius) & (P[:, 2] > P[:, 2].max() - depth)]
    return Q[np.argmin(Q[:, 2])]

def fossa_seat_grid(P, R, d, start, bound=1.0, n=21):
    "The test machine positions the ball over the fosan: centre (angled against d) within the ‘bound' mm from the fossil point\n    is chosen where the ball sits most deeply (best wedged between the cups); lowering along d without lateral movement."
    d = d / np.linalg.norm(d)
    u = np.cross(d, [1.0, 0, 0])
    if np.linalg.norm(u) < 1e-06:
        u = np.cross(d, [0, 1.0, 0])
    u /= np.linalg.norm(u)
    v = np.cross(d, u)
    best = None
    for a in np.linspace(-bound, bound, n):
        for b in np.linspace(-bound, bound, n):
            if a * a + b * b > bound * bound:
                continue
            (c, s) = seat_sphere(P, R, d, np.asarray(start, float) + a * u + b * v)
            if best is None or s > best[1]:
                best = (c, s)
    return best[0]

def load_case(model, lc, E_crown, nu_crown, F_ref=1500.0, a_min=0.25):
    "lc: dict(kind='fossa'|'directed', R, d (3,), center (3,), a_fixed (valfri: gummiskiva/film), sticking).\n    Floss: the sphere is positioned over the central fossil (deepest insertion within 1,0 mm, no lateral movement under load);\n    Directed: the sphere is taken along d from the point of the fossa. A single contact -> force along d (friction/machine stiffness).\n    Active contacts = surface points within 2 · δ_H of the rigidly set sphere, δ_H = a_H² /R = Hertz approach at F_ref\n    (elastic penetration at test load activates more contacts; ≥ 20 µm). Returns node loads for F = 1 N."
    P = model.outer_points()
    d = np.asarray(lc['d'], float)
    d /= np.linalg.norm(d)
    R = lc['R']
    fp = fossa_point(P, lc['center'])
    if lc['kind'] == 'fossa':
        C = fossa_seat_grid(P, R, d, fp, bound=1.0)
    else:
        side = d - d @ np.array([0, 0, 1.0]) * np.array([0, 0, 1.0])
        if np.linalg.norm(side) > 0.05:
            side /= np.linalg.norm(side)
            P = P[(P - fp) @ side > 0.0]
        (C, _) = seat_sphere(P, R, d, fp)
    tol = 0.15
    (K, gaps) = contacts(P, C, R, tol=tol, return_gaps=True)
    if len(K) == 0:
        raise RuntimeError('ingen kontakt')
    normals = (K - C) / np.linalg.norm(K - C, axis=1, keepdims=True)
    (share, delta) = hertz_share(gaps, normals, d, F_ref, R, E_crown, nu_crown)
    keep = share > 0.0001
    (K, gaps, normals, share) = (K[keep], gaps[keep], normals[keep], share[keep])
    resid = 0.0
    if lc.get('sticking', False):
        fmag = share / share.sum()
        dirs = [d] * len(K)
    else:
        dirs = list(normals)
        sc = float((share[:, None] * normals).sum(0) @ d)
        fmag = share / sc
        resid = float(np.linalg.norm((fmag[:, None] * normals).sum(0) - d))
    loads = {}
    patches = []
    for (i, p) in enumerate(K):
        if fmag[i] <= 1e-09:
            continue
        pen = max(delta * float(np.clip(normals[i] @ d, 0.001, 1)) - gaps[i], 0.0)
        a = lc.get('a_fixed') or max(a_min, float(np.sqrt(R * pen)))
        r = np.linalg.norm(model.cc_free - p, axis=1)
        sel = np.nonzero(r < a)[0]
        if len(sel) < 3:
            sel = np.argsort(r)[:6]
        w = np.sqrt(np.clip(1 - (r[sel] / a) ** 2, 0.05, 1)) * model.a_free[sel]
        w = w / w.sum() * fmag[i]
        for (f_i, fi) in zip(sel, w):
            for nd in model.Fc_free[f_i, 5:8]:
                loads[int(nd)] = loads.get(int(nd), 0) + dirs[i] * fi / 3.0
        patches.append({'p': p.tolist(), 'f_per_N': float(fmag[i]), 'a_mm': float(a), 'n_faces': int(len(sel)), 'dir': np.asarray(dirs[i]).tolist()})
    tot = np.sum(list(loads.values()), axis=0)
    return (loads, {'center': C.tolist(), 'fossa_point': fp.tolist(), 'n_contacts': int(len(K)), 'patches': patches, 'contact_tol_mm': float(tol), 'gaps_mm': gaps.tolist(), 'approach_mm': float(delta), 'total_force_per_N': tot.tolist(), 'lateral_resid_per_N': float(resid)})

def fossa_seat_bounded(P, R, d, start, bound=1.0, step0=0.2, step_min=0.005):
    "As fossa_seat but the position of the sphere (angled to d) is held within'bound' mm from the starting point."
    d = d / np.linalg.norm(d)
    u = np.cross(d, [1.0, 0, 0])
    if np.linalg.norm(u) < 1e-06:
        u = np.cross(d, [0, 1.0, 0])
    u /= np.linalg.norm(u)
    v = np.cross(d, u)
    x0 = np.asarray(start, float)
    x = x0.copy()
    (c, s) = seat_sphere(P, R, d, x)
    step = step0
    while step > step_min:
        moved = False
        for (du, dv) in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)):
            xn = x + step * (du * u + dv * v)
            off = xn - x0
            if np.linalg.norm(off - off @ d * d) > bound:
                continue
            (cn, sn) = seat_sphere(P, R, d, xn)
            if sn > s + 1e-07:
                (x, c, s, moved) = (xn, cn, sn, True)
                break
        if not moved:
            step /= 2
    return c

def write_deck(path, model, mat_crown, mat_die, cases_loads, F_scale=1.0, solver='PARDISO'):
    (Vc, Tc, Vd, Td, nc) = (model.Vc, model.Tc, model.Vd, model.Td, model.nc)
    L = []
    L.append('*HEADING\nDENT-DESIGN-CROWN krona+die, mm N MPa')
    L.append('*NODE, NSET=NALL')
    L += [f'{i + 1},{x:.6f},{y:.6f},{z:.6f}' for (i, (x, y, z)) in enumerate(Vc)]
    L += [f'{i + 1 + nc},{x:.6f},{y:.6f},{z:.6f}' for (i, (x, y, z)) in enumerate(Vd)]
    ne_c = len(Tc)
    L.append('*ELEMENT, TYPE=C3D10, ELSET=ECROWN')
    L += [f'{e + 1},' + ','.join((str(int(n) + 1) for n in row)) for (e, row) in enumerate(Tc)]
    L.append('*ELEMENT, TYPE=C3D10, ELSET=EDIE')
    L += [f'{e + 1 + ne_c},' + ','.join((str(int(n) + 1 + nc) for n in row)) for (e, row) in enumerate(Td)]
    L.append('*MATERIAL, NAME=MCROWN\n*ELASTIC\n%.1f,%.4f' % mat_crown)
    L.append('*MATERIAL, NAME=MDIE\n*ELASTIC\n%.1f,%.4f' % mat_die)
    L.append('*SOLID SECTION, ELSET=ECROWN, MATERIAL=MCROWN')
    L.append('*SOLID SECTION, ELSET=EDIE, MATERIAL=MDIE')
    slave = np.unique(model.Fc_tied[:, 2:8]) + 1
    L.append('*NSET, NSET=NSLAVE')
    L += [','.join((str(x) for x in slave[i:i + 12])) for i in range(0, len(slave), 12)]
    L.append('*SURFACE, NAME=SSLAVE, TYPE=NODE\nNSLAVE')
    L.append('*SURFACE, NAME=SMASTER')
    L += [f'{int(e) + 1 + ne_c},S{int(f)}' for (e, f) in model.Fd_master[:, :2]]
    L.append('*TIE, NAME=TCEM, POSITION TOLERANCE=%.3f, ADJUST=NO\nSSLAVE,SMASTER' % (2.5 * model.grid.h))
    L.append('*NSET, NSET=NFIX')
    fx = model.fix_nodes + 1 + nc
    L += [','.join((str(x) for x in fx[i:i + 12])) for i in range(0, len(fx), 12)]
    L.append('*BOUNDARY\nNFIX,1,3,0.0')
    L.append('*NSET, NSET=NCROWN, GENERATE\n1,%d,1' % nc)
    for (name, loads) in cases_loads:
        L.append(f'*STEP\n*STATIC, SOLVER={solver}')
        L.append('*CLOAD, OP=NEW')
        for (nd, v) in loads.items():
            for k in range(3):
                if abs(v[k]) > 0:
                    L.append(f'{nd + 1},{k + 1},{v[k] * F_scale:.8e}')
        L.append('*EL FILE, OUTPUT=3D, NSET=NCROWN\nS\n*END STEP')
    open(path, 'w').write('\n'.join(L) + '\n')

def read_frd_stress(path, nmax):
    "Nodala stresses per step from .frd ( ASCII ). Returns a list of (N,6) [SXX SYY SZZ SXY SYZ SZX]."
    out = []
    cur = None
    mode = None
    with open(path) as f:
        for line in f:
            if line.startswith(' -4'):
                mode = line.split()[1]
                if mode == 'STRESS':
                    cur = np.full((nmax, 6), np.nan)
                continue
            if mode == 'STRESS':
                if line.startswith(' -1'):
                    nid = int(line[3:13])
                    vals = [float(line[13 + 12 * k:25 + 12 * k]) for k in range(6)]
                    if nid <= nmax:
                        cur[nid - 1] = vals
                elif line.startswith(' -3'):
                    out.append(cur)
                    mode = None
    return out

def sigma1(S):
    "Maximum main voltage clock (N , 6 ) Voigt [xx yy zz xy yz zx]."
    M = np.zeros((len(S), 3, 3))
    (M[:, 0, 0], M[:, 1, 1], M[:, 2, 2]) = (S[:, 0], S[:, 1], S[:, 2])
    M[:, 0, 1] = M[:, 1, 0] = S[:, 3]
    M[:, 1, 2] = M[:, 2, 1] = S[:, 4]
    M[:, 0, 2] = M[:, 2, 0] = S[:, 5]
    return np.linalg.eigvalsh(M)[:, -1]

def tangential_sigma1(S, n):
    "Maximum main voltage in the facet key plane (surface cracks are opened by tangential tensile voltage)."
    M = np.zeros((len(S), 3, 3))
    (M[:, 0, 0], M[:, 1, 1], M[:, 2, 2]) = (S[:, 0], S[:, 1], S[:, 2])
    M[:, 0, 1] = M[:, 1, 0] = S[:, 3]
    M[:, 1, 2] = M[:, 2, 1] = S[:, 4]
    M[:, 0, 2] = M[:, 2, 0] = S[:, 5]
    t1 = np.cross(n, np.where(np.abs(n[:, :1]) < 0.9, [[1.0, 0, 0]], [[0, 1.0, 0]]))
    t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    t2 = np.cross(n, t1)
    Tm = np.stack([t1, t2], 1)
    M2 = np.einsum('nai,nij,nbj->nab', Tm, M, Tm)
    return np.linalg.eigvalsh(M2)[:, -1]

def face_elem_stress(model, S_nodes):
    "Voltage stones in the center of the crown's edge facets (6-node square interpolation: vertices -1/9, mitt 4/9) and i\n    elementens centroid (tet10: corner -1/8, mitt 1/4). Saves compactly so that the post-treatment can be re-made without FE."
    S = S_nodes[:model.nc]
    F = model.Fc_all
    Sf = -1 / 9 * (S[F[:, 2]] + S[F[:, 3]] + S[F[:, 4]]) + 4 / 9 * (S[F[:, 5]] + S[F[:, 6]] + S[F[:, 7]])
    T = model.Tc
    Se = -1 / 8 * S[T[:, :4]].sum(1) + 0.25 * S[T[:, 4:]].sum(1)
    return (Sf.astype(np.float32), Se.astype(np.float32))

def model_geometry_arrays(model):
    (cf, nf, af) = face_geom(model.Vc, model.Fc_all)
    X = model.Vc[model.Tc[:, :4]]
    ve = np.abs(np.einsum('ij,ij->i', np.cross(X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]), X[:, 3] - X[:, 0])) / 6
    return dict(cf=cf, nf=nf, af=af, tied=model.tied_mask, ce=X.mean(1), ve=ve)

def postprocess_arrays(ga, Sf, Se, load_centers, a_list, m_values):
    "Weibull-integraler per lastfall vid F = 1 N (σ ∝ F).\n    Yta: tangentiell σ1 i facettmitt. The ZONE · a around the loading centre is excluded from the BARA on the outer surface\n    (Hertz ring strokes/conscraps are not modelled); intakelio is always taken with (flexions during the load are exactly the same\n    brottmod som tunna kronor har). Volym: σ1 i elementcentroid, zon = a (only the proximity of the contact)."
    (cf, nf, af, tied, ce, ve) = (ga['cf'], ga['nf'], ga['af'], ga['tied'], ga['ce'], ga['ve'])
    s1f = tangential_sigma1(Sf.astype(float), nf)
    s1e = sigma1(Se.astype(float))
    near_f = np.zeros(len(cf), bool)
    near_e = np.zeros(len(ce), bool)
    for (p, a) in zip(load_centers, a_list):
        near_f |= np.linalg.norm(cf - p, axis=1) < ZONE * a
        near_e |= np.linalg.norm(ce - p, axis=1) < a
    zone_f = near_f & ~tied
    outs = ~zone_f
    res = {'s1_surface_max_per_N': float(s1f.max()), 's1_surface_max_outside_zone_per_N': float(s1f[outs].max()), 's1_intaglio_max_per_N': float(s1f[tied].max()) if tied.any() else None, 's1_outer_max_outside_zone_per_N': float(s1f[~tied & outs].max()), 's1_volume_max_per_N': float(s1e.max()), 'argmax_surface_outside_zone': cf[outs][np.argmax(s1f[outs])].tolist(), 'argmax_on_intaglio': bool(tied[outs][np.argmax(s1f[outs])]), 'area_crown_mm2': float(af.sum()), 'volume_crown_mm3': float(ve.sum())}
    for m in m_values:
        pf = np.clip(s1f, 0, None) ** m * af
        pe = np.clip(s1e, 0, None) ** m * ve
        res[f'IS_m{m}'] = float(pf.sum())
        res[f'IS_out_m{m}'] = float(pf[outs].sum())
        res[f'IS_intaglio_frac_m{m}'] = float(pf[tied].sum() / max(pf[outs].sum(), 1e-300))
        res[f'IV_m{m}'] = float(pe.sum())
        res[f'IV_out_m{m}'] = float(pe[~near_e].sum())
    return res

def prepare_run(workdir, name, model, mat_crown, mat_die, lcs, solver='ITERATIVE CHOLESKY'):
    "Phase 1 ( local): loads + tyres. Returns (name list, diagnostics); the tyre is written to workdir/name.inp ."
    os.makedirs(workdir, exist_ok=True)
    (cases, diags) = ([], [])
    for lc in lcs:
        (loads, diag) = load_case(model, lc, mat_crown[0], mat_crown[1], F_ref=lc.get('F_ref', 1500.0))
        cases.append((lc['name'], loads))
        diags.append(diag)
    write_deck(os.path.join(workdir, name + '.inp'), model, mat_crown, mat_die, cases, solver=solver)
    return ([c[0] for c in cases], diags)

def finish_run(frd_path, model, names, diags, m_values, stress_npz):
    "Phase 2 : read .frd (also .frd.gz), Weibull finishing, save compact stresses."
    import gzip
    if frd_path.endswith('.gz'):
        import shutil, tempfile
        tmp = frd_path[:-3]
        with gzip.open(frd_path, 'rb') as fi, open(tmp, 'wb') as fo:
            shutil.copyfileobj(fi, fo)
        Ss = read_frd_stress(tmp, model.nc + len(model.Vd))
        os.remove(tmp)
    else:
        Ss = read_frd_stress(frd_path, model.nc + len(model.Vd))
    assert len(Ss) == len(names), f'{len(Ss)} voltage blocks against: {len(names)} lastfall'
    out = {'n_nodes': int(model.nc + len(model.Vd)), 'n_el_crown': int(len(model.Tc)), 'n_el_die': int(len(model.Td)), 'dof': int(3 * (model.nc + len(model.Vd))), 'cases': {}}
    ga = model_geometry_arrays(model)
    save = {k: v.astype(np.float32) if v.dtype.kind == 'f' else v for (k, v) in ga.items()}
    for (nm, diag, S) in zip(names, diags, Ss):
        centers = [p['p'] for p in diag['patches']]
        alist = [p['a_mm'] for p in diag['patches']]
        (Sf, Se) = face_elem_stress(model, S)
        (save[f'Sf__{nm}'], save[f'Se__{nm}']) = (Sf, Se)
        save[f'lc__{nm}'] = np.array([list(c) + [a] for (c, a) in zip(centers, alist)])
        pp = postprocess_arrays(ga, Sf, Se, np.array(centers), alist, m_values)
        pp['load'] = diag
        out['cases'][nm] = pp
    np.savez_compressed(stress_npz, **save)
    out['stress_npz'] = stress_npz
    return out

def run_model(workdir, name, model, mat_crown, mat_die, lcs, m_values, keep_frd=False, ram_gb=8, solver='PARDISO'):
    "Build loads, write tires, run PARDISO -ccx via heavy_run , finishing. Returns result per load case."
    os.makedirs(workdir, exist_ok=True)
    (cases, diags) = ([], [])
    for lc in lcs:
        (loads, diag) = load_case(model, lc, mat_crown[0], mat_crown[1], F_ref=lc.get('F_ref', 1500.0))
        cases.append((lc['name'], loads))
        diags.append(diag)
    job = os.path.join(workdir, name)
    write_deck(job + '.inp', model, mat_crown, mat_die, cases, solver=solver)
    t0 = time.time()
    r = ccx_highres.run(job, ram_gb=ram_gb, threads=4)
    if not r['ok']:
        raise RuntimeError(f'ccx misslyckades: {r} ' + open(job + '.log').read()[-3000:])
    Ss = read_frd_stress(job + '.frd', model.nc + len(model.Vd))
    out = {'ccx': r, 'n_nodes': int(model.nc + len(model.Vd)), 'n_el_crown': int(len(model.Tc)), 'n_el_die': int(len(model.Td)), 'dof': int(3 * (model.nc + len(model.Vd))), 'cases': {}}
    ga = model_geometry_arrays(model)
    save = {k: v.astype(np.float32) if v.dtype.kind == 'f' else v for (k, v) in ga.items()}
    for ((nm, _), diag, S) in zip(cases, diags, Ss):
        centers = [p['p'] for p in diag['patches']]
        alist = [p['a_mm'] for p in diag['patches']]
        (Sf, Se) = face_elem_stress(model, S)
        (save[f'Sf__{nm}'], save[f'Se__{nm}']) = (Sf, Se)
        save[f'lc__{nm}'] = np.array([list(c) + [a] for (c, a) in zip(centers, alist)])
        pp = postprocess_arrays(ga, Sf, Se, np.array(centers), alist, m_values)
        pp['load'] = diag
        out['cases'][nm] = pp
    np.savez_compressed(job + '_stress.npz', **save)
    out['stress_npz'] = job + '_stress.npz'
    if not keep_frd:
        for ext in ('.frd', '.dat', '.sta', '.cvg', '.12d'):
            try:
                os.remove(job + ext)
            except FileNotFoundError:
                pass
        try:
            os.remove(job + '.inp')
        except FileNotFoundError:
            pass
    out['seconds'] = time.time() - t0
    return out
