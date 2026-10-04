"""Conditional quasi-static compliance contact. External engine sources are read only.

All lengths mm, forces N, rotations rad, torques N mm. This is a reduced
spring/elastic-patch model, not a resolved crown/PDL finite element calculation.
"""
from dental_release.paths import expand as _release_expand
import os, sys
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '4'
sys.dont_write_bytecode = True
import hashlib, importlib.util, json, time, datetime, zipfile, resource
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import nnls, root
P = Path(__file__).resolve().parents[1]
D = Path(os.environ.get('X54_DATA_DIR', _release_expand('@DENTAL_WORK_ROOT@/X54')))
BASE = P.parent
ZIP = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Bits2Bites/Bits2Bites_v01.zip'))
LABELS = Path(_release_expand('@DENTAL_WORK_ROOT@/X11/targets/Bits2Bites'))
ENGINE = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/3fold_release_motion_staging_20260922'))
_BUNDLE = json.loads((P / 'REVIEW_BUNDLE.json').read_text())['files'] if (P / 'REVIEW_BUNDLE.json').exists() else {}

def local_path(path):
    return P / _BUNDLE[str(path)] if str(path) in _BUNDLE else Path(path)

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, local_path(path))
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m
loader = load_module('x54_p1_loader', ENGINE / 'reproducibility/p1/data/ncp/loader.py')
ncp = load_module('x54_p2_ncp', ENGINE / 'reproducibility/p2/solver/ncp_ref.py')

def labels_adapter(lam, G, b, mu, tol_lam=1e-10, tol_rel=1e-06):
    l = np.asarray(lam).reshape(-1, 3)
    scale = max(np.max(np.abs(l)), 1e-30)
    ln = l[:, 0]
    lt = np.linalg.norm(l[:, 1:], axis=1)
    return np.where(ln <= tol_lam + tol_rel * scale, 'open', np.where(lt >= mu * ln * (1 - tol_rel) - tol_lam, 'slip', 'stick'))
ncp.classify = labels_adapter
sys.path.insert(0, str(P / 'review_bundle/imports'))
sys.path.insert(0, str(BASE / 'LANE_X2_OCCLUSION_B2B/code'))
from occlusion_operator import raster
sys.path.insert(0, str(BASE / 'LANE_X11_TOOTH_SEG/code'))
from morphology import mesh_stl

def sha(path):
    h = hashlib.sha256()
    with open(local_path(path), 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, Path):
        return str(x)
    return x

def write(path, x):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(clean(x), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(path)

def state(stage, gate, next_op, **extra):
    write(P / 'CURRENT_WORK_STATE.json', dict(lane='X54-deformable-contact', status=stage, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), latest_gate=gate, next_operation=next_op, data_directory=D, disk_limit_bytes=3000000000, **extra))

def geometry(case, grid=0.2):
    """Reuses X2's exact barycentric raster. No source pose adjustment."""
    meta = json.loads((LABELS / str(case) / 'labels+landmarks.json').read_text())
    fr = meta['arches']['lower']['frame']
    R = np.column_stack([fr['right_unit'], fr['anterior_unit'], fr['superior_unit']])
    origin3 = np.asarray(fr['center_mm'])
    if np.max(np.abs(R.T @ R - np.eye(3))) > 1e-10:
        raise ValueError('FRAME_NOT_ORTHOGONAL')
    with zipfile.ZipFile(ZIP) as z:
        names = {int(n.split('/')[-2]): n for n in z.namelist() if n.endswith('upper.stl')}
        arches = []
        sources = []
        for jaw in ('upper', 'lower'):
            member = names[case].replace('upper.stl', jaw + '.stl')
            blob = z.read(member)
            a = meta['arches'][jaw]
            if hashlib.sha256(blob).hexdigest() != a['input_sha256']:
                raise ValueError('SOURCE_HASH_MISMATCH')
            lp = Path(a['labels']['path'])
            if sha(lp) != a['labels']['sha256']:
                raise ValueError('LABEL_HASH_MISMATCH')
            (v, f) = mesh_stl(blob, fast=True)
            lab = np.load(lp)['labels']
            if len(v) != len(lab):
                raise ValueError('VERTEX_ORDER_MISMATCH')
            fl = lab[f]
            (aa, bb, cc) = fl.T
            owner = np.where((aa == bb) | (aa == cc), aa, np.where(bb == cc, bb, 0))
            allowed = [10 * q + t for q in ([1, 2] if jaw == 'upper' else [3, 4]) for t in range(1, 8)]
            keep = np.isin(owner, allowed)
            allT = ((v.astype(float) - origin3) @ R)[f]
            T = allT[keep]
            lab = owner[keep]
            ids = np.flatnonzero(keep)
            cross = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
            length = np.linalg.norm(cross, axis=1)
            normals = cross / np.maximum(length[:, None], 1e-30)
            normals *= np.where(normals[:, 2] < 0, -1.0, 1.0)[:, None]
            centers = {int(k): np.mean(T[lab == k], axis=(0, 1)) for k in np.unique(lab)}
            sources.append(dict(jaw=jaw, zip_member=member, sha256=a['input_sha256'], label_path=lp, label_sha256=a['labels']['sha256'], input_faces=len(f), retained_faces=len(T), omitted_face_fraction=1 - float(keep.mean()), mixed_vertex_label_fraction=float(np.any(fl != fl[:, 0, None], axis=1).mean())))
            arches.append(dict(T=T, lab=lab, ids=ids, normals=normals, centers=centers))
    (U, L) = arches
    lo = np.minimum(U['T'][:, :, :2].min((0, 1)), L['T'][:, :, :2].min((0, 1)))
    hi = np.maximum(U['T'][:, :, :2].max((0, 1)), L['T'][:, :, :2].max((0, 1)))
    origin = np.floor(lo / grid) * grid - 2 * grid
    shape = np.ceil((hi - origin) / grid).astype(int) + 3
    (uz, ui) = raster(U['T'], grid, origin, shape, True)
    (lz, li) = raster(L['T'], grid, origin, shape, False)
    good = np.isfinite(uz) & np.isfinite(lz)
    ix = np.flatnonzero(good)
    u = ui[good]
    l = li[good]
    xy = origin + (np.c_[ix // shape[1], ix % shape[1]] + 0.5) * grid
    g = uz[good] - lz[good]
    un = U['normals'][u]
    ln = L['normals'][l]
    norm = un + ln
    norm /= np.maximum(np.linalg.norm(norm, axis=1, keepdims=True), 1e-30)
    orient = (un[:, 2] >= 0.35) & (ln[:, 2] >= 0.35) & (np.einsum('ij,ij->i', un, ln) >= 0.5)
    candidate = (g <= 0.5) & orient
    ul = U['lab'][u]
    ll = L['lab'][l]
    area = {jaw: {str(int(k)): float(np.sum((labs == k) & (g <= 0.1)) * grid ** 2) for k in np.unique(labs)} for (jaw, labs) in [('upper', ul), ('lower', ll)]}
    nearest = {jaw: {str(int(k)): float(np.min(g[labs == k])) for k in np.unique(labs)} for (jaw, labs) in [('upper', ul), ('lower', ll)]}
    patchgroups = {}
    for i in np.flatnonzero(candidate):
        a = int(ul[i])
        b = int(ll[i])
        c = U['centers'][a]
        radial = np.array([c[0], c[1] + 15.0])
        radial /= max(np.linalg.norm(radial), 1e-12)
        mesial = np.array([-radial[1], radial[0]])
        delta = xy[i] - c[:2]
        region = 2 * int(delta @ radial >= 0) + int(delta @ mesial >= 0)
        patchgroups.setdefault((a, b, region), []).append(i)
    patches = []
    for (key, indices) in sorted(patchgroups.items()):
        (a, b, region) = key
        indices = np.asarray(indices)
        i = indices[np.argmin(g[indices])]
        n = norm[i]
        t = np.array([1.0, 0, 0])
        t -= np.dot(t, n) * n
        t /= np.linalg.norm(t)
        t2 = np.cross(n, t)
        patches.append(dict(upper_fdi=a, lower_fdi=b, region=region, gap_mm=float(g[i] * n[2]), vertical_gap_mm=float(g[i]), basis=np.vstack([n, t, t2]), upper_point_mm=np.r_[xy[i], uz[good][i]], lower_point_mm=np.r_[xy[i], lz[good][i]], upper_source_face=int(U['ids'][u[i]]), lower_source_face=int(L['ids'][l[i]]), raster_witness_index=int(ix[i]), area_mm2=float(np.sum(grid ** 2 / norm[indices, 2])), projected_near_area_mm2=float(np.sum(g[indices] <= 0.1) * grid ** 2), sample_count=len(indices)))
    if len(patches) > 192:
        raise ValueError('CONTACT_BUDGET_EXCEEDED:' + str(len(patches)))
    return clean(dict(case=case, grid_mm=grid, frame=dict(R=R, center=origin3), sources=sources, centers={jaw: a['centers'] for (jaw, a) in [('upper', U), ('lower', L)]}, patches=patches, rigid_area_mm2=area, nearest_gap_mm=nearest, census=dict(overlap_samples=len(g), near_samples=int(np.sum(g <= 0.1)), candidate_samples=int(candidate.sum()), rejected_normal_candidates=int(np.sum((g <= 0.5) & ~orient)), candidate_dropout_fraction=float(np.sum((g <= 0.5) & ~orient) / max(np.sum(g <= 0.5), 1)), negative_gap_below_minus03=int(np.sum(g < -0.03)), minimum_gap_mm=float(g.min()), unloaded_pose_compatible=bool(g.min() >= -0.03)), resolution='PER_SURFACE_REGION', scope='predictedFDI, sampled projected envelope, regional minimum witness and equivalent area; continuous3DphysicalcontactUNKNOWN'))

def operator(geom, support_scale=None, crown=True, rotation=True, mu=0.2):
    """Use P1 rebuild for incidence and for W. Minv here means K^-1, not mass."""
    patches = geom['patches']
    teeth = [(jaw, int(t)) for jaw in ('upper', 'lower') for t in sorted(geom['centers'][jaw], key=int)]
    lookup = {k: i for (i, k) in enumerate(teeth)}
    bodies = []
    for (i, (jaw, t)) in enumerate(teeth):
        scale = 1 if support_scale is None else support_scale.get((jaw, t), 1.0)
        c = np.array([1 / 75.0, 1 / 75.0, 1 / 750.0, 1 / 7500.0, 1 / 7500.0, 1 / 7500.0]) / scale
        if not rotation:
            c[3:] = 0
        bodies.append(dict(id=i, dof=6, minv_block=np.diag(c).tolist(), block_dim=6))
    contacts = []
    for (i, p) in enumerate(patches):
        a = lookup['upper', p['upper_fdi']]
        b = lookup['lower', p['lower_fdi']]
        contacts.append(dict(id=i, a=a, b=b, normal=p['basis'][0], t1=p['basis'][1], t2=p['basis'][2], lever_a=(np.asarray(p['upper_point_mm']) - geom['centers']['upper'][str(p['upper_fdi'])]).tolist(), lever_b=(np.asarray(p['lower_point_mm']) - geom['centers']['lower'][str(p['lower_fdi'])]).tolist()))
    scene = dict(model='body_contacts', n_c=len(patches), bodies=bodies, contacts=contacts, eta=[0.0] * len(patches), b=[0.0] * (3 * len(patches)), mu=[mu] * len(patches))
    (_, _, _, J, C) = loader.rebuild(scene)
    radius = np.sqrt(np.array([p['area_mm2'] for p in patches]) / np.pi)
    E = 97000.0
    nu = 0.3
    shear = E / (2 * (1 + nu))
    cn = 2 * (1 - nu ** 2) / (2 * radius * E)
    ct = 2 * (2 - nu) / (8 * radius * shear)
    ce = np.column_stack([cn, ct, ct]).ravel() if crown else np.zeros(3 * len(patches))
    JJ = np.column_stack([J, np.eye(3 * len(patches))])
    CC = np.zeros((JJ.shape[1], JJ.shape[1]))
    CC[:len(C), :len(C)] = C
    CC[len(C):, len(C):] = np.diag(ce)
    s = dict(model='operator', n_c=len(patches), operator=dict(J=JJ.tolist(), Minv=CC.tolist()), eta=[0.0] * len(patches), b=[0.0] * (3 * len(patches)), mu=[mu] * len(patches))
    (W, _, mus, _, _) = loader.rebuild(s)
    return (W, J, C, ce, mus, teeth)

def load_vector(geom, closure=0.05, lateral=0.0, offset=0.0):
    b = []
    for p in geom['patches']:
        Dd = np.asarray(p['basis'])
        disp = Dd @ np.array([-lateral, 0.0, -closure])
        disp[0] += p['gap_mm'] + offset * Dd[0, 2]
        b.extend(disp)
    return np.array(b)

def normal_nnls(W, b):
    """Independent frictionless QP: L L.T=A, min||L.T*x+L^-1*b||²."""
    inds = np.arange(0, len(b), 3)
    A = W[np.ix_(inds, inds)]
    bb = b[inds]
    L = cholesky(A, lower=True)
    rhs = solve_triangular(L, -bb, lower=True)
    (x, _) = nnls(L.T, rhs, maxiter=10000, atol=1e-12)
    lam = np.zeros(len(b))
    lam[inds] = x
    return lam

def solve(geom, closure=0.05, lateral=0.0, offset=0.0, mu=0.2, support_scale=None, crown=True, rotation=True, controls=False, iters=2000):
    st = time.perf_counter()
    cpu = time.process_time()
    (W, J, C, ce, mus, teeth) = operator(geom, support_scale, crown, rotation, mu)
    b = load_vector(geom, closure, lateral, offset)
    prep = time.perf_counter() - st
    tick = time.perf_counter()
    warm = normal_nnls(W, b)
    (lam, hist, labels) = ncp.solve_ncp_pgs(W, b, mus, iters=iters, tol=1e-07, warm=warm)
    elapsed = time.perf_counter() - tick
    L = lam.reshape(-1, 3)
    gap = (W @ lam + b).reshape(-1, 3)
    active = L[:, 0] > 1e-07
    rho = ncp._rho_of(W, len(mus))
    res = ncp.natural_residual(lam, W, b, mus, rho)
    nc = float(max(0, -gap[:, 0].min()))
    comp = float(np.max(np.abs(L[:, 0] * gap[:, 0])))
    cone = float(max(0, np.max(np.linalg.norm(L[:, 1:], axis=1) - mu * L[:, 0])))
    forces = np.array([np.asarray(p['basis']).T @ f for (p, f) in zip(geom['patches'], L)])
    normalforce = {jaw: {} for jaw in ('upper', 'lower')}
    verticalforce = {jaw: {} for jaw in ('upper', 'lower')}
    vectors = {jaw: {} for jaw in ('upper', 'lower')}
    for (jaw, t) in teeth:
        ids = np.array([p[jaw + '_fdi'] == t for p in geom['patches']])
        normalforce[jaw][str(t)] = float(L[ids, 0].sum())
        verticalforce[jaw][str(t)] = float(forces[ids, 2].sum())
        vectors[jaw][str(t)] = ((1 if jaw == 'upper' else -1) * forces[ids].sum(0)).tolist()
    shares = {jaw: {t: 100 * f / max(sum(normalforce[jaw].values()), 1e-30) for (t, f) in normalforce[jaw].items()} for jaw in ('upper', 'lower')}
    wrench = J.T @ lam
    q = C @ wrench
    rigid_q = q.reshape(-1, 6)
    rows = [dict(**p, lambda_N=f, vector_on_upper_N=force, post_gap_mm=d, state=str(s), mean_normal_traction_MPa=float(f[0] / p['area_mm2']), traction_scope='regional equivalent-area mean, not peak pressure') for (p, f, force, d, s) in zip(geom['patches'], L, forces, gap, labels)]
    equal = float(np.linalg.norm(sum((np.asarray(v) for jaw in vectors.values() for v in jaw.values()), np.zeros(3))))
    gates = dict(natural_residual=res <= 1e-06, nonpenetration=nc <= 1e-06, complementarity=comp <= 1e-05, coulomb_cone=cone <= 1e-07, equal_opposite=equal <= 1e-08)
    ctl = {}
    if controls:
        t = time.perf_counter()
        zero = np.zeros_like(mus)
        (a, h, s) = ncp.solve_ncp_pgs(W, b, zero, iters=2000, tol=1e-07, warm=warm)
        ctl['frictionless_NNLS'] = dict(max_difference_N=float(np.max(np.abs(a - warm))), natural_residual_N=ncp.natural_residual(warm, W, b, zero), wall_s=time.perf_counter() - t, pass_gate=bool(np.max(np.abs(a - warm)) <= 0.0001))

        def mapping(x):
            u = (W @ x + b).reshape(-1, 3)
            xx = x.reshape(-1, 3)
            return (xx - ncp.proj_cone(xx - rho[:, None] * (u + ncp.desaxce(u, mus)), mus)).ravel()
        t = time.perf_counter()
        guess = lam + np.random.default_rng(54).normal(0, 0.01, len(lam)) * max(1, np.max(np.abs(lam)))
        rr = root(mapping, guess, method='krylov', options=dict(fatol=1e-07, maxiter=100))
        diff = float(np.max(np.abs(rr.x - lam)))
        rrres = float(np.max(np.abs(mapping(rr.x))))
        ctl['independent_scipy_root'] = dict(max_difference_N=diff, natural_residual_N=rrres, success=bool(rr.success), wall_s=time.perf_counter() - t, pass_gate=bool(diff <= 0.0001 and rrres <= 1e-06))
        bad = lam.copy()
        bad[int(np.argmax(np.abs(bad)))] += max(1, 0.2 * np.max(np.abs(lam)))
        badres = float(np.max(np.abs(mapping(bad))))
        ctl['injected_force_error'] = dict(natural_residual_N=badres, rejected=bool(badres > 1e-06))
    out = dict(case=geom['case'], closure_mm=closure, lateral_mm=lateral, gap_offset_scenario_mm=offset, mu=mu, normal_force_N=normalforce, vertical_force_N=verticalforce, force_shares_pp=shares, net_force_on_teeth_N=vectors, total_normal_N=float(L[:, 0].sum()), total_vertical_N=float(forces[:, 2].sum()), argmax={jaw: max(shares[jaw], key=shares[jaw].get) for jaw in shares}, rows=rows, support_motion=[dict(jaw=jaw, fdi=t, translation_mm=rigid_q[i, :3], rotation_rad=rigid_q[i, 3:], wrench_N_Nmm=wrench.reshape(-1, 6)[i]) for (i, (jaw, t)) in enumerate(teeth)], numerics=dict(natural_residual_N=res, minimum_normal_gap_mm=float(gap[:, 0].min()), normal_nonpenetration_violation_mm=nc, normal_complementarity_Nmm=comp, cone_violation_N=cone, equal_opposite_N=equal, iterations=len(hist), gates=gates, all_pass=all(gates.values()), W_min_eigenvalue_mm_per_N=float(np.linalg.eigvalsh(W)[0])), controls=ctl, cost=dict(operator_s=prep, solve_s=elapsed, total_s=time.perf_counter() - st, cpu_s=time.process_time() - cpu, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), resolution='PER_SURFACE_REGION', claim_type=['information_link', 'capability'], physical_status='CONDITIONAL_SIMULATION; support/load/preload/FDI/sensor UNKNOWN', branches={str(k): int(np.sum(labels == k)) for k in np.unique(labels)}, rigorous_error_enclosure='MISSING: sampled geometry, reduced elasticity and floating NCP residual are not a continuum or parameter-box enclosure')
    return (clean(out), dict(W=W, b=b, lambda_N=lam, J=J, C=C, crown_C=ce, mu=mus, residual_history=hist))
