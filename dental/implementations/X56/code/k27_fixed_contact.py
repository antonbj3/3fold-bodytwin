"""DENT - PROC - MICROMOTION — FE with contact: threaded implant in bone / leg block directly after insertion .

Cell (MICRO_motion, 2026-09-24). Independent: in = Spec-JSON (geometry, material, preload, load, network, contact),
Output = a JSONL line with displacement at the point/cronee and relative tangential displacement along the boundary surface.
Units: mm , N , MPa , µm in outputs .

Model (3D semimodel y ≥ 0 , symmetry plane y = 0 contains the load vector):
- Implant: rotationally symmetric body with ANNULAR V threads (profile from cells/physics/implant_thread_fe.py
  thread_profile_points , K1 ; helix angle ~ 3° neglected), smooth collar top, Ti-Distance (cylinder r_ab ) to
  load height e above the crest. Load evenly distributed over the top surface of the distance (consistent node weights).
- Bone : block (width W, depth Hb) with conventional threaded holes = surface of the implant (thread shaped/cut in the leg).
  Condition: cortical crown strap 0 .. t_c ( E_c ) + cancellous /skum ( E_s ). Linear-elastic.
- Interface : the nodes are duplicated (leg side) -> surface contact ccx SURFACE TO SURFACE (one pair per bone layer, penalty
  K = kfac·E_skikt/h_if, Coulomb-friktion mu, stick-lutning λ = lamfac·K/10). Slav = benytan.
- Preload (initial state of PROC - INSERT): the implant is administered radial/tangential thermal expansion
  ( ORTHO , cylindrical orientation, α_z = 0 ) with temperature T(z) per layer so that the mean radial pressure per layer =
  p0_radial from cells/procedure/implant_insertion.py (calibrated in a pre-run, see calibrate()).
- Step 1 preload, step 2 load F in direction (sin α, 0, − cos α) (α from implant axis).
- Micromotion = relative displacement bone − implant for each duplicated node pair, step 2 minus step 1, divided into:
  tangential (sliding) and normal (opening) component against local surface normal. Maximum and area-weighted funds are reported,
  whether or not per layer. the displacement of the implant at the loading point (mean of the peak surface nodes) and at the crest (axis, z = 0);
RV: 'bottom' = blockets botten fast; 'lateral_lower' = lower portion (fraction fr) of the lateral surfaces fixed (PMC5577443).
Solution: CalculiX 2.23 PARDISO building ( cells/solvers/ccx_highres.py ), gmsh 4.15 . Runs via tasks/heavy_run.sh.
Run: python3.10 micromotion_contact_fe.py <spec.json> <ut.jsonl> <working-directory>
"""
from dental_release.paths import expand as _release_expand
import json
import math
import os
import re
import subprocess
import sys
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from implant_thread_fe import thread_profile_points
PARDISO = _release_expand('@DENTAL_EXTERNAL_ROOT@/shared-storage/scratch/dental_F2/opt/src/CalculiX/ccx_2.23/src/ccx_2.23_pardiso')
CCX = os.environ.get('CCX', PARDISO if os.path.exists(PARDISO) else _release_expand('@DENTAL_EXTERNAL_ROOT@/.local/opt/squashfs-root/usr/bin/ccx'))
SOLVER = 'PARDISO' if CCX.endswith('pardiso') else 'SPOOLES'
ALPHA = 0.001
DEFAULT = dict(D=4.1, L=10.0, h=0.35, P=0.8, rho=0.08, collar=0.4, r_ab=None, e=7.0, W=16.0, Hb=None, t_c=1.0, E_c=11500.0, nu_c=0.3, E_s=200.0, nu_s=0.3, E_ti=110000.0, nu_ti=0.34, mu=0.3, p0_c=0.0, p0_s=0.0, F=30.0, angle=90.0, bc='bottom', bc_fr=0.5, h_if=0.16, h_far=1.6, sig_inf=0.01, c0=0.001, K_pen=None, lam=None, T_c=None, T_s=None, inc0=0.5, kfac=10.0, lamfac=1.0, calib=True, Tfac_c=1.0, Tfac_s=1.0)

def build(p):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.option.setNumber('General.NumThreads', int(os.environ.get('OMP_NUM_THREADS', '4')))
    gmsh.model.add('mm')
    occ = gmsh.model.occ
    (D, L, h, P) = (p['D'], p['L'], p['h'], p['P'])
    (rM, rm) = (D / 2, D / 2 - h)
    r_ab = p['r_ab'] or min(rM, 1.75)
    e = p['e']
    Hb = p['Hb'] or L + 6.0
    W = p['W']
    n_thr = int(math.floor((L - p['collar']) / P))
    z0 = -L
    (segs, crest) = thread_profile_points(rM, rm, P, p['rho'], n_thr)
    cache = {}

    def PT(r, z):
        k = (round(r, 9), round(z, 9))
        if k not in cache:
            cache[k] = occ.addPoint(r, 0.0, z)
        return cache[k]
    cs = [occ.addLine(PT(0.0, z0), PT(rM, z0))]
    for s in segs:
        if s[0] == 'line':
            (a, b) = ((s[1][0], s[1][1] + z0), (s[2][0], s[2][1] + z0))
            if abs(a[0] - b[0]) + abs(a[1] - b[1]) < 1e-09:
                continue
            cs.append(occ.addLine(PT(*a), PT(*b)))
        else:
            c = occ.addPoint(rm + p['rho'], 0.0, s[2][1] + z0)
            cs.append(occ.addCircleArc(PT(s[1][0], s[1][1] + z0), c, PT(s[3][0], s[3][1] + z0)))
    zt = z0 + n_thr * P
    rest = [(rM, zt), (rM, 0.0)]
    if abs(r_ab - rM) > 1e-09:
        rest.append((r_ab, 0.0))
    rest += [(r_ab, e), (0.0, e), (0.0, z0)]
    for i in range(len(rest) - 1):
        (a, b) = (rest[i], rest[i + 1])
        if abs(a[0] - b[0]) + abs(a[1] - b[1]) < 1e-09:
            continue
        cs.append(occ.addLine(PT(*a), PT(*b)))
    f_imp = occ.addPlaneSurface([occ.addCurveLoop(cs)])
    occ.synchronize()
    v_imp = [t for (d, t) in occ.revolve([(2, f_imp)], 0, 0, 0, 0, 0, 1, math.pi) if d == 3]
    b_c = occ.addBox(-W / 2, 0, -p['t_c'], W, W / 2, p['t_c']) if p['t_c'] > 0 else None
    b_s = occ.addBox(-W / 2, 0, -Hb, W, W / 2, Hb - p['t_c'])
    boxes = [b for b in (b_c, b_s) if b is not None]
    (frag, fmap) = occ.fragment([(3, v) for v in v_imp], [(3, b) for b in boxes])
    occ.synchronize()
    imp_set = set((t for (d, t) in fmap[0]))
    groups = {'IMP': sorted(imp_set), 'CORT': [], 'CAN': []}
    k = len(v_imp)
    for b in boxes:
        ch = [t for (d, t) in fmap[k] if t not in imp_set]
        groups['CORT' if b == b_c else 'CAN'] += ch
        k += 1
    groups = {g: sorted(set(v)) for (g, v) in groups.items() if v}

    def bnd(vs):
        return set((abs(t) for (d, t) in gmsh.model.getBoundary([(3, v) for v in vs], oriented=False)))
    bone_v = groups.get('CORT', []) + groups['CAN']
    iface = sorted(bnd(groups['IMP']) & bnd(bone_v))
    F = gmsh.model.mesh.field
    fd = F.add('Distance')
    F.setNumbers(fd, 'SurfacesList', iface)
    F.setNumber(fd, 'Sampling', 60)
    ft = F.add('Threshold')
    F.setNumber(ft, 'InField', fd)
    F.setNumber(ft, 'SizeMin', p['h_if'])
    F.setNumber(ft, 'SizeMax', p['h_far'])
    F.setNumber(ft, 'DistMin', 0.25)
    F.setNumber(ft, 'DistMax', 5.0)
    F.setAsBackgroundMesh(ft)
    for key in ('Mesh.MeshSizeExtendFromBoundary', 'Mesh.MeshSizeFromPoints', 'Mesh.MeshSizeFromCurvature'):
        gmsh.option.setNumber(key, 0)
    gmsh.option.setNumber('Mesh.Algorithm3D', 10)
    gmsh.option.setNumber('Mesh.ElementOrder', 2)
    gmsh.option.setNumber('Mesh.SecondOrderLinear', 1)
    gmsh.option.setNumber('Mesh.HighOrderOptimize', 0)
    gmsh.model.mesh.generate(3)
    (ntags, coords, _) = gmsh.model.mesh.getNodes()
    nodes = coords.reshape(-1, 3)
    tag2i = np.full(int(ntags.max()) + 1, -1, np.int64)
    tag2i[ntags.astype(np.int64)] = np.arange(len(ntags))
    elems = {}
    for (g, vs) in groups.items():
        Ls = []
        for v in vs:
            (et, _, en) = gmsh.model.mesh.getElements(3, v)
            for (ty, n) in zip(et, en):
                if ty == 11:
                    Ls.append(tag2i[np.asarray(n, np.int64)].reshape(-1, 10))
        elems[g] = np.vstack(Ls)[:, [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]]
    gmsh.finalize()
    for g in elems:
        t = elems[g]
        v = np.einsum('ij,ij->i', np.cross(nodes[t[:, 1]] - nodes[t[:, 0]], nodes[t[:, 2]] - nodes[t[:, 0]]), nodes[t[:, 3]] - nodes[t[:, 0]])
        neg = v < 0
        t[neg] = t[neg][:, [0, 2, 1, 3, 6, 5, 4, 7, 9, 8]]
    used = np.unique(np.concatenate([e_.ravel() for e_ in elems.values()]))
    remap = -np.ones(len(nodes), np.int64)
    remap[used] = np.arange(len(used))
    nodes = nodes[used]
    for g in elems:
        elems[g] = remap[elems[g]]
    imp_n = np.unique(elems['IMP'])
    bone_n = np.unique(np.concatenate([elems[g].ravel() for g in elems if g != 'IMP']))
    shared = np.intersect1d(imp_n, bone_n)
    new_id = -np.ones(len(nodes), np.int64)
    new_id[shared] = len(nodes) + np.arange(len(shared))
    nodes = np.vstack([nodes, nodes[shared]])
    for g in elems:
        if g == 'IMP':
            continue
        t = elems[g]
        m = new_id[t] >= 0
        t[m] = new_id[t][m]
    info = dict(n_thr=n_thr, crest_flat=crest, n_nodes=int(len(nodes)), n_pairs=int(len(shared)), n_tets={g: int(len(e_)) for (g, e_) in elems.items()}, Hb=Hb, r_ab=r_ab)
    pairs = np.stack([shared, new_id[shared]], 1)
    return (nodes, elems, pairs, info)
FACES = {1: (0, 1, 2, 4, 5, 6), 2: (0, 3, 1, 7, 8, 4), 3: (1, 3, 2, 8, 9, 5), 4: (2, 3, 0, 9, 7, 6)}

def faces_on(tets, mask):
    out = []
    for (fno, idx) in FACES.items():
        m = np.all(mask[tets[:, list(idx)]], axis=1)
        out += [(int(e), fno) for e in np.where(m)[0]]
    return out

def face_geom(nodes, tets, lst):
    """Area, unit normal (outwards from the element), centroid of the corner triangle."""
    (A, N, C, IDX) = ([], [], [], [])
    for (e, fno) in lst:
        idx = tets[e, list(FACES[fno])]
        X = nodes[idx[:3]]
        n = np.cross(X[1] - X[0], X[2] - X[0])
        a = 0.5 * np.linalg.norm(n)
        c = X.mean(0)
        inner = nodes[tets[e, :4]].mean(0)
        n = n / (2 * a)
        if np.dot(n, c - inner) < 0:
            n = -n
        A.append(a)
        N.append(n)
        C.append(c)
        IDX.append(idx)
    return (np.array(A), np.array(N), np.array(C), np.array(IDX))

def write_inp(path, nodes, elems, pairs, info, p, Tz, loads=True):
    (offs, off) = ({}, 0)
    order = [g for g in ('IMP', 'CORT', 'CAN') if g in elems]
    for g in order:
        offs[g] = off
        off += len(elems[g])
    (x, y, z) = nodes.T
    (W, Hb) = (p['W'], info['Hb'])
    sym = np.where(np.abs(y) < 1e-07)[0]
    bone_nodes = np.unique(np.concatenate([elems[g].ravel() for g in order if g != 'IMP']))
    bm = np.zeros(len(nodes), bool)
    bm[bone_nodes] = True
    if p['bc'] == 'bottom':
        fix = np.where(bm & (np.abs(z + Hb) < 1e-06))[0]
    elif p['bc'] == 'lateral_lower':
        side = (np.abs(np.abs(x) - W / 2) < 1e-06) | (np.abs(y - W / 2) < 1e-06)
        fix = np.where(bm & side & (z <= -Hb * (1 - p['bc_fr']) + 1e-09))[0]
    else:
        raise ValueError(p['bc'])
    apex = np.abs(nodes[:, 2] + p['L']) < 1e-07
    imask = np.zeros(len(nodes), bool)
    imask[pairs[:, 0]] = True
    imask &= ~apex
    bmask = np.zeros(len(nodes), bool)
    bmask[pairs[:, 1]] = True
    bmask &= ~apex
    s_imp = faces_on(elems['IMP'], imask)
    s_bone = {g: faces_on(elems[g], bmask) for g in order if g != 'IMP'}
    am_i = np.zeros(len(nodes), bool)
    am_i[pairs[:, 0]] = True
    am_i &= apex
    am_b = np.zeros(len(nodes), bool)
    am_b[pairs[:, 1]] = True
    am_b &= apex
    s_imp_apex = faces_on(elems['IMP'], am_i)
    s_bone_apex = {g: faces_on(elems[g], am_b) for g in order if g != 'IMP'}
    zt = p['e']
    tmask = np.abs(z - zt) < 1e-06
    s_top = faces_on(elems['IMP'], tmask)
    ref = len(nodes) + 1
    Kd = dict(p['K_pen']) if isinstance(p['K_pen'], dict) else {'CORT': p['kfac'] * p['E_c'] / p['h_if'], 'CAN': p['kfac'] * p['E_s'] / p['h_if']}
    lamd = dict(p['lam']) if isinstance(p['lam'], dict) else {g: p['lamfac'] * K / 10.0 for (g, K) in Kd.items()}
    (Kp, lam) = (Kd['CAN'], lamd['CAN'])
    imp_nodes = np.unique(elems['IMP'])
    with open(path, 'w') as f:
        f.write('*NODE, NSET=NALL\n')
        for (i, (a, b, c)) in enumerate(nodes, 1):
            f.write(f'{i},{a:.9e},{b:.9e},{c:.9e}\n')
        for g in order:
            f.write(f'*ELEMENT, TYPE=C3D10, ELSET={g}\n')
            for (k, t) in enumerate(elems[g], offs[g] + 1):
                f.write(f'{k},' + ','.join((str(int(i) + 1) for i in t)) + '\n')

        def nset(name, idx):
            f.write(f'*NSET, NSET={name}\n')
            idx = list(idx)
            for j in range(0, len(idx), 12):
                f.write(','.join((str(int(i) + 1) for i in idx[j:j + 12])) + '\n')

        def surf(name, lst, g):
            f.write(f'*SURFACE, NAME={name}, TYPE=ELEMENT\n')
            for (e_, fno) in lst:
                f.write(f'{offs[g] + e_ + 1},S{fno}\n')
        nset('NSYM', sym)
        nset('NFIX', fix)
        surf('S_IMP', s_imp, 'IMP')
        for (g, lst) in s_bone.items():
            if lst:
                surf('S_BONE_' + g, lst, g)
        surf('S_TOP', s_top, 'IMP')
        apex_on = p.get('apex_contact', True) and s_imp_apex and any(s_bone_apex.values())
        if apex_on:
            surf('S_IMP_APEX', s_imp_apex, 'IMP')
            f.write('*SURFACE, NAME=S_BONE_APEX, TYPE=ELEMENT\n')
            for (g, lst) in s_bone_apex.items():
                for (e_, fno) in lst:
                    f.write(f'{offs[g] + e_ + 1},S{fno}\n')
        f.write('*ORIENTATION, NAME=CYL, SYSTEM=CYLINDRICAL\n0.,0.,-100.,0.,0.,100.\n')
        G_ti = p['E_ti'] / (2 * (1 + p['nu_ti']))
        f.write(f"*MATERIAL, NAME=MTI\n*ELASTIC, TYPE=ENGINEERING CONSTANTS\n{p['E_ti']},{p['E_ti']},{p['E_ti']},{p['nu_ti']},{p['nu_ti']},{p['nu_ti']},{G_ti},{G_ti}\n{G_ti}\n*EXPANSION, TYPE=ORTHO\n{ALPHA},{ALPHA},0.\n")
        f.write(f"*MATERIAL, NAME=MCORT\n*ELASTIC\n{p['E_c']},{p['nu_c']}\n*EXPANSION\n0.\n")
        f.write(f"*MATERIAL, NAME=MCAN\n*ELASTIC\n{p['E_s']},{p['nu_s']}\n*EXPANSION\n0.\n")
        f.write('*SOLID SECTION, ELSET=IMP, MATERIAL=MTI, ORIENTATION=CYL\n')
        if 'CORT' in elems:
            f.write('*SOLID SECTION, ELSET=CORT, MATERIAL=MCORT\n')
        f.write('*SOLID SECTION, ELSET=CAN, MATERIAL=MCAN\n')
        for (g, lst) in s_bone.items():
            if not lst:
                continue
            f.write(f"*SURFACE INTERACTION, NAME=SI_{g}\n*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR\n{Kd[g]:.4e},{p['sig_inf']:.4e},{p['c0']:.4e}\n*FRICTION\n{p['mu']},{lamd[g]:.4e}\n")
            f.write(f'*CONTACT PAIR, INTERACTION=SI_{g}, TYPE=SURFACE TO SURFACE\nS_BONE_{g},S_IMP\n')
        if apex_on:
            f.write(f"*SURFACE INTERACTION, NAME=SI_APEX\n*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR\n{Kd['CAN']:.4e},{p['sig_inf']:.4e},{p['c0']:.4e}\n")
            f.write('*CONTACT PAIR, INTERACTION=SI_APEX, TYPE=SURFACE TO SURFACE\nS_BONE_APEX,S_IMP_APEX\n')
        f.write('*BOUNDARY\nNFIX,1,3,0.\nNSYM,2,2,0.\n')
        top_nodes = np.unique(np.concatenate([elems['IMP'][e_, list(FACES[fn])] for (e_, fn) in s_top]))
        nset('NTOP', top_nodes)
        f.write('*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL,0.\n')
        f.write(f'*STEP, INC=500, NLGEOM=NO\n*STATIC, SOLVER={SOLVER}\n0.5,1.0,1e-5,0.5\n')
        f.write('*TEMPERATURE\n')
        for i in imp_nodes:
            if Tz[i] != 0.0:
                f.write(f'{i + 1},{Tz[i]:.6e}\n')
        f.write('*NODE FILE\nU\n*END STEP\n')
        if loads:
            a = math.radians(p['angle'])
            Fh = p['F'] / 2.0
            f.write(f"*STEP, INC=500, NLGEOM=NO\n*STATIC, SOLVER={SOLVER}\n{p['inc0']:.4f},1.0,1e-5,{p['inc0']:.4f}\n")
            f.write('*CLOAD\n')
            wn = np.zeros(len(nodes))
            (At, Nt, Ct, It) = face_geom(nodes, elems['IMP'], s_top)
            for (a_, idx) in zip(At, It):
                wn[idx[3:]] += a_ / 3.0
            wn /= wn.sum()
            for i in np.where(wn > 0)[0]:
                f.write(f'{i + 1},1,{Fh * math.sin(a) * wn[i]:.8e}\n{i + 1},3,{-Fh * math.cos(a) * wn[i]:.8e}\n')
            f.write('*NODE FILE\nU\n*END STEP\n')
    top_nodes = np.unique(np.concatenate([elems['IMP'][e_, list(FACES[fn])] for (e_, fn) in s_top]))
    return dict(ref=ref, top=top_nodes, K_pen=Kp, lam=lam, Kd=Kd, lamd=lamd, n_fix=int(len(fix)), n_s_imp=len(s_imp), n_s_bone=int(sum((len(v) for v in s_bone.values()))), s_bone=s_bone, offs=offs, order=order)

def read_frd_disp(path, n_nodes):
    """The last DISP block in each step (frd writes every increment; step number from the 1PSTEP row)."""
    (blocks, steps) = ([], [])
    step = None
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if '1PSTEP' in lines[i]:
            step = int(lines[i].split()[3])
        if lines[i].startswith(' -4') and 'DISP' in lines[i]:
            steps.append(step)
            U = np.zeros((n_nodes, 3))
            i += 1
            while not lines[i].startswith(' -3'):
                ln = lines[i]
                if ln.startswith(' -1'):
                    nid = int(ln[3:13])
                    if nid <= n_nodes:
                        U[nid - 1] = [float(ln[13 + 12 * k:25 + 12 * k]) for k in range(3)]
                i += 1
            blocks.append(U)
        i += 1
    last = {}
    for (b, st) in zip(blocks, steps):
        last[st] = b
    return [last[k] for k in sorted(last, key=lambda x: (x is None, x))]

def run_ccx(job, wd):
    env = dict(os.environ)
    t = time.time()
    lp = os.path.join(wd, job + '.log')
    with open(lp, 'w') as fo:
        r = subprocess.run([CCX, '-i', job], cwd=wd, stdout=fo, stderr=subprocess.STDOUT, text=True, env=env)
    so = open(lp).read()
    ok = r.returncode == 0 and 'Job finished' in so
    if ok and (not os.environ.get('MM_KEEP')):
        open(lp, 'w').write(so[-20000:])
    return (ok, time.time() - t, so[-3000:])

def postprocess(nodes, elems, pairs, info, p, deck, U1, U2):
    ref = deck['ref']
    dU = U2 - U1
    nn = np.zeros((len(nodes), 3))
    na = np.zeros(len(nodes))
    lay = {}
    for (g, lst) in deck['s_bone'].items():
        (A, N, C, IDX) = face_geom(nodes, elems[g], lst)
        for (a, n, idx) in zip(A, N, IDX):
            for j in idx:
                nn[j] += a * n
                na[j] += a / 6.0
            for j in idx:
                lay[j] = g
    keep = np.abs(nodes[pairs[:, 1], 2] + p['L']) > 1e-07
    (ib, bb) = (pairs[keep, 0], pairs[keep, 1])
    n = nn[bb]
    n = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    rel = dU[bb] - dU[ib]
    rn = np.einsum('ij,ij->i', rel, n)
    rt = rel - rn[:, None] * n
    ut = np.linalg.norm(rt, axis=1) * 1000.0
    un = rn * 1000.0
    w = na[bb]
    zc = nodes[bb, 2]
    uL = dU[deck['top']].mean(0)
    out = dict(u_load_um=(uL * 1000.0).tolist())
    imp_n = np.unique(elems['IMP'])
    d = np.hypot(nodes[imp_n, 0], nodes[imp_n, 1]) + 10 * np.abs(nodes[imp_n, 2])
    ic = imp_n[np.argmin(d)]
    out['u_crest_um'] = (dU[ic] * 1000.0).tolist()
    a = math.radians(p['angle'])
    ld = np.array([math.sin(a), 0.0, -math.cos(a)])
    out['u_load_along_F_um'] = float(np.dot(uL, ld) * 1000.0)
    out['u_load_lat_um'] = float(uL[0] * 1000.0)
    out['u_crest_lat_um'] = float(dU[ic][0] * 1000.0)
    out['slip_max_um'] = float(ut.max())
    out['slip_mean_um'] = float((ut * w).sum() / w.sum())
    out['slip_p95_um'] = float(np.percentile(ut, 95))
    out['gap_open_max_um'] = float(un.max())
    out['rel_max_um'] = float(np.linalg.norm(rel, axis=1).max() * 1000.0)
    for g in set(lay.values()):
        m = np.array([lay.get(j) == g for j in bb])
        if m.any():
            out['slip_max_%s_um' % g] = float(ut[m].max())
            out['slip_mean_%s_um' % g] = float((ut[m] * w[m]).sum() / w[m].sum())
    return (out, ut, zc)

def read_contact_press(frd, n_nodes):
    """Last/ all CONTACT Block (CPRESS etc.) from .frd; returns list of dict name->(n_nodes,) for each step."""
    res = []
    with open(frd) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith(' -4') and 'CONTACT' in lines[i]:
            comps = []
            i += 1
            while lines[i].startswith(' -5'):
                comps.append(lines[i].split()[1])
                i += 1
            V = {c: np.zeros(n_nodes) for c in comps}
            while not lines[i].startswith(' -3'):
                ln = lines[i]
                if ln.startswith(' -1'):
                    nid = int(ln[3:13])
                    vals = [float(ln[13 + 12 * k:25 + 12 * k]) for k in range(len(comps))]
                    if nid <= n_nodes:
                        for (c, v) in zip(comps, vals):
                            V[c][nid - 1] = v
                i += 1
            res.append(V)
        i += 1
    return res

def radial_prestress(nodes, elems, pairs, deck, U1, p):
    """Radial mean pressure per layer in steps 1 = radial force per length/ ( π r_h) (half model), from the penalty law:
    pressure = K_pen · overlap in the normal direction (overlap = − (u_ben − u_imp ) · n, n from the leg towards the implant)."""
    out = {}
    for (g, lst) in deck['s_bone'].items():
        Kp = deck['Kd'][g]
        (A, N, C, IDX) = face_geom(nodes, elems[g], lst)
        b2i = dict(zip(pairs[:, 1].tolist(), pairs[:, 0].tolist()))
        Fr = 0.0
        for (a, n, idx) in zip(A, N, IDX):
            ov = []
            for j in idx[:3]:
                rel = U1[j] - U1[b2i[j]]
                ov.append(max(0.0, np.dot(rel, n)))
            pr = Kp * np.mean(ov)
            r_hat = np.array([C[0][0], C[0][1], 0.0])
            cc = nodes[idx[:3]].mean(0)
            rh = np.array([cc[0], cc[1], 0.0])
            rh /= max(np.linalg.norm(rh), 1e-12)
            Fr += pr * a * abs(np.dot(n, rh))
        zs = C[:, 2]
        Lz = zs.max() - zs.min() + 1e-09
        rr = np.hypot(C[:, 0], C[:, 1]).mean()
        out[g] = dict(Fr=Fr, p_rad=Fr / (math.pi * rr * Lz), Lz=Lz, r=rr)
    return out

def temps(nodes, elems, p, T_c, T_s):
    Tz = np.zeros(len(nodes))
    imp_n = np.unique(elems['IMP'])
    z = nodes[imp_n, 2]
    T = np.where(z >= 0.0, 0.0, np.where(z > -p['t_c'], T_c, T_s))
    Tz[imp_n] = T
    return Tz

def solve(p, wd, tag='job'):
    os.makedirs(wd, exist_ok=True)
    t0 = time.time()
    (nodes, elems, pairs, info) = build(p)
    tmesh = time.time() - t0
    rM = p['D'] / 2

    def T0(p0, E, nu):
        return p0 * (1 + nu) / E / ALPHA
    T_c = p['T_c'] if p['T_c'] is not None else T0(p['p0_c'], p['E_c'], p['nu_c']) if p['t_c'] > 0 else 0.0
    T_s = p['T_s'] if p['T_s'] is not None else T0(p['p0_s'], p['E_s'], p['nu_s'])
    rec = dict(spec=p, info=info, t_mesh=round(tmesh, 1), solver=SOLVER)
    calib = []
    T_c *= p['Tfac_c']
    T_s *= p['Tfac_s']
    do_cal = p['calib'] and p['T_c'] is None and (p['p0_c'] > 0 or p['p0_s'] > 0)
    has_c = p['t_c'] > 0
    T_pert = max(5.0, 5.0 * abs(T_c))
    runs = ([(0.0, T_s), (T_pert, T_s)] if has_c and p['p0_s'] > 0 and (p['p0_c'] > 0) else [(T_c, T_s)]) if do_cal else []
    res_cal = []
    for (it, (tc_, ts_)) in enumerate(runs):
        Tz = temps(nodes, elems, p, tc_, ts_)
        job = f'{tag}_cal{it}'
        deck = write_inp(os.path.join(wd, job + '.inp'), nodes, elems, pairs, info, p, Tz, loads=False)
        (ok, ts, _) = run_ccx(job, wd)
        if not ok:
            rec.update(ok=False, fail=f'calib {it}')
            return rec
        U = read_frd_disp(os.path.join(wd, job + '.frd'), len(nodes))
        pr = radial_prestress(nodes, elems, pairs, deck, U[-1], p)
        prd = {g: v['p_rad'] for (g, v) in pr.items()}
        res_cal.append(prd)
        calib.append(dict(T_c=tc_, T_s=ts_, pr=prd, t=round(ts, 1)))
        for ext in ('.frd', '.dat', '.sta', '.cvg', '.12d', '.inp'):
            fp = os.path.join(wd, job + ext)
            if os.path.exists(fp) and (not os.environ.get('MM_KEEP')):
                os.remove(fp)
    if len(res_cal) == 2:
        a_cs = res_cal[0]['CORT'] / T_s
        a_ss = res_cal[0]['CAN'] / T_s
        a_cc = (res_cal[1]['CORT'] - res_cal[0]['CORT']) / T_pert
        a_sc = (res_cal[1]['CAN'] - res_cal[0]['CAN']) / T_pert
        Amat = np.array([[a_cc, a_cs], [a_sc, a_ss]])
        (T_c, T_s) = [float(x) for x in np.linalg.solve(Amat, np.array([p['p0_c'], p['p0_s']]))]
        calib.append(dict(A=Amat.tolist()))
    elif len(res_cal) == 1:
        g = 'CAN' if 'CAN' in res_cal[0] else 'CORT'
        if (not has_c or p['p0_c'] <= 0) and res_cal[0].get('CAN', 0) > 0:
            T_s *= p['p0_s'] / res_cal[0]['CAN']
        elif has_c and res_cal[0].get('CORT', 0) > 0 and (p['p0_c'] > 0):
            T_c *= p['p0_c'] / res_cal[0]['CORT']
    rec['calib'] = calib
    Tz = temps(nodes, elems, p, T_c, T_s)
    deck = write_inp(os.path.join(wd, tag + '.inp'), nodes, elems, pairs, info, p, Tz)
    (ok, ts, so) = run_ccx(tag, wd)
    rec.update(ok=ok, t_solve=round(ts, 1), K_pen=deck['Kd'], lam=deck['lamd'], n_fix=deck['n_fix'], n_s_bone=deck['n_s_bone'], T_c=T_c, T_s=T_s)
    if not ok:
        rec['fail'] = so[-800:]
        return rec
    U = read_frd_disp(os.path.join(wd, tag + '.frd'), len(nodes))
    if len(U) < 2:
        rec.update(ok=False, fail='less than two DISP blocks')
        return rec
    pr = radial_prestress(nodes, elems, pairs, deck, U[0], p)
    rec['prestress_p_rad'] = {g: v['p_rad'] for (g, v) in pr.items()}
    (out, ut, zc) = postprocess(nodes, elems, pairs, info, p, deck, U[0], U[1])
    rec.update(out)
    bins = np.arange(0.0, p['L'] + 0.5, 0.5)
    prof = []
    for (a_, b_) in zip(bins[:-1], bins[1:]):
        m = (-zc >= a_) & (-zc < b_)
        prof.append(float(ut[m].max()) if m.any() else None)
    rec['slip_profile_max_um'] = prof
    for ext in ('.frd', '.dat', '.12d'):
        fp = os.path.join(wd, tag + ext)
        if os.path.exists(fp) and (not os.environ.get('MM_KEEP')):
            os.remove(fp)
    fp = os.path.join(wd, tag + '.inp')
    if os.path.exists(fp) and (not os.environ.get('MM_KEEP')):
        os.remove(fp)
    rec['t_total'] = round(time.time() - t0, 1)
    return rec

def main():
    spec = json.load(open(sys.argv[1]))
    (out, wd) = (sys.argv[2], sys.argv[3])
    specs = spec if isinstance(spec, list) else [spec]
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)['id'] for l in open(out) if l.strip()}
    for s in specs:
        if s['id'] in done:
            continue
        p = dict(DEFAULT)
        p.update(s)
        try:
            rec = solve(p, os.path.join(wd, s['id']), tag='mm')
        except Exception:
            import traceback
            rec = dict(ok=False, fail=traceback.format_exc()[-1500:])
        rec['id'] = s['id']
        with open(out, 'a') as f:
            f.write(json.dumps(rec, default=float) + '\n')
        print(s['id'], rec.get('ok'), rec.get('t_total'), rec.get('slip_max_um'), rec.get('u_load_lat_um'), flush=True)
if __name__ == '__main__':
    main()
