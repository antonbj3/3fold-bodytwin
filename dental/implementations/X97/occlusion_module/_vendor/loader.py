"""Standalone loader for the `contact-scene-v1` contact-problem package.

Each numbered stem is one contact problem:

  scenes/NN.json   compact description (bodies/contacts or a supplied operator)
  scenes/NN.hdf5   the same operator in a compact FCLIB local-problem container

The loader reconstructs, from the JSON alone,

    W = G = J M^-1 J^T + diag(eta),    q = b,    mu

and can compare it element-by-element with the HDF5 reference.  It uses only
NumPy and h5py; it imports nothing from any project build tree and needs no
generator code.

Scene JSON layout (`schema = "contact-scene-v1"`):

  schema, index, scene, dt, n_b, n_c, dof, mu[], b[], eta[]
  model : "body_contacts" | "operator"
  provenance : public origin/construction text
  bodies   : [ {id, dof, mass, inertia, inv_mass, inv_inertia} ]
  contacts : [ {id, a, b, normal[3], t1[3], t2[3], lever_a[3], lever_b[3], mu} ]
  operator : { J_b64, J_shape, minv_diag | minv_b64, minv_shape }   (model == operator)

  body_contacts: J = per-contact blocks [D ; r x D] (rows n,t1,t2; cols v,omega),
                 M^-1 = block diagonal of per-body inverse mass / inertia.
  operator:     J and M^-1 given directly (dense random or floating-base multibody).

HDF5 (compact FCLIB local-problem layout):

  /fclib_local/W/{x,i,p,n,m,nz,nzmax} with attribute symmetric_packed=1:
      lower triangle incl. diagonal in CSC; full W = L + L^T - diag(L).
  /fclib_local/vectors/{q,mu}, /fclib_local/info/...

  This is a compact symmetric-packed container derived from the operator, not an
  untouched standard full-matrix FCLIB export; see README.md.

Command line (run from the package root):

    python loader.py            # rebuild and verify all 21 scenes
    python loader.py 0 1 12     # only the listed scene indices
    python loader.py --list     # print the scene index/name table
"""
from __future__ import annotations
import base64
import json
import os
import numpy as np
import h5py
HERE = os.path.dirname(os.path.abspath(__file__))
SCENES = os.path.join(HERE, 'scenes')
N_SCENES = 21

def b64_f8(a):
    return base64.b64encode(np.ascontiguousarray(a, dtype='<f8').tobytes()).decode()

def b64_i4(a):
    return base64.b64encode(np.ascontiguousarray(a, dtype='<i4').tobytes()).decode()

def unb64_f8(s, shape):
    return np.frombuffer(base64.b64decode(s), dtype='<f8').reshape(shape).astype(np.float64)

def unb64_i4(s, shape):
    return np.frombuffer(base64.b64decode(s), dtype='<i4').reshape(shape).astype(np.int64)

def read_hdf5_problem(path):
    """Return (G, b, mu) from a compact symmetric-packed FCLIB hdf5."""
    with h5py.File(path, 'r') as f:
        loc = f['fclib_local']
        W = loc['W']
        n = int(np.asarray(W['n']).ravel()[0])
        m = int(np.asarray(W['m']).ravel()[0])
        nz = int(np.asarray(W['nz']).ravel()[0])
        if nz != -1:
            raise ValueError(f'{path}: expected CSC nz=-1, got {nz}')
        x = np.asarray(W['x'], float)
        i = np.asarray(W['i'], np.int64)
        p = np.asarray(W['p'], np.int64)
        sym = int(W.attrs.get('symmetric_packed', 0))
        n_cols = int(p.shape[0]) - 1
        L = np.zeros((n, m))
        for c in range(n_cols):
            for k in range(p[c], p[c + 1]):
                L[i[k], c] = x[k]
        if sym:
            G = L + L.T - np.diag(np.diag(L))
        else:
            G = L
        q = np.asarray(loc['vectors']['q'], float)
        mu = np.asarray(loc['vectors']['mu'], float)
    return (G, q, mu)

def rebuild(scene):
    """Reconstruct (G, b, mu, J, Minv) from a scene dict."""
    model = scene['model']
    nc = int(scene['n_c'])
    n3 = 3 * nc
    if model == 'body_contacts':
        bodies = scene['bodies']
        contacts = scene['contacts']
        dof = int(bodies[0]['dof'])
        nb = len(bodies)
        nv = nb * dof
        Minv = np.zeros((nv, nv))
        for bd in bodies:
            b = int(bd['id'])
            if 'minv_block' in bd:
                blk = np.asarray(bd['minv_block'], float).reshape(bd.get('block_dim', dof), -1)
                Minv[dof * b:dof * b + dof, dof * b:dof * b + dof] = blk
            elif dof == 3:
                Minv[3 * b, 3 * b] = float(bd['inv_mass'])
                Minv[3 * b + 1, 3 * b + 1] = float(bd['inv_mass'])
                Minv[3 * b + 2, 3 * b + 2] = float(bd['inv_mass'])
            else:
                im = float(bd['inv_mass'])
                ii = bd['inv_inertia']
                Minv[dof * b + 0, dof * b + 0] = im
                Minv[dof * b + 1, dof * b + 1] = im
                Minv[dof * b + 2, dof * b + 2] = im
                Minv[dof * b + 3, dof * b + 3] = ii[0]
                Minv[dof * b + 4, dof * b + 4] = ii[1]
                Minv[dof * b + 5, dof * b + 5] = ii[2]
        J = np.zeros((n3, nv))
        for ct in contacts:
            c = int(ct['id'])
            D = np.vstack([ct['normal'], ct['t1'], ct['t2']]).astype(float)
            (ia, ib) = (int(ct['a']), int(ct['b']))
            if dof == 3:
                if ia >= 0:
                    J[3 * c:3 * c + 3, 3 * ia:3 * ia + 3] = D
                if ib >= 0:
                    J[3 * c:3 * c + 3, 3 * ib:3 * ib + 3] = -D
            else:
                if ia >= 0:
                    la = np.asarray(ct['lever_a'], float)
                    J[3 * c:3 * c + 3, 6 * ia:6 * ia + 3] = D
                    J[3 * c:3 * c + 3, 6 * ia + 3:6 * ia + 6] = np.cross(la[None, :], D)
                if ib >= 0:
                    lb = np.asarray(ct['lever_b'], float)
                    J[3 * c:3 * c + 3, 6 * ib:6 * ib + 3] = -D
                    J[3 * c:3 * c + 3, 6 * ib + 3:6 * ib + 6] = -np.cross(lb[None, :], D)
    elif model == 'operator':
        op = scene['operator']
        J = unb64_f8(op['J_b64'], tuple(op['J_shape'])) if 'J_b64' in op else np.asarray(op['J'], float)
        nv = J.shape[1]
        if 'minv_diag' in op:
            Minv = np.diag(np.asarray(op['minv_diag'], float))
        elif 'minv_b64' in op:
            Minv = unb64_f8(op['minv_b64'], tuple(op['minv_shape']))
        else:
            Minv = np.asarray(op['Minv'], float)
    else:
        raise ValueError(f'unknown model {model!r}')
    G = J @ Minv @ J.T
    eta = np.asarray(scene['eta'], float)
    if np.any(eta):
        G = G + np.diag(np.repeat(eta, 3))
    b = np.asarray(scene['b'], float)
    mu = np.asarray(scene['mu'], float)
    return (G, b, mu, J, Minv)

def load_scene(index):
    with open(os.path.join(SCENES, f'{index:02d}.json')) as f:
        return json.load(f)

def verify(json_path, hdf5_path):
    with open(json_path) as f:
        scene = json.load(f)
    (Gj, bj, muj, J, Minv) = rebuild(scene)
    (Gh, bh, muh) = read_hdf5_problem(hdf5_path)
    scale = float(np.abs(Gh).max()) or 1.0
    nrm = float(np.linalg.norm(Gh)) or 1.0
    rel = float(np.abs(Gj - Gh).max() / scale)
    rel_fro = float(np.linalg.norm(Gj - Gh) / nrm)
    absd = float(np.abs(Gj - Gh).max())
    bd = float(np.abs(bj - bh).max()) if bh.shape == bj.shape else float('nan')
    mdm = float(np.abs(muj - muh).max()) if muh.shape == muj.shape else float('nan')
    return dict(scene=scene['scene'], index=scene['index'], n_c=int(scene['n_c']), n_b=int(scene['n_b']), model=scene['model'], G_rel_fro=rel_fro, G_rel=rel, G_abs=absd, b_abs=bd, mu_abs=mdm, G_scale=scale, nv=int(J.shape[1]))

def _main(argv):
    if '--list' in argv:
        manifest = json.load(open(os.path.join(SCENES, 'manifest.json')))
        for key in sorted(manifest['scenes']):
            e = manifest['scenes'][key]
            print(f"{int(key):2d} {e['scene']:14s} n_c={e['n_c']:5d} model={e['model']}")
        return 0
    idx = [int(x) for x in argv if not x.startswith('-')] or list(range(N_SCENES))
    ok = True
    for i in idx:
        r = verify(os.path.join(SCENES, f'{i:02d}.json'), os.path.join(SCENES, f'{i:02d}.hdf5'))
        flag = 'OK ' if r['G_rel_fro'] <= 1e-15 and r['b_abs'] <= 1e-15 and (r['mu_abs'] == 0.0) else 'FAIL'
        ok = ok and flag == 'OK '
        print(f"{flag} {i:2d} {r['scene']:14s} n_c={r['n_c']:5d} model={r['model']:16s} G_rel_fro={r['G_rel_fro']:.3e} G_rel_max={r['G_rel']:.3e} G_abs={r['G_abs']:.3e} b_abs={r['b_abs']:.3e} mu_abs={r['mu_abs']:.3e}")
    print('ALL_OK' if ok else 'SOME_FAILED')
    return 0 if ok else 1
if __name__ == '__main__':
    import sys
    raise SystemExit(_main(sys.argv[1:]))
