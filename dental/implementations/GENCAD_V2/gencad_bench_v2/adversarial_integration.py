"""Supplemental controls on native task artifacts and full-cap representations."""
import tempfile, shutil, sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from .common import *
from .tasks import load_task
from .axial_round3 import cone, cap_planes
sys.path.insert(0, str(ROOT / 'gencad_bench_v2/vendor'))
from continuous import exact_batch, plane

def lp_gap(u, l):
    A = []
    b = []
    for t in [u, l]:
        xy = t[:, :2]
        orientation = np.sign(np.cross(xy[1] - xy[0], xy[2] - xy[0]))
        for (p, q) in zip(xy, np.roll(xy, -1, axis=0)):
            edge = q - p
            n = orientation * np.array([edge[1], -edge[0]])
            A.append(n)
            b.append(n @ p)
    c = np.linalg.solve(np.c_[u[:, :2], np.ones(3)], u[:, 2]) - np.linalg.solve(np.c_[l[:, :2], np.ones(3)], l[:, 2])
    r = linprog(c[:2], A_ub=A, b_ub=b, bounds=[(None, None)] * 2, method='highs', options={'threads': 1})
    return r.fun + c[2] if r.success else None

def run():
    controls = []

    def check(name, positive, rejected, **kw):
        row = dict(name=name, positive=bool(positive), injected_fault_rejected=bool(rejected), **kw)
        controls.append(row)
        if not positive or not rejected:
            raise ValueError('integration control failed ' + name)
    (v, f) = cone(2.0, 1.8, 0.0, 3.0)
    good = bool(cap_planes(v, f))
    (u, g) = cone(1.8, 2.0, 0.0, 3.0)
    bad = False
    try:
        cap_planes(u, g)
    except ValueError:
        bad = True
    check('conical insertion normals', good, bad)
    rng = np.random.default_rng(20261002)
    errs = []
    for _ in range(48):
        base = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]]) + rng.normal(0, 0.12, (3, 2))
        other = np.array([[0.1, 0.1], [1.8, 0.2], [0.2, 1.7]]) + rng.normal(0, 0.1, (3, 2))
        u = np.c_[base, rng.uniform(1, 2, 3)]
        l = np.c_[other, rng.uniform(0, 1, 3)]
        a = lp_gap(u, l)
        if a is None:
            continue
        e = exact_batch(u[None], l[None], plane(u[None]), plane(l[None]), np.array([[0, 0]]))[0]
        errs.append(abs(a - e))
    check('continuous triangle intersections vs independent halfplane LP', bool(errs) and max(errs) < 1e-07, all((abs(e - 0.2) > 1e-07 for e in errs)), pairs=len(errs), max_error_mm=max(errs), external_referent=dict(kind='published_code', locator='https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html', compared_quantity='minimum affine gap on a projected intersection polygon', refutes_us=True))
    p = next((p for p in sorted((DATA / 'public').glob('*.json')) if read(p)['status'] == 'READY'))
    original = load_task(p)
    with tempfile.TemporaryDirectory(dir=DATA) as td:
        dest = Path(td)
        (dest / 'scenes').mkdir()
        shutil.copy2(p, dest / p.name)
        gfile = DATA / 'public/scenes' / original['geometry_file']
        shutil.copy2(gfile, dest / 'scenes' / gfile.name)
        okay = load_task(dest / p.name)['status'] == 'READY'
        with (dest / 'scenes' / gfile.name).open('ab') as out:
            out.write(b'INJECTED_POSE_OR_GEOMETRY_CHANGE')
        caught = False
        try:
            load_task(dest / p.name)
        except ValueError:
            caught = True
        check('actual public geometry mutation', okay, caught)
        task = read(dest / p.name)
        task['geometry_file'] = '../../private/secret.npz'
        dump(dest / p.name, task)
        caught = False
        try:
            load_task(dest / p.name)
        except ValueError:
            caught = True
        check('input path traversal', True, caught)
    import trimesh
    p = next((DATA / 'axial_R3').glob('*axial_flare.npz'))
    a = np.load(p)
    mesh = trimesh.Trimesh(a['vertices'], a['faces'], process=False)
    bad = trimesh.Trimesh(a['vertices'], a['faces'][1:], process=False)
    check('actual complete crown missing facet', mesh.is_watertight, not bad.is_watertight)
    from gencad_bench.checks import milling
    pl = [[0, 0, 1, 0]]
    good = milling.verify_positive(pl, [[0, 0, 0]], 0.5, 0, [0, 0, -1], [[0, 0, -0.5]])
    bad = milling.verify_positive(pl, [[0, 0, 0]], 0.5, 0, [0, 0, -1], [[0, 0, 0]])
    check('ideal milling ball intrudes into solid', good, not bad)
    cohort = read(ROOT / 'data/COHORT_LOCK.json')['payload']['cases']
    keys = {s: {r['case_key'] for r in cohort if r['split'] == s} for s in ['train', 'dev', 'test']}
    good = not (keys['train'] & keys['test'] or keys['train'] & keys['dev'] or keys['dev'] & keys['test'])
    badkeys = keys['test'] | {next(iter(keys['train']))}
    check('case split membership', good, bool(keys['train'] & badkeys), latent_cross_dataset_identity='UNKNOWN')
    dump(ROOT / 'raw/INTEGRATION_CONTROLS.json', dict(checks=controls, count=len(controls), all_passed=True, scope='Supplemental executed software/mathematical controls; not external physical validation'))
    return controls
if __name__ == '__main__':
    print(run())
