from dental_release.paths import expand as _release_expand
import os, sys, json, hashlib, time, datetime
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT.parent / 'PROOF_LANE_GENCAD_V3'
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X42_gen_beat_standard'))
os.environ.setdefault('PYTHONDONTWRITEBYTECODE', '1')
sys.dont_write_bytecode = True
sys.path.insert(0, str(BENCH / 'code'))
from task_io import attach
from legacy.generators import submit, population, parametric, crown_loop, antagonist, optimize as original_lp
from fast_geometry import optimize as standard_lp, feasibility
from legacy.checks import all_checks
from legacy.geometry import section_weights
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, eye, hstack, vstack

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, allow_nan=False, default=lambda a: a.tolist() if hasattr(a, 'tolist') else str(a)) + '\n')

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def split():
    return json.loads((ROOT / 'SPLITS.json').read_text())

def state(phase, gate, next_operation, **kw):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X42-gen-beat-standard', claim_type='algorithm', phase=phase, latest_gate=gate, next_operation=next_operation, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW', **kw))

def tasks(key, public=None):
    return json.loads((Path(public or BENCH / 'payload/public') / 'tasks' / (key + '.json')).read_text())

def scene(t, public=None):
    with np.load(Path(public or BENCH / 'payload/public') / t['geometry_file'], allow_pickle=False) as a:
        return attach(t, dict(a))

def dev_reference(key):
    if key not in split()['dev_all']:
        raise ValueError('Only dev references permitted in development')
    with np.load(BENCH / 'payload/private/references' / (key + '.npz'), allow_pickle=False) as a:
        return dict(a)

def metric(t, z, ref):
    valid = np.isfinite(ref)
    w = t['weights']
    den = w[valid].sum()
    if den <= 0:
        return {'rmse_mm': None, 'contact_error_mm2': None, 'coverage': 0}
    ok = valid & np.isfinite(t['ceiling'])
    cg = t['ceiling'] - ref
    pg = t['ceiling'] - z
    gold = ok & (cg >= 0) & (cg <= 0.1)
    pred = ok & (pg >= 0) & (pg <= 0.1)
    return {'rmse_mm': float(np.sqrt(np.sum(w[valid] * (z[valid] - ref[valid]) ** 2) / den)), 'contact_error_mm2': float(abs(w[pred].sum() - w[gold].sum())) if ok.any() else None, 'coverage': float(den / w.sum()), 'gold_contact_mm2': float(w[gold].sum()), 'pred_contact_mm2': float(w[pred].sum())}

def structural_refusal(t):
    r = t['requirements']
    if t['family'] == 'implant_crown':
        a = np.array(t['implant_axis'])
        ang = np.degrees(np.arccos(a[2] / np.linalg.norm(a)))
        if ang > r['channel_max_deg'] + 1e-06:
            return 'fixed implant channel'
    if t['family'] == 'lattice_onlay' and t['minimum_strut_mm'] < r['strut_mm'] - 1e-06:
        return 'fixed strut width'
    return None

def project(t, target):
    if structural_refusal(t):
        return {'status': 'ABSTAIN', 'reason': structural_refusal(t), 'task_id': t['task_id']}
    n = len(target)
    r = t['requirements']
    inner = t['preparation_z'] + r['film_min_mm']
    lower = inner + r['wall_mm']
    z = np.maximum(target, lower)
    A = t['A']
    b = t['obstacle_b']
    rows = []
    rhs = []
    if A.shape[0] and np.min(b - A @ lower) < -1e-06:
        return {'status': 'ABSTAIN', 'reason': 'wall-film obstacle obstruction', 'task_id': t['task_id']}
    feasible = not A.shape[0] or np.max(A @ z - b) <= 0
    if t['family'] == 'bridge3':
        for x in t['connector_x']:
            w = section_weights(t['xy'], t['faces'], x)
            rows.append(csr_matrix(np.r_[-w, np.zeros(n)][None, :]))
            rhs.append(np.array([-r['connector_mm2'] - w @ inner]))
            feasible &= w @ (z - inner) >= r['connector_mm2']
    if feasible:
        return submit(t, z, inner, {'projection': 'separable lower-bound certificate'})
    I = eye(n, format='csr')
    rows = [hstack([A, csr_matrix((A.shape[0], n))]), hstack([I, -I]), hstack([-I, -I])] + rows
    rhs = [b, target, -target] + rhs
    sol = linprog(np.r_[np.zeros(n), t['weights']], A_ub=vstack(rows), b_ub=np.concatenate(rhs), bounds=[(v, None) for v in lower] + [(0, None)] * n, method='highs', options={'threads': 1})
    if not sol.success:
        return {'status': 'ABSTAIN', 'reason': 'LP ' + str(sol.status), 'task_id': t['task_id']}
    return submit(t, sol.x[:n], inner, {'projection': 'exact constrained target LP', 'iterations': int(sol.nit)})

def contact_r1(t):
    z = np.array(t['prior'])
    ok = np.isfinite(t['ceiling'])
    z[ok] += 0.2 * (t['ceiling'][ok] - 0.06 - z[ok])
    return project(t, z)

def summarize(rows):
    import collections
    out = {}
    for name in sorted({r['method'] for r in rows}):
        rr = [r for r in rows if r['method'] == name]
        passed = [r for r in rr if r['verdict'] == 'PASS']
        out[name] = {'counts': dict(collections.Counter((r['verdict'] for r in rr))), 'rows': len(rr), 'median_height_rmse_mm': float(np.median([r['rmse_mm'] for r in passed if r.get('rmse_mm') is not None])) if passed else None, 'median_contact_error_mm2': float(np.median([r['contact_error_mm2'] for r in passed if r.get('contact_error_mm2') is not None])) if passed else None}
    return out
