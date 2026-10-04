from dental_release.paths import expand as _release_expand
import os, sys, json, time, hashlib, resource
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', IGL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
import numpy as np
import scipy.linalg, scipy.optimize
sys.dont_write_bytecode = True
R = Path(__file__).resolve().parents[1]
B = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts'))
sys.path.insert(0, str(B / 'PROOF_LANE_FULL_CROWN_R5/code'))
from construct_a import normal_offset
from construct_b import fast_nearest, boundary_gap, compact
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', IGL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
from geometry import np, trimesh, loops, sample, closest
from threadpoolctl import threadpool_limits
threadpool_limits(1)
sys.path.insert(0, str(B / _release_expand('PROOF_LANE')))
from gencad_bench.checks.exact import triangle_distance
from fractions import Fraction as F

def read(p):
    return json.loads(Path(p).read_text())

def dump(p, x):
    Path(p).write_text(json.dumps(x, indent=2) + '\n')

def load(p):
    return dict(np.load(p))

def exact_pair(a, b):
    conv = lambda t: tuple((tuple((F(float(c)) for c in p)) for p in t))
    (d, x, y) = triangle_distance(conv(a), conv(b))
    return dict(distance_squared_mm2=str(d), outer_point=[[str(c) for c in x]][0], inner_point=[str(c) for c in y], fails_0p5mm=d < F(1, 4), triangles=[a.tolist(), b.tolist()], scope='Exact stored float coordinates; source measurement uncertainty absent')

def pair_witness(ext, inn):
    (q, d, j) = fast_nearest(inn, ext.mean(1))
    i = int(d.argmin())
    return dict(outer_triangle_index=i, inner_triangle_index=int(j[i]), centroid_distance_mm=float(d[i]), exact=exact_pair(ext[i], inn[j[i]]))

def stats(a):
    return dict(min=float(np.min(a)), median=float(np.median(a)), max=float(np.max(a)), p95=float(np.quantile(a, 0.95)))

def run():
    st = time.perf_counter()
    pr = read(R / 'PREREG_A.json')
    keys = pr['selection']['keys']
    r4 = {r['key']: r for r in read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_PREDICTIONS_C.json')['rows'] if r['method'] == 'exact_margin'}
    r5 = {r['key']: r for r in read(B / 'PROOF_LANE_FULL_CROWN_R5/FROZEN_PREDICTIONS_A.json')['rows'] if r['method'] == 'offset_only'}
    r5b = {r['key']: r for r in read(B / 'PROOF_LANE_FULL_CROWN_R5/RESULTS_B.json')['rows'] if r['method'] == 'distance_offset_only'}
    inp = {r['key']: r for r in read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']}
    rows = []
    for key in keys:
        m = load(r4[key]['mesh_path'])
        a = load(r5[key]['mesh_path'])
        p = load(inp[key]['public_path'])
        target = load(inp[key]['private_path'])['target']
        (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
        ext = v[f[roles == 0]]
        prep = v[f[roles == 1]]
        inner = a['vertices'][a['faces'][a['roles'] == 1]]
        (replay, offset) = normal_offset(v, f, roles)
        (pv, pf) = compact(v, f[roles == 1])
        rim = pv[loops(trimesh.Trimesh(pv, pf, process=False))[0]]
        nn = -np.cross(prep[:, 1] - prep[:, 0], prep[:, 2] - prep[:, 0])
        area = np.linalg.norm(nn, axis=1) / 2
        nn /= np.maximum(np.linalg.norm(nn, axis=1)[:, None], 1e-30)
        ids = np.unique(f[roles == 1])
        apex = ids[-1]
        delta = replay[apex] - v[apex]
        residual = nn @ delta - 0.05
        k = int(np.argmax(abs(residual)))
        pts = sample(inner, 4096)
        (q, d, j) = fast_nearest(prep, pts)
        desired = boundary_gap(q, rim)
        errors = abs(d - desired)
        worst = int(errors.argmax())
        (q2, d2, j2) = closest(prep, pts[::64])
        parity = float(np.max(abs(d2 - d[::64])))
        native = pair_witness(target, prep)
        base = pair_witness(ext, prep)
        new = pair_witness(ext, inner)
        native_d = fast_nearest(prep, sample(target, 8192))[1]
        r = dict(key=key, family=inp[key]['family'], resolution='PER_POINT', replay_max_coordinate_error_mm=float(abs(replay - a['vertices']).max()), source_preparation=dict(formula='rim_xy=c+0.65*(margin_xy-c); rim_z=margin_z; apex_xy=c; apex_z=max(native_z)-1.5', apex=v[apex].tolist(), rim_z_span_mm=float(np.ptp(rim[:, 2])), facets=len(prep), area_mm2=stats(area), tiny_facets_lt1e_minus12=int(np.sum(area < 1e-12)), normal_z=stats(nn[:, 2]), apex_ls_delta_mm=delta.tolist(), apex_ls_residual_mm=stats(abs(residual)), worst_apex_facet=k, worst_apex_area_mm2=float(area[k]), worst_apex_normal=nn[k].tolist()), inherited_offset=offset, gap=dict(actual_unsigned_mm=stats(d), desired_mm=stats(desired), max_error_mm=float(errors.max()), backend_parity_mm=parity, worst_probe=pts[worst].tolist(), nearest_preparation=q[worst].tolist(), nearest_triangle=int(j[worst]), desired_at_worst_mm=float(desired[worst]), actual_at_worst_mm=float(d[worst]), sample_count=len(pts)), wall_R4=base, wall_R5_A=new, native_to_preparation=native, native_clearance_mm=stats(native_d), R5_B_record=dict(gap=r5b[key]['gap_unsigned_sampled_mm'], offset=r5b[key]['offset'], gates=r5b[key]['gates']), checks=dict(replay_exact=bool(np.array_equal(replay, a['vertices'])), backends_agree=parity <= 1e-08, R4_actual_wall_failure=base['exact']['fails_0p5mm'], R5_actual_wall_failure=new['exact']['fails_0p5mm']))
        rows.append(r)
        dump(R / 'raw/A_TRACE.json', dict(rows=rows))
        print(key, 'wall', base['centroid_distance_mm'], 'native', native['centroid_distance_mm'], 'gap maxerr', errors.max(), 'apexres', abs(residual).max(), flush=True)
    suff = dict(summary_mean_mm=[0.75, 0.75], identity_error_mm=0.0, state_A=[0.625, 0.875], state_B=[0.375, 1.125], downstream_min_mm=[0.625, 0.375], downstream_difference_mm=0.25, decisions=[True, False], smallest_extension='Minimum for this fixed threshold; a local field for design edits', kind='our_own_fixture', resolution=dict(summary='PER_TOOTH', downstream='PER_SURFACE_REGION'))
    out = dict(claim_type='capability', rows=rows, sufficiency=suff, external_referent=pr['external_referent'], cost=dict(seconds=time.perf_counter() - st, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    dump(R / 'RESULTS_A.json', out)
    dump(R / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-crown-root-cause', phase='A_DIAGNOSED', latest_gate=[r['checks'] for r in rows], next_operation='Change preparation: jointly nested solids with exact film; preserve original failure and external shape checks', review_state='PENDING_INDEPENDENT_REVIEW'))
if __name__ == '__main__':
    run()
