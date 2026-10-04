"""Same-patient source-facet measurement; fixed axial projection, physical pose unknown."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import numpy as np, json, hashlib, time, resource, datetime
from scipy.optimize import linprog
from geometry_kernel import affine, boxes, buckets, extrema
P = Path(__file__).resolve().parents[1]
BASE = P.parent
DATA = P / 'local_inputs/Demo_1'
DATA_ORIGIN = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/dental_sol_night/X69/hao_demo/Demo/Demo_1'))
MAP = P / 'local_inputs/R2_RESULTS.json'
MAP_ORIGIN = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X69/raw/R2_RESULTS.json'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, v):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(v, indent=2, ensure_ascii=False, allow_nan=False, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + '\n')

def stl(p):
    b = p.read_bytes()
    n = int.from_bytes(b[80:84], 'little')
    if len(b) != 84 + 50 * n:
        raise ValueError('STL format')
    dt = np.dtype([('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attribute', '<u2')])
    return np.frombuffer(b, dt, count=n, offset=84)['vertices'].astype(float)

def lp_pair(U, L):
    c = np.r_[U[:, 2], -L[:, 2]]
    A = np.zeros((4, 6))
    A[0, :3] = 1
    A[1, 3:] = 1
    A[2, :3] = U[:, 0]
    A[2, 3:] = -L[:, 0]
    A[3, :3] = U[:, 1]
    A[3, 3:] = -L[:, 1]
    return linprog(c, A_eq=A, b_eq=[1, 1, 0, 0], bounds=(0, None), method='highs', options={'threads': 4})

def main():
    tick = time.perf_counter()
    source = json.load(open(MAP))
    arches = []
    manifest = []
    drop = []
    for j in source['jaws']:
        Ts = []
        labs = []
        faces = []
        T = np.array(j['source_to_target'])
        folder = 'L' if j['jaw'] == 'lower' else 'U'
        for p in sorted((DATA / f'IOS/teeth_seg/{folder}').glob('*.stl')):
            fdi = int(p.stem)
            tr = stl(p)
            tr = tr @ T[:3, :3].T + T[:3, 3]
            (co, valid) = affine(tr)
            Ts.append(tr[valid])
            labs.extend([fdi] * valid.sum())
            faces.extend(np.flatnonzero(valid).tolist())
            manifest.append(dict(path=str(DATA_ORIGIN / p.relative_to(DATA)), sha256=sha(p), fdi=fdi, triangles=len(tr), jaw=j['jaw']))
            drop.append(dict(fdi=fdi, faces=len(tr), rejected_degenerate_axial_projection=int((~valid).sum()), rejected_fraction=float((~valid).mean()), reason='No single-valued z over projected triangle'))
        tr = np.ascontiguousarray(np.concatenate(Ts))
        (co, _) = affine(tr)
        arches.append((j['jaw'], tr, np.array(labs, np.int32), np.array(faces, np.int32), co))
    a = {j: (tr, la, fa, co) for (j, tr, la, fa, co) in arches}
    (U, ul, uf, uc) = a['upper']
    (L, ll, lf, lc) = a['lower']
    ub = boxes(U)
    lb = boxes(L)
    origin = np.minimum(U[:, :, :2].min((0, 1)), L[:, :, :2].min((0, 1))) - 0.8
    hi = np.maximum(U[:, :, :2].max((0, 1)), L[:, :, :2].max((0, 1)))
    shape = np.ceil((hi - origin) / 0.8).astype(np.int64) + 2
    (ptr, ids) = buckets(lb, origin, 0.8, shape)
    (best, wi, xy, checked, tested) = extrema(U, L, ul, ll, uc, lc, ub, lb, origin, 0.8, shape, ptr, ids, True)
    rows = []
    controls = []
    for (u, l) in np.argwhere(np.isfinite(best)):
        (i, j) = wi[u, l]
        lp = lp_pair(U[i], L[j])
        err = abs(float(lp.fun) - best[u, l]) if lp.success else None
        controls.append(dict(upper_fdi=int(u), lower_fdi=int(l), error_mm=err, pass_gate=err is not None and err <= 1e-07, injected_001mm_rejected=lp.success and abs(float(lp.fun) - best[u, l] - 0.001) > 1e-07))
        point = xy[u, l]
        up = np.r_[point, uc[i, :2] @ point + uc[i, 2]]
        lo = np.r_[point, lc[j, :2] @ point + lc[j, 2]]
        rows.append(dict(upper_fdi=int(u), lower_fdi=int(l), minimum_axial_gap_mm=float(best[u, l]), upper_source_face=int(uf[i]), lower_source_face=int(lf[j]), upper_point_cbct_mm=up.tolist(), lower_point_cbct_mm=lo.tolist(), upper_triangle_cbct_mm=U[i].tolist(), lower_triangle_cbct_mm=L[j].tolist(), resolution='PER_POINT', timescale='SIMULTANEOUS'))
    ix = np.flatnonzero(ll == 36)
    Ls = L[ix]
    ls = ll[ix]
    lcs = lc[ix]
    lbs = lb[ix]
    (ptr2, ids2) = buckets(lbs, origin, 0.8, shape)
    (ex, _, _, _, _) = extrema(U, Ls, ul, ls, uc, lcs, ub, lbs, origin, 0.8, shape, ptr2, ids2, False)
    exerr = float(np.max(np.abs(ex[np.isfinite(ex)] - best[np.isfinite(ex)])))
    g0 = min((r['minimum_axial_gap_mm'] for r in rows))
    target = min((r['minimum_axial_gap_mm'] for r in rows if r['lower_fdi'] == 36))
    other = min((r['minimum_axial_gap_mm'] for r in rows if r['lower_fdi'] != 36))
    edits = []
    for h in [-0.05, 0, 0.05]:
        first = min(other, target - h)
        edits.append(dict(crown_fdi=36, height_edit_mm=h, first_touch_translation_mm=first, target_gap_after_virtual_first_touch_mm=target - h - first, target_becomes_first_touch=target - h <= other, force_decision='UNKNOWN', physical_height_decision='UNKNOWN_MISSING_LOADED_BITE_POSE_AND_DIRECTED_ERRORS', resolution='PER_TOOTH'))
    e36 = next((r for r in rows if r['lower_fdi'] == 36))
    elsewhere = next((r for r in rows if r['lower_fdi'] != 36 and r['upper_fdi'] != e36['upper_fdi']))
    suff = dict(kind='Identifiability construction on real patient incidence; not measured patient states', summary_total_force_N_A=100.0, summary_total_force_N_B=100.0, identity_error_N=0.0, tooth36_N_A=100.0, tooth36_N_B=0.0, downstream_difference_N=100.0, state_A_edge=[e36['upper_fdi'], 36], state_B_edge=[elsewhere['upper_fdi'], elsewhere['lower_fdi']], minimal_extension='Calibrated target-tooth force channel with registered same-patient contact allocation; total alone cannot identify per-tooth forces', note='Nonnegative balanced contact feasibility only; support/pose can differ between worlds, not a prediction at fixed measured pose', resolution='PER_TOOTH')
    force = []
    for m in manifest:
        force.append(dict(fdi=m['fdi'], jaw=m['jaw'], patient_force_interval_N=[0, None], upper_endpoint_meaning='UNBOUNDED_WITHOUT_MEASURED_TOTAL', scenario_force_interval_N=[0, 1963], scenario_total_N=[205, 1963], scenario_quantity='population-range scenario, not patient bound or confidence interval', physical_status='UNKNOWN', resolution='PER_TOOTH'))
    write(P / 'raw/SOURCE_MANIFEST.json', dict(dataset='LiuHao2023Demo1', locator='https://doi.org/10.5281/zenodo.8027553', license='CC BY4.0 per stored Zenodo metadata', stls=manifest, map_source_path=str(MAP_ORIGIN), map_source_sha256=sha(MAP)))
    out = dict(round='R1', claim_type='capability', patient_id='LiuHao2023Demo1', model='UNKNOWN_PHYSICAL;FIXED_AXIAL_PROJECTED_GEOMETRY', candidate_acquisition_pose='UNKNOWN_OPEN_OR_CLOSED', metric_scale='mm conditional STL interpretation', source_crowns=len(manifest), missing_expected_FDI=[24], missing_fdi_interpretation='No source mesh, not verified anatomical absence', source_triangles=sum((x['triangles'] for x in manifest)), dropout=drop, pairs=rows, possible_projected_pairs=len(rows), all_upper_lower_pairs=len(set(ul)) * len(set(ll)), pair_projection_dropout_fraction=1 - len(rows) / (len(set(ul)) * len(set(ll))), pair_dropout_reason='No projected triangle overlap, not tooth absence or zero load in another pose', source_pose_minimum_gap_mm=g0, source_pose_interpenetrates_projected_model=g0 < 0, virtual_first_touch_translation_mm=g0, target36_first_touch_gap_mm=target - g0, height_for_target36_to_become_first_touch_mm=target - other, edits=edits, forces=force, sufficiency=suff, controls=dict(winner_LP=controls, unpruned_target36_error_mm=exerr, all_pass=all((c['pass_gate'] and c['injected_001mm_rejected'] for c in controls)) and exerr <= 1e-07), physical_gate='UNKNOWN_NO_SAME_PATIENT_CLOSED_BITE_FORCE_MEASUREMENT', rigorous_enclosure='Exact affine polygon operation under fixed digital coordinates; formal float enclosure MISSING; physical directional pose/scale error UNKNOWN', cost=dict(wall_s=time.perf_counter() - tick, maxrss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, upstream_fit='UNKNOWN_REUSED', discovery='UNKNOWN', physical_validation='NOT_RUN', fallback='UNKNOWN'))
    write(P / 'rounds/R1/results.json', out)
    write(P / 'exports/PATIENT_LOAD_CASE_R1.json', dict(schema='patient-bound-load-v1', patient_id=out['patient_id'], claim_type='capability', status=out['physical_gate'], contact_pose='VIRTUAL_AXIAL_FIRST_TOUCH_NOT_PATIENT_MEASUREMENT', contact_witnesses=rows, tooth_force_intervals=force, crown_decisions=edits, unknown=['loaded_bite_pose', 'metric_scale', 'directed_cross_jaw_error', 'absolute_force', 'force_distribution', 'support_and_preload', 'tangential_force', 'time_history']))
    frozen = P / 'FROZEN_PREDICTIONS.json'
    if not frozen.exists():
        write(frozen, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), patient_id=out['patient_id'], prereg_sha256=sha(P / 'PREREG_R1.json'), result_sha256=sha(P / 'rounds/R1/results.json'), predictions=dict(target36_first_touch_gap_mm=target - g0, height_edits=edits), physical_prediction='UNKNOWN; no force fit before measurement', external_search_started=False))
        frozen.with_suffix('.json.sha256').write_text(sha(frozen) + '\n')
    write(P / 'CURRENT_WORK_STATE.json', dict(status='R1_COMPLETE', latest_gate=out['physical_gate'], next_operation='External closed-bite per-tooth force facit; change missing measurement representation', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    print(json.dumps({k: out[k] for k in ['source_crowns', 'possible_projected_pairs', 'source_pose_minimum_gap_mm', 'target36_first_touch_gap_mm', 'height_for_target36_to_become_first_touch_mm', 'cost']}, indent=2))
    print('controls', out['controls']['all_pass'])
if __name__ == '__main__':
    main()
