"""Exact exterior carrier, separate preparation edit, leave-one-tooth replacement."""
from common import *
from scipy import ndimage
from scipy.spatial import cKDTree
import trimesh, resource

def crop(tri, lo, hi):
    return tri[np.all(tri[:, :, :2].max(1) >= lo, axis=1) & np.all(tri[:, :, :2].min(1) <= hi, axis=1)]

def specimen(a, fdi, depth):
    ids = np.flatnonzero(a['owner'] == fdi)
    if len(ids) < 100:
        raise ValueError('Fewer than100 predicted source faces; no interpolation')
    ff = a['f'][ids]
    (vi, inv) = np.unique(ff, return_inverse=True)
    faces = inv.reshape(-1, 3)
    v = a['v'][vi].copy()
    normals = np.cross(v[faces[:, 1]] - v[faces[:, 0]], v[faces[:, 2]] - v[faces[:, 0]])
    vn = np.zeros_like(v)
    for k in range(3):
        np.add.at(vn, faces[:, k], normals)
    vn /= np.maximum(np.linalg.norm(vn, axis=1)[:, None], 1e-15)
    prep = v - depth * vn
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    ee = np.sort(edges, axis=1)
    (_, ix, counts) = np.unique(ee, axis=0, return_index=True, return_counts=True)
    boundary = edges[ix[counts == 1]]
    nv = len(v)
    collar = np.asarray([[a, b, b + nv] for (a, b) in boundary] + [[a, b + nv, a + nv] for (a, b) in boundary], int).reshape(-1, 3)
    allv = np.r_[v, prep]
    allf = np.r_[faces, faces[:, ::-1] + nv, collar]
    role = np.r_[np.full(len(faces), 0), np.full(len(faces), 1), np.full(len(collar), 2)]
    pn = np.cross(prep[faces[:, 1]] - prep[faces[:, 0]], prep[faces[:, 2]] - prep[faces[:, 0]])
    folds = np.einsum('ij,ij->i', pn, normals) <= 0
    eall = np.sort(np.concatenate([allf[:, [0, 1]], allf[:, [1, 2]], allf[:, [2, 0]]]), axis=1)
    ct = np.unique(eall, axis=0, return_counts=True)[1]
    return dict(vertices=allv, faces=allf, face_roles=role, source_face_ids=ids, source_vertex_ids=vi, external_vertices=v, external_faces=faces, preparation_vertices=prep, preparation_fold_count=int(folds.sum()), preparation_fold_fraction=float(folds.mean()), nonmanifold_edge_count=int(np.sum(ct != 2)), boundary_edges=boundary)

def read_stl(p):
    blob = Path(p).read_bytes()
    n = int.from_bytes(blob[80:84], 'little')
    dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    assert len(blob) == 84 + 50 * n
    return np.frombuffer(blob, dt, offset=84)['v'].astype(float)

def predict():
    pr = json.loads((H / 'PREREG_R1.json').read_text())
    start = time.perf_counter()
    rows = []
    files = []
    for case in pr['cases']:
        state('R1_CONSTRUCTING', 'Frozen exact-surface contract', 'Source carrier ' + str(case), completed_sites=len(rows))
        (data, man) = pair(case)
        for (jaw, fdilist) in pr['fdi_by_jaw'].items():
            for fdi in fdilist:
                try:
                    s = specimen(data[jaw], fdi, pr['preparation']['normal_offset_mm'])
                    stem = f'{case}_{fdi}'
                    path = D / (stem + '_carrier.npz')
                    np.savez_compressed(path, **s)
                    ex = D / (stem + '_external.stl')
                    inner = D / (stem + '_preparation.stl')
                    shell = D / (stem + '_specimen_shell.stl')
                    write_stl(ex, s['external_vertices'], s['external_faces'])
                    write_stl(inner, s['preparation_vertices'], s['external_faces'])
                    write_stl(shell, s['vertices'], s['faces'])
                    pf = [artifact(x) for x in [path, ex, inner, shell]]
                    files.extend(pf)
                    rows.append(dict(case=case, jaw=jaw, fdi=fdi, status='FROZEN', carrier_file=str(path), external_stl=str(ex), preparation_stl=str(inner), specimen_stl=str(shell), source=man, source_faces=len(s['source_face_ids']), prep_fold_count=s['preparation_fold_count'], prep_fold_fraction=s['preparation_fold_fraction'], nonmanifold_edge_count=s['nonmanifold_edge_count']))
                except ValueError as e:
                    rows.append(dict(case=case, jaw=jaw, fdi=fdi, status='REJECTED', reason=str(e)))
        dump(H / 'raw/PREDICTIONS_R1_checkpoint.json', rows)
        print('Carrier frozen', case, len(rows), flush=True)
        data_budget()
    dump(H / 'raw/PREDICTIONS_R1.json', rows)
    files.append(artifact(H / 'raw/PREDICTIONS_R1.json'))
    dump(H / 'FROZEN_PREDICTIONS_R1.json', dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R1.json'), files=files, code_sha256={str(p): sha(p) for p in [H / 'code/surface.py', H / 'code/common.py']}, reference_contact_maps_read=False, target_exterior_used_as_authorized_preoperative_input=True, physical_preparation_measurement_status='NOT_RUN'))
    dump(H / 'raw/R1_PREPARATION_COST.json', dict(wall_s=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, data_bytes=data_budget()))
    state('R1_PREDICTIONS_FROZEN', 'Exact carriers hashed before reference contact evaluation', 'Query untouched original arch and exported exteriors')

def grid(tri, step):
    lo = tri[:, :, :2].min((0, 1))
    hi = tri[:, :, :2].max((0, 1))
    axes = [np.arange(lo[k], hi[k] + step * 0.25, step) for k in range(2)]
    (x, y) = np.meshgrid(*axes, indexing='ij')
    xy = np.c_[x.ravel(), y.ravel()]
    index = np.arange(len(xy)).reshape(x.shape)
    return (xy, index, np.full(len(xy), step ** 2))

def minimum(u, l):
    up = u[np.abs(cg.cross(u[:, 1, :2] - u[:, 0, :2], u[:, 2, :2] - u[:, 0, :2])) > 1e-12]
    lp = l[np.abs(cg.cross(l[:, 1, :2] - l[:, 0, :2], l[:, 2, :2] - l[:, 0, :2])) > 1e-12]
    pairs = cg.broad_phase(up, lp)
    if not len(pairs):
        return (dict(minimum_gap_mm=None), dict(minimum_gap_mm=None), 0)
    return (cg.extremum(up, lp, pairs, True), cg.extremum(up, lp, pairs, False), len(pairs))

def evaluate():
    verify_freeze('R1')
    pr = json.loads((H / 'PREREG_R1.json').read_text())
    preds = json.loads((H / 'raw/PREDICTIONS_R1.json').read_text())
    rows = []
    cur = None
    start = time.perf_counter()
    metric = pr['metrics']
    for r in preds:
        if r['status'] != 'FROZEN':
            rows.append(r)
            continue
        if cur != r['case']:
            (data, man) = pair(r['case'])
            cur = r['case']
        t = time.perf_counter()
        a = data[r['jaw']]
        opp = data['lower' if r['jaw'] == 'upper' else 'upper']
        src = a['tri'][a['owner'] == r['fdi']]
        candidate = read_stl(r['external_stl'])
        (xy, index, weights) = grid(src, pr['grid_step_mm'])
        opposing = crop(opp['tri'], xy.min(0), xy.max(0))
        isupper = r['jaw'] == 'upper'
        (reference_z, _) = query_height(src, xy, isupper)
        (candidate_z, _) = query_height(candidate, xy, isupper)
        (oz, _) = query_height(opposing, xy, not isupper)
        refgap = reference_z - oz if isupper else oz - reference_z
        gap = candidate_z - oz if isupper else oz - candidate_z
        fault = candidate.copy()
        fault[:, :, 2] += 0.2
        (fz, _) = query_height(fault, xy, isupper)
        fgap = fz - oz if isupper else oz - fz
        (U, L) = (candidate, opposing) if isupper else (opposing, candidate)
        (fast, full, n_pairs) = minimum(U, L)
        roundtrip = float(np.max(np.abs(src - candidate)))
        maps = {}
        array = dict(xy=xy, index=index, weights=weights, source_z=reference_z, candidate_z=candidate_z, opposing_z=oz, reference_gap=refgap)
        finite = np.isfinite(reference_z) & np.isfinite(oz)
        pose_areas = {}
        for band in [pr['contact_band_mm']] + pr['sensitivity_bands_mm']:
            truth = contact(refgap, reference_z, xy, index, weights, band, pr['minimum_patch_area_mm2'])
            pred = contact(gap, candidate_z, xy, index, weights, band, pr['minimum_patch_area_mm2'])
            bad = contact(fgap, fz, xy, index, weights, band, pr['minimum_patch_area_mm2'])
            err = score(pred, truth, xy)
            ferr = score(bad, truth, xy)
            key = str(band)
            array[key + '_reference_contact'] = truth.pop('mask')
            array[key + '_candidate_contact'] = pred.pop('mask')
            bad.pop('mask')
            maps[key] = dict(reference=truth, candidate=pred, errors=err, injected_shift_errors=ferr)
        for delta in [-0.03, 0, 0.03]:
            pose_areas[str(delta)] = contact(refgap + delta, reference_z, xy, index, weights, 0.1, 0.1)['area_mm2']
        parity = abs(fast['minimum_gap_mm'] - full['minimum_gap_mm']) if fast['minimum_gap_mm'] is not None else 0
        (original_fast, original_full, _) = minimum(*((src, opposing) if isupper else (opposing, src)))
        original_parity = abs(fast['minimum_gap_mm'] - original_full['minimum_gap_mm']) if fast['minimum_gap_mm'] is not None else 0
        co = maps['0.1']['errors']
        gates = dict(surface=roundtrip <= metric['surface_roundtrip_max_mm'], gap=parity <= metric['continuous_gap_parity_mm'] and original_parity <= metric['continuous_gap_parity_mm'], count=co['absolute_patch_count_error'] <= metric['patch_count_absolute_error_max'], centroid=co['centroid_error_mm'] <= metric['contact_centroid_max_mm'], IoU=co['IoU'] >= metric['contact_IoU_min'])
        fault_gates = dict(surface_shift=float(np.max(np.abs(src - fault))) > metric['surface_roundtrip_max_mm'], gap_shift=0.2 > metric['continuous_gap_parity_mm'], contact_shift=maps['0.1']['injected_shift_errors']['IoU'] < metric['contact_IoU_min'] if maps['0.1']['reference']['patch_count'] else None)
        arrpath = D / f"{r['case']}_{r['fdi']}_contact_reference.npz"
        np.savez_compressed(arrpath, **array)
        rows.append(dict(**{k: v for (k, v) in r.items() if k != 'status'}, status='SCORED', resolution='PER_TOOTH', contact_resolution='PER_SURFACE_REGION', surface_max_error_mm=roundtrip, finite_opposing_support_fraction=float(finite.mean()), continuous=fast, full_enumeration=full, continuous_parity_mm=parity, original_parity_mm=original_parity, broad_phase_pairs=n_pairs, contacts=maps, pose_shift_area_mm2=pose_areas, gates=gates, fault_rejections=fault_gates, physical_pose_admissible=fast['minimum_gap_mm'] is not None and fast['minimum_gap_mm'] >= -metric['physical_penetration_max_mm'], array_file=artifact(arrpath), query_wall_s=time.perf_counter() - t))
        dump(H / 'raw/RESULTS_R1_checkpoint.json', rows)
        print('Surface evaluated', r['case'], r['fdi'], 'patches', maps['0.1']['reference']['patch_count'], 'min gap', fast['minimum_gap_mm'], flush=True)
    good = [r for r in rows if r['status'] == 'SCORED']
    gates = dict(coverage=len(good) / len(preds) >= metric['coverage_min'], surface=all((all(r['gates'].values()) for r in good)), injected_shift=all((all((v for v in r['fault_rejections'].values() if v is not None)) for r in good)))
    out = dict(round='R1', claim_type='information_link', external_referent=pr['external_referent'], decision='SURFACE_PASS' if all(gates.values()) else 'SURFACE_FAIL', gates=gates, attempted=len(preds), retained=len(good), rejected=len(preds) - len(good), rejection_fraction=(len(preds) - len(good)) / len(preds), rows=rows, physical_contact_pose_admissible=sum((r['physical_pose_admissible'] for r in good)), physical_contact_pose_rejected=sum((not r['physical_pose_admissible'] for r in good)), preparation_orientation_pass_count=sum((r['prep_fold_count'] == 0 for r in good)), preparation_physical_validation='UNKNOWN; normal-offset specimen may fold and margins/self-intersections/tool accessibility are not externally validated', cost=dict(wall_s=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), scope='Lossless preoperative surface roundtrip only. Leave-one-tooth replacement; original target exterior was an input, no held-out anatomy inference.')
    dump(H / 'rounds/R1.json', out)
    (H / 'HANDOFF_R1.md').write_text(f"R1 {out['decision']}: {len(good)}/{len(preds)} source-indexed exteriors scored. Physical-pose rejections {out['physical_contact_pose_rejected']}. Preparation folds retained in rows. No de novo reconstruction or validated clinical preparation. Next construction: jaw-specific tooth-force intervals and pressure-simplex stress propagation (PREREG_R2).\n")
    state('R1_DECIDED', out['decision'], 'Run jaw-matched pressure intervals; preserve source-pose and preparation failures')
if __name__ == '__main__':
    {'predict': predict, 'evaluate': evaluate}[sys.argv[1]]()
