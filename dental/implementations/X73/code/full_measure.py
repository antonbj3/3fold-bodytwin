"""One patient in RAM. Fresh measurements; cache replay is deliberately absent."""
from dental_release.paths import expand as _release_expand
from common import *
import gc, resource, traceback
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
sys.path.insert(0, str(X58 / 'code'))
from sections import ray_exit
from measure import region_map
sys.path.insert(0, str(X8))
from full_geometry import voxel_cylinder_bracket
SEG = {0: 'mental_endpoint_5mm_proxy', 1: 'molar_nearest_FDI', 2: 'posterior_endpoint_10mm_proxy', -1: 'UNRESOLVED'}

def stations(case):
    p = DENT / 'results/NV1_canals/per_case' / f'{case}.json'
    if not p.exists():
        return None
    nv = json.loads(p.read_text())
    C = []
    U = []
    SG = []
    FD = []
    ST = []
    SIDE = []
    for (side, c) in nv.get('canals', {}).items():
        if not c.get('present') or not c.get('C'):
            continue
        centers = np.array(c['C'])
        tangent = np.gradient(centers, axis=0)
        norm = np.linalg.norm(tangent, axis=1)
        if np.any(norm == 0):
            raise ValueError('Degenerate centerline tangent')
        tangent /= norm[:, None]
        s = np.array(c['s'])
        front = c.get('anterior_end')
        ai = int(np.argmin(np.linalg.norm(centers - np.array(front['point']), axis=1))) if front else None
        ds = np.abs(s - s[ai]) if ai is not None else None
        for (i, (cc, t)) in enumerate(zip(centers, tangent)):
            v = np.cross(t, [1.0, 0, 0] if abs(t[0]) < 0.9 else [0, 1.0, 0])
            v /= np.linalg.norm(v)
            w = np.cross(t, v)
            phi = np.arange(36) * 2 * np.pi / 36
            dirs = np.cos(phi)[:, None] * v + np.sin(phi)[:, None] * w
            fd = -1 if c['region'][i] is None else int(c['region'][i])
            sg = -1
            if ds is not None and ds[i] <= 5:
                sg = 0
            elif fd in (36, 37, 38, 46, 47, 48):
                sg = 1
            elif ds is not None and ds[i] >= ds.max() - 10:
                sg = 2
            C.extend([cc] * 36)
            U.extend(dirs)
            SG.extend([sg] * 36)
            FD.extend([fd] * 36)
            ST.extend([i] * 36)
            SIDE.extend([side] * 36)
    if not C:
        return None
    return dict(center_zyx_mm=np.array(C), direction_zyx=np.array(U), segment=np.array(SG, dtype=np.int8), fdi=np.array(FD, dtype=np.int16), station=np.array(ST, dtype=np.int16), side=np.array(SIDE), nv1_path=str(p), nv1_sha256=sha(p))

def bracket(points, s):
    assert s['radius_mm'] <= 2 and abs(np.linalg.norm(s['axis_zyx']) - 1) < 1e-12
    x = voxel_cylinder_bracket(points * SP, np.full(3, SP), np.array(s['entry_zyx_mm']), np.array(s['axis_zyx']), s['length_mm'], s['radius_mm'])
    return {k: x[k] for k in ('lower_mm', 'upper_mm', 'gap_mm', 'active_voxel_boxes', 'witness')}

def boundary(m):
    return np.argwhere(m & ~ndi.binary_erosion(m, structure=ndi.generate_binary_structure(3, 1), border_value=0))

def main():
    check_frozen()
    diskcheck()
    start = time.perf_counter()
    census = json.loads((ROOT / 'raw/PATIENT_CENSUS.json').read_text())
    patients = [p for p in census if p['dense_pair_candidate']]
    certpath = Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk/R5_FULL_GEOMETRY.jsonl'))
    certsha = sha(certpath)
    cert = [json.loads(l) for l in certpath.read_text().splitlines()]
    groups = {}
    for s in cert:
        groups.setdefault(s['case'], []).append(s)
    out = []
    sites = []
    fails = []
    pairs = []
    artifacts = []
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with zipfile.ZipFile(TF1) as a, zipfile.ZipFile(MAX) as b, zipfile.ZipFile(TF2) as c:
        ais = {i.filename.split('/')[-2]: i.filename for i in a.infolist() if i.filename.endswith('data.npy')}
        bis = {i.filename.split('/')[-2]: i.filename for i in b.infolist() if i.filename.endswith('data.npy')}
        for (ii, p) in enumerate(patients, 1):
            k = p['patient']
            case = p['case']
            t = time.perf_counter()
            row = dict(patient=k, case=case, resolution='PER_POINT', time_scale='SIMULTANEOUS', independent_annotations=False, regions=[], source_bindings=[])
            masks = {}
            pointsets = {}
            images = {}
            st = None
            try:
                (image, ih) = npy(a, ais[k])
                images['tf1'] = dict(member=ais[k], sha256=ih, shape=list(image.shape))
                shape = image.shape
                if any((w['dataset'] == 'maxillo' for w in p['walls'])):
                    ib = member(b, bis[k])
                    maxsha = hashlib.sha256(ib).hexdigest()
                    del ib
                    assert maxsha == ih, 'Original image SHA mismatch'
                    images['maxillo'] = dict(member=bis[k], sha256=maxsha, exact_npy_byte_equal=True)
                if any((w['dataset'] == 'tf2' for w in p['walls'])):
                    (im2, h, h2) = mha(c, f'Dataset112_ToothFairy2/imagesTr/{case}_0000.mha')
                    sp = list(map(float, h['ElementSpacing'].split()[::-1]))
                    assert all((x == SP for x in sp))
                    assert h.get('TransformMatrix') == '1 0 0 0 1 0 0 0 1' and h.get('Offset') == '0 0 0'
                    eq = image.shape == im2.shape and np.array_equal(image, im2)
                    images['tf2'] = dict(member=f'Dataset112_ToothFairy2/imagesTr/{case}_0000.mha', sha256=h2, exact_numeric_equal=bool(eq), spacing_mm=sp, transform=h.get('TransformMatrix'), origin=h.get('Offset'))
                    if not eq:
                        raise ValueError('TF1-TF2 image array mismatch')
                    del im2
                del image
                row['image_identity'] = dict(images)
                row['shape'] = list(shape)
                for w in p['walls']:
                    key = w['dataset']
                    z = {'tf1': a, 'maxillo': b, 'tf2': c}[key]
                    if key == 'tf2':
                        (lab, h, hh) = mha(z, w['member'])
                        m = np.isin(lab, [3, 4])
                        del lab
                    else:
                        (lab, hh) = npy(z, w['member'])
                        m = lab > 0
                        del lab
                    assert m.shape == shape and m.any(), 'Invalid wall occupancy'
                    masks[key] = m
                    pointsets[key] = np.argwhere(m).astype(np.int16)
                    row['source_bindings'].append(dict(dataset=key, archive=str(z.filename), member=w['member'], member_sha256=hh, mask_numeric_sha256=arrsha(m), voxels=int(m.sum())))
                bounds = {key: boundary(m) for (key, m) in masks.items()}
                pairrows = []
                for (x, y) in [('maxillo', 'tf1'), ('tf1', 'tf2'), ('maxillo', 'tf2')]:
                    if x not in masks or y not in masks:
                        continue
                    same = np.array_equal(masks[x], masks[y])
                    d1 = cKDTree(bounds[y] * SP).query(bounds[x] * SP, workers=1)[0]
                    d2 = cKDTree(bounds[x] * SP).query(bounds[y] * SP, workers=1)[0]
                    inter = int((masks[x] & masks[y]).sum())
                    n1 = int(masks[x].sum())
                    n2 = int(masks[y].sum())
                    pr = dict(patient=k, pair=x + '__' + y, equal=bool(same), dice=2 * inter / (n1 + n2), added_percent=100 * (n2 - inter) / n1, resolution='PER_SURFACE_REGION', **summary(np.r_[d1, d2]))
                    pairrows.append(pr)
                    pairs.append(pr)
                row['pairs'] = pairrows
                del bounds
                if 'tf2' in masks:
                    st = stations(case)
                data = {key + '_voxels_zyx': v for (key, v) in pointsets.items()}
                if st:
                    row['nv1_sha256'] = st['nv1_sha256']
                    rad = {key: ray_exit(m, st['center_zyx_mm'], st['direction_zyx']) for (key, m) in masks.items()}
                    vals = np.array(list(rad.values()))
                    valid = np.isfinite(vals)
                    count = valid.sum(0)
                    complete = valid.all(0)
                    low = np.min(np.where(valid, vals, np.inf), axis=0)
                    high = np.max(np.where(valid, vals, -np.inf), axis=0)
                    spread = high - low
                    spread[count < 2] = np.nan
                    row['stations'] = len(st['segment']) // 36
                    row['rays'] = len(spread)
                    row['all_masks_valid_rays'] = int(complete.sum())
                    row['at_least_two_valid_rays'] = int((count >= 2).sum())
                    for (sg, name) in SEG.items():
                        chosen = st['segment'] == sg
                        strict = bool(chosen.any() and complete[chosen].all())
                        row['regions'].append(dict(segment=name, resolution='PER_SURFACE_REGION', all_rays=int(chosen.sum()), all_masks_valid_rays=int((chosen & complete).sum()), at_least_two_valid_rays=int((chosen & (count >= 2)).sum()), strict_complete=strict, strict_score_mm=float(spread[chosen].max()) if strict else 'INF', **summary(spread[chosen])))
                    for (key, v) in st.items():
                        if isinstance(v, np.ndarray):
                            data[key] = v
                    data.update({key + '_radius_mm': v for (key, v) in rad.items()})
                    data['wall_range_mm'] = spread
                    data['all_masks_supported'] = complete
                else:
                    row['regional_status'] = 'NO_VERIFIED_TF2_CENTERLINE; not imputed'
                ownrows = []
                if 'tf2' in masks:
                    for s in groups.get(case, []):
                        br = {key: bracket(points, s) for (key, points) in pointsets.items()}
                        assert all((v['gap_mm'] <= 0.0001 for v in br.values())), 'Cylinder bracket wider than frozen tolerance'
                        ref = br['tf2']
                        referr = max((abs(ref[q] - s[q]) for q in ['lower_mm', 'upper_mm']))
                        assert referr <= 0.0001, 'X8 TF2 replay mismatch'
                        lo = min((x['lower_mm'] for x in br.values()))
                        hi = min((x['upper_mm'] for x in br.values()))
                        cls = {key: classify(v['lower_mm'], v['upper_mm']) for (key, v) in br.items()}
                        (fd, sg, info) = region_map(case, np.array([ref['witness']['witness_zyx_mm']]) / SP)
                        rr = dict(patient=k, case=case, fdi=s['fdi'], resolution='PER_TOOTH', time_scale='SIMULTANEOUS', segment=SEG[int(sg[0])], witness_fdi=int(fd[0]), pose={q: s[q] for q in ['entry_zyx_mm', 'axis_zyx', 'length_mm', 'radius_mm']}, distances=br, classes=cls, union=dict(lower_mm=lo, upper_mm=hi, class_2mm=classify(lo, hi)), loss_lower_mm=max(0.0, ref['lower_mm'] - hi), loss_upper_mm=max(0.0, ref['upper_mm'] - lo), tf2_replay_max_error_mm=referr, reference_certificate_sha256=certsha, source_bindings=row['source_bindings'])
                        ownrows.append(rr)
                        sites.append(rr)
                row['site_count'] = len(ownrows)
                artifact = DATA / f'{k}_points.npz'
                np.savez_compressed(artifact, **data)
                binding = dict(path=str(artifact), sha256=sha(artifact), bytes=artifact.stat().st_size, resolution='PER_POINT')
                artifacts.append(binding)
                row['point_artifact'] = binding
                for rr in ownrows:
                    rr['point_artifact'] = binding
                dump(ROOT / 'raw/sites' / f'{k}.json', ownrows)
                row.update(status='COMPLETE', wall_seconds=time.perf_counter() - t)
                out.append(row)
                dump(ROOT / 'raw/patients' / f'{k}.json', row)
                print(f"{ii}/{len(patients)} {k} masks={len(masks)} sites={len(ownrows)} sec={row['wall_seconds']:.2f}", flush=True)
            except Exception as e:
                failure = dict(patient=k, error=repr(e), traceback=traceback.format_exc(), row=row)
                fails.append(failure)
                dump(ROOT / 'raw/MEASURE_FAILURES.json', fails)
                print(k, repr(e), flush=True)
            finally:
                masks.clear()
                pointsets.clear()
                images.clear()
                st = None
                if 'data' in locals():
                    del data
                if 'rad' in locals():
                    del rad, vals, valid, count, complete, low, high, spread
                gc.collect()
            dump(ROOT / 'raw/PATIENTS.json', out)
            dump(ROOT / 'raw/SITES.json', sites)
            dump(ROOT / 'raw/POINT_MANIFEST.json', artifacts)
            state('R2_R3_MEASURING', 'R1_INDEPENDENCE_UNKNOWN', 'Continue full dense patient stream', patients_complete=len(out), patients_candidates=len(patients), sites_complete=len(sites), failures=len(fails))
            diskcheck()
    dump(ROOT / 'raw/PAIRS.json', pairs)
    dump(ROOT / 'raw/MEASURE_FAILURES.json', fails)
    n = len(out)
    result = dict(claim_type='information_link', started_utc=started, completed_patients=n, dense_pair_candidates=len(patients), exact_tf2_image_pairs=sum(('tf2' in r['image_identity'] for r in out)), total_sites=len(sites), x8_total_sites=len(cert), x8_paired_patients=len(set((r['patient'] for r in sites))), all_candidate_patients_completed=n == len(patients), failures=len(fails), dropout=dict(dense_candidate_rejections=len(patients) - n, site_unmatched=len(cert) - len(sites), site_unmatched_fraction=1 - len(sites) / len(cert), all_masks_supported_ray_fraction=sum((r.get('all_masks_valid_rays', 0) for r in out)) / sum((r.get('rays', 0) for r in out))), cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, threads=1, gpu=False, scratch_bytes=diskcheck()), formal_floating_point_enclosure='MISSING; numerical brackets do not certify outward IEEE arithmetic')
    dump(ROOT / 'raw/MEASURE_RESULT.json', result)
    print(json.dumps(result, indent=2), flush=True)
    state('R2_R3_GEOMETRY_COMPLETE', dict(all_candidates_completed=result['all_candidate_patients_completed'], independent_spread=False), 'Patient-disjoint calibration and fault controls', sites_complete=len(sites))
if __name__ == '__main__':
    main()
