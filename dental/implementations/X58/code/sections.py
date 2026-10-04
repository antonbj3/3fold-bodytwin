from measure import *

def ray_exit(mask, origin, direction, max_length=6.0):
    """First contiguous exit from union of closed occupied voxel cubes centered at i*SP.
 Vectorized rays; ray starting outside yields NaN. No resampling interpolation.
 """
    origin = np.asarray(origin, float)
    direction = np.asarray(direction, float)
    p = origin / SP + 0.5
    i = np.floor(p).astype(int)
    shape = np.array(mask.shape)
    step = np.sign(direction).astype(int)
    inside = np.all((i >= 0) & (i < shape), axis=1)
    ok = np.zeros(len(i), bool)
    ok[inside] = mask[tuple(i[inside].T)]
    out = np.full(len(i), np.nan)
    active = ok.copy()
    face = i + (step > 0)
    safe = np.abs(direction) > 1e-14
    tmax = np.divide((face - p) * SP, direction, out=np.full_like(p, np.inf), where=safe)
    tdelta = np.divide(SP, np.abs(direction), out=np.full_like(p, np.inf), where=safe)
    for it in range(100):
        ids = np.flatnonzero(active)
        if len(ids) == 0:
            break
        tm = np.min(tmax[ids], axis=1)
        too = tm > max_length
        active[ids[too]] = False
        ids = ids[~too]
        tm = tm[~too]
        if not len(ids):
            continue
        crossed = np.abs(tmax[ids] - tm[:, None]) <= 1e-12
        i[ids] += crossed * step[ids]
        tmax[ids] = np.where(crossed, tmax[ids] + tdelta[ids], tmax[ids])
        in_bounds = np.all((i[ids] >= 0) & (i[ids] < shape), axis=1)
        occ = np.zeros(len(ids), bool)
        occ[in_bounds] = mask[tuple(i[ids[in_bounds]].T)]
        exiting = ~occ
        out[ids[exiting]] = tm[exiting]
        active[ids[exiting]] = False
    return out

def summarize(a):
    a = np.asarray(a)
    a = a[np.isfinite(a)]
    return {'n': len(a), 'median_mm': float(np.median(np.abs(a))) if len(a) else None, 'p95_mm': float(np.quantile(np.abs(a), 0.95)) if len(a) else None, 'signed_median_mm': float(np.median(a)) if len(a) else None, 'max_mm': float(np.max(np.abs(a))) if len(a) else None}

def main():
    start = time.perf_counter()
    pat = json.loads((ROOT / 'raw/PAIRED_PATIENTS.json').read_text())
    reports = []
    artifacts = []
    names = {0: 'mental_endpoint_5mm_proxy', 1: 'molar_nearest_FDI', 2: 'posterior_endpoint_10mm_proxy', -1: 'UNRESOLVED'}
    for patient in pat:
        if not patient['image_identity']['exact_numeric_equal']:
            continue
        case = patient['case_tf2']
        nvpath = DENT / 'results/NV1_canals/per_case' / f'{case}.json'
        if not nvpath.exists():
            reports.append({'patient': patient['patient'], 'status': 'NO_NV1_CENTERLINE'})
            continue
        nv = json.loads(nvpath.read_text())
        C = []
        U = []
        sgs = []
        fd = []
        station = []
        sides = []
        for (sd, c) in nv.get('canals', {}).items():
            if not c.get('present') or not c.get('C'):
                continue
            centers = np.array(c['C'])
            tang = np.gradient(centers, axis=0)
            tang /= np.linalg.norm(tang, axis=1)[:, None]
            s = np.array(c['s'])
            front = c.get('anterior_end')
            ai = int(np.argmin(np.linalg.norm(centers - np.array(front['point']), axis=1))) if front else None
            ds = np.abs(s - s[ai]) if ai is not None else None
            for (k, (cc, t)) in enumerate(zip(centers, tang)):
                v = np.cross(t, [1.0, 0, 0] if abs(t[0]) < 0.9 else [0, 1.0, 0])
                v /= np.linalg.norm(v)
                w = np.cross(t, v)
                th = np.arange(36) * 2 * np.pi / 36
                us = np.cos(th)[:, None] * v + np.sin(th)[:, None] * w
                fdi = -1 if c['region'][k] is None else c['region'][k]
                sg = -1
                if ds is not None and ds[k] <= 5:
                    sg = 0
                elif fdi in (36, 37, 38, 46, 47, 48):
                    sg = 1
                elif ds is not None and ds[k] >= ds.max() - 10:
                    sg = 2
                C.extend([cc] * 36)
                U.extend(us)
                sgs.extend([sg] * 36)
                fd.extend([fdi] * 36)
                station.extend([k] * 36)
                sides.extend([0 if sd == 'L' else 1] * 36)
        if not C:
            reports.append({'patient': patient['patient'], 'status': 'NO_NV1_CENTERLINE'})
            continue
        C = np.array(C)
        U = np.array(U)
        sgs = np.array(sgs)
        shape = patient['shape']
        dat = np.load(patient['point_artifact']['path'])
        radius = []
        center_support = []
        for name in ('old_voxels_zyx', 'new_voxels_zyx', 'tf2_voxels_zyx'):
            mask = np.zeros(shape, bool)
            mask[tuple(dat[name].T)] = True
            inds = np.floor(C / SP + 0.5).astype(int)
            inside = np.all((inds >= 0) & (inds < np.array(shape)), 1)
            supported = np.zeros(len(C), bool)
            supported[inside] = mask[tuple(inds[inside].T)]
            radius.append(ray_exit(mask, C, U))
            center_support.append(supported)
            del mask
        (old, new, tf) = radius
        (os, ns, ts) = center_support
        common = os & ns & np.isfinite(old) & np.isfinite(new)
        delta = new - old
        report = {'patient': patient['patient'], 'case': case, 'revised': not patient['numeric_mask_equal'], 'resolution': 'PER_POINT', 'status': 'COMPLETE', 'stations': len(C) // 36, 'rays': len(C), 'stations_inside_old': int(os[::36].sum()), 'stations_inside_new': int(ns[::36].sum()), 'stations_new_only': int((ns & ~os)[::36].sum()), 'stations_old_only': int((os & ~ns)[::36].sum()), 'stations_neither': int((~os & ~ns)[::36].sum()), 'common_valid_rays': int(common.sum()), 'common_ray_fraction': float(common.mean()), 'old_new_supported_radial': summarize(delta[common]), 'tf2_old_supported_radial': summarize((tf - old)[os & ts & np.isfinite(old) & np.isfinite(tf)]), 'tf2_new_supported_radial': summarize((tf - new)[ns & ts & np.isfinite(new) & np.isfinite(tf)]), 'regions': [], 'global_surface_p95_mm': patient['symmetric']['p95_mm'], 'limitations': 'Radial first exits are conditional on shared TF2-derived centers; missing centers retained as coverage loss. Neither masks nor stations are independent anatomical truth.'}
        for (sg, name) in names.items():
            report['regions'].append({'segment': name, 'resolution': 'PER_SURFACE_REGION', **summarize(delta[common & (sgs == sg)]), 'common_valid_rays': int((common & (sgs == sg)).sum()), 'all_rays': int((sgs == sg).sum())})
        file = DATA / f"{patient['patient']}_sections.npz"
        np.savez_compressed(file, center_zyx_mm=C, direction_zyx=U, old_radius_mm=old, new_radius_mm=new, tf2_radius_mm=tf, old_center_supported=os, new_center_supported=ns, tf2_center_supported=ts, segment=sgs, fdi=fd, station=station, side=sides)
        artifact = {'path': str(file), 'sha256': sha_file(file), 'bytes': file.stat().st_size, 'resolution': 'PER_POINT'}
        report['artifact'] = artifact
        artifacts.append(artifact)
        reports.append(report)
        print(patient['patient'], report['stations_new_only'], report['old_new_supported_radial']['p95_mm'], flush=True)
    dump(ROOT / 'raw/PER_PATIENT_SECTIONS.json', reports)
    revised = [r for r in reports if r.get('revised') and r.get('status') == 'COMPLETE' and r['old_new_supported_radial']['n']]
    global95 = float(np.median([r['global_surface_p95_mm'] for r in revised]))
    radial95 = float(np.median([r['old_new_supported_radial']['p95_mm'] for r in revised]))
    ratio = global95 / radial95 if radial95 else None
    cube = np.zeros((15, 15, 15), bool)
    cube[3:12, 4:11, 5:10] = True
    rng = np.random.default_rng(58)
    dirs = rng.normal(size=(127, 3))
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    orig = np.tile(np.array([7, 7, 7]) * SP, (127, 1))
    got = ray_exit(cube, orig, dirs, max_length=20)
    lo = (np.array([3, 4, 5]) - 0.5) * SP
    hi = (np.array([11, 10, 9]) + 0.5) * SP
    exact = np.min(np.where(dirs > 0, (hi - orig) / dirs, (lo - orig) / dirs), axis=1)
    error = float(np.max(np.abs(got - exact)))
    out = {'claim_type': 'information_link', 'patients_complete': sum((r.get('status') == 'COMPLETE' for r in reports)), 'revised_patients_with_shared_wall': len(revised), 'median_global_patient_p95_mm': global95, 'median_shared_wall_patient_radial_p95_mm': radial95, 'ratio_global_to_supported': ratio, 'stations_new_only': sum((r.get('stations_new_only', 0) for r in reports)), 'stations_old_only': sum((r.get('stations_old_only', 0) for r in reports)), 'stations_total': sum((r.get('stations', 0) for r in reports)), 'common_valid_rays': sum((r.get('common_valid_rays', 0) for r in reports)), 'total_rays': sum((r.get('rays', 0) for r in reports)), 'gates': {'coverage_confounding': ratio is not None and ratio > 2, 'analytic_first_exit_control': error <= 1e-10, 'injected_0p3mm_exit_rejected': float(np.max(np.abs(got + 0.3 - exact))) > 1e-10, 'independent_noise_floor': False}, 'first_exit_control': {'cuboid_bounds_mm': [lo.tolist(), hi.tolist()], 'rays': len(dirs), 'max_abs_error_mm': error, 'formal_floating_point_enclosure': 'NOT_CERTIFIED'}, 'dropout': {'patient_candidates': len(pat), 'no_exact_tf2_image_or_centerline': len(pat) - sum((r.get('status') == 'COMPLETE' for r in reports)), 'unsupported_or_censored_ray_fraction': 1 - sum((r.get('common_valid_rays', 0) for r in reports)) / sum((r.get('rays', 0) for r in reports))}, 'cost': {'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'fit_seconds': 0, 'questions': 0}, 'artifacts': artifacts}
    dump(ROOT / 'raw/R4_RESULT.json', out)
    assert all((out['gates'][k] for k in ('analytic_first_exit_control', 'injected_0p3mm_exit_rejected')))
    state('R4_COMPLETE', out['gates'], 'Package one -command demo, figure, external reference imitations and graph feedback')
    print(json.dumps({k: v for (k, v) in out.items() if k != 'artifacts'}, indent=2))
if __name__ == '__main__':
    main()
