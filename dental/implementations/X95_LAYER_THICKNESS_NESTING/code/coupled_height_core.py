"""Source-safe basefootprint and pose are jointly chosen; virtual geometry only."""
from height_core import *

def construct(allow, shape, fullidx, lo, h):
    occ = np.zeros(shape, bool)
    occ.ravel()[fullidx[allow]] = True
    if not occ.any():
        return None
    base = base_footprint(occ)
    caps = vertex_caps(cap_prefix(allow, shape, fullidx, base))
    assert np.all(caps > 0)
    return (base, caps, formula_volume(caps, h))

def independent_rectangle_area(a):
    best = 0
    for i in range(a.shape[0]):
        intersection = np.ones(a.shape[1], bool)
        for j in range(i, a.shape[0]):
            intersection &= a[j]
            longest = 0
            current = 0
            for yes in intersection:
                current = current + 1 if yes else 0
                longest = max(longest, current)
            best = max(best, (j - i + 1) * longest)
    return best

def run():
    st = time.perf_counter()
    pr = json.loads((ROOT / 'PREREG_R4.json').read_text())
    r1 = json.loads((ROOT / 'PREREG_R1.json').read_text())
    assert source_gate(json.loads((ROOT / 'RULE_CONTRACT.json').read_text()))
    rows = []
    h = pr['model_contract']['h_mm']
    rv = math.sqrt(3) * h / 2
    for e in json.loads((ROOT / 'inputs/R4_FROZEN_EXPORTS.json').read_text())['exports']:
        t0 = time.perf_counter()
        raw = dict(np.load(e['mesh_path'], allow_pickle=False))
        a = dict(np.load(ROOT / 'exports' / f"{e['family']}_virtual_core.npz", allow_pickle=False))
        y = a['cube_centers_mm']
        fullidx = a['eligible_original_core_index']
        shape = tuple(a['grid_shape'])
        lo = a['grid_origin_mm']
        outer = raw['vertices'][raw['faces'][raw['roles'] == 0]]
        c = outer.mean(1)
        rho = np.linalg.norm(outer - c[:, None], axis=2).max(1)
        (high, low) = (0.8, 0.4) if e['family'] == 'anterior' else (1.0, 0.5)
        (iy, it, dist) = incidences(y, c, rho, high, rv)
        blocked = np.zeros(len(y), bool)
        blocked[iy[dist < low + rho[it] + rv]] = True
        uni_block = np.zeros(len(y), bool)
        uni_block[iy[dist < 0.5 + rho[it] + rv]] = True
        baseline = ~uni_block
        (ub, uc, uv) = construct(baseline, shape, fullidx, lo, h)
        grid = [p for p in poses(raw['vertices'], r1) if p['stock_ok']]
        signatures = {}
        cache = {}
        best = None
        bestallow = None
        bestbase = None
        bestcaps = None
        records = []
        for (ix, pose) in enumerate(grid):
            n = axis(pose['tilt_deg'], pose['azimuth_deg'])
            key = (pose['tilt_deg'], pose['azimuth_deg'])
            top = np.max(raw['vertices'] @ n)
            if key not in signatures:
                s = np.full(len(y), -np.inf)
                np.maximum.at(s, iy, np.max(outer @ n, axis=1)[it])
                signatures[key] = s
            allow = ~blocked & (pose['top_depth_mm'] + top - signatures[key] > pose['height_mm'] / 2)
            sig = np.packbits(allow).tobytes()
            if sig not in cache:
                cache[sig] = construct(allow, shape, fullidx, lo, h)
            result = cache[sig]
            vol = result[2] if result else 0.0
            rec = dict(pose, pose_index=ix, height_core_feasible=result is not None, retained_volume_mm3=vol)
            records.append(rec)
            if result and (best is None or vol > best['retained_volume_mm3']):
                best = rec
                (bestbase, bestcaps, _) = result
                bestallow = allow.copy()
        assert best is not None
        n = axis(best['tilt_deg'], best['azimuth_deg'])
        ztri = best['top_depth_mm'] + np.max(raw['vertices'] @ n) - outer @ n
        rr = np.where(ztri.min(1) <= best['height_mm'] / 2, high, low)
        fresh = np.sqrt(np.sum((y[iy] - c[it]) ** 2, 1))
        b = np.zeros(len(y), bool)
        b[iy[fresh < rr[it] + rho[it] + rv]] = True
        control = sequential_control(~b, shape, fullidx, bestbase)
        assert np.array_equal(control, bestcaps)
        uc_check = sequential_control(baseline, shape, fullidx, ub)
        assert np.array_equal(uc_check, uc)
        occ = np.zeros(shape, bool)
        occ.ravel()[fullidx[~b]] = True
        independent_area = max((independent_rectangle_area(occ[:, :, k]) for k in range(shape[2])))
        assert independent_area == bestbase[0]
        injected = bestcaps.copy()
        injected[0, 0] += 1
        assert np.any(injected > control)
        m = height_mesh(bestcaps, bestbase, lo, h)
        um = height_mesh(uc, ub, lo, h)
        assert m.is_watertight and m.is_winding_consistent and (m.volume > 0) and (len(m.split(only_watertight=False)) == 1)
        assert um.is_watertight and um.is_winding_consistent and (um.volume > 0)
        volumeerr = abs(float(m.volume) - best['retained_volume_mm3'])
        assert volumeerr <= 1e-09
        (_, k, (x0, x1, y0, y1)) = bestbase
        cnt = cap_prefix(bestallow, shape, fullidx, bestbase)
        used = np.zeros(shape, bool)
        for i in range(cnt.shape[0]):
            for j in range(cnt.shape[1]):
                used[x0 + i, y0 + j, k:k + cnt[i, j]] = True
        usedmask = used.ravel()[fullidx]
        guards = np.full(len(y), np.inf)
        np.minimum.at(guards, iy, dist - rr[it] - rho[it] - rv)
        guardmin = float(guards[usedmask].min())
        assert guardmin >= 0 and (not np.any(usedmask & b))
        alt = usedmask.copy()
        safeidx = int(np.flatnonzero(usedmask)[0])
        badidx = int(np.flatnonzero(b)[0])
        alt[safeidx] = False
        alt[badidx] = True
        identity = (int(usedmask.sum()) - int(alt.sum())) * h ** 3
        assert identity == 0 and guards[alt].min() < 0
        dest = ROOT / 'exports' / f"{e['family']}_source_coupled_preparation.stl"
        m.export(dest)
        um.export(ROOT / 'exports' / f"{e['family']}_source_coupled_universal0_5.stl")
        path = ROOT / 'exports' / f"{e['family']}_source_coupled_core.npz"
        np.savez_compressed(path, vertices=m.vertices, faces=m.faces, heights_cubes=bestcaps, base_plane_grid=k, base_rectangle_xy=np.array(bestbase[2]), origin_mm=lo, h_mm=h, source_used_cubes=usedmask, universal_vertices=um.vertices, universal_faces=um.faces, locked_crown_vertices=raw['vertices'], locked_crown_faces=raw['faces'], locked_crown_roles=raw['roles'])
        reread = dict(np.load(path, allow_pickle=False))
        outererr = float(abs(reread['locked_crown_vertices'] - raw['vertices']).max())
        assert outererr == 0
        with (ROOT / 'raw' / f"{e['family']}_R4_POSES.csv").open('w') as f:
            w = csv.DictWriter(f, fieldnames=list(records[0]))
            w.writeheader()
            w.writerows(records)
        row = {'family': e['family'], 'key': e['key'], 'resolution': 'PER_POINT', 'source_best_pose': best, 'source_retained_volume_mm3': best['retained_volume_mm3'], 'universal0_5_retained_volume_mm3': uv, 'source_minus_universal_retained_mm3': best['retained_volume_mm3'] - uv, 'source_base_xy_grid': [int(q) for q in bestbase[2]], 'source_base_z_grid': int(k), 'source_base_area_mm2': bestbase[0] * h * h, 'universal_base_xy_grid': [int(q) for q in ub[2]], 'universal_base_z_grid': int(ub[1]), 'watertight': bool(m.is_watertight), 'winding_consistent': bool(m.is_winding_consistent), 'components': 1, 'positive_volume': bool(m.volume > 0), 'outer_identity_error_mm': outererr, 'volume_formula_mesh_error_mm3': volumeerr, 'guard_min_used_cubes_mm': guardmin, 'control_height_error_mm': 0.0, 'control_prefix_mismatches': 0, 'independent_max_rectangle_area_cells': int(independent_area), 'injected_height_over_cap_rejected': True, 'same_volume_sufficiency': {'identity_error_mm3': identity, 'valid_guard_min_mm': guardmin, 'invalid_guard_min_mm': float(guards[alt].min()), 'downstream_guard_difference_mm': float(guards[alt].min() - guardmin), 'minimal_extension': 'spatial core membership plus source/pose guard', 'origin': 'our_own_fixture, not external measurement'}, 'unique_guard_masks_evaluated': len(cache), 'stock_fitting_poses': len(grid), 'height_core_rejected_poses': sum((not p['height_core_feasible'] for p in records)), 'stl_path': str(dest.relative_to(ROOT)), 'stl_sha256': sha(dest), 'array_path': str(path.relative_to(ROOT)), 'array_sha256': sha(path), 'qualification': 'CONDITIONAL_VIRTUAL_HEIGHT_FAMILY_ONLY_NOT_FABRICATION_QUALIFIED', 'rigorous_float_enclosure': 'MISSING', 'seconds': time.perf_counter() - t0}
        rows.append(row)
        dump(ROOT / 'raw/R4_PARTIAL.json', rows)
        print(e['family'], 'closed', row['watertight'], 'source', row['source_retained_volume_mm3'], 'universal', uv, flush=True)
    out = {'round': 'R4', 'claim_type': 'capability', 'primary_gate': 'PASS_SCOPED_CONNECTED_VIRTUAL_GEOMETRY', 'rows': rows, 'seconds': time.perf_counter() - st, 'physical_measurement': 'NOT_PERFORMED', 'global_minimum_anatomical_preparation': 'UNKNOWN', 'rigorous_float_enclosure': 'MISSING', 'timescale': 'HANDOVER'}
    dump(ROOT / 'raw/R4_RESULTS.json', out)
    return out
if __name__ == '__main__':
    run()
