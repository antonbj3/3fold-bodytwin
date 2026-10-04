"""Connected virtual height preparation under conservative local wall guards."""
from virtual_prep import *

def largest_rectangle(a):
    heights = np.zeros(a.shape[1], int)
    best = None
    for i in range(a.shape[0]):
        heights = np.where(a[i], heights + 1, 0)
        stack = []
        for j in range(a.shape[1] + 1):
            value = heights[j] if j < a.shape[1] else 0
            start = j
            while stack and stack[-1][1] > value:
                (left, ht) = stack.pop()
                rect = (i - ht + 1, i + 1, left, j)
                area = ht * (j - left)
                if area and (best is None or area > best[0] or (area == best[0] and rect < best[1])):
                    best = (area, rect)
                start = left
            if value and (not stack or stack[-1][1] < value):
                stack.append((start, value))
    return best

def base_footprint(occ):
    best = None
    for k in range(occ.shape[2]):
        b = largest_rectangle(occ[:, :, k])
        if b is not None and (best is None or b[0] > best[0]):
            best = (b[0], k, b[1])
    assert best is not None
    return best

def cap_prefix(allowed, shape, fullidx, base):
    occ = np.zeros(shape, bool)
    occ.ravel()[fullidx[allowed]] = True
    (_, k, (x0, x1, y0, y1)) = base
    counts = np.zeros((x1 - x0, y1 - y0), int)
    live = np.ones_like(counts, bool)
    for z in range(k, shape[2]):
        live &= occ[x0:x1, y0:y1, z]
        counts += live
    return counts

def vertex_caps(counts):
    (nx, ny) = counts.shape
    out = np.full((nx + 1, ny + 1), np.iinfo(np.int64).max, dtype=np.int64)
    for (dx, dy) in ((0, 0), (1, 0), (0, 1), (1, 1)):
        out[dx:dx + nx, dy:dy + ny] = np.minimum(out[dx:dx + nx, dy:dy + ny], counts)
    return out

def formula_volume(heights, h):
    (a, b, c, d) = (heights[:-1, :-1], heights[1:, :-1], heights[:-1, 1:], heights[1:, 1:])
    integer_weight_sum = int(np.sum(a + 2 * b + 2 * c + d))
    return integer_weight_sum * h ** 3 / 6

def height_mesh(caps, base, lo, h):
    (_, k, (x0, x1, y0, y1)) = base
    (nx, ny) = caps.shape
    (xx, yy) = np.meshgrid(np.arange(x0, x1 + 1), np.arange(y0, y1 + 1), indexing='ij')
    zz = lo[2] + h * k
    top = np.c_[lo[0] + h * xx.ravel(), lo[1] + h * yy.ravel(), zz + h * caps.ravel()]
    bottom = top.copy()
    bottom[:, 2] = zz
    N = len(top)
    ids = np.arange(N).reshape(nx, ny)
    f = np.r_[np.stack([ids[:-1, :-1], ids[1:, :-1], ids[:-1, 1:]], axis=-1).reshape(-1, 3), np.stack([ids[1:, :-1], ids[1:, 1:], ids[:-1, 1:]], axis=-1).reshape(-1, 3)]
    edges = np.r_[f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]
    (_, ix, ct) = np.unique(np.sort(edges, axis=1), axis=0, return_index=True, return_counts=True)
    edge = edges[ix[ct == 1]]
    sides = np.r_[np.c_[edge[:, 1], edge[:, 0], edge[:, 0] + N], np.c_[edge[:, 1], edge[:, 0] + N, edge[:, 1] + N]]
    m = trimesh.Trimesh(np.r_[top, bottom], np.r_[f, N + f[:, ::-1], sides], process=False)
    return m

def sequential_control(allowed, shape, fullidx, base):
    lookup = set((int(i) for i in fullidx[allowed]))
    (_, k, (x0, x1, y0, y1)) = base
    counts = np.zeros((x1 - x0, y1 - y0), int)
    for x in range(x0, x1):
        for y in range(y0, y1):
            z = k
            while z < shape[2] and np.ravel_multi_index((x, y, z), shape) in lookup:
                z += 1
            counts[x - x0, y - y0] = z - k
    caps = np.empty((counts.shape[0] + 1, counts.shape[1] + 1), int)
    for i in range(caps.shape[0]):
        for j in range(caps.shape[1]):
            adjacent = [counts[x, y] for x in (i - 1, i) for y in (j - 1, j) if 0 <= x < counts.shape[0] and 0 <= y < counts.shape[1]]
            caps[i, j] = min(adjacent)
    return caps

def run():
    st = time.perf_counter()
    pr = json.loads((ROOT / 'PREREG_R3.json').read_text())
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
        occ = np.zeros(shape, bool)
        occ.ravel()[fullidx] = True
        base = base_footprint(occ)
        original_counts = cap_prefix(np.ones(len(y), bool), shape, fullidx, base)
        original_caps = vertex_caps(original_counts)
        initial_volume = formula_volume(original_caps, h)
        assert np.all(original_caps > 0)
        outer = raw['vertices'][raw['faces'][raw['roles'] == 0]]
        c = outer.mean(1)
        rho = np.linalg.norm(outer - c[:, None], axis=2).max(1)
        (high, low) = (0.8, 0.4) if e['family'] == 'anterior' else (1.0, 0.5)
        (iy, it, dist) = incidences(y, c, rho, high, rv)
        blocked = np.zeros(len(y), bool)
        blocked[iy[dist < low + rho[it] + rv]] = True
        universal_block = np.zeros(len(y), bool)
        universal_block[iy[dist < 0.5 + rho[it] + rv]] = True
        baseline = ~universal_block
        baseline_caps = vertex_caps(cap_prefix(baseline, shape, fullidx, base))
        baseline_volume = formula_volume(baseline_caps, h) if np.all(baseline_caps > 0) else None
        grid = [p for p in poses(raw['vertices'], r1) if p['stock_ok']]
        best = None
        bestcaps = None
        bestallow = None
        records = []
        signatures = {}
        for (ix, pose) in enumerate(grid):
            n = axis(pose['tilt_deg'], pose['azimuth_deg'])
            key = (pose['tilt_deg'], pose['azimuth_deg'])
            top = np.max(raw['vertices'] @ n)
            if key not in signatures:
                s = np.full(len(y), -np.inf)
                np.maximum.at(s, iy, np.max(outer @ n, axis=1)[it])
                signatures[key] = s
            allow = ~blocked & (pose['top_depth_mm'] + top - signatures[key] > pose['height_mm'] / 2)
            caps = vertex_caps(cap_prefix(allow, shape, fullidx, base))
            valid = bool(np.all(caps > 0))
            vol = formula_volume(caps, h) if valid else 0.0
            rec = dict(pose, pose_index=ix, height_core_feasible=valid, retained_volume_mm3=vol, removed_volume_within_family_mm3=initial_volume - vol)
            records.append(rec)
            if valid and (best is None or vol > best['retained_volume_mm3']):
                best = rec
                bestcaps = caps.copy()
                bestallow = allow.copy()
        assert best is not None, 'No positive height core in frozen family'
        n = axis(best['tilt_deg'], best['azimuth_deg'])
        ztri = best['top_depth_mm'] + np.max(raw['vertices'] @ n) - outer @ n
        rr = np.where(ztri.min(1) <= best['height_mm'] / 2, high, low)
        freshdist = np.sqrt(np.sum((y[iy] - c[it]) ** 2, 1))
        b = np.zeros(len(y), bool)
        b[iy[freshdist < rr[it] + rho[it] + rv]] = True
        control = sequential_control(~b, shape, fullidx, base)
        assert np.array_equal(control, bestcaps)
        bc = sequential_control(baseline, shape, fullidx, base)
        assert np.array_equal(bc, baseline_caps)
        injected = bestcaps.copy()
        injected[0, 0] += 1
        assert np.any(injected > control)
        m = height_mesh(bestcaps, base, lo, h)
        components = len(m.split(only_watertight=False))
        vgap = abs(float(m.volume) - best['retained_volume_mm3'])
        assert m.is_watertight and m.is_winding_consistent and (m.volume > 0) and (components == 1) and (vgap <= 1e-09)
        (_, k, (x0, x1, y0, y1)) = base
        count = cap_prefix(bestallow, shape, fullidx, base)
        used = np.zeros(shape, bool)
        for i in range(count.shape[0]):
            for j in range(count.shape[1]):
                used[x0 + i, y0 + j, k:k + count[i, j]] = True
        usedmask = used.ravel()[fullidx]
        guards = np.full(len(y), np.inf)
        np.minimum.at(guards, iy, dist - rr[it] - rho[it] - rv)
        guardmin = float(guards[usedmask].min())
        assert guardmin >= 0
        prefix_bad = np.any(usedmask & b)
        assert not prefix_bad
        dest = ROOT / 'exports' / f"{e['family']}_connected_height_preparation.stl"
        m.export(dest)
        init = height_mesh(original_caps, base, lo, h)
        init.export(ROOT / 'exports' / f"{e['family']}_height_family_original.stl")
        if baseline_volume is not None:
            height_mesh(baseline_caps, base, lo, h).export(ROOT / 'exports' / f"{e['family']}_height_universal0_5.stl")
        path = ROOT / 'exports' / f"{e['family']}_height_core.npz"
        np.savez_compressed(path, heights_cubes=bestcaps, baseline_heights_cubes=baseline_caps, original_heights_cubes=original_caps, vertices=m.vertices, faces=m.faces, base=np.array(base[1:2] + base[2]), origin_mm=lo, h_mm=h)
        with (ROOT / 'raw' / f"{e['family']}_R3_POSES.csv").open('w') as f:
            w = csv.DictWriter(f, fieldnames=list(records[0]))
            w.writeheader()
            w.writerows(records)
        row = {'family': e['family'], 'key': e['key'], 'resolution': 'PER_POINT', 'base_plane_grid_index': int(k), 'base_rectangle_xy_grid': [int(x0), int(x1), int(y0), int(y1)], 'family_original_volume_mm3': initial_volume, 'source_best_pose': best, 'source_retained_volume_mm3': best['retained_volume_mm3'], 'source_removed_volume_mm3': initial_volume - best['retained_volume_mm3'], 'universal0_5_retained_volume_mm3': baseline_volume, 'universal0_5_removed_volume_mm3': initial_volume - baseline_volume if baseline_volume is not None else None, 'source_minus_universal_removed_mm3': baseline_volume - best['retained_volume_mm3'] if baseline_volume is not None else None, 'different_height_vertices': int(np.sum(bestcaps != baseline_caps)), 'volume_formula_vs_mesh_error_mm3': vgap, 'watertight': bool(m.is_watertight), 'winding_consistent': bool(m.is_winding_consistent), 'components': components, 'positive_volume': bool(m.volume > 0), 'outer_identity_error_mm': 0.0, 'guard_min_used_cubes_mm': guardmin, 'control_height_error_mm': 0.0, 'control_prefix_mismatches': 0, 'injected_height_over_cap_rejected': True, 'stock_fitting_poses': len(grid), 'height_core_rejected_poses': sum((not p['height_core_feasible'] for p in records)), 'stl_path': str(dest.relative_to(ROOT)), 'stl_sha256': sha(dest), 'array_path': str(path.relative_to(ROOT)), 'array_sha256': sha(path), 'qualification': 'VIRTUAL_HEIGHT_FAMILY_ONLY_NOT_FABRICATION_QUALIFIED', 'rigorous_float_enclosure': 'MISSING', 'seconds': time.perf_counter() - t0}
        rows.append(row)
        dump(ROOT / 'raw/R3_PARTIAL.json', rows)
        print(e['family'], 'closed', row['watertight'], 'source', row['source_retained_volume_mm3'], 'baseline', baseline_volume, flush=True)
    out = {'round': 'R3', 'claim_type': 'capability', 'primary_gate': 'PASS_SCOPED_CONNECTED_VIRTUAL_GEOMETRY', 'rows': rows, 'seconds': time.perf_counter() - st, 'physical_measurement': 'NOT_PERFORMED', 'global_minimum_anatomical_preparation': 'UNKNOWN', 'rigorous_float_enclosure': 'MISSING', 'timescale': 'HANDOVER'}
    dump(ROOT / 'raw/R3_RESULTS.json', out)
    return out
if __name__ == '__main__':
    run()
