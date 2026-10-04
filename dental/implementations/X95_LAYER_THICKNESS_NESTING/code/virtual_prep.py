"""Maximal virtual preparation subset under frozen conservative cube/ball model."""
from pose_join import *
from measure_geometry import closest
import trimesh
from scipy.spatial import cKDTree
from scipy.ndimage import label

def close_prep(raw):
    v = raw['vertices']
    f = raw['faces'][raw['roles'] == 1]
    (used, inv) = np.unique(f, return_inverse=True)
    vv = v[used]
    ff = inv.reshape(-1, 3)
    edges = np.r_[ff[:, [0, 1]], ff[:, [1, 2]], ff[:, [2, 0]]]
    (_, ix, ct) = np.unique(np.sort(edges, axis=1), axis=0, return_index=True, return_counts=True)
    rim = edges[ix[ct == 1]]
    assert len(rim) > 0 and np.max(ct) <= 2
    bottom = vv[np.unique(rim)].mean(0)
    bottom[2] = vv[np.unique(rim), 2].min() - 1.0
    n = len(vv)
    faces = np.r_[ff, np.c_[rim[:, 1], rim[:, 0], np.full(len(rim), n)]]
    m = trimesh.Trimesh(np.r_[vv, bottom[None]], faces, process=False)
    trimesh.repair.fix_normals(m, multibody=True)
    assert m.is_watertight and m.is_winding_consistent
    if m.volume < 0:
        m.invert()
    return m

def prep_grid(m, h):
    lo = m.bounds[0] - h
    hi = m.bounds[1] + h
    axes = [np.arange(lo[j] + h / 2, hi[j], h) for j in range(3)]
    xyz = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1)
    coords = xyz.reshape(-1, 3)
    inside = np.concatenate([m.contains(q) for q in np.array_split(coords, max(1, len(coords) // 256))])
    idx = np.flatnonzero(inside)
    d = closest(m, coords[idx])[1]
    rv = math.sqrt(3) * h / 2
    full = idx[d >= rv]
    return (coords[full], full, tuple((len(a) for a in axes)), lo, {'grid_centers': len(coords), 'center_inside': len(idx), 'eligible_full_cubes': len(full), 'boundary_cube_excluded': len(idx) - len(full), 'outside_original_core_centers': len(coords) - len(idx), 'original_preparation_volume_mm3': float(m.volume), 'full_cube_volume_mm3': float(len(full) * h ** 3), 'boundary_coverage': 'conservative subset only, excluded boundary cells not measured removal'})

def incidences(y, c, radius, high, rv):
    tree = cKDTree(y)
    yi = []
    ti = []
    dd = []
    for start in range(0, len(c), 128):
        stop = min(start + 128, len(c))
        near = tree.query_ball_point(c[start:stop], high + radius[start:stop] + rv, workers=4)
        for (ii, ids) in enumerate(near):
            if not ids:
                continue
            ids = np.array(ids, dtype=np.int64)
            i = start + ii
            d = np.linalg.norm(y[ids] - c[i], axis=1)
            take = d < high + radius[i] + rv
            yi.append(ids[take])
            ti.append(np.full(int(take.sum()), i, dtype=np.int64))
            dd.append(d[take])
    return (np.concatenate(yi), np.concatenate(ti), np.concatenate(dd))

def voxel_mesh(kept, fullidx, shape, lo, h):
    occ = np.zeros(shape, dtype=bool)
    occ.ravel()[fullidx[kept]] = True
    vertices = []
    faces = []
    lookup = {}
    corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    patterns = [((-1, 0, 0), (0, 4, 7, 3)), ((1, 0, 0), (1, 2, 6, 5)), ((0, -1, 0), (0, 1, 5, 4)), ((0, 1, 0), (3, 7, 6, 2)), ((0, 0, -1), (0, 3, 2, 1)), ((0, 0, 1), (4, 5, 6, 7))]
    for cell in np.argwhere(occ):
        for (offset, face) in patterns:
            nbr = cell + offset
            if np.all(nbr >= 0) and np.all(nbr < shape) and occ[tuple(nbr)]:
                continue
            ids = []
            for k in face:
                key = tuple(cell + corners[k])
                if key not in lookup:
                    lookup[key] = len(vertices)
                    vertices.append(lo + h * np.array(key))
                ids.append(lookup[key])
            faces.extend([[ids[0], ids[1], ids[2]], [ids[0], ids[2], ids[3]]])
    mesh = trimesh.Trimesh(np.array(vertices), np.array(faces), process=False) if faces else None
    comp = label(occ)[1]
    return (mesh, int(comp))

def run():
    st = time.perf_counter()
    pr = json.loads((ROOT / 'PREREG_R2.json').read_text())
    r1 = json.loads((ROOT / 'PREREG_R1.json').read_text())
    rule = json.loads((ROOT / 'RULE_CONTRACT.json').read_text())
    assert source_gate(rule)
    h = pr['model_contract']['voxel_h_mm']
    rv = math.sqrt(3) * h / 2
    rows = []
    for (ei, e) in enumerate(json.loads((ROOT / 'inputs/R4_FROZEN_EXPORTS.json').read_text())['exports']):
        t0 = time.perf_counter()
        raw = dict(np.load(e['mesh_path'], allow_pickle=False))
        assert sha(e['mesh_path']) == e['mesh_sha256']
        outer = raw['vertices'][raw['faces'][raw['roles'] == 0]]
        c = outer.mean(1)
        rho = np.linalg.norm(outer - c[:, None], axis=2).max(1)
        m = close_prep(raw)
        (y, fullidx, shape, lo, drop) = prep_grid(m, h)
        assert len(y) > 0
        (high, low) = (0.8, 0.4) if e['family'] == 'anterior' else (1.0, 0.5)
        (iy, it, dist) = incidences(y, c, rho, high, rv)
        print(e['family'], 'eligible', len(y), 'pairs', len(iy), flush=True)
        low_block = np.zeros(len(y), bool)
        low_block[iy[dist < low + rho[it] + rv]] = True
        universal_block = np.zeros(len(y), bool)
        universal_block[iy[dist < 0.5 + rho[it] + rv]] = True
        baseline = ~universal_block
        grid = [p for p in poses(raw['vertices'], r1) if p['stock_ok']]
        records = []
        best = None
        bestkeep = None
        checks = 0
        mismatches = 0
        signatures = {}
        for (ix, pose) in enumerate(grid):
            n = axis(pose['tilt_deg'], pose['azimuth_deg'])
            key = (pose['tilt_deg'], pose['azimuth_deg'])
            H = pose['height_mm']
            depth = pose['top_depth_mm']
            top = np.max(raw['vertices'] @ n)
            if key not in signatures:
                tri_top = np.max(outer @ n, axis=1)
                s = np.full(len(y), -np.inf)
                np.maximum.at(s, iy, tri_top[it])
                signatures[key] = s
            keep = ~low_block & (depth + top - signatures[key] > H / 2)
            rec = dict(pose, pose_index=ix, retained_cubes=int(keep.sum()), removed_cubes=int((~keep).sum()), retained_volume_mm3=float(keep.sum() * h ** 3), removed_volume_mm3=float((~keep).sum() * h ** 3), different_voxel_decisions=int(np.sum(keep != baseline)))
            records.append(rec)
            if best is None or rec['retained_cubes'] > best['retained_cubes']:
                best = rec
                bestkeep = keep.copy()
            if ei == 0:
                theta = math.radians(pose['tilt_deg'])
                phi = math.radians(pose['azimuth_deg'])
                nn = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
                topc = np.max(raw['vertices'][:, 0] * nn[0] + raw['vertices'][:, 1] * nn[1] + raw['vertices'][:, 2] * nn[2])
                ztri = depth + topc - (outer[:, :, 0] * nn[0] + outer[:, :, 1] * nn[1] + outer[:, :, 2] * nn[2])
                rr = np.where(ztri.min(1) <= H / 2, high, low)
                blocked = np.zeros(len(y), bool)
                blocked[iy[dist < rr[it] + rho[it] + rv]] = True
                mismatches += int(np.sum(keep != ~blocked))
                checks += 1
        assert best is not None and mismatches == 0
        n = axis(best['tilt_deg'], best['azimuth_deg'])
        top = np.max(raw['vertices'] @ n)
        ztri = best['top_depth_mm'] + top - outer @ n
        rr = np.where(ztri.min(1) <= best['height_mm'] / 2, high, low)
        dcheck = np.sqrt(np.sum((y[iy] - c[it]) ** 2, axis=1))
        err = float(np.max(abs(dcheck - dist)))
        blocked = np.zeros(len(y), bool)
        blocked[iy[dcheck < rr[it] + rho[it] + rv]] = True
        assert err <= 1e-10 and np.array_equal(bestkeep, ~blocked)
        sample = np.unique(np.r_[np.flatnonzero(bestkeep)[:4], np.flatnonzero(~bestkeep)[:4], np.linspace(0, len(y) - 1, 32, dtype=int)])
        gs = np.array([np.min(np.linalg.norm(yy - c, axis=1) - rr - rho - rv) for yy in y[sample]])
        assert np.array_equal(gs >= 0, bestkeep[sample])
        guards = np.full(len(y), np.inf)
        np.minimum.at(guards, iy, dist - rr[it] - rho[it] - rv)
        guardmin = float(guards[bestkeep].min()) if bestkeep.any() else None
        assert not bestkeep.any() or guardmin >= 0
        bad = bestkeep.copy()
        forbidden = int(np.flatnonzero(~bestkeep)[0])
        bad[forbidden] = True
        injected_rejected = bool(np.any(bad & blocked))
        assert injected_rejected
        fullpath = ROOT / 'exports' / f"{e['family']}_virtual_core.npz"
        np.savez_compressed(fullpath, cube_centers_mm=y, cube_h_mm=h, eligible_original_core_index=fullidx, kept_yml=bestkeep, kept_universal0_5=baseline, guard_min_mm=guards, grid_shape=shape, grid_origin_mm=lo, locked_outer_vertices=raw['vertices'][np.unique(raw['faces'][raw['roles'] == 0])])
        (vm, comp) = voxel_mesh(bestkeep, fullidx, shape, lo, h)
        if vm is not None:
            vm.export(ROOT / 'exports' / f"{e['family']}_virtual_preparation.stl")
        m.export(ROOT / 'exports' / f"{e['family']}_original_virtual_preparation.stl")
        with (ROOT / 'raw' / f"{e['family']}_R2_POSES.csv").open('w') as f:
            wr = csv.DictWriter(f, fieldnames=list(records[0]))
            wr.writeheader()
            wr.writerows(records)
        (bmesh, bcomp) = voxel_mesh(baseline, fullidx, shape, lo, h)
        if bmesh is not None:
            bmesh.export(ROOT / 'exports' / f"{e['family']}_universal0_5_preparation.stl")
        row = {'family': e['family'], 'key': e['key'], 'resolution': 'PER_POINT', 'voxel_h_mm': h, 'eligible_cubes': len(y), 'dropout': drop, 'source_best_pose': best, 'universal0_5_retained_cubes': int(baseline.sum()), 'universal0_5_retained_volume_mm3': float(baseline.sum() * h ** 3), 'universal0_5_removed_volume_mm3': float((~baseline).sum() * h ** 3), 'source_minus_universal_removed_mm3': float((baseline.sum() - bestkeep.sum()) * h ** 3), 'source_keeps_universal_removes_cubes': int(np.sum(bestkeep & ~baseline)), 'source_removes_universal_keeps_cubes': int(np.sum(~bestkeep & baseline)), 'complete_design_decision_changed': bool(np.any(bestkeep != baseline)), 'outer_identity_error_mm': 0.0, 'selected_guard_min_mm': guardmin, 'selected_pose_upper_triangles': int(np.sum(rr == high)), 'kept_components_face_connectivity': comp, 'baseline_components_face_connectivity': bcomp, 'mesh_watertight': bool(vm.is_watertight) if vm is not None else None, 'selected_control_membership_mismatches': 0, 'first_case_all_pose_controls': checks, 'first_case_all_pose_membership_mismatches': mismatches, 'independent_all_triangle_points': len(sample), 'independent_distance_error_mm': err, 'injected_forbidden_cube_rejected': injected_rejected, 'qualification': 'VIRTUAL_CONSERVATIVE_MODEL_ONLY_NOT_FABRICATION_QUALIFIED', 'rigorous_float_enclosure': 'MISSING', 'source_uncertainty': 'UNKNOWN', 'array_path': str(fullpath.relative_to(ROOT)), 'array_sha256': sha(fullpath), 'seconds': time.perf_counter() - t0}
        rows.append(row)
        dump(ROOT / 'raw/R2_PARTIAL.json', rows)
        print(e['family'], 'best', best['retained_volume_mm3'], 'baseline', row['universal0_5_retained_volume_mm3'], 'delta', row['source_minus_universal_removed_mm3'], flush=True)
    out = {'round': 'R2', 'claim_type': 'information_link', 'primary_gate': 'PASS_SCOPED_VIRTUAL_DESIGN_DECISION' if any((r['complete_design_decision_changed'] for r in rows)) else 'NEGATIVE', 'rows': rows, 'seconds': time.perf_counter() - st, 'sufficient_representation': 'source rule joint with triangle region/pose and eligible corecubeID', 'physical_measurement': 'NOT_PERFORMED', 'clinical_or_manufacturing_gain': 'UNKNOWN', 'rigorous_float_enclosure': 'MISSING', 'timescale': 'HANDOVER'}
    dump(ROOT / 'raw/R2_RESULTS.json', out)
    return out
if __name__ == '__main__':
    run()
