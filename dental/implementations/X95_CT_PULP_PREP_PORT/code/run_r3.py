import json, resource, time, warnings
import numpy as np
from scipy import sparse
from scipy.optimize import linprog, OptimizeWarning
from identity import ROOT, DATA, sha, write
from geometry import cap, rectangle_distance, measure, validate_pulp_payload

def minimum_deficit_lp(coords, valid, baseline, target):
    groups = {}
    for (i, (_, y, x)) in enumerate(coords):
        groups.setdefault((int(y), int(x)), []).append(i)
    rr = []
    cc = []
    vv = []
    rhs = []
    for group in groups.values():
        group.sort(key=lambda i: -coords[i, 0])
        for (a, b) in zip(group, group[1:]):
            row = len(rhs)
            rr.extend([row, row])
            cc.extend([b, a])
            vv.extend([1.0, -1.0])
            rhs.append(0.0)
    A = sparse.coo_matrix((vv, (rr, cc)), shape=(len(rhs), len(coords))).tocsr()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', OptimizeWarning)
        solution = linprog(-baseline.astype(float), A_ub=A, b_ub=np.array(rhs), A_eq=sparse.csr_matrix(np.ones((1, len(coords)))), b_eq=[float(target)], bounds=[(0.0, 1.0) if v else (0.0, 0.0) for v in valid], method='highs', options={'threads': 1, 'primal_feasibility_tolerance': 1e-07, 'dual_feasibility_tolerance': 1e-07})
    assert solution.success
    return float(baseline.sum() + solution.fun)

def main():
    started = time.monotonic()
    rows = []
    controls = []
    manifest = []
    frozen = json.loads((ROOT / 'FROZEN_PREDICTIONS_R3.json').read_text())
    assert sha(ROOT / 'FROZEN_PREDICTIONS_R3.json') == (ROOT / 'FROZEN_PREDICTIONS_R3.sha256').read_text().strip()
    previous = json.loads((ROOT / 'raw/R1_ROWS.json').read_text())
    cohort = json.loads((ROOT / 'FROZEN_COHORT.json').read_text())['cohort']
    for pred in frozen['predictions']:
        fdi = pred['fdi']
        old = next((r for r in previous if r['fdi'] == fdi and r['operation'] == 'X31_CAP_1.5mm'))
        a = np.load(old['spatial_export'])
        t = a['tooth']
        u = a['pulp']
        h = float(a['spacing'][0])
        base = int(old['virtual_finish_z_voxel'])
        item = next((c for c in cohort if c['canonical_fdi'] == fdi))
        validate_pulp_payload(u, item, json.loads(str(a['frame_json']))['crown_axis_sign'])
        coords = np.argwhere(t & ~cap(t, 7, base))
        pc = np.argwhere(u)
        gaps = rectangle_distance(coords, pc, np.full((len(pc), 3), 0.5), h, cube=True)
        valid = (gaps >= 1.0 - 1e-12) & ~u[tuple(coords.T)]
        baseline = (t & ~cap(t, 4, base))[tuple(coords.T)]
        groups = {}
        for (i, (_, y, x)) in enumerate(coords):
            groups.setdefault((int(y), int(x)), []).append(i)
        allowed = np.zeros(len(coords), bool)
        for group in groups.values():
            group.sort(key=lambda i: -coords[i, 0])
            for i in group:
                if not valid[i]:
                    break
                allowed[i] = True
        target = pred['target_removed_voxels']
        deficit = int((baseline & ~allowed).sum())
        assert int(allowed.sum()) == pred['max_admissible_removed_voxels'] and deficit == pred['baseline_forbidden_removed_cells']
        selected = baseline & allowed
        remaining = target - int(selected.sum())
        extras = []
        for key in sorted(groups):
            extras.extend((i for i in groups[key] if allowed[i] and (not selected[i])))
        assert len(extras) >= remaining
        selected[np.array(extras[:remaining], int)] = True
        p = t.copy()
        p[tuple(coords[selected].T)] = False
        lp = minimum_deficit_lp(coords, valid, baseline, target)
        (summary, fields) = measure(t, u, p, h)
        actual_deficit = int((baseline & ~selected).sum())
        assert actual_deficit == deficit and abs(lp - deficit) <= 1e-07
        assert summary['removed_voxels'] == target and summary['pulp_removed_voxels'] == 0
        assert summary['cut_exact_pulp_cube_min_mm'] >= 1.0 - 1e-12
        frame = json.loads(str(a['frame_json'])) | {'operation': 'R3_MINIMUM_LOCAL_CAP_DEFICIT', 'research_clearance_mm': 1.0}
        deficitmask = np.zeros(t.shape, bool)
        deficitmask[tuple(coords[baseline & ~selected].T)] = True
        path = DATA / f'R3_P48_FDI{fdi}.npz'
        np.savez_compressed(path, tooth=t, pulp=u, prepared=p, removed=t & ~p, spacing=a['spacing'], origin_zyx=a['origin_zyx'], frame_json=json.dumps(frame), uniform_cap_deficit_mask=deficitmask, candidate_removed_coords_zyx=coords, candidate_cell_pulp_cube_gap_mm=gaps, allowed_prefix_cells=allowed, **fields)
        manifest.append({'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size, 'fdi': fdi})
        row = {'case': 'P48', 'fdi': fdi, 'status': 'EXACT_DIGITAL_RESEARCH_RESERVE_WITH_LOCAL_CAP_DEBT', **summary, 'research_clearance_mm': 1.0, 'uniform_cap_deficit_voxels': deficit, 'uniform_cap_deficit_mm3': float(deficit * h ** 3), 'LP_minimum_deficit_voxels': lp, 'LP_error_voxels': abs(lp - deficit), 'uniform1p5_cut_clearance_mm': old['cut_exact_pulp_cube_min_mm'], 'same_budget_volume_hex': summary['removed_volume_mm3'].hex(), 'volume_identity_error_mm3': abs(summary['removed_volume_mm3'] - old['removed_volume_mm3']), 'downstream_difference_from_uniform_mm': summary['cut_exact_pulp_cube_min_mm'] - old['cut_exact_pulp_cube_min_mm'], 'spatial_export': str(path), 'minimum_material_space': 'UNKNOWN; locally waived4voxelcaprequirement', 'physical_anatomy': 'UNKNOWN; onlyannotatedpulp protected'}
        assert row['volume_identity_error_mm3'] == 0 and row['same_budget_volume_hex'] == old['removed_volume_mm3'].hex()
        controls.append({'check': 'minimum_deficit_plus1', 'fdi': fdi, 'valid_pass': row['LP_error_voxels'] <= 1e-07, 'injected_rejected': abs(lp - (deficit + 1)) > 1e-07})
        reread = np.load(path)
        assert np.array_equal(reread['prepared'], p) and np.array_equal(reread['uniform_cap_deficit_mask'], deficitmask)
        rows.append(row)
        print('R3 FDI', fdi, 'minimum waived cells', deficit, 'cut clearance', row['cut_exact_pulp_cube_min_mm'], flush=True)
    assert all((c['valid_pass'] and c['injected_rejected'] for c in controls))
    write(ROOT / 'raw/R3_ROWS.json', rows)
    write(ROOT / 'raw/R3_CONTROLS.json', controls)
    write(ROOT / 'raw/R3_MANIFEST.json', manifest)
    out = {'claim_type': 'capability', 'exact_budget_reserve_solids': len(rows), 'locally_waived_uniform_caps': sum((r['uniform_cap_deficit_voxels'] > 0 for r in rows)), 'mandatory_retained_cap_cells': sum((r['uniform_cap_deficit_voxels'] for r in rows)), 'minimum_cut_cube_clearance_mm': min((r['cut_exact_pulp_cube_min_mm'] for r in rows)), 'max_LP_error_voxels': max((r['LP_error_voxels'] for r in rows)), 'all_volume_identity_errors_zero': all((r['volume_identity_error_mm3'] == 0 for r in rows)), 'all_faults_rejected': all((c['injected_rejected'] for c in controls)), 'wall_s': time.monotonic() - started, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'scope': 'Virtualresearchgeometry; uniformcapruledebtretained; no material/clinical feasibility'}
    write(ROOT / 'raw/R3_OUTCOME.json', out)
    print(json.dumps(out), flush=True)
if __name__ == '__main__':
    main()
