"""Exact inverse prefix front for declared research-only removed-cell reserves."""
import json, resource, time, warnings
import numpy as np
from scipy import sparse
from scipy.optimize import linprog, OptimizeWarning
from identity import ROOT, DATA, sha, write, source_path
from geometry import cap, rectangle_distance, measure, validate_port, validate_pulp_payload

def independent_lp(coords, valid, baseline, target):
    n = len(coords)
    groups = {}
    for (i, (z, y, x)) in enumerate(coords):
        groups.setdefault((int(y), int(x)), []).append(i)
    rr = []
    cc = []
    vv = []
    rhs = []
    for group in groups.values():
        group.sort(key=lambda i: -coords[i, 0])
        for (a, b) in zip(group, group[1:]):
            k = len(rhs)
            rr.extend([k, k])
            cc.extend([b, a])
            vv.extend([1.0, -1.0])
            rhs.append(0.0)
    for i in np.flatnonzero(~valid):
        k = len(rhs)
        rr.append(k)
        cc.append(int(i))
        vv.append(1.0)
        rhs.append(0.0)
    A = sparse.coo_matrix((vv, (rr, cc)), shape=(len(rhs), n)).tocsr()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', OptimizeWarning)
        opts = {'threads': 1, 'primal_feasibility_tolerance': 1e-07, 'dual_feasibility_tolerance': 1e-07}
        maximum = linprog(-np.ones(n), A_ub=A, b_ub=np.array(rhs), bounds=(0.0, 1.0), method='highs', options=opts)
        bounds = [(1.0, 1.0) if b else (0.0, 1.0) for b in baseline]
        feasible = linprog(np.zeros(n), A_ub=A, b_ub=np.array(rhs), A_eq=sparse.csr_matrix(np.ones((1, n))), b_eq=[float(target)], bounds=bounds, method='highs', options=opts)
    assert maximum.success and feasible.status in [0, 2]
    return (-float(maximum.fun), bool(feasible.success), float(np.max(np.abs(maximum.x - np.rint(maximum.x)))))

def main():
    started = time.monotonic()
    rows = []
    controls = []
    manifest = []
    assert sha(ROOT / 'FROZEN_PREDICTIONS_R2.json') == (ROOT / 'FROZEN_PREDICTIONS_R2.sha256').read_text().strip()
    r1 = json.loads((ROOT / 'raw/R1_ROWS.json').read_text())
    cohort = json.loads((ROOT / 'FROZEN_COHORT.json').read_text())
    for item in cohort['cohort']:
        fdi = item['canonical_fdi']
        if item['status'] != 'AVAILABLE':
            continue
        previous = next((r for r in r1 if r['fdi'] == fdi and r['operation'] == 'X31_CAP_1.5mm'))
        a = np.load(previous['spatial_export'])
        t = a['tooth']
        u = a['pulp']
        h = float(a['spacing'][0])
        frame = json.loads(str(a['frame_json']))
        validate_port('P48', sha(source_path(item['source_path'])), fdi, item['source_pulpy_fdi'], a['spacing'], item)
        validate_pulp_payload(u, item, frame['crown_axis_sign'])
        base = int(previous['virtual_finish_z_voxel'])
        baseline_p = cap(t, 4, base)
        max_p = cap(t, 7, base)
        coords = np.argwhere(t & ~max_p)
        pcoords = np.argwhere(u)
        gaps = rectangle_distance(coords, pcoords, np.full((len(pcoords), 3), 0.5), h, cube=True)
        baseline = (t & ~baseline_p)[tuple(coords.T)]
        target = int(previous['removed_voxels'])
        groups = {}
        for (i, (z, y, x)) in enumerate(coords):
            groups.setdefault((int(y), int(x)), []).append(i)
        for group in groups.values():
            group.sort(key=lambda i: -coords[i, 0])
        for b in [0.0, 0.25, 1.0]:
            valid = (gaps >= b - 1e-12) & ~u[tuple(coords.T)]
            allowed = np.zeros(len(coords), bool)
            cap_yx = np.zeros(t.shape[1:], np.int16)
            for (key, group) in groups.items():
                count = 0
                for i in group:
                    if not valid[i]:
                        break
                    allowed[i] = True
                    count += 1
                cap_yx[key] = count
            capacity = int(allowed.sum())
            conflict = bool(np.any(baseline & ~allowed))
            candidate_feasible = not conflict and capacity >= target
            (lpmax, lpfeasible, integrality) = independent_lp(coords, valid, baseline, target)
            capacity_error = abs(lpmax - capacity)
            assert capacity_error <= 1e-07 and candidate_feasible == lpfeasible and (integrality <= 1e-07)
            status = 'BASELINE_PREFIX_CONFLICT' if conflict else 'EXACT_BUDGET_FEASIBLE' if candidate_feasible else 'FINITE_PREFIX_CAPACITY_OBSTRUCTION'
            row = {'case': 'P48', 'fdi': fdi, 'research_clearance_mm': b, 'status': status, 'target_removed_voxels': target, 'max_admissible_removed_voxels': capacity, 'baseline_forbidden_removed_cells': int((baseline & ~allowed).sum()), 'LP_capacity_voxels': lpmax, 'LP_capacity_error_voxels': capacity_error, 'LP_target_feasible': lpfeasible, 'LP_maximum_integrality_error': integrality, 'candidate_feasible': candidate_feasible, 'physical_reserve': 'UNKNOWN', 'dentin': 'UNKNOWN_NO_DEJ', 'scope': 'exact digital removed-cell-union reserve in this frozen column-prefix family; no general crown impossibility'}
            if candidate_feasible:
                selected = baseline.copy()
                remaining = target - int(selected.sum())
                extras = []
                for key in sorted(groups):
                    extras.extend((i for i in groups[key] if allowed[i] and (not selected[i])))
                selected[np.array(extras[:remaining], int)] = True
                p = t.copy()
                p[tuple(coords[selected].T)] = False
                (summary, fields) = measure(t, u, p, h)
                minimum = float(gaps[selected].min())
                assert int((t & ~p).sum()) == target and (not (u & ~p).any())
                assert minimum >= b - 1e-12 and summary['cut_exact_pulp_cube_min_mm'] >= b - 1e-12
                row.update(summary)
                row['removed_cell_union_exact_gap_mm'] = minimum
                row['count_identity_error'] = abs(row['removed_voxels'] - target)
                row['volume_identity_error_mm3'] = abs(row['removed_volume_mm3'] - previous['removed_volume_mm3'])
                row['volume_bit_identity'] = row['removed_volume_mm3'].hex() == previous['removed_volume_mm3'].hex()
                row['cut_clearance_change_from_uniform1p5_mm'] = summary['cut_exact_pulp_cube_min_mm'] - previous['cut_exact_pulp_cube_min_mm']
                assert row['volume_bit_identity'] and row['volume_identity_error_mm3'] == 0
                frame = frame | {'operation': 'R2_INVERSE_PREFIX', 'research_clearance_mm': b}
                path = DATA / f'R2_P48_FDI{fdi}_b{b:g}.npz'
                np.savez_compressed(path, tooth=t, pulp=u, prepared=p, removed=t & ~p, spacing=np.array([h] * 3), origin_zyx=a['origin_zyx'], frame_json=json.dumps(frame), capacity_yx=cap_yx, candidate_removed_coords_zyx=coords, candidate_cell_pulp_cube_gap_mm=gaps, allowed_prefix_cells=allowed, **fields)
                reread = np.load(path)
                assert np.array_equal(reread['prepared'], p) and json.loads(str(reread['frame_json'])) == frame
                row['spatial_export'] = str(path)
                manifest.append({'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size, 'fdi': fdi, 'research_clearance_mm': b})
            controls.append({'check': 'capacity_plus1', 'fdi': fdi, 'b_mm': b, 'valid_pass': capacity_error <= 1e-07, 'injected_rejected': abs(lpmax - (capacity + 1)) > 1e-07})
            forbidden = np.flatnonzero(~allowed)
            if len(forbidden):
                bad = forbidden[0]
                trial = allowed.copy()
                trial[bad] = True
                intrinsic = not valid[bad]
                precedent = False
                group = groups[int(coords[bad, 1]), int(coords[bad, 2])]
                pos = group.index(int(bad))
                if pos:
                    precedent = not trial[group[pos - 1]]
                controls.append({'check': 'forbidden_prefix_cell', 'fdi': fdi, 'b_mm': b, 'valid_pass': True, 'injected_rejected': intrinsic or precedent, 'bad_cell_zyx': coords[bad].tolist(), 'cell_gap_mm': float(gaps[bad]), 'predecessor_missing': precedent})
            rows.append(row)
        print('inverse FDI', fdi, [(r['research_clearance_mm'], r['status']) for r in rows if r['fdi'] == fdi], flush=True)
        write(ROOT / 'raw/R2_ROWS.json', rows)
        write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': 'X95-ct-pulp-prep-port', 'milestone': 'INVERSE_FRONT_RUNNING', 'latest_fdi': fdi, 'latest_gate': 'prefix capacity matches LP; actual feasible solids checked', 'next_operation': 'finish frozen inverse cohort and STL export'})
    assert all((c['valid_pass'] and c['injected_rejected'] for c in controls))
    write(ROOT / 'raw/R2_CONTROLS.json', controls)
    write(ROOT / 'raw/R2_MANIFEST.json', manifest)
    result = {'claim_type': 'capability', 'queries': len(rows), 'feasible': sum((r['candidate_feasible'] for r in rows)), 'research1mm_feasible': sum((r['candidate_feasible'] for r in rows if r['research_clearance_mm'] == 1)), 'baseline_conflicts': sum((r['status'] == 'BASELINE_PREFIX_CONFLICT' for r in rows)), 'capacity_obstructions': sum((r['status'] == 'FINITE_PREFIX_CAPACITY_OBSTRUCTION' for r in rows)), 'max_LP_capacity_error_voxels': max((r['LP_capacity_error_voxels'] for r in rows)), 'max_integrality_error': max((r['LP_maximum_integrality_error'] for r in rows)), 'all_faults_rejected': all((c['injected_rejected'] for c in controls)), 'wall_s': time.monotonic() - started, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'physical_clinical_scope': 'UNKNOWN; research b is not a clinical threshold'}
    write(ROOT / 'raw/R2_OUTCOME.json', result)
    print(json.dumps(result), flush=True)
if __name__ == '__main__':
    main()
