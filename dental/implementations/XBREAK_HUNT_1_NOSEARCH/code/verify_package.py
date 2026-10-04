"""Recheck raw artifact contracts, independent R2 quadrature control, and hashes."""
from dental_release.paths import expand as _release_expand
import hashlib, json, sys, time
from pathlib import Path
import numpy as np
from run_r1 import ROOT, write, sha

def main():
    start = time.perf_counter()
    checks = []
    rounds = {}
    for number in range(1, 7):
        name = f'R{number}'
        p = ROOT / f'PREREG_{name}.json'
        r = ROOT / f'rounds/{name}/results.json'
        frozen = json.loads((ROOT / f'raw/PREREG_{name}_FROZEN.json').read_text())
        data = json.loads(r.read_text())
        rounds[name] = data
        assert sha(p) == frozen['sha256'] == data['prereg_sha256']
        rows = data['records']
        assert all((x['lower_N'] - 1e-05 <= x['facit_N'] <= x['upper_N'] + 1e-05 for x in rows))
        assert all((x['lower_N'] <= x['upper_N'] + 1e-05 for x in rows))
        checks.append({'check': f'{name}_frozen_hash_and_raw_coverage', 'n_queries': len(rows), 'status': 'PASS'})
    ramps = json.loads((ROOT / 'raw/R2_RAMPS.json').read_text())
    index = {(r['case'], r['scenario'], r['ramp']): r for r in ramps}
    maxwork = 0.0
    maxforce = 0.0
    for row in rounds['R2']['records']:
        n = row['intervals']
        eta = row['error_mm']
        T = 205.0
        ix = np.arange(0, 513, 512 // n)
        candidate = []
        keys = [f"{row['FDI']}_{row['h_mm']}_-1", 'base', f"{row['FDI']}_{row['h_mm']}_1"]
        for (j, key) in enumerate(keys):
            r = index[row['case'], row['scenario'], key]
            ts = np.array(r['force_N'])[ix]
            dd = np.array(r['d_mm'])[ix]
            obs = dd + eta * (0.6 * np.sin(np.arange(len(ix)) * 0.73 + j) + 0.4 * np.cos(j))
            dt = np.diff(ts)
            center = float(dt @ ((obs[:-1] + obs[1:]) / 2))
            radius = float(0.5 * dt @ np.diff(obs)) + T * (eta + 1e-07)
            conventional = [center - radius, center + radius]
            candidate.append(conventional)
            maxwork = max(maxwork, max((abs(a - b) for (a, b) in zip(conventional, row['work_bounds_Nmm'][j]))))
        lo = max(0.0, (candidate[2][0] - candidate[1][1]) / row['h_mm'])
        hi = min(T, (candidate[1][1] - candidate[0][0]) / row['h_mm'])
        maxforce = max(maxforce, abs(lo - row['lower_N']), abs(hi - row['upper_N']))
    assert maxwork <= 1e-07 and maxforce <= 1e-07
    checks.append({'check': 'R2_independent_conventional_quadrature_control', 'status': 'PASS', 'n': len(rounds['R2']['records']), 'max_work_difference_Nmm': maxwork, 'max_force_difference_N': maxforce})
    assert rounds['R3']['LP_parity_gate']
    assert max((r['same_information_work_control_max_difference_Nmm'] for r in rounds['R4']['records'])) <= 1e-07
    assert rounds['R4']['max_LP_difference_N'] <= 1e-07
    assert rounds['R5']['max_same_information_LP_difference_N'] <= 1e-07
    assert rounds['R6']['same_information_control_gate'] and rounds['R6']['regime_guard_gate']
    checks.append({'check': 'R3_R6_actual_strong_controls_and_UNKNOWN_guards', 'status': 'PASS'})
    sources = {}
    for w in json.loads((ROOT / 'raw/R1_WORLDS.json').read_text()):
        assert sha(w['source_path']) == w['source_sha256']
        sources[w['source_path']] = w['source_sha256']
    for path in [_release_expand('@DENTAL_IMPLEMENTATIONS@/FALT_TANDLAST/code/model.py'), _release_expand('@DENTAL_IMPLEMENTATIONS@/FALT_TANDLAST/engine_snapshot/src/motion_engine/ncp/gap_box.py'), _release_expand('@DENTAL_IMPLEMENTATIONS@/X2/raw/EXTERNAL_REFERENTS.json'), _release_expand('@DENTAL_INPUT_ROOT@/artifacts/LANE_XREVIEW_BATCH18/REVIEW_SOL_FALT_TANDLAST_20261003.json')]:
        if Path(path).is_file():
            sources[path] = sha(path)
        else:
            raise FileNotFoundError(path)
    data = rounds['R4']['data_manifest']
    for a in data:
        assert Path(a['path']).stat().st_size == a['bytes'] and sha(a['path']) == a['sha256']
        with np.load(a['path']) as z:
            assert all((z[k].dtype == np.float64 and np.all(np.isfinite(z[k])) for k in z.files))
    data_bytes = sum((a['bytes'] for a in data))
    lane_bytes = sum((p.stat().st_size for p in ROOT.rglob('*') if p.is_file()))
    assert data_bytes + lane_bytes < 3000000000
    checks.append({'check': 'source_and_data_hashes_disk_budget', 'status': 'PASS', 'external_array_bytes': data_bytes, 'lane_bytes': lane_bytes})
    result = {'status': 'PASS', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'checks': checks, 'sources_postrun_hashes': sources, 'wall_seconds': time.perf_counter() - start, 'limits': 'Producer artifact verification; not independent scientific review or empirical validation'}
    write('raw/PACKAGE_VERIFICATION.json', result)
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
