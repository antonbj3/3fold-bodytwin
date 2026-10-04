"""Verify delivered numeric outputs against retained controls; no physical inference."""
import json, math, pathlib
import numpy as np

def validate(rows, records, result, prereg):
    by = {(v['case'], v['site']): v for v in rows}
    rays = [v for v in records if v['kind'] == 'ray_box']
    assert len(rays) == len(by), 'ray control count'
    assert len({(v['case'], v['site']) for v in rays}) == len(rays), 'duplicate ray controls'
    for v in rays:
        (candidate, control) = (v['candidate'], v['control'])
        error = 0.0 if candidate is None and control is None else max((abs(a - b) for (a, b) in zip(candidate, control))) if candidate is not None and control is not None else 1000.0
        assert math.isfinite(error) and error <= 1e-08, 'actual candidate/control mismatch'
        assert v['error'] == error and v['pass'] == (error <= 1e-08), 'stale cached verdict'
        row = by[v['case'], v['site']]
        q = next((q for q in row['rays'] if q['depth_mm'] == 4.0))
        expected = [None, None, None] if candidate is None else [candidate[1] - candidate[0], candidate[1], -candidate[0]]
        for (key, want) in zip(['envelope_chord_mm', 'buccal_extent_mm', 'lingual_extent_mm'], expected):
            assert q[key] == want, 'unbound reported geometry: ' + key
        for q in row['rays']:
            if q['envelope_chord_mm'] is not None:
                assert q['envelope_chord_mm'] == q['buccal_extent_mm'] + q['lingual_extent_mm'], 'chord/endpoints'
        a = row['axial_run_limits_from_pose_mm']
        assert row['axial_envelope_chord_mm'] == (None if a is None else a[1] - a[0]), 'axial endpoints'
    train_ids = set(prereg['selection']['train_cases'])
    test_ids = set(prereg['selection']['test_cases'])
    assert not train_ids & test_ids, 'case leakage'
    assert {r['case'] for r in rows} <= train_ids | test_ids, 'unassigned case'
    train = [r['rays'][1]['envelope_chord_mm'] for r in rows if r['case'] in train_ids]
    test = [r['rays'][1]['envelope_chord_mm'] for r in rows if r['case'] in test_ids]
    train = [v for v in train if v is not None]
    test = [v for v in test if v is not None]
    mean = float(np.mean(train))
    rmse = float(np.sqrt(np.mean((np.asarray(test) - mean) ** 2)))
    g = result['annotation_geometry']
    assert g['train_count'] == len(train) and g['test_count'] == len(test), 'reported split counts'
    assert g['training_population_mean_chord_mm'] == mean and g['test_RMSE_population_proxy_mm'] == rmse, 'reported fit/test result'
    assert np.array_equal(np.quantile(train + test, [0, 0.25, 0.5, 0.75, 1]), g['all_chord_quantiles_mm']), 'reported chord distribution'
    return {'sites': len(rows), 'train_sites': len(train), 'test_sites': len(test), 'case_overlap': 0, 'RMSE_mm': rmse, 'biological_patient_independence': 'UNKNOWN across source families', 'physical_validation': False}

def run(root):
    root = pathlib.Path(root)
    out = {}
    for tag in ['R1', 'R2']:
        rows = [json.loads(v) for v in (root / f'raw/PER_SITE_{tag}.jsonl').read_text().splitlines()]
        out[tag] = validate(rows, json.loads((root / f'raw/CONTROLS_{tag}.json').read_text()), json.loads((root / f'results_{tag}.json').read_text()), json.loads((root / f'PREREG_{tag}.json').read_text()))
    return out
if __name__ == '__main__':
    import sys
    print(json.dumps(run(pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[1]), indent=2))
