"""Compare independently supplied measurements to an already frozen prediction file.

No values are supplied by this program. Default CSV must describe existing frozen
TF2 query coordinates; a new specimen requires predict_local_case.py first.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from measure import HERE, dump, sha

def compare(predictions, rows):
    known = {(r['case'], int(r['fdi'])): r for r in predictions['predictions']}
    result = []
    for r in rows:
        key = (r['case'], int(r['fdi']))
        if key not in known:
            raise ValueError(f'unfrozen query: {key}')
        if r.get('reference_kind') != 'independent_measurement' or not r.get('locator'):
            raise ValueError('Independent source provenance is required; fixtures cannot be physical truth')
        if r.get('matched_fiducial_registration') != 'true' or r.get('units') != 'mm':
            raise ValueError('Need matched query frame/fiducial registration and millimetres')
        p = known[key]
        h = float(r['height_mm'])
        w = float(r['width1_mm'])
        (lo, hi) = p['height_interval_mm']
        (wl, wh) = p['width1_interval_mm']
        result.append({'case': key[0], 'fdi': key[1], 'height_measured_mm': h, 'width1_measured_mm': w, 'height_error_mm': h - p['height_mm'], 'width1_error_mm': w - p['width1_mm'], 'height_box_pass': lo <= h <= hi, 'width_box_pass': wl <= w <= wh, 'source_locator': r['locator']})
    return {'n_measured': len(result), 'gate': 'UNKNOWN_NO_VALUES' if not result else 'PASS_SCOPED_BOX_ONLY' if all((r['height_box_pass'] and r['width_box_pass'] for r in result)) else 'FAIL', 'raw_comparisons': result, 'clinical_validation': 'UNKNOWN'}
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--csv', required=True, type=Path)
    p.add_argument('--predictions', type=Path, default=HERE / 'FROZEN_LAB_PREDICTIONS.json')
    p.add_argument('--out', type=Path, default=HERE / 'PHYSICAL_COMPARISON.json')
    a = p.parse_args()
    expected = a.predictions.with_suffix('.sha256').read_text().strip()
    assert sha(a.predictions) == expected, 'Frozen predictions hash drift'
    predictions = json.loads(a.predictions.read_text())
    rows = list(csv.DictReader(a.csv.open()))
    out = compare(predictions, rows)
    out['predictions_sha256'] = expected
    out['measurement_csv_sha256'] = sha(a.csv)
    dump(a.out, out)
    print(json.dumps(out, indent=2))
