"""Compare future independently measured bottlenecks with frozen predictions."""
import argparse, json, hashlib
from pathlib import Path
if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Requires real independent measurements; no synthetic fallback')
    p.add_argument('--measurement', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    here = Path(__file__).resolve().parent
    f = here / 'FROZEN_MEASUREMENT_TARGETS.json'
    expected = (here / 'FROZEN_MEASUREMENT_TARGETS.sha256').read_text().split()[0]
    assert hashlib.sha256(f.read_bytes()).hexdigest() == expected, 'Frozen measurement targets changed'
    frozen = json.load(open(f))
    m = json.load(open(a.measurement))
    assert m['external_referent']['kind'] == 'independent_measurement' and m['external_referent']['locator'], 'Independent source required'
    assert m['compared_quantity'] == 'directional_outer_jaw_margin_mm', 'Quantity mismatch'
    assert m['independent_of_TF2_annotation'] is True, 'Annotation reuse cannot calibrate anatomy'
    assert m['registered_same_specimen'] is True, 'Other patients do not validate local edges'
    rows = []
    for v in m['observations']:
        i = v['frozen_target_index']
        pred = frozen['targets'][i]
        assert v['case'] == pred['case'] and v['tooth'] == pred['tooth']
        error = float(v['measured_margin_mm']) - pred['predicted_annotation_margin_mm']
        rows.append({'target_index': i, 'signed_prediction_error_mm': error, 'measurement_error_bound_mm': v['measurement_error_bound_mm'], 'conditional_scenario_survives': float(v['measured_margin_mm']) - float(v['measurement_error_bound_mm']) >= 0.608, 'whole_root_and_biological_status': 'UNDETERMINED'})
    a.output.write_text(json.dumps({'frozen_predictions_sha256': expected, 'external_referent': m['external_referent'], 'rows': rows, 'fit_performed': False, 'clinical_accuracy': 'UNKNOWN'}, indent=2) + '\n')
