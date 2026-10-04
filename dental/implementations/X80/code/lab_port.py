"""Future measured-lab intake: no fixture is accepted as independent validation.
CLI: /usr/bin/python3 -s code/lab_port.py path/to/measurement.json
"""
from pathlib import Path
import json, sys, hashlib, math
P = Path(__file__).resolve().parents[1]
REQUIRED = ['measurement_id', 'measured_at', 'kind', 'material_batch', 'region', 'species', 'unit', 'value', 'uncertainty', 'observation_operator', 'duration', 'medium', 'source_locator', 'independent_measurement']

def intake(record):
    missing = [k for k in REQUIRED if k not in record or record[k] in [None, '']]
    if missing:
        return {'status': 'UNKNOWN', 'reason': 'missing measured ports', 'missing': missing}
    if not isinstance(record['value'], (int, float)) or isinstance(record['value'], bool) or (not math.isfinite(record['value'])):
        return {'status': 'REJECTED', 'reason': 'Finite measured numeric value required'}
    if record['source_locator'].lower().startswith(('fixture', 'simulation', 'our_own_fixture')):
        return {'status': 'REJECTED', 'reason': 'Fixture provenance cannot validate the frozen measurement'}
    if record['independent_measurement'] is not True:
        return {'status': 'REJECTED', 'reason': 'Independent measured fact required; simulation/fixture cannot validate biology'}
    if record['kind'] == 'resin_DC_30min':
        f = json.loads((P / 'FROZEN_PREDICTIONS.json').read_text())
        match = dict(species='FTIR_C=C_conversion', unit='%', region='ATR_FTIR_same_as_PMC10892052', observation_operator='1637cm-1/1525cm-1 cured_to_uncured_peak_ratio', duration='30min', medium='not_applicable_FTIR')
        mismatch = [k for (k, v) in match.items() if record[k] != v]
        if record.get('oven') != 'Triad2000' or record.get('material_identity') != 'PMC10892052_model_resin':
            mismatch.append('material_or_oven')
        if mismatch:
            return {'status': 'UNKNOWN', 'reason': 'Frozen recipe observation is not matched', 'mismatch': mismatch}
        return {'status': 'MATCHED_MODEL_CONDITIONAL_COMPARISON', 'observed': record['value'], 'predicted': f['predicted_value'], 'absolute_error_pp': abs(record['value'] - f['predicted_value']), 'source_prediction_sha256': hashlib.sha256((P / 'FROZEN_PREDICTIONS.json').read_bytes()).hexdigest(), 'decision': 'No clinical or viability decision; this validates only the declared DC observation'}
    if record['kind'] == 'ISO10993_5_viability':
        required = ['cells', 'extraction_ratio', 'extraction_duration', 'dilution', 'negative_control', 'positive_control', 'assay']
        if any((k not in record for k in required)):
            return {'status': 'UNKNOWN', 'reason': 'ISO assay ports incomplete', 'missing': [k for k in required if k not in record]}
        if record['unit'] != '%' or record['species'] != 'relative_cell_viability':
            return {'status': 'REJECTED', 'reason': 'Viability percent required; DC, current and Ti ppm are not viability'}
        return {'status': 'LAB_ASSAY_ENDPOINT_ONLY', 'relative_viability_percent': record['value'], 'below_70pct_viability': record['value'] < 70, 'exactly70pct_boundary': record['value'] == 70, 'uncertainty': record['uncertainty'], 'cells': record['cells'], 'clinical_risk': 'UNKNOWN', 'caveat': 'A laboratory endpoint under the supplied assay; no patient safe dose or material recommendation'}
    return {'status': 'UNKNOWN', 'reason': 'No calibrated biological consumer for this measured source; preserved input, no inferred tissue dose'}
if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: lab_port.py measured_record.json')
    print(json.dumps(intake(json.loads(Path(sys.argv[1]).read_text())), indent=2))
