"""Reproduce the retained negative R3 source gate, without generating bad geometry."""
from common import ROOT, check_frozen
import json
p = check_frozen('PREREG_R3.json')
ref = json.loads((ROOT / 'rounds/R3.json').read_text())
for (k, v) in ref['external_values'].items():
    assert abs(p['source_geometry'][k] - v) == ref['absolute_errors_mm'][k]
    assert abs(p['source_geometry'][k] - v) > p['metrics']['STL_dimension_absolute_error_mm_max']
assert ref['gate'] == 'FAIL_SOURCE_GEOMETRY' and ref['did_numerical_export_run'] is False
print('PASS: retained R3 failure is reproducible; bad capsule geometry not exported.')
