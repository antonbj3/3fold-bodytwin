"""Read the actual constructed cavity/margin result with all missing inputs intact."""
from common import *
import argparse
p = argparse.ArgumentParser()
p.add_argument('--case', required=True)
args = p.parse_args()
r3 = json.loads((ROOT / 'raw/R3.json').read_text())
r4 = json.loads((ROOT / 'raw/R4.json').read_text())
row = next((q for q in r3['rows'] if q['case'] == args.case), None)
if row is None:
    raise SystemExit('case absent from R3 pilots; no answer inferred')
rim = next((q for q in r4['rows'] if q['case'] == args.case))
final = row['grids'][-1]
if sha(final['stl']) != final['stl_sha256']:
    raise SystemExit('cavity export hash drift')
out = {'case': args.case, 'digital_scenario_gate': rim['linked_gate'], 'gap_bounds_source_units': [final['gap']['lower_max'], final['gap']['upper_max']], 'gap_semantics': 'Euclidean distance to conservative cubical cavity boundary; no designed cement spacer added', 'opening_expansion': rim['opening'], 'direction': row['direction'], 'cavity_stl': final['stl'], 'source_finish_line': str(DATA / (args.case + '_SOURCE_FINISH_LINE.npz')), 'unit': 'SOURCE_NATIVE_UNKNOWN', 'consumer_status': 'GEOMETRIC_TOOL_AVAILABLE', 'complete_crown': 'UNKNOWN_NOT_CONSTRUCTED', 'physical_manufacturing_decision': 'UNKNOWN_MISSING_SAME_OBJECT_MEASUREMENTS', 'required_inputs': ['Explicit physical source unit/scale calibration', 'Calibrated source scanner and margin uncertainty', 'Original same-object crown exterior and inside surfaces', 'Independent same-object dry fit/cemented intaglio and lab acceptance bands']}
print(json.dumps(out, indent=2))
