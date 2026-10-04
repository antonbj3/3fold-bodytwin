from pathlib import Path
import json, hashlib, datetime
from series import judge
R = Path(__file__).resolve().parents[1]

def wr(n, x):
    (R / n).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
x = json.load(open(R / 'raw/PREDICTIONS_R3.json'))
f = json.load(open(R / 'FROZEN_PREDICTIONS_R3.json'))
assert hashlib.sha256((R / f['prediction_file']).read_bytes()).hexdigest() == f['prediction_sha256']
y = json.load(open(R / 'inputs/R3_FORCE_REFERENTS.json'))
groups = {}
for (name, v) in y.items():
    if name not in ['within_study', 'between_study']:
        continue
    rows = []
    for r in v['rows']:
        p = next((a for a in x['rows'] if a['sheet_mm'] == r['sheet_mm']))
        rows.append({**r, 'predicted_N': p['predicted_absolute_palatal_force_N'], 'relative_error': abs(p['predicted_absolute_palatal_force_N'] - r['force_N']) / r['force_N'], 'implied_intercept_mm': 0.25 - r['force_N'] / r.get('source_K_N_per_mm', p['predicted_K_N_per_mm'])})
    gate = judge([r['predicted_N'] for r in rows], [r['force_N'] for r in rows], 0.25)
    groups[name] = {'rows': rows, 'gate': gate, 'locator': v['locator'], 'doi': v['doi'], 'condition': v['condition']}
g = [r['implied_intercept_mm'] for r in groups['within_study']['rows']]
gr = max(g) - min(g)
truth = [r['force_N'] for r in groups['between_study']['rows']]
injection = judge([10 * t for t in truth], truth, 0.25)
out = {'claim_type': 'capability', 'decision': 'WITHIN_STUDY_CONVERSION_PASS_BETWEEN_STUDY_FORCE_TRANSFER_FAIL', 'groups': groups, 'implied_intercept_range_mm': gr, 'constant_effective_intercept_gate': 'PASS' if gr <= 0.02 else 'FAIL', 'fault_injection': {'truth_gate': judge(truth, truth, 0.25)['gate'], 'ten_times_wrong_force_gate': injection['gate']}, 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.1007/s00056-015-0307-3;PMID26446503;Abstract Results', 'compared_quantity': 'Duran palatalforce .25mm tooth21 for .5/.75sheet', 'refutes_us': True}, 'negative_result': True, 'interpretation': 'Source-conditioned interpolation is executable; source-agnostic force transfer is rejected. Carry insertion/securing condition and measured contact intercept alongside finished thickness. No single between-study factor causally identified.', 'next_construction': 'Measureactualformed regionalh and six-axisforce-displacement on same selectedrealarch replica under declared seating load; fit passive compliance+effective intercept per movementdirection; freeze.625formed-thickness predictions before independentreplicates.'}
wr('raw/RESULTS_R3.json', out)
txt = '# R3 handoff\n\nAbsoluteforceconversion within2016 sensorstudy PASS. Independent2015 study transfer FAIL: predictedvsmeasured3.01/4.49N ' + str([r['predicted_N'] for r in groups['between_study']['rows']]) + '. Relativeerrors ' + str(groups['between_study']['gate']['relative_errors']) + '.\n\nSamebrand,activation andsheetthickness do not define a transferable physicalforceport. Differentcontralateraltooth,formingandretention/temperature/reportingconditions must remain sourceconditioned; no attributing difference to30N alone.\n\nRealarch outputsfromR1 remain simulationswithunknownphysicalerror, including256N artificialprediction causedbyinvalidpredicted-toothgeometry. Next measurementcontract must includeFDI/geometryquality before numericalforceclaims.\n\nNext: same-specimen finishedthickness+zero/intermediate/target6-axisforce atknownseatingload, at leasttwo thicknessstates andoneindependentholdout; preserve independentforceandmoment originsandload-specificrelaxation.\n'
(R / 'HANDOFF_R3.md').write_text(txt)
(R / 'HANDOFF.md').write_text(txt)
wr('CURRENT_WORK_STATE.json', {'lane': 'X19-aligner-force', 'phase': 'R3_adjudicated_package_and_measurement_port', 'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'last_gate': out['decision'], 'next_operation': out['next_construction'], 'empirical_patient_force': 'UNKNOWN'})
print(json.dumps({k: out[k] for k in ['decision', 'groups', 'implied_intercept_range_mm', 'fault_injection']}, indent=2))
