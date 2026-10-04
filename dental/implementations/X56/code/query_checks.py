"""Exercise the frozen R3 held-load fault control through the public JSON CLI."""
from pathlib import Path
import copy, json, subprocess, sys
LANE = Path(__file__).resolve().parents[1]
base = json.loads((LANE / 'demo_inputs/spread_support.json').read_text())
cases = {'valid': copy.deepcopy(base), 'bad_10um': copy.deepcopy(base), 'missing_bone': copy.deepcopy(base)}
cases['bad_10um']['implant_marker_readings_um'][1][1] += 10
cases['missing_bone']['bone_marker_readings_um'] = None
outputs = {}
directory = LANE / 'raw/query_cases'
directory.mkdir(exist_ok=True)
for (name, value) in cases.items():
    file = directory / (name + '.json')
    file.write_text(json.dumps(value, indent=2) + '\n')
    output = subprocess.run([sys.executable, str(LANE / 'code/query_motion.py'), str(file)], capture_output=True, text=True, check=True)
    outputs[name] = json.loads(output.stdout)
assert outputs['valid']['held_wrench_validation']['status'] == 'PASS_CONSISTENT_WITH_DECLARED_BOXES'
assert outputs['valid']['sampled_regional_peak_interval_um'] == [27.25, 35.25]
assert outputs['bad_10um']['status'] == 'REJECT_HELD_WRENCH_OUTSIDE_READING_BOXES'
assert outputs['missing_bone']['status'] == 'UNKNOWN_MISSING_LOCAL_BONE_REFERENCE'
(LANE / 'raw/QUERY_DEMO.json').write_text(json.dumps(outputs['valid'], indent=2) + '\n')
out = {'round': 'R3A', 'claim_type': 'capability', 'status': 'COMPLETE', 'outcome': 'JSON_HELD_WRENCH_GUARD_EXECUTED', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'PREREG_R3.json and demo_inputs/spread_support.json', 'compared_quantity': 'query interval and held-wrench reading-box consistency', 'refutes_us': False}, 'cases': outputs, 'full_cost': {'fixed_CLI_cases': 3, 'new_empirical_measurements': 0}, 'resolution': 'PER_SURFACE_REGION; declared error PHENOMENOLOGICAL', 'time_scale': 'SIMULTANEOUS'}
(LANE / 'rounds/R3A').mkdir(exist_ok=True)
(LANE / 'rounds/R3A/results.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({'valid_pass': True, 'bad_10um_rejected': True, 'missing_bone_unknown': True, 'empirical_calibration': 'UNKNOWN'}))
