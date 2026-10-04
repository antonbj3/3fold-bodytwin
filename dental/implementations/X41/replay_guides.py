"""Retain inherited guide predicates and independently check source moments."""
import hashlib
import importlib.util
import json
import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from scipy.stats import gamma
ROOT = Path(__file__).resolve().parent

def main():
    start = time.perf_counter()
    p = ROOT / 'PREREG_GUIDE_REPLAY_SUPPLEMENT.json'
    if hashlib.sha256(p.read_bytes()).hexdigest() != p.with_suffix('.sha256').read_text().strip():
        raise ValueError('Guide prereg hash drift')
    for r in json.loads((ROOT / 'SOURCE_MANIFEST_GUIDE_AND_DATASET.json').read_text()):
        if hashlib.sha256((ROOT / r['local']).read_bytes()).hexdigest() != r['sha256']:
            raise ValueError('Guide source drift')
    spec = importlib.util.spec_from_file_location('guide_op', ROOT / 'inputs/guide_risk_operator.py')
    op = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(op)
    profiles = json.loads((ROOT / 'inputs/guide_profiles.json').read_text())['profiles']
    original = json.loads((ROOT / 'inputs/LANE_X31_CLINICAL_ANSWERS/results.json').read_text())
    rows = []
    checks = []
    text = (ROOT / 'references/Varga2020.txt').read_text()
    source_table = text.split('MANDIBLE')[-1].split('MAXILLA')[0]
    vals = {}
    for line in source_table.splitlines():
        m = re.match('\\s*(Fr|Pi|Fu)\\s+(37|29)\\s+([\\d. ]+)\\s*$', line)
        if m:
            vals.setdefault(m[1], []).append([float(x) for x in m[3].split()])
    wu = ' '.join(ET.parse(ROOT / 'references/Wu2020.xml').getroot().itertext())

    def check(name, observed, expected, tol, bad):
        ok = bool(abs(observed - expected) <= tol)
        reject = bool(abs(bad - expected) > tol)
        checks.append({'gate': name, 'observed': observed, 'expected': expected, 'tolerance': tol, 'pass': ok, 'injected_value': bad, 'injected_rejected': reject})
        if not ok or not reject:
            raise ValueError(name)
    for profile in profiles:
        gid = profile['id']
        r = original['guide_parameters'][gid]
        if gid == 'dynamic_navigation':
            for (mu, sd) in [profile[k] for k in ['entry_mean_sd_mm', 'apex_mean_sd_mm', 'angle_mean_sd_deg']]:
                if not re.search(f'{mu:.2f}\\s*±\\s*{sd:.2f}', wu):
                    raise ValueError('Wu moments locator failed')
        else:
            arm = {'freehand': 'Fr', 'pilot_guided': 'Pi', 'fully_guided': 'Fu'}[gid]
            (angle, apex) = vals[arm]
            for (i, (a, b)) in enumerate(zip(profile['angle_mean_sd_deg'] + profile['entry_mean_sd_mm'] + profile['apex_mean_sd_mm'], angle[:2] + angle[5:7] + apex[:2])):
                check(gid + ':source_moment' + str(i), a, b, 1e-12, a + 0.1)
        apex = op.scenario_margin(profile, x=0.05)
        whole = op.bound_margin(profile, x=0.05, boundary_budget=0.3)
        check(gid + ':apex_gap', apex, r['scenario_apex_required_gap_mm'], 1e-09, apex + 1)
        check(gid + ':whole_body_gap', whole, r['whole_body_sufficient_gap_B0p3_mm'], 1e-09, whole + 1)
        check(gid + ':independent_scalar_control', op.scalar_control(whole, profile, B=0.3), 0.05, 1e-07, op.scalar_control(whole + 1, profile, B=0.3))
        (mu, sd) = profile['apex_mean_sd_mm']
        toward = 1.3 + gamma.ppf(0.95, (mu / sd) ** 2, scale=sd * sd / mu)
        old = next((x for x in original['rounds']['R1']['orientation_counterexample'] if x['guide'] == gid))
        check(gid + ':direction_counterexample', toward, old['toward_plane_required_apex_gap_mm'], 1e-09, toward + 0.1)
        rows.append({'guide': gid, 'apex_gap_mm': apex, 'whole_body_gap_mm': whole, 'toward_plane_gap_mm': float(toward), 'physical_margin': 'UNKNOWN', 'numerical_reference': 'inherited operator plus independent scalar moment-bound control'})
    out = {'claim_type': 'capability', 'residual_clearance_mm': 1, 'nominal_comparator_mm': 2, 'geometric_tail_target': 0.05, 'guide_rows': rows, 'controls': checks, 'physical_certificates': 0, 'wall_seconds': time.perf_counter() - start}
    (ROOT / 'raw/GUIDE_PARAMETER_REPLAY.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'guide_replay': 'PASS', 'controls': len(checks), 'physical': 'UNKNOWN', 'wall_seconds': out['wall_seconds']}))
if __name__ == '__main__':
    main()
