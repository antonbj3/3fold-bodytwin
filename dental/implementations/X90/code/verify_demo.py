"""Executed falsifiers and exact-summary witnesses; no fabricated dental truth."""
import argparse
import copy
import hashlib
import subprocess
import sys
from pathlib import Path
from common import ROOT, read, dump, encoded, validate_manifest, validate_decisions, validate_force_binding, validate_pulp_observation, sha
from render import verify_report

def sufficiency():
    a = [(0, 2), (2, 0)]
    b = [(0, 0), (2, 2)]

    def summary(states):
        return [(min(q), max(q), sum(q) / len(q)) for q in zip(*states)]
    (sa, sb) = (summary(a), summary(b))
    assert encoded(sa) == encoded(sb)
    ta = [sum(v) for v in a]
    tb = [sum(v) for v in b]
    d1 = {'material_wall': 'FAIL', 'milling': 'UNKNOWN'}
    d2 = {'material_wall': 'UNKNOWN', 'milling': 'FAIL'}
    assert sum((v == 'FAIL' for v in d1.values())) == sum((v == 'FAIL' for v in d2.values()))
    qa = [(0, 0, 0.25), (1, 0, 0.25), (1, 1, 0.25), (0, 1, 0.25)]
    qb = [(x + 3, y, z) for (x, y, z) in qa]

    def plane_summary(q):
        distances = sorted((abs(z) for (x, y, z) in q))
        return dict(p95_mm=distances[-1], rms_mm=(sum((d * d for d in distances)) / len(distances)) ** 0.5, maximum_mm=distances[-1], projected_area_mm2=1)
    (ma, mb) = (plane_summary(qa), plane_summary(qb))
    assert encoded(ma) == encoded(mb)
    centroid_a = sum((x for (x, y, z) in qa)) / len(qa)
    centroid_b = sum((x for (x, y, z) in qb)) / len(qb)
    return dict(claim_type='capability', external_referent=dict(kind='our_own_fixture', locator='code/verify_demo.py:sufficiency', compared_quantity='Exact-summary sufficiency for laboratory action selection, not dental prediction', refutes_us=True), force_witness=dict(summary=sa, identity_error=0, summary_bytes_identical=True, states_A=a, states_B=b, total_force_A_N=ta, total_force_B_N=tb, downstream_max_force_difference_N=max(tb) - max(ta), resolution='PER_TOOTH -> PER_ARCH', minimum_extension='For this query, joint total-load constraint; for general consumers retain joint region-force state and wrench, not independent marginal intervals'), decision_witness=dict(summary_failed_count=1, identity_error=0, states=[d1, d2], next_measurements=['registered wall/preparation geometry', 'tool/approach/CAM accessibility'], minimum_extension='Identity of unresolved edge and geometry/state binding'), surface_witness=dict(summary=ma, identity_error=0, summary_bytes_identical=True, centroid_x_A_mm=centroid_a, centroid_x_B_mm=centroid_b, downstream_centroid_difference_mm=centroid_b - centroid_a, resolution='PER_SURFACE_REGION', minimum_extension='Contact location/region identity for this query; retained surface/contact field for general downstream questions'), enclosure='Exact enumeration of finite integer states; no affine approximation')

def faults(out):
    out = Path(out)
    case = read(out / 'case.json')
    m = read(ROOT / 'inputs' / (case['patient'].replace(':', '_') + '.json'))
    checks = []

    def rejects(name, fn):
        try:
            fn()
        except (ValueError, KeyError, AssertionError):
            checks.append(dict(name=name, rejected=True))
        else:
            raise AssertionError('Undetected fault: ' + name)
    for (key, value, name) in [('units', 'um', 'wrong_unit'), ('source_pose', 'UNREGISTERED', 'missing_common_frame')]:
        q = dict(m)
        q[key] = value
        rejects(name, lambda q=q: validate_manifest(q))
    q = copy.deepcopy(m)
    q['cbct'] = {'subject': 'another_patient', 'units': 'mm', 'frame': 'CBCT_NATIVE_MM'}
    rejects('other_patient_CBCT', lambda : validate_manifest(q))
    valid = dict(case=m['case'], geometry_sha256=case['source_geometry_sha256'], gap_convention='MEASURED_LEVELED_REFERENCE', acquisition_kind='PHYSICAL_MEASUREMENT')
    validate_force_binding(valid, m, case['source_geometry_sha256'])
    for (key, value, name) in [('case', 'another_patient', 'other_patient_force'), ('geometry_sha256', '0' * 64, 'other_geometry_force'), ('gap_convention', 'NATIVE_X21_GAPS', 'unloaded_force_transfer')]:
        q = dict(valid)
        q[key] = value
        rejects(name, lambda q=q: validate_force_binding(q, m, case['source_geometry_sha256']))
    rows = copy.deepcopy(case['decisions'])
    del rows[0]['uncertainty']
    rejects('missing_uncertainty', lambda : validate_decisions(rows))
    rows = copy.deepcopy(case['decisions'])
    rows[3]['evidence'] = 'kalibrerat'
    rejects('unknown_force_promoted_to_calibrated', lambda : validate_decisions(rows))
    h = out / 'report.html'
    original = h.read_bytes()
    key = 'p95_mm'
    old = repr(case['external_referent']['measured_original'][key])
    h.write_bytes(original.replace(('data-metric="' + key + '">' + old).encode(), ('data-metric="' + key + '">999').encode(), 1))
    try:
        rejects('reported_p95_999mm', lambda : verify_report(out))
    finally:
        h.write_bytes(original)
    f = out / 'FROZEN_PREDICTIONS.json'
    original = f.read_bytes()
    f.write_bytes(original + b' ')
    try:
        rejects('prediction_freeze_mutated', lambda : verify_report(out))
    finally:
        f.write_bytes(original)
    c = out / 'design_contract.json'
    q = read(c)
    q['input_sha256']['crown'] = '0' * 64
    dump(out / 'FAULT_DESIGN_CONTRACT.json', q)
    cfg = read(ROOT / 'SITE_CONFIG.json')
    pkg = (ROOT / Path(cfg['package'])).resolve() if not Path(cfg['package']).is_absolute() else Path(cfg['package'])
    script = "import sys;from pathlib import Path;sys.path.insert(0,sys.argv[1]);from designgate.gate import check;d=Path(sys.argv[2]);check({'crown':d/'crown_local.stl','prep':d/'virtual_prep.stl','antagonist':d/'antagonist_local.stl'},'3Y',d/'FAULT_DESIGN_CONTRACT.json',round_version='R2')"
    env = __import__('os').environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    cmd = ['bwrap', '--ro-bind', '/', '/', '--dev', '/dev', '--proc', '/proc', '--tmpfs', '/tmp', '--unshare-net', '--', cfg['python'], '-c', script, str(pkg / 'batch9/demos/PATIENT360/dependencies/results/LANE_X34_STL_DESIGN_GATE/code'), str(out)]
    p = subprocess.run(cmd, env=env, capture_output=True, text=True)
    (out / 'FAULT_X34.stderr').write_text(p.stderr)
    if p.returncode == 0 or 'Hash-bound contract mismatch: crown' not in p.stderr:
        raise AssertionError('Actual X34 source-hash mutation not rejected')
    checks.append(dict(name='actual_X34_source_hash', rejected=True, command=cmd))
    verify_report(out)
    return dict(patient=case['patient'], count=len(checks), checks=checks)

def optional_cbct_contract():
    m = read(ROOT / 'inputs/Bits2Bites_13.json')
    geometry = ROOT / 'raw/CBCT_SYNTHETIC_GEOMETRY.txt'
    geometry.write_text('Synthetic non-anatomical API binding fixture\n')
    obs = dict(patient=m['patient'], unit='mm', quantity='total_hard_tissue_to_pulp', fdi=36, resolution='PER_POINT', locator='our_own_fixture:CBCT_SYNTHETIC_GEOMETRY.txt', input_sha256=sha(geometry), distance_mm=2, digital_two_surface_radius_mm=0.2)
    validate_pulp_observation(obs, m, geometry)
    mutations = [('patient', 'wrong'), ('fdi', 46), ('resolution', 'POPULATION'), ('unit', 'um'), ('input_sha256', '0' * 64), ('distance_mm', -1), ('digital_two_surface_radius_mm', -1)]
    checks = []
    for (key, value) in mutations:
        bad = dict(obs)
        bad[key] = value
        try:
            validate_pulp_observation(bad, m, geometry)
        except ValueError:
            checks.append(dict(field=key, rejected=True))
        else:
            raise AssertionError('Optional CBCT field error accepted:' + key)
    return dict(status='PASS', controls=checks, physical_measurement='NOT_RUN', external_referent=dict(kind='our_own_fixture', locator='code/verify_demo.py:optional_cbct_contract', compared_quantity='Subject/FDI/source/resolution binding only; not pulp anatomy', refutes_us=False))

def compare(a, b):
    (a, b) = (Path(a), Path(b))
    names = ['report.html', 'case.json', 'FROZEN_PREDICTIONS.json', 'crown.stl', 'crown_local.npz', 'case_figure.png']
    rows = []
    for name in names:
        (first, second) = ((a / name).read_bytes(), (b / name).read_bytes())
        rows.append(dict(artifact=name, bytes=len(first), sha256=hashlib.sha256(first).hexdigest(), bitwise_equal=first == second, byte_count_difference=len(first) - len(second)))
    return rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--report', action='append', default=[])
    p.add_argument('--replay-root')
    a = p.parse_args()
    fault_rows = [faults(q) for q in a.report]
    replay = []
    if a.replay_root:
        for q in a.report:
            rows = compare(q, Path(a.replay_root) / Path(q).name)
            replay.append(dict(patient=read(Path(q) / 'case.json')['patient'], artifacts=rows))
        if not all((x['bitwise_equal'] for r in replay for x in r['artifacts'])):
            dump(ROOT / 'raw/REPLAY_FAILED.json', replay)
            raise AssertionError('Nonidentical fresh process output; preserved raw/REPLAY_FAILED.json')
    result = dict(faults=fault_rows, sufficiency=sufficiency(), optional_cbct_contract=optional_cbct_contract(), replay=replay, status='PASS', independent_review='NOT_PERFORMED', physical_measurements='NOT_RUN')
    dump(ROOT / 'raw/VERIFICATION.json', result)
    print('PASS', sum((r['count'] for r in fault_rows)), 'fault injections')
if __name__ == '__main__':
    main()
