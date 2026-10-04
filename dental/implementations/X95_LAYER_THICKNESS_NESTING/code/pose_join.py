"""YML rule join. Digital source contract only; no physical qualification."""
import os, sys
os.environ.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', MPLBACKEND='Agg')
sys.dont_write_bytecode = True
from pathlib import Path
import json, hashlib, math, csv, time, copy, xml.etree.ElementTree as ET
import numpy as np
ROOT = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, x):
    Path(p).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def source_gate(rule):
    expected = [('anterior_upper', 'single_anterior_crown', 'upper', 0.8), ('anterior_lower', 'single_anterior_crown', 'lower', 0.4), ('posterior_upper', 'single_posterior_crown', 'upper', 1.0), ('posterior_lower', 'single_posterior_crown', 'lower', 0.5)]
    actual = [(r['id'], r['indication'], r['half'], r['min_mm']) for r in rule['crown_rules']]
    return rule['product'] == 'KATANA Zirconia YML' and actual == expected and all((r['footnotes'] == ['*1', '*2'] and r['arrow_locator'] and (r['status'] == 'IMAGE_BOUND_RULE') for r in rule['crown_rules'])) and (sha(ROOT / rule['source_path']) == rule['source_sha256'])

def source_only_strength():
    root = ET.parse(ROOT / 'inputs/S21_article.xml').getroot()
    table = root.find('.//table-wrap[@id="bioengineering-13-00462-t003"]')
    rows = []
    for tr in table.findall('.//tbody/tr'):
        cells = [''.join(td.itertext()).strip() for td in tr.findall('td')]
        rows.append({'slice': int(cells[0]), 'Prime': {'sigma0_MPa': float(cells[1].split()[0]), 'm': float(cells[2].split()[0])}, 'YML': {'sigma0_MPa': float(cells[3].split()[0]), 'm': float(cells[4].split()[0])}})
    profiles = json.loads((ROOT / 'inputs/MATERIAL_PROFILES.json').read_text())['profiles']

    def compare(pp):
        return all((row[brand][k] == pp[row['slice'] - 1][brand][k] for row in rows for brand in ['Prime', 'YML'] for k in ['sigma0_MPa', 'm']))
    bad = copy.deepcopy(profiles)
    bad[0]['YML']['sigma0_MPa'] += 1.0
    out = {'status': 'PASS_SOURCE_ONLY' if compare(profiles) else 'FAIL', 'rows_compared': 10, 'values_compared': 40, 'all_exact': compare(profiles), 'injected_sigma_plus1_rejected': not compare(bad), 'unit_sigma0': 'MPa', 'unit_m': 'dimensionless', 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.3390/bioengineering13040462 Table3 XML bioengineering-13-00462-t003', 'compared_quantity': 'B3B coupon sigma0 MPa and m dimensionless', 'refutes_us': True}, 'crown_strength_transfer': 'UNKNOWN_NOT_USED', 'resolution': 'PER_SURFACE_REGION', 'xml_sha256': sha(ROOT / 'inputs/S21_article.xml'), 'rows': rows}
    dump(ROOT / 'raw/S21_SOURCE_CONTROL.json', out)
    return out

def point_triangle_distance(p, t):
    (a, b, c) = (t[:, 0], t[:, 1], t[:, 2])
    ab = b - a
    ac = c - a
    n = np.cross(ab, ac)
    n2 = np.sum(n * n, 1)
    height = np.sum((p - a) * n, 1) / np.maximum(n2, 1e-300)
    q = p - height[:, None] * n
    aq = q - a
    aa = np.sum(ab * ab, 1)
    bb = np.sum(ab * ac, 1)
    cc = np.sum(ac * ac, 1)
    dd = np.sum(aq * ab, 1)
    ee = np.sum(aq * ac, 1)
    det = aa * cc - bb * bb
    u = (cc * dd - bb * ee) / np.maximum(det, 1e-300)
    v = (aa * ee - bb * dd) / np.maximum(det, 1e-300)
    d = np.where((det > 1e-25) & (u >= 0) & (v >= 0) & (u + v <= 1), np.linalg.norm(p - q, axis=1), np.inf)
    for (x, y) in ((a, b), (b, c), (c, a)):
        edge = y - x
        w = np.clip(np.sum((p - x) * edge, 1) / np.maximum(np.sum(edge * edge, 1), 1e-300), 0, 1)
        d = np.minimum(d, np.linalg.norm(p - (x + w[:, None] * edge), axis=1))
    return float(d.min())

def sufficiency():
    a = np.array([1.1, 0.6])
    b = a[::-1].copy()
    req = np.array([1.0, 0.5])
    aa = bool(np.all(a >= req))
    bb = bool(np.all(b >= req))
    out = {'histogram_identity_error_mm': float(np.max(abs(np.sort(a) - np.sort(b)))), 'volume_identity_error_mm3': math.fsum(a.tolist()) - math.fsum(b.tolist()), 'minimum_identity_error_mm': float(a.min() - b.min()), 'state_A_thickness_mm': a.tolist(), 'state_B_thickness_mm': b.tolist(), 'fixed_locations': ['upper_half', 'lower_half'], 'source_requirement_mm': req.tolist(), 'rule_outcome_A': aa, 'rule_outcome_B': bb, 'downstream_pass_count_difference': int(np.sum(a >= req) - np.sum(b >= req)), 'minimum_slack_difference_mm': float((a - req).min() - (b - req).min()), 'minimal_extension': 'joint thickness/location-or-rule-ID pairing; histogram alone loses the relation', 'witness_origin': 'our_own_fixture , not empirical reference', 'resolution': 'PER_POINT'}
    assert out['histogram_identity_error_mm'] == out['volume_identity_error_mm3'] == 0 and aa != bb
    dump(ROOT / 'raw/SUFFICIENCY_WITNESS.json', out)
    return out

def axis(tilt, az):
    (t, a) = np.deg2rad([tilt, az])
    return np.array([np.sin(t) * np.cos(a), np.sin(t) * np.sin(a), np.cos(t)])

def poses(v, pr):
    g = pr['pose_grid']
    c = g['stock_clearance_mm']
    out = []
    for H in g['height_discs_mm']:
        ds = np.unique(np.r_[np.arange(c, H - c + 1e-12, g['step_mm']), H - c])
        for tilt in g['tilt_degrees']:
            for az in [0] if tilt == 0 else g['azimuth_degrees']:
                n = axis(tilt, az)
                mx = float(np.max(v @ n))
                span = float(np.ptp(v @ n))
                for depth in ds:
                    out.append({'height_mm': H, 'tilt_deg': tilt, 'azimuth_deg': az, 'top_depth_mm': float(depth), 'stock_ok': bool(depth + span <= H - c), 'span_mm': span})
    return out

def run():
    st = time.perf_counter()
    pr = json.loads((ROOT / 'PREREG_R1.json').read_text())
    rule = json.loads((ROOT / 'RULE_CONTRACT.json').read_text())
    assert sha(ROOT / 'RULE_CONTRACT.json') == pr['rule_contract_sha256'] and source_gate(rule)
    bad = copy.deepcopy(rule)
    bad['crown_rules'][0]['footnotes'] = ['*1', '*3']
    assert not source_gate(bad)
    strength = source_only_strength()
    witness = sufficiency()
    rows = []
    controls = []
    for e in json.loads((ROOT / 'inputs/R4_FROZEN_EXPORTS.json').read_text())['exports']:
        a = dict(np.load(ROOT / 'raw' / f"{e['family']}_wall.npz", allow_pickle=False))
        (p, d, area, v) = (a['points'], a['distance'], a['area'], a['vertices'])
        assert sha(e['mesh_path']) == e['mesh_sha256']
        assert len(p) > 0
        raw = dict(np.load(e['mesh_path'], allow_pickle=False))
        inner = raw['vertices'][raw['faces'][raw['roles'] == 1]]
        ids = np.unique(np.r_[int(d.argmin()), np.linspace(0, len(p) - 1, 32, dtype=int)])
        indep = np.array([point_triangle_distance(p[i], inner) for i in ids])
        err = float(np.max(abs(indep - d[ids])))
        assert err <= 1e-10
        inject = d[ids].copy()
        inject[0] += 0.01
        assert np.max(abs(indep - inject)) > 1e-10
        (high, low) = (0.8, 0.4) if e['family'] == 'anterior' else (1.0, 0.5)
        pp = poses(v, pr)
        records = []
        maxerr = 0.0
        mismatch = 0
        best = None
        for (ix, pose) in enumerate(pp):
            H = pose['height_mm']
            n = axis(pose['tilt_deg'], pose['azimuth_deg'])
            depth = pose['top_depth_mm']
            z = depth + np.max(v @ n) - p @ n
            required = np.where(z <= H / 2, high, low)
            slack = d - required
            newpass = bool(np.all(slack >= 0))
            oldpass = bool(np.all(d >= 0.5))
            viol = slack < 0
            theta = math.radians(pose['tilt_deg'])
            phi = math.radians(pose['azimuth_deg'])
            nx = math.sin(theta) * math.cos(phi)
            ny = math.sin(theta) * math.sin(phi)
            nz = math.cos(theta)
            top = max(v[:, 0] * nx + v[:, 1] * ny + v[:, 2] * nz)
            zc = depth + top - (p[:, 0] * nx + p[:, 1] * ny + p[:, 2] * nz)
            rc = np.full(len(p), low)
            rc[zc <= H * 0.5] = high
            sc = d - rc
            maxerr = max(maxerr, float(abs(slack - sc).max()))
            mismatch += int(newpass != bool(np.all(sc >= 0)))
            k = int(np.argmin(slack))
            record = dict(pose, pose_index=ix, yml_sampled_pass=newpass, universal0_5_sampled_pass=oldpass, min_slack_mm=float(slack[k]), violating_points=int(viol.sum()), point_count=len(p), violating_area_mm2=float(area[viol].sum()), area_weighted_deficit_mm3=float(np.sum(area * np.maximum(-slack, 0))), different_point_decisions=int(np.sum((slack >= 0) != (d >= 0.5))), witness_point_index=k)
            records.append(record)
            if pose['stock_ok'] and (best is None or record['area_weighted_deficit_mm3'] < best['area_weighted_deficit_mm3']):
                best = record
        assert mismatch == 0 and maxerr <= 1e-10
        dest = ROOT / 'raw' / f"{e['family']}_POSE_GRID.csv"
        with dest.open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        stock = [r for r in records if r['stock_ok']]
        allowed = [r for r in stock if r['yml_sampled_pass']]
        base = [r for r in stock if r['universal0_5_sampled_pass']]
        assert best is not None
        n = axis(best['tilt_deg'], best['azimuth_deg'])
        z = best['top_depth_mm'] + np.max(v @ n) - p @ n
        req = np.where(z <= best['height_mm'] / 2, high, low)
        np.savez_compressed(ROOT / 'exports' / f"{e['family']}_point_slack.npz", points=p, thickness_mm=d, required_mm=req, slack_mm=d - req, disc_depth_mm=z, area_mm2=area)
        row = {'family': e['family'], 'key': e['key'], 'resolution': 'PER_POINT', 'requested_poses': len(records), 'stock_fitting_poses': len(stock), 'stock_rejected_poses': len(records) - len(stock), 'allowed_yml_poses': len(allowed), 'allowed_universal0_5_poses': len(base), 'complete_decisions_changed': sum((r['yml_sampled_pass'] != r['universal0_5_sampled_pass'] for r in stock)), 'best_diagnostic_pose': best, 'sampled_min_wall_mm': float(d.min()), 'all_pose_impossibility_witness': {'point_index': int(d.argmin()), 'distance_mm': float(d.min()), 'minimum_requirement_any_pose_mm': low, 'necessary_slack_mm': float(d.min() - low)}, 'qualification': 'NOT_FABRICATION_QUALIFIED', 'unknowns': ['continuous wall proof', 'source uncertainty', 'batch IFU', 'actual prep', 'optics', 'CAM'], 'pose_grid_path': str(dest.relative_to(ROOT))}
        rows.append(row)
        controls.append({'family': e['family'], 'all_pose_points_checked': len(records) * len(p), 'max_slack_difference_mm': maxerr, 'decision_mismatches': mismatch, 'independent_distance_points': len(ids), 'independent_distance_max_error_mm': err, 'injected_distance_plus0_01_rejected': True})
        print(e['family'], len(stock), len(allowed), best['min_slack_mm'], flush=True)
    out = {'round': 'R1', 'claim_type': 'information_link', 'primary_gate': 'PASS' if any((r['complete_decisions_changed'] for r in rows)) else 'NEGATIVE_NO_COMPLETE_POSE_DECISION_CHANGED', 'rows': rows, 'controls': controls, 'source_gate': {'good_contract': True, 'bad_footnote_blocks_all': True}, 'source_only_strength_control': strength, 'sufficiency_witness': witness, 'seconds': time.perf_counter() - st, 'dropout': {'requested': sum((r['requested_poses'] for r in rows)), 'rejected': sum((r['stock_rejected_poses'] for r in rows)), 'reason': 'outside disc axial stock; rule-rejected but stock-fitting counted separately'}, 'resolution': 'PER_POINT', 'timescale': 'HANDOVER', 'physical_measurement': 'NOT_PERFORMED'}
    out['dropout']['fraction'] = out['dropout']['rejected'] / out['dropout']['requested']
    dump(ROOT / 'raw/R1_RESULTS.json', out)
    return out
if __name__ == '__main__':
    run()
