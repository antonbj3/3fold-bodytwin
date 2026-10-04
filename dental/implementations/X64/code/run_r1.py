"""Offline, two spatial manufacturing consumers on the same published D1 crown."""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
from pathlib import Path
import csv, json, math, sys, time, resource
import numpy as np
import trimesh
from lxml import etree
from freeze import ROOT, W, sha, dump
OUT = Path(os.environ.get('X64_REPLAY_DIR', str(ROOT)))
(OUT / 'raw').mkdir(parents=True, exist_ok=True)

def source_path(manifest, key):
    m = manifest[key]
    return ROOT / m['local_path'] if m['local_path'] else Path(m['original_path'])

def validate_inputs():
    p = json.loads((ROOT / 'PREREG_R1.json').read_text())
    assert sha(ROOT / 'PREREG_R1.json') == (ROOT / 'PREREG_R1.sha256').read_text().strip()
    assert sha(ROOT / 'SOURCE_MANIFEST.json') == p['source_manifest_sha256']
    manifest = json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text())
    for (key, m) in manifest.items():
        assert sha(source_path(manifest, key)) == m['sha256'], key
    return (p, manifest)

def parse_tables(manifest):
    import re

    def mean_sd(s):
        (a, b) = re.match('\\s*([\\d.]+)\\s*±\\s*([\\d.]+)', s).groups()
        return (float(a), float(b))
    (optical, cure) = ([], [])
    root = etree.parse(str(source_path(manifest, 'optical_source')))
    trs = root.xpath('//table-wrap[@id="T3"]//tr')
    times = {'ST': 420.0, 'SP': 105.0}
    for ri in (2, 3):
        tr = trs[ri]
        route = ''.join(tr[0].itertext())
        for (ci, h) in ((1, 0.5), (2, 1.0)):
            text = ''.join(tr[ci].itertext())
            (mu, sd) = mean_sd(text)
            ranges = re.findall('\\(([\\d.]+)[-–]([\\d.]+)\\)', text)
            (lo, hi) = map(float, ranges[-1])
            optical.append({'route': route, 'h_mm': h, 'mean_TP': mu, 'SD_TP': sd, 'observed_min_TP': lo, 'observed_max_TP': hi, 'duration_min': times[route], 'source': 'PMC10088447', 'locator': f'//table-wrap[@id="T3"]//tr[{ri + 1}]/*[{ci + 1}]', 'resolution_level': 'POPULATION', 'time_scale': 'HANDOVER'})
    root = etree.parse(str(source_path(manifest, 'cure_source')))
    for (ri, tr) in enumerate(root.xpath('//table-wrap[@id="materials-17-01496-t001"]//tr')[1:], 2):
        name = ''.join(tr[0].itertext())
        (mu, sd) = mean_sd(''.join(tr[1].itertext()))
        minutes = 0.0 if name == 'Control' else float(name.split()[0])
        cure.append({'cure_min': minutes, 'mean_MPa': mu, 'SD_MPa': sd, 'n': 9, 'source': 'PMC11012777', 'locator': f'//table-wrap[@id="materials-17-01496-t001"]//tr[{ri}]/*[2]', 'resolution_level': 'POPULATION', 'time_scale': 'HANDOVER'})
    facit_o = [(16.46, 1.09, 15.17, 18.38), (12.7, 0.62, 11.93, 13.87), (16.61, 1.35, 14.03, 18.37), (11.84, 0.45, 11.17, 12.6)]
    facit_c = [(15.9, 3.8), (80.5, 3.2), (76.5, 1.6), (83.2, 2.2)]
    observed = [(r['mean_TP'], r['SD_TP'], r['observed_min_TP'], r['observed_max_TP']) for r in optical]
    observed += [(r['mean_MPa'], r['SD_MPa']) for r in cure]
    ref = facit_o + facit_c
    err = max((abs(a - b) for (row, exp) in zip(observed, ref) for (a, b) in zip(row, exp)))
    injected = [list(row) for row in observed]
    injected[0][0] += 1.0
    injected_err = max((abs(a - b) for (row, exp) in zip(injected, ref) for (a, b) in zip(row, exp)))
    return (optical, cure, {'max_error': err, 'injected_plus1_rejected': injected_err > 1e-10, 'locator': 'Source XML table cells compared with independently transcribed published table', 'quantity': 'TP and flexural strength mean/SD', 'resolution_level': 'POPULATION'})

def brute_ray(vertices, faces, origin, direction):
    """Independent vectorized Moller-Trumbore check, every original triangle."""
    triangles = vertices[faces]
    e1 = triangles[:, 1] - triangles[:, 0]
    e2 = triangles[:, 2] - triangles[:, 0]
    q = np.cross(np.broadcast_to(direction, e2.shape), e2)
    det = np.einsum('ij,ij->i', e1, q)
    valid = abs(det) > 1e-12
    inv = np.zeros_like(det)
    inv[valid] = 1.0 / det[valid]
    s = origin - triangles[:, 0]
    u = np.einsum('ij,ij->i', s, q) * inv
    r = np.cross(s, e1)
    v = r @ direction * inv
    t = np.einsum('ij,ij->i', e2, r) * inv
    valid &= (u >= -1e-10) & (v >= -1e-10) & (u + v <= 1 + 1e-10) & (t > 1e-09)
    return float(t[valid].min()) if valid.any() else None

def optical_geometry(p, manifest):
    mesh = trimesh.load_mesh(source_path(manifest, 'crown'), process=True)
    g = json.loads(source_path(manifest, 'geometry').read_text())
    (cf, nf) = (mesh.triangles_center, mesh.face_normals)
    radial = cf - np.array([0.0, 0.0, g['z_m']])
    outer = (np.einsum('ij,ij->i', radial, nf) > 0) & (cf[:, 2] > g['z_m'] + 0.3)
    ids = np.flatnonzero(outer)
    ids = ids[np.linspace(0, len(ids) - 1, min(p['TM02']['ray_sample_max'], len(ids))).astype(int)]
    offset = p['TM02']['ray_offset_mm']
    origins = cf[ids] - offset * nf[ids]
    dirs = -nf[ids]
    (hits, rays, _) = mesh.ray.intersects_location(origins, dirs, multiple_hits=True)
    distance = np.einsum('ij,ij->i', hits - origins[rays], dirs[rays])
    ray_len = np.full(len(ids), np.inf)
    valid = distance > 1e-09
    np.minimum.at(ray_len, rays[valid], distance[valid] + offset)
    ray_len[~np.isfinite(ray_len)] = np.nan
    visible = (cf[ids, 2] > g['z_m'] + 1.0) & (nf[ids, 1] < -0.2)
    checks = []
    for k in np.linspace(0, len(ids) - 1, 12).astype(int):
        control = brute_ray(mesh.vertices, mesh.faces, origins[k], dirs[k])
        control = None if control is None else control + offset
        err = abs(ray_len[k] - control) if control is not None and np.isfinite(ray_len[k]) else None
        checks.append({'ray_id': int(k), 'error_mm': err})
    maxerr = max((r['error_mm'] for r in checks if r['error_mm'] is not None))
    np.savez_compressed(OUT / 'raw/OPTICAL_RAYS.npz', points=cf[ids], normals=nf[ids], thickness_mm=ray_len, visible=visible, face_ids=ids)
    with (OUT / 'raw/OPTICAL_RAYS.csv').open('w') as f:
        writer = csv.writer(f)
        writer.writerow(['ray', 'face', 'x_mm', 'y_mm', 'z_mm', 'h_mm', 'visible', 'resolution_level', 'time_scale'])
        for k in range(len(ids)):
            writer.writerow([k, int(ids[k]), *cf[ids[k]], ray_len[k], int(visible[k]), 'PER_POINT', 'HANDOVER'])
    return (mesh, ray_len, visible, {'checks': checks, 'max_error_mm': maxerr, 'injected_plus0_01mm_rejected': bool(maxerr + 0.01 > p['metrics']['ray_distance_control_error_mm_max'])})

def tp_mean(h, rows):
    (r0, r1) = rows
    w = (h - 0.5) / 0.5
    return np.exp((1 - w) * math.log(r0['mean_TP']) + w * math.log(r1['mean_TP']))

def get_stresses(manifest):
    d = np.load(source_path(manifest, 'stress'))
    nf = d['nf'].astype(float)
    tensors = {}
    sigmas = {}
    controls = {}
    for name in ('axial', 'offaxis30'):
        a = d['Sf__' + name].astype(float)
        S = np.zeros((len(a), 3, 3))
        S[:, 0, 0] = a[:, 0]
        S[:, 1, 1] = a[:, 1]
        S[:, 2, 2] = a[:, 2]
        S[:, 0, 1] = S[:, 1, 0] = a[:, 3]
        S[:, 1, 2] = S[:, 2, 1] = a[:, 4]
        S[:, 0, 2] = S[:, 2, 0] = a[:, 5]
        P = np.eye(3)[None] - nf[:, :, None] * nf[:, None, :]
        tangent = P @ S @ P
        sig = np.maximum(np.linalg.eigvalsh(tangent)[:, -1], 0.0)
        basis = np.zeros_like(nf)
        basis[:, 0] = 1
        basis[abs(nf[:, 0]) > 0.8] = [0, 1, 0]
        e1 = np.cross(nf, basis)
        e1 /= np.linalg.norm(e1, axis=1)[:, None]
        e2 = np.cross(nf, e1)
        aa = np.einsum('ij,ijk,ik->i', e1, S, e1)
        bb = np.einsum('ij,ijk,ik->i', e2, S, e2)
        ab = np.einsum('ij,ijk,ik->i', e1, S, e2)
        sig_control = np.maximum(0.5 * (aa + bb + np.sqrt((aa - bb) ** 2 + 4 * ab * ab)), 0.0)
        controls[name] = {'direct_2x2_error_MPa_per_N': float(abs(sig - sig_control).max())}
        sigmas[name] = sig
        tensors[name] = S
    return (d, sigmas, controls)

def main():
    tic = time.perf_counter()
    (p, manifest) = validate_inputs()
    (optical, cure, table_check) = parse_tables(manifest)
    dump(OUT / 'raw/SOURCE_ROWS.json', {'optical': optical, 'cure': cure})
    (mesh, h, vis, ray_check) = optical_geometry(p, manifest)
    source_support = np.isfinite(h) & (h >= 0.5) & (h <= 1.0)
    op = []
    for route in ('ST', 'SP'):
        rows = [r for r in optical if r['route'] == route]
        pred = tp_mean(h, rows)
        supported = vis & source_support
        op.append({'route': route, 'duration_min': rows[0]['duration_min'], 'endpoint_nominal_TP': [r['mean_TP'] for r in rows], 'coupon_endpoint_acceptance': [r['mean_TP'] >= p['TM02']['TP_floor'] for r in rows], 'visible_rays': int(vis.sum()), 'visible_supported_rays': int(supported.sum()), 'visible_support_fraction': float(supported.sum() / vis.sum()), 'nominal_pass_fraction_among_supported': float((pred[supported] >= p['TM02']['TP_floor']).mean()) if supported.any() else None, 'physical_D1_release': 'UNKNOWN', 'resolution_level': 'PHENOMENOLOGICAL mapped to PER_POINT geometry; source POPULATION'})
    i = int(np.flatnonzero(vis)[0])
    j = int(np.flatnonzero(~vis)[0])
    ha = np.full(len(h), 0.5)
    hb = ha.copy()
    ha[i] = 1.0
    hb[j] = 1.0
    sp = [r for r in optical if r['route'] == 'SP']
    ta = tp_mean(ha, sp)
    tb = tp_mean(hb, sp)
    oa = float(ta[vis].min())
    ob = float(tb[vis].min())
    optical_witness = {'summary': 'Complete thickness histogram and arithmetic mean thickness; canonical mean TP also retained', 'mean_h_A_mm': float(ha.mean()), 'mean_h_B_mm': float(hb.mean()), 'identity_error_mm': float(abs(ha.mean() - hb.mean())), 'mean_TP_identity_error': abs(math.fsum(ta) / len(ta) - math.fsum(tb) / len(tb)), 'same_histogram': bool(np.array_equal(np.sort(ha), np.sort(hb))), 'downstream_min_visible_TP_A': oa, 'downstream_min_visible_TP_B': ob, 'downstream_difference_TP': abs(oa - ob), 'accept_A': oa >= 12.0, 'accept_B': ob >= 12.0, 'minimal_extension_for_fixed_question': 'TP minimum on the declared visible subset; to change view retain thickness/location pairing', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'raw/SUFFICIENCY_STATES.npz', 'compared_quantity': 'Synthetic rearrangement on real D1 locations using published endpoint TP', 'refutes_us': True}, 'resolution_level': 'PER_POINT', 'time_scale': 'HANDOVER'}
    (d, sigmas, stress_check) = get_stresses(manifest)
    fe = json.loads(source_path(manifest, 'resin_fe').read_text())
    F = p['TM05']['load_N']
    mechanical = []
    for (name, sig) in sigmas.items():
        peak = float(sig.max())
        expected = fe['cases'][name]['s1_surface_max_per_N']
        stress_check[name].update({'summary_error_MPa_per_N': abs(peak - expected), 'injected_summary_plus1percent_rejected': abs(peak * 1.01 - expected) > p['metrics']['stress_summary_error_MPa_per_N_max']})
        choices = [r['cure_min'] for r in cure if r['mean_MPa'] >= F * peak]
        mechanical.append({'loadcase': name, 'force_N': F, 'peak_demand_MPa': F * peak, 'nominal_shortest_cure_min': min(choices) if choices else None, 'rows': [dict(r, nominal_ratio=r['mean_MPa'] / (F * peak), nominal_accept=r['mean_MPa'] >= F * peak, conditional_force_limit_N=r['mean_MPa'] / peak) for r in cure], 'physical_D1_release': 'UNKNOWN', 'resolution_level': 'PER_SURFACE_REGION'})
    demand = F * sigmas['offaxis30']
    hi = int(np.argmax(demand))
    lo = int(np.argmin(demand))
    sa = np.full(len(demand), 100.0)
    sb = sa.copy()
    sa[hi] = 60.0
    sb[lo] = 60.0
    positive = demand > 1e-12
    ra = float((sa[positive] / demand[positive]).min())
    rb = float((sb[positive] / demand[positive]).min())
    strength_witness = {'summary': 'Identical mean, min, max and complete strength histogram of all D1 facets', 'mean_A_MPa': float(sa.mean()), 'mean_B_MPa': float(sb.mean()), 'identity_error_MPa': abs(float(sa.mean() - sb.mean())), 'same_histogram': bool(np.array_equal(np.sort(sa), np.sort(sb))), 'downstream_min_ratio_A': ra, 'downstream_min_ratio_B': rb, 'downstream_difference_ratio': abs(ra - rb), 'accept_A': ra >= 1.0, 'accept_B': rb >= 1.0, 'minimal_extension_for_fixed_question': 'min_i(S_i/sigma_i) for each specified loadcase; changing loads requires retaining local strength/stress association', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'raw/SUFFICIENCY_STATES.npz', 'compared_quantity': 'Synthetic equal-histogram strength fields on actual D1 stress map', 'refutes_us': True}, 'resolution_level': 'PER_POINT', 'time_scale': 'HANDOVER'}
    np.savez_compressed(OUT / 'raw/SUFFICIENCY_STATES.npz', optical_h_A=ha, optical_h_B=hb, visible=vis, strength_A=sa, strength_B=sb, demand_MPa=demand)
    np.savez_compressed(OUT / 'raw/DEMANDS.npz', cf=d['cf'], nf=d['nf'], af=d['af'], tied=d['tied'], axial_MPa=F * sigmas['axial'], offaxis30_MPa=demand)
    source_checks = table_check['max_error'] <= p['metrics']['table_absolute_error_max'] and table_check['injected_plus1_rejected']
    numeric_checks = ray_check['max_error_mm'] <= p['metrics']['ray_distance_control_error_mm_max'] and all((max(v['direct_2x2_error_MPa_per_N'], v['summary_error_MPa_per_N']) <= p['metrics']['stress_summary_error_MPa_per_N_max'] and v['injected_summary_plus1percent_rejected'] for v in stress_check.values()))
    witness_checks = optical_witness['identity_error_mm'] == 0 and optical_witness['mean_TP_identity_error'] == 0 and (strength_witness['identity_error_MPa'] == 0)
    result = {'round': 'R1', 'claim_type': 'capability', 'PREREG_sha256': sha(ROOT / 'PREREG_R1.json'), 'TM02': op, 'TM05': mechanical, 'sufficiency': {'optical': optical_witness, 'strength': strength_witness}, 'verification': {'source': table_check, 'rays': ray_check, 'stress': stress_check, 'execution_gate_pass': bool(source_checks and numeric_checks and witness_checks)}, 'summary_sufficient': False, 'physical_D1_validation': 'UNKNOWN', 'external_referent': [p['TM02']['external_referent'], p['TM05']['external_referent']], 'cost': {'seconds': time.perf_counter() - tic, 'peak_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'new_FE_solves': 0, 'new_physical_measurements': 0, 'fit': 0, 'questions': 0, 'inherited_resin_FE_seconds': fe['seconds_total'], 'inherited_geometry_seconds': json.loads(source_path(manifest, 'geometry').read_text())['seconds'], 'all_source_discovery_and_prior_campaign_cost': 'UNMEASURED; no total-cost superiority claim'}, 'next_operation': 'R2 replace scalar transfer by explicit interval/inverse local assay contract; avoid extrapolation and keep physical release UNKNOWN'}
    dump(OUT / 'raw/RESULTS_R1.json', result)
    dump(OUT / 'CURRENT_WORK_STATE.json', {'phase': 'R1_COMPLETE', 'last_gate': result['verification'], 'next_operation': result['next_operation'], 'lane': 'X64-manufacturing'})
    (OUT / 'HANDOFF_R1.md').write_text('R1 executed on D1. Both global material summaries fail exact sufficiency under spatial rearrangement; see raw/RESULTS_R1.json and raw/SUFFICIENCY_STATES.npz. Published rows are external facit only for coupon TP/flexural strength. Crown transfer and physical release UNKNOWN. Next construction: R2 rigorous finite-box model enclosure and inverse matched assay targets; no new generic source adapter.\n')
    print(json.dumps({'execution_gate': result['verification']['execution_gate_pass'], 'optical': op, 'mechanical': [{k: r[k] for k in ('loadcase', 'peak_demand_MPa', 'nominal_shortest_cure_min')} for r in mechanical], 'witnesses': [optical_witness['downstream_difference_TP'], strength_witness['downstream_difference_ratio']], 'cost': result['cost']}, indent=2))
    assert result['verification']['execution_gate_pass'], 'Preserve failed R1 artifact before any changed construction'
if __name__ == '__main__':
    main()
