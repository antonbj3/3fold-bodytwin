"""Independent STL triple gate. Digital scopes and physical debts are explicit."""
import argparse
import datetime as dt
import json
import math
import resource
import sys
import time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from parents import HERE, original_gate, geometry as G, distance as D, transport as T, milling, inherited, sha, verify_dependencies
RULES = ['mesh_health', 'scale', 'insertion', 'material_wall', 'film_min', 'film_max', 'cement_physical', 'ball_milling', 'occlusal_contact', 'export']

def result(status='UNKNOWN', reason='Required information missing', **kw):
    value = dict(status=status, reason=reason, resolution='PER_SURFACE_REGION')
    value.update(kw)
    return value

def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def raw_normals(path):
    b = Path(path).read_bytes()
    n = int.from_bytes(b[80:84], 'little') if len(b) >= 84 else 0
    if len(b) == 84 + 50 * n:
        dtype = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
        return np.frombuffer(b, dtype=dtype, count=n, offset=84)['n'].astype(float)
    lines = b.decode('ascii').splitlines()
    return np.array([[float(v) for v in line.split()[2:]] for line in lines if line.strip().startswith('facet normal')])

def mesh_health(data):
    issues = []
    rows = {}
    for (name, d) in data.items():
        diag = dict(d['diagnostics'])
        m = d['mesh']
        stored = raw_normals(d['path'])
        if stored.shape != (len(m.faces), 3) or not np.isfinite(stored).all():
            issues.append(dict(mesh=name, fault='invalid stored normals', face_id=0))
            continue
        norms = np.linalg.norm(stored, axis=1)
        nonzero = norms > 1e-12
        agreement = np.ones(len(stored))
        agreement[nonzero] = np.sum(stored[nonzero] * m.face_normals[nonzero], axis=1) / norms[nonzero]
        wrong = np.flatnonzero(agreement < 1 - 1e-05)
        diag.update(stored_normal_disagreements=len(wrong), stored_zero_normals=int((~nonzero).sum()), oriented_volume_mm3=float(m.volume), self_intersection_status='UNKNOWN')
        rows[name] = diag
        if len(wrong):
            i = int(wrong[0])
            issues.append(dict(mesh=name, fault='stored normal disagrees with facet winding', face_id=i, point_mm=m.triangles_center[i].tolist(), stored_normal=stored[i].tolist(), geometric_normal=m.face_normals[i].tolist()))
        if diag['degenerate_faces'] or diag['duplicate_faces'] or (not diag['winding_consistent']):
            tri = m.triangles
            deg = np.flatnonzero(np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) <= 1e-12)
            sf = np.sort(m.faces, axis=1)
            (_, inv, counts) = np.unique(sf, axis=0, return_inverse=True, return_counts=True)
            dup = np.flatnonzero(counts[inv] > 1)
            i = int(deg[0] if len(deg) else dup[0] if len(dup) else 0)
            issues.append(dict(mesh=name, fault='degenerate, duplicate or inconsistent winding', face_id=i, point_mm=tri[i].mean(0).tolist(), triangle_mm=tri[i].tolist(), degeneracy_cross_norm_mm2=1e-12, numerical_policy_not_manufactured_failure=True))
        if name == 'crown' and (not diag['watertight'] or m.volume <= 0):
            edges = m.edges_unique[np.flatnonzero(np.bincount(m.edges_unique_inverse) != 2)]
            issues.append(dict(mesh=name, fault='crown open or globally inverted', face_id=0, boundary_edge_mm=m.vertices[edges[0]].tolist() if len(edges) else None))
    return result('FAIL' if issues else 'PASS', 'Digital topology and normal consistency; self-intersection remains UNKNOWN', diagnostics=rows, witnesses=issues, scope='Closed oriented crown; prep/antagonist may be scan patches')

def scale_rule(crown, units):
    extent = float(max(crown['diagnostics']['extent_mm']))
    guess = 'um likely' if extent > 60 else 'mm plausible' if extent >= 0.5 else 'meter or undersized patch possible'
    if units not in ('mm', 'um'):
        return result(reason='STL has no unit; heuristic is not unit evidence', guessed_unit=guess, resolution='PER_TOOTH')
    return result('PASS' if 0.5 <= extent <= 60 else 'FAIL', 'Frozen single-crown extent policy [0.5,60] mm; not anatomy validation', extent_mm=extent, policy_mm=[0.5, 60], guessed_unit=guess, resolution='PER_TOOTH', witness=dict(point_mm=crown['vertices'][int(np.argmax(crown['vertices'][:, 0]))].tolist()), policy_resolution='PHENOMENOLOGICAL', replacement_measurement='calibrated scanner/unit export metadata')

def classify_interval(interval, threshold, minimum=True):
    (lo, hi) = interval
    if minimum:
        return 'PASS' if lo >= threshold else 'FAIL' if hi < threshold else 'UNKNOWN'
    return 'PASS' if hi <= threshold else 'FAIL' if lo > threshold else 'UNKNOWN'

def robust_unsigned(source, target, threshold, errors, minimum=True, face_ids=None):
    """Exact set perturbation inequality; numerical enclosure is conditional."""
    e = sum(errors)
    d = D.enclosure(source, target, minimum_required=threshold + e if minimum else None, maximum_allowed=threshold - e if not minimum else None, source_face_ids=face_ids)
    which = 'minimum' if minimum else 'maximum'
    nominal = d[which + '_interval_mm']
    interval = [nominal[0] - e, nominal[1] + e]
    return result(classify_interval(interval, threshold, minimum), 'Hausdorff set-error envelope propagated by triangle inequality', nominal_interval_mm=nominal, robust_interval_mm=interval, threshold_mm=threshold, supplied_errors_mm=errors, witness=d[which + '_witness'], measurement=d, uncertainty_scope='unsigned labelled surface separation under supplied Hausdorff bounds', rigorous_geometry_perturbation_bound=True, rigorous_floating_point_enclosure=False, physical_validity='UNKNOWN until supplied bounds and region semantics independently measured')

def cavity_planes(crown, ids):
    """Exact unnormalised facet planes, with convex consistency check."""
    (triangles, ids) = G.faces(crown, ids)
    planes = []
    for tri in triangles:
        p = [[Q(float(x)) for x in v] for v in tri]
        u = [p[1][k] - p[0][k] for k in range(3)]
        v = [p[2][k] - p[0][k] for k in range(3)]
        a = [-(u[1] * v[2] - u[2] * v[1]), -(u[2] * v[0] - u[0] * v[2]), -(u[0] * v[1] - u[1] * v[0])]
        if sum((x * x for x in a)) == 0:
            raise ValueError('degenerate cavity facet')
        b = sum((a[k] * p[0][k] for k in range(3)))
        planes.append(a + [b])
    vertices = [[Q(float(x)) for x in v] for v in np.unique(triangles.reshape(-1, 3), axis=0)]
    if any((sum((a[k] * p[k] for k in range(3))) > a[3] for a in planes for p in vertices)):
        raise ValueError('intaglio is not a convex polyhedral cavity in exact stored coordinates')
    unique = list(dict.fromkeys((tuple(p) for p in planes)))
    return ([[str(x) for x in p] for p in unique], triangles.reshape(-1, 3).tolist())

def convex_controls(data, meta, rules):
    if meta.get('convex_cavity_complete') is not True:
        return
    ids = meta.get('regions', {}).get('intaglio')
    if not ids:
        return
    try:
        (planes, vertices) = cavity_planes(data['crown'], ids)
    except ValueError as e:
        rules['insertion'] = result(reason=str(e))
        rules['ball_milling'] = result(reason=str(e))
        return
    normal = [[float(Q(x)) for x in p[:3]] for p in planes]
    cone = inherited.insertion_cone(normal)
    direction = np.asarray(meta.get('extraction_axis', [0, 0, 1]), float)
    if direction.shape != (3,) or not np.isfinite(direction).all() or np.linalg.norm(direction) == 0:
        raise ValueError('invalid extraction direction')
    dq = [Q(float(v)) for v in direction]
    pq = [[Q(float(v)) for v in p] for p in data['prep']['vertices']]
    aq = [[Q(x) for x in a] for a in planes]
    contained = all((sum((a[k] * p[k] for k in range(3))) <= a[3] for a in aq for p in pq))
    recession = all((sum((a[k] * dq[k] for k in range(3))) >= 0 for a in aq))
    if contained and recession:
        rules['insertion'] = result('PASS', 'Exact rational cavity containment and recession; full straight extraction in declared convex cavity', planes=planes, extraction_direction=direction.tolist(), reused_finite_cone=cone, scope='ideal rigid STL geometry, complete convex cavity labeling asserted; no other scene obstacles', physical_uncertainty='UNKNOWN; no normal-orientation error envelope')
    else:
        rules['insertion']['reused_finite_cone'] = cone
        rules['insertion']['convex_containment'] = contained
        rules['insertion']['convex_recession'] = recession
    burs = meta.get('bur_radii_final_mm', [])
    if not burs:
        return
    allowances = meta.get('milling_allowance_mm', 0.0)
    reports = []
    for radius in burs:
        if not isinstance(radius, (int, float)) or not math.isfinite(radius) or radius <= 0:
            raise ValueError('invalid bur radius')
        r = milling.check(planes, vertices, radius, allowances, (-direction).tolist())
        if r['status'] == 'PASS':
            cert = r['certificate']
            verified = milling.verify_positive(planes, vertices, radius, allowances, cert['direction'], cert['centres'])
        elif r['status'] == 'FAIL':
            w = r['witness']
            verified = milling.verify_negative(planes, w['point'], radius, allowances, w['weights'])
        else:
            verified = False
        if r['status'] in ('PASS', 'FAIL') and (not verified):
            raise ValueError('unchanged GenCAD certificate failed independent replay')
        reports.append(dict(radius_final_mm=radius, certificate_rechecked=verified, **r))
    statuses = [r['status'] for r in reports]
    rules['ball_milling'] = result('PASS' if 'PASS' in statuses else 'FAIL' if all((x == 'FAIL' for x in statuses)) else 'UNKNOWN', 'Unchanged GenCAD ideal-ball convex-access check per bur', bur_reports=reports, planes=planes, vertices=vertices, scope='At least one listed ball meets declared allowance; shaft/fixture/vendor machine UNKNOWN')

def check(paths, material='3Y', contract=None, units=None, sinter_factor=None):
    verify_dependencies()
    start = time.perf_counter()
    cpu = time.process_time()
    meta = json.loads(Path(contract).read_text()) if contract else {}
    if not isinstance(meta, dict):
        raise ValueError('contract must be JSON object')
    unit = units or meta.get('units')
    if units and meta.get('units') not in (None, units):
        raise ValueError('conflicting units')
    scale = {'mm': 1.0, 'um': 0.001}.get(unit, 1.0)
    data = {k: G.load(p, scale) for (k, p) in paths.items()}
    health = mesh_health(data)
    sc = scale_rule(data['crown'], unit)
    rules = {k: result() for k in RULES}
    rules.update(mesh_health=health, scale=sc)
    out = dict(schema='x49-design-gate-v1', claim_type='capability', created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), inputs={k: dict(path=d['path'], sha256=d['sha256'], diagnostics=d['diagnostics']) for (k, d) in data.items()}, rules=rules, material=material, sinter_factor=sinter_factor, clinical_eligibility='UNKNOWN', declared_units=unit, preparation_s=time.perf_counter() - start)
    if sinter_factor is not None and (not math.isfinite(sinter_factor) or sinter_factor <= 0):
        raise ValueError('invalid sinter factor')
    if health['status'] == 'FAIL' or sc['status'] != 'PASS':
        for key in RULES[2:]:
            rules[key] = result(reason='Input health/explicit scale gate blocks geometric interpretation')
        return finish(out, start, cpu)
    parent = original_gate.check(paths, material=material, contract=contract, units=units, round_version='R2')
    out['predecessor_report'] = parent
    for name in ('insertion', 'material_wall', 'film_min', 'film_max', 'ball_milling', 'occlusal_contact'):
        rules[name] = parent['rules'][name]
    rules['cement_physical'] = result(reason='X13: nominal dry geometric gap does not identify seated spatial cement film, viscosity/fill/venting', source=str(HERE.parents[1] / 'results/LANE_X13_CEMENT_GAP/CEMENT_FIELD_PORT.json'), replacement_measurement='same-object seated regional film, dry/wet state, viscosity and boundary/fill conditions')
    errs = meta.get('surface_error_mm')
    regions = meta.get('regions', {})
    if errs is not None:
        if not isinstance(errs, dict) or any((k not in errs for k in paths)) or any((type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in errs.values())):
            raise ValueError('surface_error_mm requires finite nonnegative bound per mesh')
        out['uncertainty'] = dict(surface_error_mm=errs, physical_calibration=meta.get('error_source', 'UNKNOWN'), rigorous_numeric_enclosure=False, registration_error='must be included in relevant surface bounds')
        if all((k in regions for k in ('intaglio', 'exterior'))):
            (inner, ii) = G.faces(data['crown'], regions['intaglio'])
            (outer, oi) = G.faces(data['crown'], regions['exterior'])
            target = T.physical_mesh(outer)
            wall = parent['ifu_profile'].get('minimum_wall_mm', {}).get(meta.get('indication'))
            if wall is not None:
                rules['material_wall'] = robust_unsigned(inner, target, wall, [errs['crown'], errs['crown']], face_ids=ii)
                rules['material_wall']['source'] = parent['ifu_profile']['source']
                rules['material_wall']['product_profile'] = parent['ifu_profile']
                rules['material_wall']['witness']['target_STL_face'] = int(oi[rules['material_wall']['witness']['target_region_face']])
            if 'preparation' in regions and meta.get('film_limits_mm'):
                (pt, pi) = G.faces(data['prep'], regions['preparation'])
                target = T.physical_mesh(pt)
                if D.planar_exterior(inner, target, meta.get('preparation_outward_axis', [0, 0, 0])):
                    (low, high) = meta['film_limits_mm']
                    axis = np.asarray(meta['preparation_outward_axis'], float)
                    axis /= np.linalg.norm(axis)
                    side_clearance = float(np.min(inner.reshape(-1, 3) @ axis) - np.max(pt.reshape(-1, 3) @ axis))
                    for (name, limit, minimum) in [('film_min', low, True), ('film_max', high, False)]:
                        r = robust_unsigned(inner, target, limit, [errs['crown'], errs['prep']], minimum, ii)
                        r['exterior_slab_clearance_mm'] = side_clearance
                        r['exterior_slab_error_mm'] = errs['crown'] + errs['prep']
                        if minimum:
                            r['robust_interval_mm'][0] = max(r['robust_interval_mm'][0], side_clearance - r['exterior_slab_error_mm'] - D.EPS)
                            r['status'] = classify_interval(r['robust_interval_mm'], limit, minimum)
                        if side_clearance <= errs['crown'] + errs['prep'] + D.EPS:
                            r['status'] = 'UNKNOWN'
                            r['reason'] = 'Surface-error envelope can cross preparation; signed topology unsupported'
                        r['witness']['target_STL_face'] = int(pi[r['witness']['target_region_face']])
                        rules[name] = r
                        r['source'] = meta.get('film_limits_source', 'Caller laboratory protocol; PHENOMENOLOGICAL')
                elif errs['crown'] + errs['prep'] > 0:
                    for name in ('film_min', 'film_max'):
                        rules[name] = result(reason='Signed-film perturbation near arbitrary solid topology unsupported')
        if errs['antagonist'] + errs['crown'] > 0:
            rules['occlusal_contact'] = result(reason='Projected gap cannot use unsigned Hausdorff bound when XY overlap changes; specify measured axial error model')
    else:
        out['uncertainty'] = dict(physical_calibration='UNKNOWN; nominal digital decisions only')
    axial = meta.get('axial_error_mm')
    if axial is not None:
        if not isinstance(axial, dict) or any((type(axial.get(k)) not in (int, float) or not math.isfinite(axial[k]) or axial[k] < 0 for k in ('crown', 'antagonist'))):
            raise ValueError('axial_error_mm needs crown and antagonist finite nonnegative bounds')
        contact = parent['rules']['occlusal_contact']
        measure = contact.get('measurement', {})
        nominal = measure.get('minimum_gap_mm')
        band = meta.get('occlusal_gap_limits_mm')
        if meta.get('fixed_projected_domain_confirmed') is True and nominal is not None and band and (not measure.get('excluded_projected_faces')):
            error = axial['crown'] + axial['antagonist'] + G.TOL
            interval = [nominal - error, nominal + error]
            status = 'PASS' if interval[0] >= band[0] and interval[1] <= band[1] else 'FAIL' if interval[1] < band[0] or interval[0] > band[1] else 'UNKNOWN'
            rules['occlusal_contact'] = result(status, 'Barycentric axial-error envelope at fixed projected overlap', nominal_minimum_gap_mm=nominal, robust_interval_mm=interval, limits_mm=band, measurement=measure, witness=measure['witness'], axial_error_mm=axial, scope='Fixed XY vertices/triangulations and labelled overlap; axial displacements bounded at every vertex', rigorous_analytic_enclosure=True, rigorous_floating_point_enclosure=False, error_source=meta.get('axial_error_source', 'UNKNOWN physical calibration'), equation='|min(g_perturbed)-min(g_nominal)| <= e_upper+e_lower by barycentric convexity; EPS numeric allowance conditional')
        else:
            rules['occlusal_contact'] = result(reason='Axial bound requires confirmed fixed XY projection/overlap and nondegenerate labelled facets')
    if meta.get('common_frame_confirmed') is True:
        convex_controls(data, meta, rules)
    rules['export'] = result(reason='Request --export to run report-bound 3MF round-trip')
    return finish(out, start, cpu)

def finish(out, start, cpu):
    statuses = [x['status'] for x in out['rules'].values()]
    out['verdict'] = 'FAIL' if 'FAIL' in statuses else 'UNKNOWN' if 'UNKNOWN' in statuses else 'PASS'
    out['cost'] = dict(wall_s=time.perf_counter() - start, cpu_s=time.process_time() - cpu, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, fit_s=0, discovery_s=None, validation_s=None, questions_s=None, fallback_s=None)
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ('prep', 'antagonist', 'crown'):
        ap.add_argument(name)
    ap.add_argument('--material', default='3Y')
    ap.add_argument('--sinter-factor', type=float)
    ap.add_argument('--units', choices=['mm', 'um'])
    ap.add_argument('--contract')
    ap.add_argument('--output')
    ap.add_argument('--export', help='new directory for source-bound 3MF export')
    a = ap.parse_args()
    paths = {k: getattr(a, k) for k in ('prep', 'antagonist', 'crown')}
    try:
        r = check(paths, a.material, a.contract, a.units, a.sinter_factor)
        if a.export:
            from export_bridge import export
            r['rules']['export'] = export(paths, r, a.export, a.contract, a.sinter_factor)
            r = finish(r, time.perf_counter() - r['cost']['wall_s'], time.process_time() - r['cost']['cpu_s'])
        if a.output:
            dump(a.output, r)
        print(json.dumps(r, indent=2, ensure_ascii=False, allow_nan=False))
        return {'PASS': 0, 'FAIL': 2, 'UNKNOWN': 3}[r['verdict']]
    except (ValueError, KeyError, OSError, TypeError, IndexError) as e:
        r = dict(verdict='UNKNOWN', input_rejected=True, reason=str(e), claim_type='capability')
        if a.output:
            dump(a.output, r)
        print(json.dumps(r, ensure_ascii=False, allow_nan=False))
        return 4
if __name__ == '__main__':
    sys.exit(main())
