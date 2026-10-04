"""Bite + crown -> contact, force bounds and research actuation proposals.

Operator bodies are preserved in _vendor. This file performs validation,
identity binding, unit conversion-free composition and result formatting only.
No contact area is converted to force. Caller arrays are never mutated.
"""
import hashlib, json
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.sparse import csr_matrix
from ._vendor import contact, mesh, optimize, reachability, tandlast, height_force, regional, deformable, force_port, local_boundary, height_uncertainty
from ._solver import linprog
optimize.linprog = linprog
regional.linprog = linprog

class ContractError(ValueError):
    """Invalid or mismatched input. Missing measurements are reported as UNKNOWN."""

def clean(x):
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, Fraction):
        return str(x)
    if isinstance(x, float) and (not np.isfinite(x)):
        return None
    return x

def _number(x, name, low=0.0):
    if isinstance(x, bool):
        raise ContractError(name + ' must be numeric')
    try:
        y = float(x)
    except (TypeError, ValueError):
        raise ContractError(name + ' must be numeric') from None
    if not np.isfinite(y) or y < low:
        raise ContractError(name + ' outside finite domain')
    return y

def _array(x, name, cols=None, finite=True):
    try:
        a = np.asarray(x, dtype=float)
    except (TypeError, ValueError):
        raise ContractError(name + ' invalid array') from None
    if cols is not None and (a.ndim != 2 or a.shape[1] != cols):
        raise ContractError(name + ' invalid shape')
    if finite and (not np.isfinite(a).all()):
        raise ContractError(name + ' nonfinite')
    return a.copy()

def _faces(x, n, name):
    a = _array(x, name, 3)
    if not len(a) or np.any(a != np.floor(a)) or a.min() < 0 or (a.max() >= n):
        raise ContractError(name + ' invalid indices')
    return a.astype(int)

def _fdi(t):
    if isinstance(t, bool) or not isinstance(t, (int, np.integer)) or t // 10 not in (1, 2, 3, 4) or (not 1 <= t % 10 <= 8):
        raise ContractError('Invalid FDI')
    return int(t)

def geometry_sha256(crown):
    """Digest geometry AND identities; lengths are mm, no tolerance-based identity."""
    h = hashlib.sha256()
    for k in ['case_id', 'frame_id', 'pose_id', 'unit', 'tooth_fdi', 'representation', 'margin_z_mm', 'wall_lower_mm', 'relief_cap_mm']:
        h.update(json.dumps(crown.get(k), sort_keys=True).encode() + b'\x00')
    for k in ['vertices_mm', 'faces', 'face_roles', 'xy_mm', 'outer_z_mm', 'inner_z_mm', 'actuation_basis']:
        if k in crown:
            a = np.asarray(crown[k], dtype='<f8')
            h.update(k.encode() + str(a.shape).encode() + a.tobytes())
    return h.hexdigest()

def bite_sha256(bite):
    h = hashlib.sha256()
    for k in ['case_id', 'frame_id', 'pose_id', 'unit', 'source_locator', 'evidence_kind']:
        h.update(json.dumps(bite.get(k), sort_keys=True).encode() + b'\x00')
    for k in ['xy_mm', 'faces', 'antagonist_z_mm', 'reference_gap_mm']:
        a = np.asarray(bite[k], dtype='<f8')
        a = np.where(np.isnan(a), np.nan, a)
        h.update(k.encode() + str(a.shape).encode() + a.tobytes())
    return h.hexdigest()

def _constraint_sha256(cache):
    h = hashlib.sha256()
    for k in ['data', 'indices', 'indptr', 'shape', 'b']:
        a = np.asarray(cache[k], dtype='<f8')
        h.update(k.encode() + str(a.shape).encode() + a.tobytes())
    h.update(json.dumps(cache['statistics'], sort_keys=True).encode())
    return h.hexdigest()

def _binding(bite, crown):
    if bite.get('schema') != 'registered-bite-v1' or crown.get('schema') != 'research-crown-v1':
        raise ContractError('Unsupported input schema')
    for k in ['case_id', 'frame_id', 'pose_id', 'unit']:
        if not bite.get(k) or bite[k] != crown.get(k):
            raise ContractError('Bite/crown identity mismatch: ' + k)
    if bite['unit'] != 'mm':
        raise ContractError('Length unit must be mm')
    _fdi(crown['tooth_fdi'])
    if not bite.get('source_locator'):
        raise ContractError('Bite source locator required')
    if bite.get('evidence_kind') not in ['REGISTERED_SCAN', 'SIMULATION']:
        raise ContractError('Bite evidence kind must be explicit')

def _prepare(bite, crown):
    """Dictionary adaptation of X95 prepare; inherited all-hit stencil reused."""
    xy = _array(bite['xy_mm'], 'bite.xy_mm', 2)
    ff = _faces(bite['faces'], len(xy), 'bite.faces')
    ceiling = _array(bite['antagonist_z_mm'], 'antagonist_z_mm', finite=False)
    rg = _array(bite['reference_gap_mm'], 'reference_gap_mm', finite=False)
    if ceiling.shape != (len(xy),) or rg.shape != ceiling.shape or np.isinf(ceiling).any() or np.isinf(rg).any():
        raise ContractError('Gap field shape or infinity')
    v = _array(crown['vertices_mm'], 'vertices_mm', 3)
    f = _faces(crown['faces'], len(v), 'crown.faces')
    roles = np.asarray(crown['face_roles'])
    if roles.shape != (len(f),) or not np.isin(roles, [0, 1, 2]).all():
        raise ContractError('face_roles: 0 exterior, 1 intaglio, 2 margin required')
    extf = f[roles == 0]
    if not len(extf):
        raise ContractError('Exterior faces absent')
    area2 = np.cross(xy[ff[:, 1]] - xy[ff[:, 0]], xy[ff[:, 2]] - xy[ff[:, 0]])
    if np.any(abs(area2) <= 1e-12):
        raise ContractError('Degenerate query triangle')
    (H, qr, z, owner, top) = mesh.stencil(v, extf, xy)
    g0 = ceiling - top
    protected = np.unique(f[roles != 0])
    nz = np.cross(v[extf[:, 1], :2] - v[extf[:, 0], :2], v[extf[:, 2], :2] - v[extf[:, 0], :2])
    nonupper = np.unique(extf[nz <= 1e-12])
    editable = np.setdiff1d(np.unique(extf), np.union1d(protected, nonupper))
    margin = _number(crown['margin_z_mm'], 'margin_z_mm', -np.inf)
    taper = np.clip((v[:, 2] - margin) / 0.6, 0, 1)
    taper[np.setdiff1d(np.arange(len(v)), editable)] = 0.0
    area = np.linalg.norm(np.cross(v[extf[:, 1]] - v[extf[:, 0]], v[extf[:, 2]] - v[extf[:, 0]]), axis=1) / 2
    vw = np.bincount(extf.ravel(), weights=np.repeat(area / 3, 3), minlength=len(v))
    common = ff[np.isfinite(g0[ff]).all(1) & np.isfinite(rg[ff]).all(1)]
    good = np.zeros(len(xy), bool)
    good[np.unique(common)] = True
    weights = np.bincount(common.ravel(), weights=np.repeat(abs(np.cross(xy[common[:, 1]] - xy[common[:, 0]], xy[common[:, 2]] - xy[common[:, 0]])) / 6, 3), minlength=len(xy))
    maxdown = float(v[extf[nz < -1e-12], 2].max()) if np.any(nz < -1e-12) else -1e+30
    return dict(vertices=v, faces=f, face_roles=roles, xy=xy, grid_faces=ff, ceiling=ceiling, reference_gap=rg, original_gap=g0, margin_z=margin, protected=protected, editable=editable, taper=taper, vertex_area_weights=vw, query_weights=weights, H=H, hit_queries=qr, hit_z=z, owner=owner, query_good=good, max_downward_exterior_z=maxdown)

def _force_bounds(bite, crown):
    teeth = sorted({_fdi(t) for t in bite.get('teeth_fdi', [crown['tooth_fdi']])} | {crown['tooth_fdi']})
    pairs = bite.get('tooth_pairs')
    if pairs is None:
        return {'status': 'UNKNOWN', 'reason': 'Declared tooth-pair incidence/gap enclosures absent', 'per_tooth': [dict(fdi=t, classification='UNKNOWN', force_interval_N=None, resolution='PER_TOOTH') for t in teeth]}
    if not pairs:
        raise ContractError('Empty declared tooth-pair graph')
    if bite.get('pair_model') != 'AXIAL_POSITIVE_SUPPORT_V1':
        raise ContractError('Unsupported force-hull model; Coulomb and wrench models need separate ports')
    edges = []
    for p in pairs:
        (u, v) = (_fdi(p['upper_fdi']), _fdi(p['lower_fdi']))
        if u // 10 not in (1, 2) or v // 10 not in (3, 4):
            raise ContractError('Upper/lower FDI orientation invalid')
        try:
            (lo, hi) = map(Fraction, map(str, p['gap_interval_mm']))
        except (ValueError, TypeError, ZeroDivisionError):
            raise ContractError('Malformed gap enclosure') from None
        if lo > hi:
            raise ContractError('Reversed gap enclosure')
        edges.append((u, v, lo, hi))
        teeth.extend([u, v])
    if len({e[:2] for e in edges}) != len(edges):
        raise ContractError('Duplicate tooth pair')
    total = bite.get('total_force_interval_N')
    if total is not None:
        ev = bite.get('total_force_evidence', {})
        if ev.get('kind') not in ['MEASURED', 'SIMULATED'] or not ev.get('locator'):
            raise ContractError('Total force requires explicit evidence and locator')
        if ev['kind'] == 'MEASURED' and (not ev.get('calibration_locator')):
            raise ContractError('Absolute total force requires calibration locator')
        if len(total) != 2:
            raise ContractError('Total-force interval requires two endpoints')
        total = [_number(x, 'total_force_interval_N') for x in total]
        if total[0] <= 0 or total[0] > total[1]:
            raise ContractError('Force-controlled model requires strictly positive total force interval')
    rows = []
    for t in sorted(set(teeth)):
        a = tandlast.classify(edges, [t])
        share = list(map(float, a['share_hull_pp']))
        hull = [0.0, 0.0] if a['classification'] == 'NEVER' else ([total[0], total[1]] if share == [100.0, 100.0] else [0.0, total[1]]) if total else [0.0, None]
        rows.append(dict(fdi=t, classification=a['classification'], force_interval_N=hull, force_unit='N', share_hull_pp=share, strict_positive_in_all_model_states=a['must_bear'], infimum_attained=a['infimum_attained'], necessity_margin_mm=None if a['necessity_margin_mm'] is None else str(a['necessity_margin_mm']), resolution='PER_TOOTH', scope='Declared pairs; normal-only positive-support axial model', physical_classification='UNKNOWN'))
    return dict(status='CONDITIONAL_MODEL_HULL', per_tooth=rows, total_force_interval_N=total, upper_bound_missing=total is None, rigor='Exact rational gap classification; force hull follows conserved total, no measured stiffness. NEVER refers only to declared-pair model.', timescale='SIMULTANEOUS')

def query_height(binding, state, request):
    """X82 whole-tooth height query; never transports it to arbitrary local edits."""
    for (k, source) in [('case_id', 'case'), ('response_geometry_sha256', 'geometry_sha256')]:
        if binding.get(k) != state.get(source):
            raise ContractError('Height response identity mismatch: ' + k)
    if not binding.get('pose_id') or binding.get('pose_id') != binding.get('calibration_pose_id'):
        raise ContractError('Height response pose mismatch')
    if not binding.get('matched_load_rate_history'):
        return dict(status='UNKNOWN', reason='Load/rate/history not matched')
    if request.get('operation') != 'uniform_tooth_height':
        return dict(status='UNKNOWN', reason='Tooth response does not identify regional actuation')
    if 'target_N' in request:
        return height_force.inverse(state, _fdi(request['fdi']), request['target_N'], request.get('target_halfwidth_N', 1.0), request.get('load_N'))
    i = state['predicted_fdi'].index(_fdi(request['fdi']))
    edit = height_force.make_edit(state, i, float(request['height_change_mm']))
    if 'load_N' in request:
        edit['load_N'] = request['load_N']
    return height_force.predict(state, edit)

def _missing_measurements():
    return [dict(quantity='absolute force per tooth', resolution='PER_TOOTH', measurement='Same specimen and registered MIP pose: calibrated per-tooth axial reaction channels with gain/error bounds, synchronized total force and denominator. T-Scan percentages alone do not supply absolute Newton values.'), dict(quantity='regional height to force', resolution='PER_SURFACE_REGION', measurement='Baseline and +/- known height for each named actuation basis direction at matched wrench/rate/history; bounded force and height errors, positive-contact corridor and independent fourth held-out height. Three probes do not bound arbitrary nonlinear tissue.'), dict(quantity='loaded contact and preload', resolution='PER_POINT', measurement='Independently registered loaded antagonist/crown geometry at known load, scanner/pose uncertainty and support initial state; repeat after manufacture.'), dict(quantity='lower-arch regional force from upper observations', resolution='PER_SURFACE_REGION', measurement='Synchronized upper-region force plus signed net cross-boundary force on the named pair graph, calibrated to a common total (X94).')]

def _local_constraints(crown, t, r, enabled):
    if not enabled:
        return (None, None)
    cache = crown.get('local_boundary_certificate')
    if cache is None:
        (A, b, stats, _) = local_boundary.build(contact, t, r)
    else:
        if cache.get('geometry_sha256') != geometry_sha256(crown) or cache.get('relief_cap_mm') != r['relief_cap_mm']:
            raise ContractError('Local-boundary certificate identity or cap mismatch')
        if not cache.get('source_locator'):
            raise ContractError('Local-boundary certificate source locator absent')
        if cache.get('payload_sha256') != _constraint_sha256(cache):
            raise ContractError('Local-boundary certificate payload drift')
        shape = tuple(cache['shape'])
        A = csr_matrix((cache['data'], cache['indices'], cache['indptr']), shape=shape)
        b = np.asarray(cache['b'], float)
        stats = cache['statistics']
        if shape != (len(b), len(t['vertices'])) or not np.isfinite(A.data).all() or (not np.isfinite(b).all()):
            raise ContractError('Malformed local-boundary certificate')
    return ((A, b), stats)

def _local_gates(gates, d, extra, stats):
    if extra is None:
        return gates
    (A, b) = extra
    slack = b - A @ d
    lo = float(slack.min()) if len(slack) else None
    good = stats['original_mixed_order_overlap_pairs'] == 0 and (lo is None or lo >= -1e-08) and (gates['protected_identity_error_mm'] == 0.0)
    return dict(gates, local_vertical_interval_order_pass=good, local_interval_slack_lower_mm=lo, cap_implied_pruned_slack_lower_mm=stats['cap_implied_minimum_slack_lower_mm'], removal_status_R3='CONDITIONAL_ON_ORIGINAL_VALID_SOLID' if good else 'UNKNOWN_OR_ORIGINAL_ORDER_FAIL')

def _regional_full_mesh(bite, crown, t, r, extra, stats):
    packet = bite.get('regional_actuation')
    if packet is None:
        return dict(status='UNKNOWN', reason='Regional probes and actuation basis absent')
    p = packet['probes']
    B = _array(packet['basis'], 'regional basis')
    if B.ndim != 2 or B.shape[0] != len(t['vertices']) or B.shape[1] != 2 or (B.min() < 0) or np.any(B[t['protected']] != 0):
        raise ContractError('Two nonnegative regional basis columns with protected vertices fixed required')
    if np.linalg.matrix_rank(B) != 2:
        raise ContractError('Regional basis must contain two independent actual mesh operations')
    digest = hashlib.sha256(B.astype('<f8').tobytes()).hexdigest()
    for (k, val) in [('case', bite['case_id']), ('frame', bite['frame_id']), ('pose_id', bite['pose_id']), ('geometry_sha256', geometry_sha256(crown)), ('bite_sha256', bite_sha256(bite)), ('basis_sha256', digest)]:
        if p.get(k) != val:
            raise ContractError('Regional probe identity mismatch: ' + k)
    if p.get('kind') not in ['SIMULATION', 'MEASURED']:
        raise ContractError('Regional evidence kind required')
    if p['kind'] == 'MEASURED' and (not p.get('calibration_locator') or not p.get('closure_locator')):
        return dict(status='UNKNOWN', reason='Sensor-to-reaction calibration or linear-closure evidence absent')
    regions = packet.get('reaction_region_ids')
    mapping = packet.get('reaction_tooth_fdi')
    if regions is None or mapping is None or len(regions) != len(p['baseline_force_N']) or (len(mapping) != len(regions)) or (len(set(regions)) != len(regions)):
        raise ContractError('Unique named reaction regions and tooth map required')
    for tooth in mapping:
        _fdi(tooth)
    if len(regions) < 2:
        raise ContractError('X85 port needs two target regions')
    model = height_uncertainty.fit_bounded_heights(p, packet.get('final_coefficient_error_mm'))
    if model.get('status') != 'CONDITIONAL_CALIBRATION':
        return model
    if np.asarray(model['J_N_per_mm']).shape != (len(regions), 2):
        raise ContractError('Regional reaction dimensions mismatch')
    u = np.asarray(model['final_coefficient_error_mm'])
    allow = np.minimum(r['wall_lower_mm'] - 0.5, t['taper'] * r['relief_cap_mm']) - B @ u
    G = B
    weights = t['vertex_area_weights']
    if extra is not None:
        (A, b) = extra
        AB = np.asarray(A @ B)
        G = np.vstack([B, AB])
        allow = np.r_[allow, b - abs(AB) @ u]
        weights = np.r_[weights, np.zeros(len(b))]
    a = regional.inverse(model, G, weights, allow)
    a['evidence_kind'] = p['kind']
    a['resolution'] = 'PER_SURFACE_REGION'
    a['reaction_region_ids'] = regions
    a['reaction_tooth_fdi'] = mapping
    a['physical_validity'] = 'UNKNOWN; declared sensor/linear closure and held-out validation must be independently checked'
    if a.get('status') != 'CONDITIONAL_TARGET_DESIGN':
        return a
    coefficients = np.asarray(a['coefficients_mm'])
    d = -B @ coefficients
    (new, gates) = mesh.wall_and_removal(t, d, r)
    gates = _local_gates(gates, d, extra, stats)
    a.update(height_change_per_vertex_mm=-d, edited_vertices_mm=new, geometry_gates=gates, source_contact=contact.compare(t['xy'], t['grid_faces'], mesh.gap(t, d), t['reference_gap']), error_contract=model, manufacturing_approval='UNKNOWN', positive_height_requires_new_manufacture=False, scope='Two target regions at inherited7..9 N experimental band; downward fixed-XY full-mesh operation, affine positive-contact closure only. Newton values are not inferred from geometry.')
    return a

def analyze(bite, crown, request=None):
    """Single public dictionary API. Raises ContractError for malformed identities.

    Missing observations return UNKNOWN in their own output branch. Geometry
    proposals are conditional research candidates, not manufacturing approval.
    """
    request = {} if request is None else dict(request)
    _binding(bite, crown)
    if crown.get('representation') != 'full_mesh':
        raise ContractError('This API accepts full_mesh; use regional roof calibration through the explicit regional_design port')
    t = _prepare(bite, crown)
    before = contact.compare(t['xy'], t['grid_faces'], t['original_gap'], t['reference_gap'])
    forces = _force_bounds(bite, crown)
    cap = _number(crown['relief_cap_mm'], 'relief_cap_mm')
    wall = _number(crown['wall_lower_mm'], 'wall_lower_mm')
    if cap > max(0, wall - 0.5):
        raise ContractError('Relief cap exceeds inherited 0.5 mm wall allowance')
    r = dict(relief_cap_mm=cap, wall_lower_mm=wall)
    (extra, local_stats) = _local_constraints(crown, t, r, request.get('local_boundary', False))
    hi = mesh.gap(t, t['taper'] * cap)
    possible = reachability.reachable(contact, t, t['original_gap'], hi, t['reference_gap'], True) if before['status'] == 'SCORED' else {'status': 'UNKNOWN_NO_COMMON_SUPPORT'}
    unbounded = reachability.reachable(contact, t, t['original_gap'], np.full(len(t['xy']), 1e+30), t['reference_gap'], False) if before['status'] == 'SCORED' else None
    status = 'NO_PROPOSAL_REQUESTED'
    adjustment = None
    if request.get('propose_adjustment', True) and before['status'] == 'SCORED':
        (d, opt) = optimize.solve(t, r, regional=True, extra=extra)
        (new, gates) = mesh.wall_and_removal(t, d, r)
        gates = _local_gates(gates, d, extra, local_stats)
        pg = mesh.gap(t, d)
        after = contact.compare(t['xy'], t['grid_faces'], pg, t['reference_gap'])
        ids = np.flatnonzero((t['face_roles'] == 0) & np.any(d[t['faces']] > 0, axis=1))
        faces = t['faces'][ids]
        adjustment = dict(status='CONDITIONAL_GEOMETRIC_PROPOSAL', objective='X95 nodal band surrogate then integrated displacement; not global minimum binary-mask repair', regional_height_adjustments=[dict(region_id=int(fi), vertex_ids=f.tolist(), height_change_mm=(-d[f]).tolist(), resolution='PER_SURFACE_REGION') for (fi, f) in zip(ids, faces)], height_change_per_vertex_mm=-d, edited_vertices_mm=new, faces=t['faces'], contact=after, gates=gates, optimization=opt, physical_force='UNKNOWN', force_response_applied=False, manufacturing_approval='UNKNOWN', roundoff_enclosure='MISSING')
        status = 'CONDITIONAL_GEOMETRY_RESULT'
    height = dict(status='UNKNOWN', reason='Same-specimen full response absent; regional edit not identifiable from tooth total')
    if 'height_query' in request:
        if 'height_response' not in bite:
            height = dict(status='UNKNOWN', reason='Height response absent')
        else:
            ctx = bite.get('force_context', {})
            if ctx.get('case_id') != bite['case_id'] or ctx.get('pose_id') != bite['pose_id'] or ctx.get('reference_crown_sha256') != geometry_sha256(crown) or (ctx.get('reference_bite_sha256') != bite_sha256(bite)):
                raise ContractError('Height response is not bound to this bite and crown')
            height = query_height(ctx, bite['height_response'], request['height_query'])
    out = dict(schema='occlusion-result-v1', claim_type='capability', status=status, case_id=bite['case_id'], frame_id=bite['frame_id'], pose_id=bite['pose_id'], geometry_sha256=geometry_sha256(crown), contact_map=dict(metrics=before, xy_mm=t['xy'], faces=t['grid_faces'], gap_mm=t['original_gap'], polygons_barycentric=contact.polygons(t['xy'], t['grid_faces'], t['original_gap'])[0], resolution='PER_POINT', definition='PL vertical gap band 0..0.1 mm; geometric proximity, not pressure', uncertainty=contact.enclosure(t['xy'], t['grid_faces'], t['original_gap'], t['reference_gap'], _number(request.get('pose_error_mm', 0.01), 'pose_error_mm')), rigorous_machine_enclosure='MISSING'), force_intervals=forces, proposed_adjustment=adjustment, contacts_requiring_buildup=dict(within_removal_cap=possible, without_depth_cap=unbounded, resolution='PER_SURFACE_REGION', interpretation='Positive-height witness: removal cannot restore an already-low source patch; new manufacture or buildup is required to reproduce the operational source band.'), height_force_query=height, deformable_contact=dict(status='UNKNOWN', reason='No bound same-specimen support/load/preload/friction scenario provided'), cannot_determine=_missing_measurements(), physical_validation='NOT_RUN', review_state='PENDING_INDEPENDENT_REVIEW', timescale='SIMULTANEOUS', coverage=dict(query_triangles=len(t['grid_faces']), retained_triangles=int(np.sum(np.isfinite(t['original_gap'][t['grid_faces']]).all(1) & np.isfinite(t['reference_gap'][t['grid_faces']]).all(1))), reason='Missing exterior/antagonist/source field support'))
    if 'deformable_scenario' in bite:
        sc = bite['deformable_scenario']
        for k in ['case_id', 'frame_id', 'pose_id']:
            if sc.get(k) != bite[k]:
                raise ContractError('Deformable scenario identity mismatch: ' + k)
        if sc.get('geometry_sha256') != geometry_sha256(crown):
            raise ContractError('Deformable geometry mismatch')
        if sc.get('evidence_kind') != 'SIMULATION':
            raise ContractError('X54 port is a declared conditional scenario; measured model admission requires independent validation')
        if sc['geometry'].get('case') != bite['case_id']:
            raise ContractError('Deformable patch case mismatch')
        pars = sc.get('parameters', {})
        if 'closure' not in pars or 'mu' not in pars:
            raise ContractError('Declare closure and friction for the conditional contact scenario')
        _number(pars['mu'], 'mu')
        _number(pars['closure'], 'closure', -np.inf)
        (ans, _) = deformable.solve(sc['geometry'], **pars)
        ans['default_support_assumed'] = pars.get('support_scale') is None
        ans['evidence_kind'] = 'SIMULATION'
        out['deformable_contact'] = ans
    if 'cross_boundary_observation' in bite:
        out['cross_boundary_force'] = force_port.project(bite['cross_boundary_observation'])
    out['coverage']['reference_supported_triangles'] = int(np.isfinite(t['reference_gap'][t['grid_faces']]).all(1).sum())
    out['coverage']['rejected_triangle_fraction'] = 1.0 - out['coverage']['retained_triangles'] / out['coverage']['query_triangles']
    out['coverage']['reference_common_area_fraction'] = before.get('coverage', 0.0)
    out['quantity_resolution'] = dict(gap_mm='PER_POINT', height_change_mm='PER_POINT', polygon_area_mm2='PER_SURFACE_REGION', force_interval_N='PER_TOOTH', total_force_interval_N='PER_ARCH', case_counts='POPULATION', constitutive_and_sensor_errors='PHENOMENOLOGICAL')
    out['local_boundary'] = dict(status='NOT_REQUESTED' if extra is None else 'CONDITIONAL_INTERVAL_ORDER', statistics=local_stats, rigorous_machine_enclosure='MISSING', original_solid_validity='UNKNOWN')
    out['regional_force_design'] = dict(status='UNKNOWN', reason='Same-specimen regional actuation probes absent')
    if request.get('regional_force_design'):
        out['regional_force_design'] = _regional_full_mesh(bite, crown, t, r, extra, local_stats)
    return clean(out)

def load_demo_input(root, row):
    """Read a hash-bound X95 small geometry derivative, no sibling dependencies."""
    p = Path(root) / row['package_input']
    if hashlib.sha256(p.read_bytes()).hexdigest() != row['input_sha256']:
        raise ContractError('Demo input hash drift')
    with np.load(p, allow_pickle=False) as a:
        t = dict(a)
    meta = dict(case_id=row['case_key'], frame_id=row['frame'], pose_id='source_registered', unit='mm')
    bite = dict(schema='registered-bite-v1', **meta, xy_mm=t['xy'], faces=t['grid_faces'], antagonist_z_mm=t['ceiling'], reference_gap_mm=t['reference_gap'], source_locator=row['license']['locator'], evidence_kind='REGISTERED_SCAN', teeth_fdi=[row['source_fdi']])
    crown = dict(schema='research-crown-v1', **meta, representation='full_mesh', tooth_fdi=row['source_fdi'], vertices_mm=t['vertices'], faces=t['faces'], face_roles=t['face_roles'], margin_z_mm=float(t['margin_z']), wall_lower_mm=row['wall_lower_mm'], relief_cap_mm=row['relief_cap_mm'])
    if 'local_constraints' in row:
        cc = row['local_constraints']
        path = Path(root) / cc['path']
        if hashlib.sha256(path.read_bytes()).hexdigest() != cc['sha256']:
            raise ContractError('Local constraint source hash drift')
        with np.load(path, allow_pickle=False) as a:
            cal = {k: a[k].copy() for k in a.files}
        crown['local_boundary_certificate'] = dict(geometry_sha256=geometry_sha256(crown), relief_cap_mm=crown['relief_cap_mm'], data=cal['A_data'], indices=cal['A_indices'], indptr=cal['A_indptr'], shape=cal['A_shape'], b=cal['b'], statistics=cc['statistics'], source_locator=cc['source_locator'])
        crown['local_boundary_certificate']['payload_sha256'] = _constraint_sha256(crown['local_boundary_certificate'])
    return (bite, crown)

def regional_design(binding, probes, B, weights, wall_allowance, band=(7.0, 9.0)):
    """Explicit X85 roof port: two target regions, positive affine contact only."""
    for k in ['case', 'frame', 'geometry_sha256', 'basis_sha256']:
        if not binding.get(k) or binding[k] != probes.get(k):
            raise ContractError('Regional calibration binding mismatch: ' + k)
    B = _array(B, 'actuation_basis')
    actual = hashlib.sha256(B.astype('<f8').tobytes()).hexdigest()
    if actual != binding['basis_sha256']:
        raise ContractError('Actual actuation basis differs')
    model = regional.fit_response(probes)
    if model.get('status') != 'CONDITIONAL_CALIBRATION':
        return model
    ans = regional.inverse(model, B, np.asarray(weights), np.asarray(wall_allowance), band=band)
    ans['scope'] = 'X85 PL roof port, two target regions; not a whole-crown force calibration'
    return clean(ans)
