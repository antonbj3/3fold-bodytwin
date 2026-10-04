"""Capability contract: geometric rule decisions carry scope and refusal reasons."""
import datetime as dt
import json
import resource
import time
from pathlib import Path
import numpy as np
from .geometry import load, faces, continuous, graph_task, TOL
from .vendor.checks import all_checks
RULES = ['material_wall', 'film_min', 'film_max', 'ball_milling', 'insertion', 'occlusal_contact', 'approximal_contact']
IFU_URL = 'https://www.kuraraynoritake.eu/media/pdfs/KATANA_Zirconia_HT_TI-010-ver.19_WEB250825.pdf'
PROFILES = {'katana-ht': dict(product='KATANA Zirconia HT', source=IFU_URL, locator='TI-010-EN ver.019 08/2025, PDF page 2, thickness table', minimum_wall_mm={'anterior': 0.4, 'posterior': 0.5}, scope='named product, crown; no material-class generalization'), '3Y': dict(product=None, source=None, minimum_wall_mm={}, scope='Class alone does not select a product IFU'), '4Y': dict(product=None, source=None, minimum_wall_mm={}, scope='No product profile supplied'), '5Y': dict(product=None, source=None, minimum_wall_mm={}, scope='No product profile supplied'), 'LS2': dict(product=None, source=None, minimum_wall_mm={}, scope='No product profile supplied')}

def rule(status='UNKNOWN', reason='Input not supplied', **kw):
    return dict(status=status, reason=reason, resolution='PER_SURFACE_REGION', manual_adjustment_minutes=None, manual_adjustment_source='UNKNOWN: no timed lab adjustment records', **kw)

def check(paths, material='3Y', contract=None, units=None, axis=None, common_frame=False, round_version='R1'):
    start = time.perf_counter()
    cpu = time.process_time()
    meta = {} if contract is None else json.loads(Path(contract).read_text())
    if not isinstance(meta, dict):
        raise ValueError('Contract must be an object')
    if units is not None and meta.get('units') not in (None, units):
        raise ValueError('CLI/contract units conflict')
    unit = units or meta.get('units')
    scale = {'mm': 1.0, 'um': 0.001}.get(unit, 1.0)
    declared_unit = unit in ('mm', 'um')
    frame = common_frame or meta.get('common_frame_confirmed') is True
    axis = axis if axis is not None else meta.get('contact_axis')
    out = dict(schema='designgate-0.1', claim_type='capability', operator_round=round_version, created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), inputs={}, material=material, ifu_profile=PROFILES.get(meta.get('ifu_profile', material), {}), provenance=dict(contract_sha256=None if contract is None else __import__('hashlib').sha256(Path(contract).read_bytes()).hexdigest(), units=unit or 'UNKNOWN', unit_source='caller declaration; STL not authoritative', frame_declared=frame, contact_axis=axis, input_geometry_error_mm=meta.get('geometry_error_mm'), physical_geometry_status=meta.get('physical_geometry_status', 'UNKNOWN')), rules={k: rule() for k in RULES}, complete_restoration_eligibility='UNKNOWN: full restoration, manufacturing and physical validation not certified')
    data = {}
    for (k, p) in paths.items():
        if p is None:
            continue
        data[k] = load(p, scale)
        out['inputs'][k] = {t: data[k][t] for t in ('path', 'sha256', 'diagnostics')}
    bind = meta.get('input_sha256', {})
    if not isinstance(bind, dict):
        raise ValueError('input_sha256 must be an object')
    for (k, v) in bind.items():
        if k not in data or data[k]['sha256'] != v:
            raise ValueError('Hash-bound contract mismatch: ' + k)
    regions = meta.get('regions', {})
    if regions and any((k not in bind for k in data)):
        raise ValueError('Region metadata requires hashes for every input')
    if not declared_unit or not frame:
        why = 'Physical units UNKNOWN' if not declared_unit else 'Common registration frame UNKNOWN'
        out['rules'] = {k: rule(reason=why) for k in RULES}
        return finish(out, start, cpu)
    if any((x['diagnostics']['degenerate_faces'] or x['diagnostics']['duplicate_faces'] for x in data.values())):
        out['rules'] = {k: rule(reason='Degenerate/duplicate faces require explicit source repair') for k in RULES}
        return finish(out, start, cpu)
    profile = out['ifu_profile']
    indication = meta.get('indication')
    wall = profile.get('minimum_wall_mm', {}).get(indication)
    regions = meta.get('regions', {})
    c = data['crown']
    p = data['prep']
    a = data['antagonist']

    def region(name, d):
        return faces(d, regions[name])[0]
    gap_bounds = meta.get('film_limits_mm')
    if gap_bounds is not None:
        if not isinstance(gap_bounds, list) or len(gap_bounds) != 2 or any((type(x) not in (float, int) or not np.isfinite(x) for x in gap_bounds)) or (gap_bounds[0] < 0) or (gap_bounds[1] < gap_bounds[0]):
            raise ValueError('Invalid film acceptance interval')
    if wall is None:
        out['rules']['material_wall'] = rule(reason='Exact IFU product/indication required', source=profile.get('source'))
    graph = None
    if axis is not None and all((k in regions for k in ('intaglio', 'exterior', 'preparation'))):
        try:
            (task, design, wit) = graph_task(region('preparation', p), region('intaglio', c), region('exterior', c), axis, *(gap_bounds or [0, 1000.0]), wall or 0)
            graph = all_checks(task, design)
            out['inherited_graph_checks'] = graph
            if not graph.get('checks'):
                raise ValueError('Inherited graph validity refusal: ' + graph.get('reason', 'UNKNOWN'))
            if wall is not None:
                r = graph['checks']['wall']
                out['rules']['material_wall'] = rule(r['status'], 'Inherited entire graph bound', minimum_interval_mm=[r['lower_mm'], r['upper_mm']], witness=wit, source=profile['source'], threshold_mm=wall, scope=r['scope'])
            if gap_bounds is not None:
                for name in ('film_min', 'film_max'):
                    r = graph['checks'][name]
                    delta = np.asarray(design['inner_vertices'])[:, 2] - np.asarray(task['preparation_z'])
                    i = int(np.argmin(delta) if name == 'film_min' else np.argmax(delta))
                    from .geometry import basis
                    R = basis(axis)
                    point = np.asarray(design['inner_vertices'])[i] @ R
                    prep_point = np.r_[task['xy'][i], task['preparation_z'][i]] @ R
                    out['rules'][name] = rule(r['status'], 'Axial nominal film against caller lab limits', value_mm=r['value'], threshold_mm=gap_bounds[0 if name == 'film_min' else 1], scope=r['scope'], witness=dict(intaglio_point_mm=point.tolist(), preparation_point_mm=prep_point.tolist()), source=meta.get('film_limits_source', 'caller lab protocol; PHENOMENOLOGICAL'))
            for (key, orig) in [('insertion', 'insertion'), ('ball_milling', 'milling')]:
                r = graph['checks'][orig]
                if meta.get('geometry_scope') == 'matched_height_graph_patches':
                    out['rules'][key] = rule(r['status'], 'Inherited restricted patch proof', scope=r['scope'])
                else:
                    out['rules'][key] = rule(reason='Full restoration exceeds inherited patch proof', scope=r['scope'])
        except ValueError as e:
            out['graph_adapter_refusal'] = str(e)
    if graph is None:
        for k in ('material_wall', 'film_min', 'film_max', 'ball_milling', 'insertion'):
            if k == 'material_wall' and wall is None:
                continue
            out['rules'][k] = rule(reason=out.get('graph_adapter_refusal', 'Missing axis/regions for inherited graph contract'))
    if gap_bounds is None:
        for k in ('film_min', 'film_max'):
            out['rules'][k] = rule(reason='Lab film min/max acceptance limits required')
    band = meta.get('occlusal_gap_limits_mm')
    if band is not None and (not isinstance(band, list) or len(band) != 2 or any((type(x) not in (float, int) or not np.isfinite(x) for x in band)) or (band[1] < band[0])):
        raise ValueError('Invalid occlusal gap interval')
    if axis is None:
        out['rules']['occlusal_contact'] = rule(reason='Declared contact axis required')
    else:
        (u, uid) = faces(a, regions.get('antagonist'))
        (l, lid) = faces(c, regions.get('occlusal', regions.get('exterior')))
        g = continuous(u, l, axis)
        status = 'UNKNOWN'
        reason = 'Caller contact gap limits required'
        if band is not None and g['minimum_gap_mm'] is not None and (not g['excluded_projected_faces']):
            x = g['minimum_gap_mm']
            status = 'PASS' if band[0] - TOL <= x <= band[1] + TOL else 'FAIL'
            reason = 'Continuous minimum against caller gap limits'
        if g['minimum_gap_mm'] is not None:
            g['witness']['upper_STL_face'] = int(uid[g['witness']['upper_region_face']])
            g['witness']['lower_STL_face'] = int(lid[g['witness']['lower_region_face']])
        out['rules']['occlusal_contact'] = rule(status, reason, measurement=g, limits_mm=band, source=meta.get('occlusal_limits_source', 'caller lab protocol; PHENOMENOLOGICAL'))
    if 'neighbor' not in data:
        out['rules']['approximal_contact'] = rule(reason='Optional neighbor mesh/region and local contact policy required')
    if round_version == 'R2':
        out = general_rules(out, data, meta, regions, wall, gap_bounds)
    return finish(out, start, cpu)

def general_rules(out, data, meta, regions, wall, gap_bounds):
    from .distance import enclosure, planar_exterior, tangent_ball_obstruction, insertion_obstruction
    import trimesh
    c = data['crown']
    p = data['prep']
    if all((x in regions for x in ('intaglio', 'exterior'))):
        (inner, ii) = faces(c, regions['intaglio'])
        (outer, oi) = faces(c, regions['exterior'])
        target = trimesh.Trimesh(outer.reshape(-1, 3), np.arange(3 * len(outer)).reshape(-1, 3), process=False)
        if wall is not None:
            d = enclosure(inner, target, minimum_required=wall, source_face_ids=ii)
            d['minimum_witness']['target_STL_face'] = int(oi[d['minimum_witness']['target_region_face']])
            out['rules']['material_wall'] = rule(d['minimum_status'], '3D Euclidean region separation against exact product IFU', measurement=d, witness=d['minimum_witness'], threshold_mm=wall, source=out['ifu_profile']['source'], scope='Labelled inner/exterior regions only; rim excluded if caller did not label it; complete crown remains UNKNOWN')
        if gap_bounds is not None and 'preparation' in regions:
            (pt, pi) = faces(p, regions['preparation'])
            target = trimesh.Trimesh(pt.reshape(-1, 3), np.arange(3 * len(pt)).reshape(-1, 3), process=False)
            solid = p['diagnostics']['watertight'] and p['diagnostics']['winding_consistent'] and (meta.get('prep_solid_integrity_asserted') is True)
            planar = planar_exterior(inner, target, meta.get('preparation_outward_axis', [0, 0, 0]))
            if solid:
                target = p['mesh']
            if solid or planar:
                d = enclosure(inner, target, minimum_required=gap_bounds[0], maximum_allowed=gap_bounds[1], signed=solid, source_face_ids=ii)
                for (key, which) in [('film_min', 'minimum'), ('film_max', 'maximum')]:
                    w = d[which + '_witness'].copy()
                    w['target_STL_face'] = w['target_region_face'] if solid else int(pi[w['target_region_face']])
                    out['rules'][key] = rule(d[which + '_status'], '3D shortest signed/declared-exterior nominal gap against caller lab limits', measurement=d, witness=w, threshold_mm=gap_bounds[0 if key == 'film_min' else 1], source=meta.get('film_limits_source', 'caller protocol; PHENOMENOLOGICAL'), scope='Nominal shortest distance over labelled intaglio; no measured seated cement film', sign_basis='closed prep material' if solid else 'declared planar outward halfspace')
            else:
                for key in ('film_min', 'film_max'):
                    out['rules'][key] = rule(reason='Signed film needs asserted prep solid integrity or declared planar outward halfspace')
    solid_crown = c['diagnostics']['watertight'] and c['diagnostics']['winding_consistent'] and (c['mesh'].volume > 0) and (meta.get('crown_solid_integrity_asserted') is True)
    if meta.get('geometry_scope') != 'matched_height_graph_patches':
        out['rules']['ball_milling'] = rule(reason='Closed oriented crown integrity, intaglio labels and final-design ball radius required')
        out['rules']['insertion'] = rule(reason='Closed crown integrity and declared extraction axis/travel required')
    if solid_crown and 'intaglio' in regions and (meta.get('ball_radius_final_mm') is not None):
        d = tangent_ball_obstruction(c, regions['intaglio'], meta['ball_radius_final_mm'])
        out['rules']['ball_milling'] = rule(**d)
    if solid_crown and meta.get('extraction_axis') is not None and (meta.get('extraction_travel_mm') is not None):
        d = insertion_obstruction(p, c, meta['extraction_axis'], meta['extraction_travel_mm'])
        out['rules']['insertion'] = rule(**d)
    if 'neighbor' in data:
        axis = meta.get('neighbor_contact_axis')
        band = meta.get('approximal_gap_limits_mm')
        if axis is None or band is None or (not all((k in regions for k in ('approximal', 'neighbor')))):
            out['rules']['approximal_contact'] = rule(reason='Neighbor axis, region labels and lab min/max gap limits required')
        else:
            if not isinstance(band, list) or len(band) != 2 or any((type(x) not in (int, float) or not np.isfinite(x) for x in band)) or (band[0] > band[1]):
                raise ValueError('Invalid proximal gap limits')
            (nt, ni) = faces(data['neighbor'], regions['neighbor'])
            (ct, ci) = faces(c, regions['approximal'])
            g = continuous(nt, ct, axis)
            status = 'UNKNOWN'
            if g['minimum_gap_mm'] is not None and (not g['excluded_projected_faces']):
                status = 'PASS' if band[0] - TOL <= g['minimum_gap_mm'] <= band[1] + TOL else 'FAIL'
                g['witness']['upper_STL_face'] = int(ni[g['witness']['upper_region_face']])
                g['witness']['lower_STL_face'] = int(ci[g['witness']['lower_region_face']])
            out['rules']['approximal_contact'] = rule(status, 'Continuous local projected proximal gap against caller limits', measurement=g, limits_mm=band, source=meta.get('approximal_limits_source', 'caller lab protocol; PHENOMENOLOGICAL'))
    return out

def finish(out, start, cpu):
    statuses = [x['status'] for x in out['rules'].values()]
    out['verdict'] = 'FAIL' if 'FAIL' in statuses else 'UNKNOWN' if 'UNKNOWN' in statuses else 'PASS'
    out['cost'] = dict(wall_s=time.perf_counter() - start, cpu_s=time.process_time() - cpu, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, fit_s=0, manual_adjustment_minutes=None, discovery_s=None, fallback_s=None)
    return out
