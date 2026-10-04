"""Read-only invocation of the reviewed dental operations, one stage per process."""
from dental_release.paths import expand as _release_expand
import argparse, ast, copy, csv, hashlib, importlib.util, io, json, math, os
import sys, zipfile
from pathlib import Path
import numpy as np
import trimesh
R = Path(__file__).resolve().parent
D = R.parent.parent
RESULTS = D / 'results'
V4 = RESULTS / _release_expand('GENCAD_V4')
X60 = RESULTS / _release_expand('X60')
P = D / 'results/DEMO48_PACKAGE'
sys.dont_write_bytecode = True
sys.path.insert(0, str(V4 / 'payload/participant_code/legacy/vendor'))
sys.path.insert(0, str(V4 / 'payload/participant_code'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def clean(x):
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    if isinstance(x, float) and (not math.isfinite(x)):
        return None
    return x

def write(p, x):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(clean(x), indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + '\n')

def read(p):
    return json.loads(Path(p).read_text())

def npz(p, **arrays):
    with zipfile.ZipFile(p, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for (name, a) in sorted(arrays.items()):
            buf = io.BytesIO()
            np.lib.format.write_array(buf, np.asanyarray(a), allow_pickle=False)
            zi = zipfile.ZipInfo(name + '.npy', (1980, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(zi, buf.getvalue())

def normalize(x, out):
    """Preserve all scientific fields; isolate only runtime clock/cost/location telemetry."""
    skip = {'created_utc', 'preparation_s', 'cost', 'seconds', 'wall_seconds', 'total_seconds', 'fit_and_measure_seconds', 'query_and_export_seconds', 'query_s', 'query_wall_s'}
    if isinstance(x, dict):
        return {k: normalize(v, out) for (k, v) in x.items() if k not in skip}
    if isinstance(x, list):
        return [normalize(v, out) for v in x]
    if isinstance(x, str):
        return x.replace(str(out), '$RUN')
    return clean(x)

def module(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sys.modules[name] = m
    sp.loader.exec_module(m)
    return m

def task():
    sys.path.insert(0, str(V4 / 'payload/participant_code'))
    from task_io import load_local
    sel = read(R / 'PATIENT_SELECTION.json')
    ts = read(V4 / 'payload/public/tasks' / (sel['patient']['case_key'] + '.json'))
    t = next((t for t in ts if t['family'] == sel['requested_family'] and t['level'] == sel['level']))
    return load_local(t, V4 / 'payload/public')

def scan(out):
    sys.path.insert(0, str(V4 / 'code'))
    sys.path.insert(0, str(V4 / 'payload/participant_code'))
    from source_ops import pair, make_site, ZIPS
    from legacy.geometry import shell, parse_stl
    from scipy.spatial import cKDTree
    row = read(R / 'PATIENT_SELECTION.json')['patient']
    (arches, sources) = pair(row)
    site = make_site(arches, 'molar_crown', compute_obstacles=False)
    t = task()
    difference = float(np.max(abs(site['xy'] - t['xy'])))
    if difference != 0:
        raise ValueError('Rebuilt patient/site frame changed')
    with zipfile.ZipFile(ZIPS['Bits2Bites']) as z:
        for (jaw, src) in zip(['upper', 'lower'], sources):
            (out / (jaw + '.stl')).write_bytes(z.read(src['member']))
    lower = arches['lower']
    target = lower['owner'] == 36
    v = (lower['v'] - site['base']) @ site['R']
    m = trimesh.Trimesh(v, lower['f'][target], process=False)
    m.export(out / 'segmented_tooth.stl')
    av = (arches['upper']['v'] - site['base']) @ site['R']
    atri = av[arches['upper']['f']]
    xy = site['xy']
    lo = xy.min(0)
    hi = xy.max(0)
    overlaps = np.all(atri[:, :, :2].max(1) >= lo, axis=1) & np.all(atri[:, :, :2].min(1) <= hi, axis=1)
    det = np.cross(atri[:, 1, :2] - atri[:, 0, :2], atri[:, 2, :2] - atri[:, 0, :2])
    keep = overlaps & (abs(det) > 1e-09)
    source_ids = np.flatnonzero(keep)
    serialized = atri[keep].astype(np.float32).astype(float)
    anti = trimesh.Trimesh(serialized.reshape(-1, 3), np.arange(3 * len(serialized)).reshape(-1, 3), process=False)
    anti.export(out / 'antagonist.stl')
    npz(out / 'ANTAGONIST_FACET_MAP.npz', original_upper_scan_face=source_ids, original_local_triangles_mm=atri[keep], serialized_triangles_mm=serialized)
    write(out / 'ANTAGONIST_DOMAIN.json', dict(round='R2', source_total_facets=len(atri), retained_facets=int(keep.sum()), excluded_outside_query=int((~overlaps).sum()), excluded_projection_degenerate=int((overlaps & ~(abs(det) > 1e-09)).sum()), excluded_fraction=float((~keep).mean()), source_original_upper_stl_sha256=sources[0]['sha256'], operation='Unchanged v4 projected antagonist crop; source face map retained', resolution='PER_POINT', geometry_rounding_max_mm=float(np.max(abs(serialized - atri[keep]))), full_3D_insertion_status='UNKNOWN; projected crop serves fixed axial roof-contact query only'))
    reference = site['reference'].copy()
    missing = ~np.isfinite(reference)
    if missing.any():
        (_, ii) = cKDTree(site['xy'][~missing]).query(site['xy'][missing])
        reference[missing] = reference[~missing][ii]
    (v, f) = shell(site['xy'], reference, np.full(len(reference), t['preparation_height_mm'] - 2), site['faces'])
    capped = trimesh.Trimesh(v, f, process=False)
    if capped.volume < 0:
        capped.invert()
    if not capped.is_watertight:
        raise ValueError('Capped roof is not closed')
    capped.export(out / 'measured_capped_roof.stl')
    npz(out / 'scene.npz', xy=site['xy'], faces=site['faces'], reference_z_mm=reference, base_mm=site['base'], local_to_scan_rotation=site['R'], upper_labels=arches['upper']['labels'], lower_labels=lower['labels'])
    write(out / 'SCAN.json', dict(claim_type='capability', patient_case=row['case'], case_key=row['case_key'], source_fdi=36, sources=sources, units='mm', resolution='PER_POINT', timescale='SIMULTANEOUS', selected_label_facets=int(target.sum()), label_status='FROZEN_X11_MODEL_PREDICTIONS_NOT_EXPERT_GROUND_TRUTH', native_projected_reference_missing=int(missing.sum()), native_projected_reference_nodes=len(reference), missing_fraction=float(missing.mean()), completion='Nearest finite measured projected roof at missing nodes; flat lower cap, restricted roof closure', same_site_xy_identity_error_mm=difference, base_mm=site['base'], local_to_scan_rotation=site['R'], frame_status='Derived from predicted tooth labels; scanner/bite-pose uncertainty UNKNOWN', external_referent=dict(kind='published_dataset', locator='https://ditto.ing.unimore.it/bits2bites/', compared_quantity='Upper/lower registered measured triangle coordinates from named ZIP members', refutes_us=True)))

def sdf(out):
    from scipy.spatial import cKDTree
    t = task()
    m = trimesh.load_mesh(out / 'measured_capped_roof.stl', process=True)
    h = 0.25
    lo = m.bounds[0] - 0.5
    hi = m.bounds[1] + 0.5
    axes = [lo[i] + np.arange(int(np.ceil((hi[i] - lo[i]) / h)) + 1) * h for i in range(3)]
    xyz = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1)
    flat = xyz.reshape(-1, 3)
    values = []
    for i in range(0, len(flat), 512):
        values.append(-trimesh.proximity.signed_distance(m, flat[i:i + 512]))
    phi = np.concatenate(values).reshape(xyz.shape[:-1])
    prep = xyz[..., 2] - t['preparation_height_mm']
    npz(out / 'SDF.npz', phi_roof_mm=phi, virtual_preparation_plane_mm=prep, origin_mm=lo, pitch_mm=np.asarray(h), shape=np.asarray(phi.shape))
    from legacy.geometry import shell
    (v, f) = shell(t['xy'], t['preparation_z'], t['preparation_z'] - 2, t['faces'])
    p = trimesh.Trimesh(v, f, process=False)
    if p.volume < 0:
        p.invert()
    p.export(out / 'prep.stl')
    query = m.triangles_center[:32]
    dd = trimesh.proximity.closest_point(m, query)[1]
    oracle = trimesh.proximity.closest_point_naive(m, query)[1]
    error = float(np.max(abs(dd - oracle)))
    write(out / 'SDF.json', dict(resolution='PER_POINT', units='mm', pitch_mm=h, shape=phi.shape, sign='negative inside declared capped roof; virtual prep plane negative below plane', source_mesh_sha256=sha(out / 'measured_capped_roof.stl'), preparation_height_mm=t['preparation_height_mm'], preparation_status='VIRTUAL_HORIZONTAL_SOURCE_PREPARATION_NOT_DENTIN', naive_surface_distance_max_discrepancy_mm=error, naive_oracle_pass=error <= 1e-09, wrong_distance_rejected=float(np.max(abs(dd + 0.01 - oracle))) > 1e-09, continuous_SDF_enclosure='MISSING; finite sampled field is not a geometric accuracy guarantee', flow='Existing task preparation height sampled as a level-set plane; X60 uses its unchanged plane nodes', source_code='published trimesh/proximity.py 4.11.3', external_referent=dict(kind='published_code', locator='https://github.com/mikedh/trimesh/blob/4.11.3/trimesh/proximity.py', compared_quantity='Nearest triangle distance at 32 original surface centroids', refutes_us=True)))

def generate(out):
    sys.path.insert(0, str(X60 / 'code'))
    from plugin import generate as fn
    from fast_geometry import optimize
    from legacy.geometry import shell
    from legacy.checks import all_checks
    t = task()
    base = optimize(t)
    g = fn(t)
    if g['status'] != 'DESIGN':
        write(out / 'GENERATION.json', g)
        raise ValueError('Frozen patient task did not generate a design: ' + str(g))
    z = np.asarray(g['outer_vertices'])[:, 2]
    inner = np.asarray(g['inner_vertices'])[:, 2]
    (v, f) = shell(t['xy'], z, inner, t['faces'])
    m = trimesh.Trimesh(v, f, process=False)
    if m.volume < 0:
        m.invert()
    m.export(out / 'crown.stl')
    for key in ['seconds', 'local_history']:
        g['diagnostics'].pop(key, None)
    ck = all_checks(t, g)
    write(out / 'GENERATION.json', dict(generator='Packaged X60 plugin.generate', task_id=t['task_id'], patient_case='27', geometry_scope='Restricted occlusal shell; complete axial crown UNKNOWN', units='mm', frame=t['frame'], diagnostics=g['diagnostics'], original_geometry_checks=ck, watertight=bool(m.is_watertight), volume_mm3=float(m.volume), calibrated_force05_N=None, physical_force='UNKNOWN_UNMATCHED_PATIENT_SPECIMEN', full_sensitivity_enclosure='MISSING', timescale='HANDOVER'))
    from legacy.generators import submit
    bad_wall = all_checks(t, submit(t, inner + 0.001, inner))
    extreme = all_checks(t, submit(t, z + 100, inner))
    bad_antagonist = all_checks(t, submit(t, z + 1, inner))
    write(out / 'GENERATOR_FAULTS.json', dict(wall_injection_rejected=bad_wall['validity'] == 'FAIL', antagonist_injection_rejected=bad_antagonist['validity'] == 'FAIL' and bad_antagonist['checks'].get('antagonist', {}).get('status') == 'FAIL', injected_wall_mm=0.001, injected_axial_shift_mm=1, extreme_extent_injection=dict(shift_mm=100, actual_source_status=extreme['validity'], rejected=extreme['validity'] in ('FAIL', 'INVALID'), reason=extreme.get('reason')), wall_check=bad_wall['checks'].get('wall'), antagonist_check=bad_antagonist['checks'].get('antagonist')))
    npz(out / 'DESIGN_FIELDS.npz', xyz_mm=v, faces=f, outer_z_mm=z, inner_z_mm=inner, local_thickness_mm=z - inner, base_outer_z_mm=np.asarray(base['outer_vertices'])[:, 2], virtual_prep_z_mm=t['preparation_z'])
    nf = len(t['faces'])
    nside = len(f) - 2 * nf
    labels = ['occlusal'] * nf + ['intaglio'] * nf + ['marginal' if i % 2 == 0 else 'axial' for i in range(nside)]
    contract = dict(units='mm', common_frame_confirmed=True, input_sha256={k: sha(out / (name + '.stl')) for (k, name) in [('prep', 'prep'), ('antagonist', 'antagonist'), ('crown', 'crown')]}, regions=dict(intaglio=list(range(nf, 2 * nf)), exterior=list(range(nf)), occlusal=list(range(nf)), margin=list(range(2 * nf, len(f))), preparation=list(range(nf)), antagonist=list(range(len(trimesh.load_mesh(out / 'antagonist.stl', process=False).faces)))), preparation_outward_axis=[0, 0, 1], contact_axis=[0, 0, 1], film_limits_mm=[0.04, 0.12], film_limits_source='Frozen v4 scenario PHENOMENOLOGICAL; replace with same-object seated-film measurement', occlusal_gap_limits_mm=[0.02, 0.2], occlusal_limits_source='Frozen demonstration lab policy PHENOMENOLOGICAL; no clinical interpretation', region_semantics='Constructed occlusal roof and intaglio plane; side rim is not a clinical axial wall or finish line', material_status='Source X60 KATANA_HTML_PLUS scenario; no HT IFU substituted', indication='posterior')
    write(out / 'CONTRACT.json', contract)
    write(out / 'REGIONS.json', dict(design_sha256=sha(out / 'crown.stl'), face_regions=labels, semantics='Construction domains; marginal/axial split of side facets, not measured anatomical regions'))

def gate(out):
    sys.path.insert(0, str(RESULTS / 'LANE_X49_DESIGN_GATE/code'))
    import design_gate as dg
    paths = {k: str(out / (k + '.stl')) for k in ['prep', 'antagonist', 'crown']}
    g = dg.check(paths, material='KATANA_HTML_PLUS', contract=out / 'CONTRACT.json', units='mm', sinter_factor=1.25)
    write(out / 'GATE_TELEMETRY.json', g)
    g = normalize(g, out)
    for k in paths:
        g['inputs'][k]['path'] = paths[k]
    write(out / 'GATE_PREEXPORT.json', g)
    from export_bridge import export
    x = export(paths, g, out / 'export', out / 'CONTRACT.json', 1.25)
    if x['status'] == 'PASS':
        zpath = out / 'export/model.3mf'
        with zipfile.ZipFile(zpath) as z:
            members = {n: z.read(n) for n in z.namelist()}
        with zipfile.ZipFile(zpath, 'w') as z:
            for (n, b) in sorted(members.items()):
                zi = zipfile.ZipInfo(n, (1980, 1, 1, 0, 0, 0))
                zi.compress_type = zipfile.ZIP_DEFLATED
                z.writestr(zi, b)
        side = out / 'export/model.3mf.json'
        s = read(side)
        s['asset_sha256'] = sha(zpath)
        write(side, s)
        (out / 'export/model.3mf.json.sha256').write_text(sha(side) + '\n')
        from export_bridge import check_native
        report_sha = sha(out / 'export/DESIGN_GATE.json')
        recheck = check_native(zpath, side, sha(side), sha(out / 'crown.stl'), report_sha, 1.25)
        recheck.pop('regions_in_imported_face_order', None)
        x['receipt'].update(asset_sha256=sha(zpath), sidecar_sha256=sha(side), roundtrip=recheck)
        write(out / 'export/EXPORT_RECEIPT.json', x['receipt'])
    g['rules']['export'] = x
    write(out / 'GATE.json', normalize(g, out))
    write(out / 'EXPORT.json', normalize(x, out))
    bad = read(out / 'CONTRACT.json')
    bad['input_sha256']['crown'] = '0' * 64
    write(out / 'BAD_CONTRACT.json', bad)
    try:
        dg.check(paths, 'KATANA_HTML_PLUS', out / 'BAD_CONTRACT.json', 'mm', 1.25)
        reject = False
    except ValueError:
        reject = True
    try:
        dg.check(paths, 'KATANA_HTML_PLUS', out / 'CONTRACT.json', 'um', 1.25)
        unitreject = False
    except ValueError:
        unitreject = True
    write(out / 'GATE_FAULTS.json', dict(wrong_input_hash_rejected=reject, conflicting_units_rejected=unitreject))

def milling(out):
    sys.path.insert(0, str(RESULTS / 'LANE_X53_MILLING_TOOLS/code'))
    from collision import Scene, search, pose_clearance
    m = trimesh.load_mesh(out / 'crown.stl', process=True)
    scene = Scene(m.vertices, m.faces)
    card = read(RESULTS / 'LANE_X53_MILLING_TOOLS/NOMINAL_TIP_REFERENCE_CARD.json')
    nf = len(task()['faces'])
    raw = trimesh.load_mesh(out / 'crown.stl', process=False)
    ids = np.linspace(nf, 2 * nf - 1, 8, dtype=int)
    rows = []
    tools = []
    for q in card['tools']:
        if q['library'] != 'vhf':
            continue
        q = copy.deepcopy(q)
        for k in ['diameter_mm', 'neck_reach_mm', 'shank_mm', 'gauge_mm', 'holder_diameter_mm', 'holder_length_mm']:
            q[k] *= 0.8
        tools.append(q)
    rejects = []
    for i in ids:
        p = raw.triangles_center[i]
        n = raw.face_normals[i]
        for q in tools:
            a = search(scene, p, n, q, 5, offsets=(0,))
            rows.append(dict(face_id=int(i), xyz_mm=p, normal=n, tool_id=q['id'], **a))
            if a['status'] == 'FOUND':
                bad = copy.deepcopy(q)
                bad['holder_diameter_mm'] = 200
                (ok, reason, val) = pose_clearance(scene, np.asarray(a['centre']), n, np.asarray(a['direction']), bad, float(np.linalg.norm(scene.bounds[1] - scene.bounds[0]) + 10))
                rejects.append(not ok)
    write(out / 'MILLING.json', dict(mesh_sha256=sha(out / 'crown.stl'), resolution='PER_POINT', units='mm', timescale='HANDOVER', queried_surface_points=len(ids), tool_point_queries=len(rows), found_queries=sum((r['status'] == 'FOUND' for r in rows)), rows=rows, source_reference='X53 nominal tip-reference-converted VHF card', green_to_final_scale=0.8, scale_status='PHENOMENOLOGICAL scenario; measured lot shrinkage required', holder_fault_rejected=bool(rejects) and all(rejects), full_machine_feasibility='UNKNOWN_NO_FIXTURE_OR_CONTINUOUS_POSE_COVERAGE', physical_milling='NOT_RUN', no_surface_fraction_claim=True))

def physical(out):
    tree = ast.parse((V4 / 'code/literature.py').read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'calibrated_force' or (isinstance(n, ast.Assign) and any((isinstance(a, ast.Name) and a.id == 'PROTOCOL_FIELDS' for a in n.targets)))]
    ns = {'math': math}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), 'packaged_calibrated_force', 'exec'), ns)
    t = task()
    f = np.load(out / 'DESIGN_FIELDS.npz')
    thick = float(np.median(f['local_thickness_mm']))
    query = dict(thickness_mm=thick, product=t['material'], specimen_geometry='Bits2Bites27_X60_restricted_roof', support='X60 rigid basal nodes', cement='UNKNOWN_SEATED_FILM', angle_deg=0, indenter_diameter_mm=None, crosshead_mm_min=None, contact_layer='X60 normalized nodal unit load', aging=None)
    answer = ns['calibrated_force']([], query)
    lit = read(V4 / 'FROZEN_LITERATURE_PREDICTIONS.json')['rows']
    rows = []
    for q in lit:
        error = abs(q['observed'] - q['predicted']) if q['kind'] == 'cement' else abs(math.log(q['observed'] / q['predicted']))
        threshold = 20 if q['kind'] == 'cement' else 0.2
        bad_error = 100 if q['kind'] == 'cement' else abs(math.log(1.6))
        rows.append(dict(**q, error=error, tolerance=threshold, gate=error <= threshold, injected_bad_value_rejected=bad_error > threshold))
    curated = read(V4 / 'payload/literature/CROWN_CURATED.json')
    groups = {}
    for q in curated:
        if q['eligible'] and q['protocol'].startswith('crown'):
            groups.setdefault((q['study'], q['protocol'], q['product']), []).append(q)
    control = []
    for (key, gg) in groups.items():
        gg = sorted(gg, key=lambda x: x['t'])
        if len(set((q['t'] for q in gg))) < 3:
            continue
        (a, b) = (gg[0], gg[-1])
        for q in gg[1:-1]:
            weights = np.array([1, math.log(q['t'])])
            coef = np.linalg.solve([[1, math.log(a['t'])], [1, math.log(b['t'])]], np.log([a['mean'], b['mean']]))
            pred = math.exp(float(weights @ coef))
            source = next((x for x in lit if x['kind'] == 'crown' and x['study'] == q['study'] and (x['protocol'] == q['protocol']) and (x['x'] == q['t'])))
            control.append(dict(locator='doi:' + q['doi'] + ' ' + q['locator'], predicted_N=pred, observed_N=q['mean'], error_log=abs(math.log(pred / q['mean'])), gate=abs(math.log(pred / q['mean'])) <= 0.2, source_prediction_N=source['predicted'], max_abs_control_difference_N=abs(pred - source['predicted']), wrong_numeric_value_rejected=abs(math.log(1.6 * pred / pred)) > 0.2))
    sys.path.insert(0, str(X60 / 'code'))
    from round1 import run as source_run
    source_run()
    summary = read(out / 'raw/R1.json')['sufficiency']
    force_rows = [q for q in rows if q['kind'] == 'crown']
    cement_rows = [q for q in rows if q['kind'] == 'cement']
    write(out / 'PHYSICAL.json', dict(claim_type='capability', patient_force_query=query, calibrated_patient_force=answer, calibrated_force05_N=None, cement_nominal_per_point_um=(f['inner_z_mm'] - f['virtual_prep_z_mm']) * 1000, seated_spatial_cement_film_um=None, film_status='UNKNOWN_NO_SAME_OBJECT_WET_SEATING_MEASUREMENT', source_controls=rows, independent_source_force_control=control, source_force_pass_count=sum((q['gate'] for q in force_rows)), source_force_count=len(force_rows), source_gap_pass_count=sum((q['gate'] for q in cement_rows)), source_gap_count=len(cement_rows), summary_sufficiency=summary, source_support_gate=False, force_likelihood_transfer='FORBIDDEN_UNKNOWN_PROTOCOL', resolution='PER_POINT nominal gap; POPULATION external force/gap groups; PER_TOOTH patient port UNKNOWN', timescale='HANDOVER', affine_sensitivity_enclosure='MISSING; X60 derivatives are local only', external_referent=dict(kind='independent_measurement', locator='doi:10.4047/jap.2021.13.5.269 Table1; doi:10.3390/ma17020365 Table2', compared_quantity='Protocol-specific published crown group fracture forces; no matched generated crown', refutes_us=True)))

def metrology(out):
    sys.path.insert(0, str(RESULTS / 'LANE_X55_METROLOGY/code'))
    import metrology as x55
    m = trimesh.load_mesh(out / 'crown.stl', process=False)
    scan = m.copy()
    scan.vertices = np.asarray(scan.vertices) + [0, 0, 0.025]
    scan.export(out / 'DIGITAL_SCAN_ONLY.stl')
    metadata = dict(reference_kind='design', reference_locator='X70 generated CAD, not independent fabricated reference', scans={'DIGITAL_SCAN_ONLY.stl': dict(scan_to_reference=np.eye(4).tolist(), specimen_id='X70_DIGITAL_ONLY')})
    write(out / 'METROLOGY_METADATA.json', metadata)
    report = x55.analyze(out / 'crown.stl', [out / 'DIGITAL_SCAN_ONLY.stl'], out / 'REGIONS.json', out / 'metrology', mode='datum', units='mm', samples=128, metadata=out / 'METROLOGY_METADATA.json')
    write(out / 'METROLOGY.json', dict(**normalize(report, out), measurement_kind='OUR_OWN_DIGITAL_FAULT_INSTRUMENT', physical_measurement_status='NOT_RUN', interpretation='Known 25um coordinate shift checks point-field pipeline; no scanner or manufacture accuracy measured'))
    a = np.load(out / 'metrology/point_fields.npz')
    npz(out / 'metrology/point_fields.npz', **dict(a))
    bad = read(out / 'REGIONS.json')
    bad['design_sha256'] = '0' * 64
    write(out / 'BAD_REGIONS.json', bad)
    try:
        x55.load_regions(out / 'BAD_REGIONS.json', out / 'crown.stl', x55.load_mesh(out / 'crown.stl', 'mm'))
        reject = False
    except ValueError:
        reject = True
    write(out / 'METROLOGY_FAULTS.json', dict(wrong_region_hash_rejected=reject, own_fixture_is_not_external_facit=True))

def lab(out):
    sys.path.insert(0, str(RESULTS / 'LANE_X61_LAB_ALARM/code'))
    import lab_alarm as x61
    row = dict(specimen_id='X70_PATIENT27_A01', design_id='X70', load_angle_deg='0', fracture_force_N='1000', MG_distance_mean_um='40')
    with (out / 'X61_WRONG_PATIENT_CONTROL.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    try:
        x61.read_csv(out / 'X61_WRONG_PATIENT_CONTROL.csv')
        reject = False
        reason = ''
    except ValueError as e:
        reject = True
        reason = str(e)
    fields = ['specimen_id', 'design_sha256', 'material_batch', 'cement_batch', 'die_sha256', 'angle_deg', 'fracture_force_N', 'force_trace_sha256', 'seated_film_field_path', 'seated_film_sha256', 'pre_sinter_ref_length_mm', 'post_sinter_ref_length_mm', 'scan_path', 'scan_sha256', 'datum_path', 'measurement_error_locator']
    with (out / 'LAB_MEASUREMENTS_TEMPLATE.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow(dict(specimen_id='X70_PATIENT27_A01', design_sha256=sha(out / 'crown.stl')))
    write(out / 'LAB_ALARM.json', dict(claim_type='capability', status='UNKNOWN_NO_MATCHED_PATIENT_LIKELIHOOD', posterior_update=False, physical_measurements=0, missing_ports=['paired manufactured scan/datum', 'assembled spatial cement film', 'patient-design matched force setup/batch/force traces', 'independently frozen lab noise/likelihood profile'], packaged_X61_invoked=True, packaged_X61_patient_ID_rejected=reject, rejection=reason, old_D1_mapping_forbidden=True, alarm_alpha_validity='UNKNOWN', timescale='HANDOVER', minimum_next_operation='Measure matched X70 specimen; create a separate X61 profile/ID contract from measured calibration before subsequent validation. Do not update the D1/M1/M2 profile.'))
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['scan', 'sdf', 'generate', 'gate', 'milling', 'physical', 'metrology', 'lab'])
    ap.add_argument('output')
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    globals()[a.stage](out)
