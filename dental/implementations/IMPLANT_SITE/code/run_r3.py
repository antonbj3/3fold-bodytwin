from common_local import *
from control_geometry import control_box
import time, resource
start = time.perf_counter()
check_frozen()
pr = load(P / 'PREREG_R3.json')
geo = module('x8_geometry3', X8 / 'full_geometry.py')
top = module('x8_printable', X8 / 'printable_geometry.py')
register_inputs([X8 / 'printable_geometry.py'])
choices = load(P / 'raw/RESEARCH_CHOICES_R2.json')
choices = [s for s in choices if s['guide'] == 'fully_guided' and s['target'] == 0.95]
ss = {(s['case'], s['fdi']): s for s in load(X73 / 'raw/LINEAGE_SITES.json')}
out = []
roots = []
ctrl = []
raw = []
for s in choices:
    k = (s['case'], s['fdi'])
    src = ss[k]
    pose = src['pose']
    choice = s['fixed_2mm_research_choice']
    assert choice is not None
    L = choice['L_mm']
    r = choice['D_mm'] / 2
    e = np.array(pose['entry_zyx_mm'])
    a = np.array(pose['axis_zyx'])
    with np.load(src['point_artifact']['path']) as ar:
        points = {key: ar[key].astype(float) * 0.3 for key in ar.files if key.endswith('_voxels_zyx')}
    cache = {}

    def measure(y, control=False):
        if y in cache:
            return cache[y]
        v = {}
        for (name, pt) in points.items():
            rr = geo.voxel_cylinder_bracket(pt, np.full(3, 0.3), e, a, L + y, r)
            v[name] = {q: rr[q] for q in ['lower_mm', 'upper_mm', 'gap_mm', 'witness']}
            if control:
                rel = pt - e
                t = rel @ a
                rad = np.sqrt(np.maximum(0, (rel * rel).sum(1) - t * t))
                dc = np.hypot(np.maximum(0, np.maximum(-t, t - (L + y))), np.maximum(0, rad - r))
                ix = np.where(dc - 0.3 * np.sqrt(3) / 2 <= rr['upper_mm'] + 1e-06)[0]
                ref = min((control_box(pt[j], e, a, L + y, r, np.full(3, 0.3))[0] for j in ix))
                ctrl.append(dict(case=k[0], fdi=k[1], advance_mm=y, mask=name, error_mm=max(0, rr['lower_mm'] - ref, ref - rr['upper_mm']), injected_plus1mm_rejected=abs(rr['upper_mm'] + 1 - ref) > 1e-06))
        lo = min((q['lower_mm'] for q in v.values()))
        hi = min((q['upper_mm'] for q in v.values()))
        cache[y] = (lo, hi)
        raw.append(dict(case=k[0], fdi=k[1], L_mm=L, D_mm=2 * r, advance_mm=y, distances=v, union_lower_mm=lo, union_upper_mm=hi))
        return (lo, hi)
    (lo0, hi0) = measure(0.0, True)
    for y in pr['selection']['advance_mm']:
        (lo, hi) = measure(y, y == 2.0)
        assert lo >= max(0, lo0 - y) - 1e-07 and hi <= hi0 + 1e-07
        out.append(dict(case=k[0], fdi=k[1], L_mm=L, D_mm=2 * r, advance_mm=y, gap_lower_mm=lo, gap_upper_mm=hi, isotropic_lower_mm=max(0, lo0 - y), geometry_gain_over_isotropic_mm=lo - max(0, lo0 - y), digital_2mm_pass=lo >= 2.0, physical_injury_probability='UNKNOWN', resolution='PER_TOOTH', tool='coaxial full-radius cylindrical sweep; bit taper and actual path not measured'))
    if hi0 < 2.0:
        interval = [None, None]
        status = 'BASELINE_ALREADY_BELOW'
    elif lo0 < 2.0:
        interval = [None, None]
        status = 'BASELINE_UNRESOLVED'
    elif measure(2.0)[0] >= 2.0:
        interval = [2.0, None]
        status = 'RIGHT_CENSORED_AT_2MM_ADVANCE'
    else:
        low = 0.0
        high = 2.0
        while high - low > 0.01:
            y = (low + high) / 2
            (lo, hi) = measure(y)
            if lo >= 2:
                low = y
            elif hi < 2:
                high = y
            else:
                break
        interval = [low, high]
        status = 'ROOT_BRACKETED' if high - low <= 0.01 else 'NUMERIC_UNRESOLVED'
    roots.append(dict(case=k[0], fdi=k[1], L_mm=L, D_mm=2 * r, allowable_advance_interval_mm=interval, status=status, criterion='>=2mm to observed mask union, at fixed pose', resolution='PER_TOOTH', rigorous_floating_point_enclosure='MISSING', physical_tolerance='UNKNOWN_UNTIL_ANATOMY_AND_TOOL_ENVELOPE_VALIDATED'))
s = choices[0]
src = ss[s['case'], s['fdi']]
c = s['fixed_2mm_research_choice']
pose = src['pose']
tri = geo.cylinder_triangles(np.array(pose['entry_zyx_mm'])[::-1], np.array(pose['axis_zyx'])[::-1], c['L_mm'], c['D_mm'] / 2, 96).astype(np.float32)
audit = top.topology(tri)
assert all((audit[k] == 0 for k in ['boundary_edges', 'edges_incidence_not2', 'orientation_inconsistent_edges', 'bad_vertex_links', 'duplicate_triangles', 'zero_area_triangles']))
(P / 'export').mkdir(exist_ok=True)
geo.write_stl(P / 'export/research_cylinder.stl', tri)
with np.load(src['point_artifact']['path']) as z:
    allpt = np.unique(np.vstack([z[k] for k in z.files if k.endswith('_voxels_zyx')]), axis=0) * 0.3
np.savetxt(P / 'export/canal_voxel_centers_xyz_mm.csv', allpt[:, ::-1], delimiter=',', header='x_mm,y_mm,z_mm', comments='', fmt='%.4f')
export = dict(case=s['case'], fdi=s['fdi'], coordinate_frame='TF2 index frame mm, xyz exported; source points centered at index*0.3, voxel halfwidth 0.15mm', length_mm=c['L_mm'], diameter_mm=c['D_mm'], pose=pose, topology=audit, facets=96, inscribed_radial_sag_mm=c['D_mm'] / 2 * (1 - np.cos(np.pi / 96)), floating_point_and_manufacturing_error='Not bounded by polygon sag', claim='Research geometry only; not a clinical implant or validated manufactured part', source=src['point_artifact'])
dump(P / 'export/GEOMETRY.json', export)
csvout(P / 'raw/DRILL_ADVANCE_R3.csv', out)
dump(P / 'raw/DRILL_LIMITS_R3.json', roots)
dump(P / 'raw/DRILL_CONTROLS_R3.json', ctrl)
dump(DATA / 'DRILL_GEOMETRY_R3.json', raw)
res = dict(claim_type=['information_link', 'capability'], status='DIRECTIONAL_DIGITAL_DRILL_LIMITS_COMPUTED_PHYSICAL_UNKNOWN', sites=len(roots), queries=len(raw), root_counts={st: sum((r['status'] == st for r in roots)) for st in sorted({r['status'] for r in roots})}, max_gain_over_isotropic_mm=max((r['geometry_gain_over_isotropic_mm'] for r in out)), digital_pass_at_1p5mm=sum((r['digital_2mm_pass'] for r in out if r['advance_mm'] == 1.5)), digital_fail_at_1p5mm=sum((not r['digital_2mm_pass'] for r in out if r['advance_mm'] == 1.5)), max_control_error_mm=max((r['error_mm'] for r in ctrl)), gates={'geometry_controls': all((r['error_mm'] <= 1e-06 for r in ctrl)), 'decisive_directional_gain': max((r['geometry_gain_over_isotropic_mm'] for r in out)) > 0.1}, external_referent=pr['external_referent'], raw_artifact=dict(path=str(DATA / 'DRILL_GEOMETRY_R3.json'), sha256=sha(DATA / 'DRILL_GEOMETRY_R3.json'), bytes=(DATA / 'DRILL_GEOMETRY_R3.json').stat().st_size), cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
dump(P / 'rounds/R3.json', res)
f = P / 'FROZEN_PREDICTIONS.json'
payload = {'schema': 'prospective-lab-metrology-targets-v1', 'source': 'R3 retrospective digital computation; future fabrication/metrology not performed', 'targets': out, 'geometry_sha256': sha(P / 'export/research_cylinder.stl'), 'raw_source_sha256': sha(DATA / 'DRILL_GEOMETRY_R3.json'), 'claim': 'Expected digital reference distances, no human injury/temperature prediction', 'absolute_acceptance_tolerance_mm': 0.3, 'tolerance_basis': 'Prescribed laboratory pilot criterion equal to one source voxel; not clinical accuracy'}
if f.exists():
    previous = load(f)
    assert all((previous[k] == v for (k, v) in payload.items())), 'Frozen future targets changed'
else:
    payload['frozen_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    dump(f, payload)
    f.with_suffix('.json.sha256').write_text(sha(f) + '\n')
state('R3_COMPLETE', res['gates'], 'Compile external clinical referent, source-control faults, demo and pending feedback')
(P / 'HANDOFF_R3.md').write_text("# R3 : Digital boundary of the drilling track\n\n" + json.dumps(res, indent=2) + "\n\nFrozen metrology targets exist. Minimum new information: achieved shaft/deep and independent canal wall in the same recorded frame; local drilling torque, temperature history and separate material calibration remain.\n")
print(json.dumps(res, indent=2))
