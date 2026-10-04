from dental_release.paths import expand as _release_expand
import sys, time, resource
import numpy as np
from shared import *
run = Path(sys.argv[1])
out = run / _release_expand('@DENTAL_CASE_C@')
d = DATA / run.name / _release_expand('@DENTAL_CASE_C@')
start = time.perf_counter()
x18 = RESULTS / _release_expand('X18')
sys.path.insert(0, str(x18 / 'code'))
import geometry as g
from experiment import local_tri, contact, score
from fe import solve_cases, normalize
meta = load(d / 'labels+landmarks.json')
frame = load(out / 'FRAME.json')
rot = np.array(frame['local_to_native_rotation'])
center = np.array(frame['native_origin_xyz_mm'])
arches = {}
for jaw in ['upper', 'lower']:
    b = (d / (jaw + '.stl')).read_bytes()
    n = int.from_bytes(b[80:84], 'little')
    dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    native = np.frombuffer(b, dt, offset=84, count=n)['v'].astype(float)
    (v, inv) = np.unique(native.reshape(-1, 3), axis=0, return_inverse=True)
    f = inv.reshape(-1, 3)
    v = (v - center) @ rot
    labs = np.load(meta['arches'][jaw]['labels']['path'])['labels']
    ll = labs[f]
    owner = np.where(ll[:, 0] == ll[:, 1], ll[:, 0], np.where(ll[:, 0] == ll[:, 2], ll[:, 0], np.where(ll[:, 1] == ll[:, 2], ll[:, 1], 0)))
    arches[jaw] = dict(v=v, f=f, tri=v[f], labels=labs, owner=owner)
pr = load(R / 'PREREG_R1.json')
ss = pr['numerical_scenarios']
rejected = []
selected = None
for fdi in ss['IOS_fdi_priority']:
    try:
        s = g.site(arches['lower'], fdi, ss['grid_n'])
        selected = fdi
        break
    except ValueError as e:
        rejected.append(dict(fdi=fdi, reason=str(e)))
if selected is None:
    dump(out / 'CROWN_RESULTS.json', dict(status='ABSTAIN_NO_SITE', rejected=rejected))
    raise SystemExit(0)
fdi = selected
U = local_tri(arches['upper']['tri'], s['xy'])
(ceiling, _) = g.query_height(U, s['xy'], True)
donor = Path(_release_expand('@DENTAL_WORK_ROOT@/X18_crown_antagonist/donor_template.npz'))
template = np.load(donor)[str(fdi % 10)]
assert len(template) == len(s['xy'])
z0 = s['z_top'] + template
z = np.minimum(z0, ceiling - ss['roof_clearance_mm'])
(z, repair) = g.repair(U, s, z, ss['roof_clearance_mm'])
ctrl = z0 - np.maximum(z0 - ceiling + ss['roof_clearance_mm'], 0)
(ctrl, ctrlrepair) = g.repair(U, s, ctrl, ss['roof_clearance_mm'])
np.savez_compressed(d / 'roof.npz', xy=s['xy'], faces=s['faces'], uv=s['uv'], z=z, z_prior=z0, ceiling=ceiling, weights=s['area_weight'], index=s['index'])
(v, fs) = g.shell(s['xy'], z, s['faces'], ss['roof_axial_thickness_mm'])
g.write_stl(d / 'crown.stl', v, fs)
plane = np.c_[s['xy'], np.full(len(z), float(z.min() - ss['roof_axial_thickness_mm'] - 0.05))]
g.write_stl(d / 'prep.stl', plane, s['faces'])
keep = np.abs(np.cross(U[:, 1, :2] - U[:, 0, :2], U[:, 2, :2] - U[:, 0, :2])) > 1e-12
ant = U[keep]
g.write_stl(d / 'antagonist.stl', ant.reshape(-1, 3), np.arange(3 * len(ant)).reshape(-1, 3))
nf = len(s['faces'])
paths = {k: d / (k + '.stl') for k in ['crown', 'prep', 'antagonist']}
contract = dict(units='mm', common_frame_confirmed=True, contact_axis=[0, 0, 1], input_sha256={k: sha(p) for (k, p) in paths.items()}, regions=dict(intaglio=list(range(nf, 2 * nf)), exterior=list(range(nf)), occlusal=list(range(nf)), preparation=list(range(nf)), antagonist=list(range(len(ant)))), ifu_profile='katana-ht', indication='posterior', film_limits_mm=[0.02, 0.12], film_limits_source='X34 illustrative PHENOMENOLOGICAL protocol', occlusal_gap_limits_mm=[-0.03, 0.1], preparation_outward_axis=[0, 0, 1], physical_geometry_status='Patient registered antagonist, generated roof and VIRTUAL plane; not measured preparation', patient_id=_release_expand('@DENTAL_CASE_C@'), frame=frame)
dump(d / 'design_contract.json', contract)
co = contact(ceiling - z, z, s['xy'], s['index'], s['area_weight'], 0.1, 0.1)
force_loads = {str(F): normalize(co['mask'], s['area_weight'], F) for F in ss['FE_force_N']}
(ans, cost) = solve_cases(s['xy'], z, s['faces'], s['uv'], force_loads, dict(thickness_mm=1.2, E_MPa=210000, nu=0.3))
freeze(out / 'FROZEN_CROWN_PREDICTIONS.json', dict(patient_id=_release_expand('@DENTAL_CASE_C@'), FDI=fdi, frame=frame, files=[artifact(d / 'roof.npz'), artifact(d / 'crown.stl'), artifact(d / 'design_contract.json')], FE=ans, FE_source=str(x18 / 'code/fe.py'), FE_scope='Conditional roof proxy; no pressure, support or strength measurement', donor_prior=artifact(donor), original_surface_scored=False, physical_measurement='NOT_RUN', theta_refs=['IOS_pose', 'force', 'support', 'material', 'preparation']))
(orig, _) = g.query_height(arches['lower']['tri'][arches['lower']['owner'] == fdi], s['xy'], False)
truth = contact(ceiling - orig, orig, s['xy'], s['index'], s['area_weight'], 0.1, 0.1)
prior = contact(ceiling - z0, z0, s['xy'], s['index'], s['area_weight'], 0.1, 0.1)
scored = {a: score(q, truth, s['xy']) for (a, q) in [('generated', co), ('donor_prior', prior)]}
x18b = module('surface', RESULTS / 'LANE_X18B_PREP_SURFACE/code/surface.py')
carrier = x18b.specimen(arches['lower'], fdi, 1.2)
g.write_stl(d / 'measured_exterior.stl', carrier['external_vertices'], carrier['external_faces'])
back = x18b.read_stl(d / 'measured_exterior.stl')
source = arches['lower']['tri'][carrier['source_face_ids']]
err = float(np.max(np.abs(source - back)))
rawg = module('designgate', RESULTS / 'LANE_X34_STL_DESIGN_GATE/code/designgate/__init__.py') if False else None
sys.path.insert(0, str(RESULTS / 'LANE_X34_STL_DESIGN_GATE/code'))
from designgate.gate import check
report = check(paths, '3Y', d / 'design_contract.json', round_version='R2')
dump(out / 'DESIGN_GATE.json', report)
bad = dict(contract)
bad['input_sha256'] = {**contract['input_sha256'], 'crown': '0' * 64}
dump(d / 'bad_design_contract.json', bad)
try:
    check(paths, '3Y', d / 'bad_design_contract.json', round_version='R2')
    bad_rejected = False
except ValueError:
    bad_rejected = True
export = module('exportgate', RESULTS / 'LANE_X38_EXPORT_GATE/code/exportgate.py')
dump(d / 'regions.json', dict(labels=['exterior'] * nf + ['intaglio'] * nf + ['virtual_closure'] * (len(fs) - 2 * nf), semantics='Generated roof roles; not anatomical crown margin/CEJ'))
try:
    ex = export.export(d / 'crown.stl', d / 'export', 'mm', d / 'regions.json', '3Y proxy; batch UNKNOWN', 1.0, 'final_sintered', out / 'FROZEN_CROWN_PREDICTIONS.json', out / 'DESIGN_GATE.json')
except ValueError as e:
    ex = dict(status='REFUSED', reason=str(e))
exportfault = None
if ex['status'] == 'PASS':
    side = d / 'export/model.stl.json'
    asset = d / 'export/model.stl'
    try:
        export.check(asset, side, sha(side), 'um', True)
        exportfault = False
    except ValueError:
        exportfault = True
lin = None
if ans['50']['status'] == ans['100']['status'] == 'SIMULATED':
    lin = abs(ans['100']['maximum_tensile_MPa'] / (2 * ans['50']['maximum_tensile_MPa']) - 1)
x1b = load(RESULTS / 'LANE_X1B_CROWN_LOOP/FROZEN_PREDICTIONS.json')
x1bstatus = dict(status='ABSTAIN_DIFFERENT_GEOMETRY', source=artifact(RESULTS / 'LANE_X1B_CROWN_LOOP/FROZEN_PREDICTIONS.json'), reason=_release_expand('Frozen STS L005 geometry cannot bind to Bite2Text @DENTAL_SURFACE_ID_C@; X18 same-patient roof FE run separately; no absolute fracture/lifetime claim'))
result = dict(patient_id=_release_expand('@DENTAL_CASE_C@'), claim_type='capability', FDI=fdi, rejected_sites=rejected, model='X18 generated conditional roof', measured_preparation=False, continuous_gap=g.signed_gap(U, s['xy'], z, s['faces']), repair=repair, control_max_height_difference_mm=float(np.max(np.abs(z - ctrl))), original_contact_area_mm2=truth['area_mm2'], generated_contact_area_mm2=co['area_mm2'], contact_patch_count=co['patch_count'], scores=scored, original_finite_fraction=float(np.isfinite(orig).mean()), original_exterior=dict(export_error_mm=err, source_faces=len(source), normal_offset_fold_count=carrier['preparation_fold_count'], nonmanifold_edges=carrier['nonmanifold_edge_count'], scope='Exact measured exterior; normal offset is diagnostic only'), FE=ans, FE_cost=cost, FE_linearity_error=lin, X1b=x1bstatus, design_gate_verdict=report['verdict'], design_gate_rules={k: v['status'] for (k, v) in report['rules'].items()}, export=ex, controls=dict(design_hash_fault_rejected=bad_rejected, wrong_export_units_rejected=exportfault, external_displacement_rejected=bool(np.max(np.abs(back + [0, 0, 0.2] - source)) > 1e-05), FE_doubled_force_would_reject_original=bool(ans['100'].get('total_load_N', 0) != ans['50'].get('total_load_N', 0))), external_referent=dict(kind='published_dataset', locator=load(out / 'IOS_RESULTS.json')['source']['lower']['member'], compared_quantity='Original source tooth roof versus generated roof against same measured antagonist; contact geometry, not force', refutes_us=scored['generated']['IoU'] < 0.5), runtime_s=time.perf_counter() - start)
dump(out / 'CROWN_RESULTS.json', result)
print('CROWN', report['verdict'], ex['status'], flush=True)
