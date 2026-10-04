from dental_release.paths import expand as _release_expand
import sys, time, copy
import numpy as np
from scipy import ndimage as ndi
from shared import *
run = Path(sys.argv[1])
out = run / 'R3'
out.mkdir(exist_ok=True)
di = DATA / run.name / _release_expand('@DENTAL_CASE_C@')
d = DATA / run.name / 'R3'
d.mkdir(parents=True, exist_ok=True)
start = time.perf_counter()
x18 = RESULTS / _release_expand('X18')
sys.path.insert(0, str(x18 / 'code'))
import geometry as g
from experiment import contact
from fe import solve_cases, normalize
ex = module('mesh_transport', RESULTS / 'LANE_X38_EXPORT_GATE/code/mesh_transport.py')
(U, _) = ex.load_mesh(di / 'antagonist.stl', 'mm')
z = np.load(di / 'roof.npz')
map0 = np.load(di / 'contact_map.npz')
rows = []
for delta in [-0.05, 0.0, 0.05]:
    tid = 'pose_' + str(delta)
    whole = map0['gap'] + delta
    c = contact(z['ceiling'] + delta - z['z'], z['z'], z['xy'], z['index'], z['weights'], 0.1, 0.1)
    loads = {tid: normalize(c['mask'], z['weights'], 100.0)}
    (fe, cost) = solve_cases(z['xy'], z['z'], z['faces'], z['uv'], loads, dict(thickness_mm=1.2, E_MPa=210000, nu=0.3))
    base = g.signed_gap(U, z['xy'], z['z'], z['faces'])
    moved = g.signed_gap(U + [0, 0, delta], z['xy'], z['z'], z['faces'])
    diff = abs(moved['minimum_gap_mm'] - (base['minimum_gap_mm'] + delta))
    rows.append(dict(patient_id=_release_expand('@DENTAL_CASE_C@'), theta_row_id=tid, shared_theta=dict(IOS_pose_delta_mm=delta, force_N=100.0, material_E_MPa=210000.0, material_nu=0.3, axial_wall_mm=1.2), whole_bite=dict(absolute_near_area_mm2=float(np.sum(np.abs(whole) <= 0.1) * float(map0['grid_mm']) ** 2), resolution='PER_ARCH'), crown=dict(contact_area_mm2=c['area_mm2'], patch_count=c['patch_count'], minimum_gap_mm=moved['minimum_gap_mm'], resolution='PER_SURFACE_REGION'), FE=fe[tid], direct_triangle_control_error_mm=diff, scope='One shared deterministic scenario; no joint probability or physical robustness claim'))

def bind(producer, consumer):
    for k in ['patient_id', 'frame_sha256', 'unit', 'source_sha256']:
        if producer[k] != consumer[k]:
            raise ValueError('PORT_MISMATCH:' + k)
    return True
framefile = run / _release_expand('@DENTAL_CASE_C@/FRAME.json')
base = dict(patient_id=_release_expand('@DENTAL_CASE_C@'), frame_sha256=sha(framefile), unit='mm', source_sha256=load(run / _release_expand('@DENTAL_CASE_C@/IOS_RESULTS.json'))['source']['upper']['sha256'])
checks = {'valid_port_pass': bind(base, base)}
for (k, value) in [('patient_id', _release_expand('@DENTAL_CASE_A@')), ('frame_sha256', 'wrong'), ('unit', 'um'), ('source_sha256', '0' * 64)]:
    bad = {**base, k: value}
    try:
        bind(base, bad)
        checks[k + '_fault_rejected'] = False
    except ValueError:
        checks[k + '_fault_rejected'] = True
fields = []
for case in [_release_expand('@DENTAL_CASE_A@'), _release_expand('@DENTAL_CASE_B@')]:
    p = DATA / run.name / case / '36_paired.npz'
    q = np.load(p)
    phi = {}
    for name in ['tooth', 'pulp']:
        mask = q[name]
        phi[name + '_grid_center_sdf_mm'] = (ndi.distance_transform_edt(~mask, sampling=q['spacing']) - ndi.distance_transform_edt(mask, sampling=q['spacing'])).astype(np.float32)
    owner = np.where(q['pulp'], 2, np.where(q['tooth'], 1, 0)).astype(np.uint8)
    file = d / (case + '_region_field.npz')
    np.savez_compressed(file, **phi, region_owner=owner, spacing=q['spacing'], crop_origin_zyx=q['crop_origin_zyx'])
    fields.append(dict(patient_id=case, FDI=36, crop_axes='zyx; z reversed when crown_axis_sign<0', crown_axis_sign=load(run / case / 'FRAME.json')['crops']['36']['crown_axis_sign'], file=artifact(file), frame=artifact(run / case / 'FRAME.json'), source_crop=artifact(p), owner_ids={'0': 'outside', '1': 'total_hard_tissue', '2': 'pulp'}, resolution='PER_POINT', sign='negative inside; discrete center EDT, not exact continuous SDF', time_scale='SIMULTANEOUS', consumer='X33 finite-volume tooth/pulp masks; E/heat properties are separate shared closures', UNKNOWN=['DEJ', 'PDL', 'cement', 'IOS_registration']))
freeze(out / 'FROZEN_PREDICTIONS.json', dict(round='R3', prereg_R3_sha256=sha(R / 'PREREG_R3.json'), patient_id=_release_expand('@DENTAL_CASE_C@'), theta_rows=rows, frame=load(framefile), clinical_release='ABSTAIN', physical_measurement='NOT_RUN'))
exporter = module('exportgate', RESULTS / 'LANE_X38_EXPORT_GATE/code/exportgate.py')
export = exporter.export(DATA / run.name / 'R2/crown_oriented.stl', d / 'export', 'mm', di / 'regions.json', '3Y proxy; patient material UNKNOWN', 1.0, 'final_sintered', out / 'FROZEN_PREDICTIONS.json', run / 'R2/DESIGN_GATE.json')
exported = load(d / 'export/FROZEN_PREDICTIONS.json')
export_row_identity = exported['theta_rows'] == clean(rows) and exported['patient_id'] == _release_expand('@DENTAL_CASE_C@') and (exported['frame'] == load(framefile))
checks['export_frozen_theta_and_frame_identity'] = export_row_identity
for field in fields:
    case = field['patient_id']
    freeze(out / (case + '_FROZEN_FIELD.json'), dict(patient_id=case, region_field=field, clinical_release='ABSTAIN'))
dump(out / 'RESULTS_R3.json', dict(final_patient_export=export, claim_type='capability', theta_rows=rows, ports=checks, region_fields=fields, maximum_direct_control_error_mm=max((r['direct_triangle_control_error_mm'] for r in rows)), runtime_s=time.perf_counter() - start, external_referent=load(R / 'PREREG_R3.json')['external_referent']))
print('R3 done', [(r['theta_row_id'], r['crown']['contact_area_mm2'], r['FE']['status']) for r in rows], flush=True)
