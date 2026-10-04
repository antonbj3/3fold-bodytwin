from dental_release.paths import expand as _release_expand
import sys, time, resource
import numpy as np
import trimesh
from shared import *
run = Path(sys.argv[1])
out = run / 'R2'
out.mkdir(exist_ok=True)
d = DATA / run.name / 'R2'
d.mkdir(parents=True, exist_ok=True)
start = time.perf_counter()
ios = run / _release_expand('@DENTAL_CASE_C@')
di = DATA / run.name / _release_expand('@DENTAL_CASE_C@')
checks = []
ex = module('exportgate', RESULTS / 'LANE_X38_EXPORT_GATE/code/exportgate.py')
mt = sys.modules['mesh_transport']
(tri, _) = mt.load_mesh(di / 'crown.stl', 'mm')
mesh = mt.physical_mesh(tri)
trimesh.repair.fix_winding(mesh)
trimesh.repair.fix_inversion(mesh)
newtri = mesh.triangles

def canonical_vertices(t):
    return np.asarray([v[np.lexsort((v[:, 2], v[:, 1], v[:, 0]))] for v in t])
a = canonical_vertices(tri)
b = canonical_vertices(newtri)
err = float(np.max(np.abs(a - b)))
assert err <= 1e-05
bad = newtri.copy()
i = next((i for (i, v) in enumerate(bad) if len(np.unique(v[:, 0])) > 1 and len(np.unique(v[:, 1:3], axis=0)) > 1))
bad[i, :, 0] = np.roll(bad[i, :, 0], 1)
assert np.max(np.abs(a - canonical_vertices(bad))) > 1e-05, 'Facet identity control cannot reject a changed facet'
mt.write_stl(d / 'crown_oriented.stl', newtri)
contract = load(di / 'design_contract.json')
contract['input_sha256']['crown'] = sha(d / 'crown_oriented.stl')
dump(d / 'design_contract.json', contract)
sys.path.insert(0, str(RESULTS / 'LANE_X34_STL_DESIGN_GATE/code'))
from designgate.gate import check
paths = dict(crown=d / 'crown_oriented.stl', prep=di / 'prep.stl', antagonist=di / 'antagonist.stl')
dg = check(paths, '3Y', d / 'design_contract.json', round_version='R2')
dump(out / 'DESIGN_GATE.json', dg)
freeze(out / 'FROZEN_PREDICTIONS.json', dict(patient_id=_release_expand('@DENTAL_CASE_C@'), frame=load(ios / 'FRAME.json'), round='R2', prereg_R2_sha256=sha(R / 'PREREG_R2.json'), inputs=[artifact(di / 'crown.stl'), artifact(d / 'crown_oriented.stl'), artifact(d / 'design_contract.json'), artifact(out / 'DESIGN_GATE.json')], expected='Transport may pass; wall/film fail must remain', original_FE=artifact(ios / 'FROZEN_CROWN_PREDICTIONS.json'), physical_observation='NOT_RUN'))
export = ex.export(d / 'crown_oriented.stl', d / 'export', 'mm', di / 'regions.json', '3Y proxy; measured material batch UNKNOWN', 1.0, 'final_sintered', out / 'FROZEN_PREDICTIONS.json', out / 'DESIGN_GATE.json')
side = d / 'export/model.stl.json'
asset = d / 'export/model.stl'
try:
    ex.check(asset, side, sha(side), 'um', True)
    unitfault = False
except ValueError:
    unitfault = True
poison = newtri.copy()
poison[0, 0, 0] += 0.2
mt.write_stl(d / 'bad_vertex.stl', poison)
try:
    ex.check(d / 'bad_vertex.stl', side, sha(side), 'mm', True)
    vertexfault = False
except ValueError:
    vertexfault = True
r = dict(patient_id=_release_expand('@DENTAL_CASE_C@'), facet_coordinate_change_mm=err, faces=len(newtri), winding_consistent=bool(mesh.is_winding_consistent), volume_mm3=float(mesh.volume), design_gate_verdict=dg['verdict'], design_gate_rules={k: v['status'] for (k, v) in dg['rules'].items()}, export=export, controls=dict(unit_fault_rejected=unitfault, changed_vertex_rejected=vertexfault), scope='Orientation-only repair; no clinical or manufacturing release')
dump(out / 'EXPORT_RESULTS.json', r)
x33 = module('laser_heat', RESULTS / 'LANE_X33_LASER_PULP/laser_heat.py')
base = load(RESULTS / 'LANE_X33_LASER_PULP/PREREG_R1.json')
thermal = []
for case in [_release_expand('@DENTAL_CASE_A@'), _release_expand('@DENTAL_CASE_B@')]:
    o = out / case
    o.mkdir(exist_ok=True)
    (o / 'raw').mkdir(exist_ok=True)
    x33.ROOT = o
    crop = DATA / run.name / case / '36_paired.npz'
    source = load(run / case / 'X12_MEASUREMENTS.json')
    assert source['patient_id'] == case and source['crops']['36']['file']['sha256'] == sha(crop)
    pr = {**base, 'sources': {**base['sources'], 'geometry': artifact(crop)}}
    q = np.load(crop)
    prepared = q['tooth'] & (__import__('scipy').ndimage.distance_transform_edt(q['tooth'], sampling=0.3) >= 0.65)
    exposed = int(np.sum(q['pulp'] & ~prepared))
    try:
        geo = x33.geometry(1, pr)
    except AssertionError:
        rec = dict(patient_id=case, FDI=36, status='ABSTAIN_PREPARATION_EXPOSES_PULP', exposed_pulp_voxels=exposed, crop=artifact(crop), preparation_mm=0.65, simulation='NOT_RUN', physical='UNKNOWN')
        thermal.append(rec)
        dump(o / 'THERMAL_RESULTS.json', rec)
        continue
    (neighbors, g, C, pm, pi, pw, wet, k, h, prep, pulp) = geo
    tests = x33.verify(pr)
    scenarios = []
    for hw in [0.0, 500.0, 2000.0]:
        cool = wet * h * h / (1 / max(hw, 1e-100) + h / (2 * k)) if hw else np.zeros(len(C))
        dt = 0.8 / np.max((g.sum(1) + cool) / C)
        (hist, peak, hard, stored, loss, ein, steps, cem) = x33.simulate(neighbors, g, C, pm, pi, pw, cool, 0.25, 0.02, -17.0, 10.0, dt)
        r = dict(h_W_m2K=hw, retained_heat_fraction=0.02, peak_pulp_delta_C=peak, peak_hard_delta_C=hard, stored_J=stored, removed_J=loss, input_J=ein, relative_energy_balance=abs(stored + loss - ein) / ein, steps=steps, status='CONDITIONAL_SIMULATION', resolution='PER_POINT', closure_resolution='PHENOMENOLOGICAL')
        scenarios.append(r)
        np.savetxt(o / 'raw' / ('history_h' + str(int(hw)) + '.csv'), hist, delimiter=',', header='t_s,max_pulp_rise_C,mean_pulp_rise_C,max_hard_rise_C', comments='')
    rec = dict(patient_id=case, FDI=36, status='SIMULATED_CONDITIONAL', crop=artifact(crop), nodes=len(C), exposed_pulp_voxels=exposed, scenarios=scenarios, numerical_controls=tests, physical_validation='UNKNOWN_NO_SAME_SPECIMEN_TEMPERATURE', preparation_mm=0.65, uncertainty='Unmeasured retained heat and cooling; no calibrated prediction interval')
    freeze(o / 'FROZEN_PREDICTIONS.json', dict(**rec, physical_measurement='NOT_RUN', prereg_R2_sha256=sha(R / 'PREREG_R2.json')))
    thermal.append(rec)
    dump(o / 'THERMAL_RESULTS.json', rec)
dump(out / 'RESULTS_R2.json', dict(claim_type='capability', export=load(out / 'EXPORT_RESULTS.json'), thermal=thermal, runtime_s=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, external_referent=load(R / 'PREREG_R2.json')['external_referent']))
print('R2 done', [(x['patient_id'], x['status']) for x in thermal], flush=True)
