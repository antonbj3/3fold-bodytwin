from dental_release.paths import expand as _release_expand
import sys, zipfile, time, resource, subprocess, os
import numpy as np
from shared import *
run = Path(sys.argv[1])
out = run / _release_expand('@DENTAL_CASE_C@')
out.mkdir(exist_ok=True)
d = DATA / run.name / _release_expand('@DENTAL_CASE_C@')
d.mkdir(parents=True, exist_ok=True)
x7 = RESULTS / _release_expand('X7')
x11 = RESULTS / 'LANE_X11_TOOTH_SEG'
x21 = RESULTS / _release_expand('X21')
zpath = load(x7 / 'raw/DATA_MANIFEST.json')['zip']
case = _release_expand('@DENTAL_SURFACE_ID_C@')
start = time.perf_counter()
src = {}
with zipfile.ZipFile(zpath) as z:
    names = z.namelist()
    for jaw in ['upper', 'lower']:
        member = case + '/ios/ios_' + jaw + '.stl'
        p = d / (jaw + '.stl')
        p.write_bytes(z.read(member))
        src[jaw] = dict(archive=zpath, member=member, **artifact(p))
    modalities = {'ios': True, 'report_members': [n for n in names if n.startswith(case + '/') and '/reports_ios_en/' in n and n.endswith('.txt')], 'image_members': [n for n in names if n.startswith(case + '/') and n.lower().endswith(('.jpg', '.png', '.jpeg'))], 'cbct_members': [n for n in names if n.startswith(case + '/') and n.lower().endswith(('.dcm', '.mha', '.nii.gz'))]}
cmd = [sys.executable, str(x11 / 'segment_ios.py'), str(d / 'upper.stl'), str(d / 'lower.stl'), '--out', str(d / 'labels+landmarks.json')]
with (out / 'X11.log').open('w') as log:
    p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if p.returncode:
        raise RuntimeError('X11 failed; see ' + str(out / 'X11.log'))
meta = load(d / 'labels+landmarks.json')
ct = module('contact', x21 / 'code/contact.py')
ct.D = d
fr = meta['arches']['lower']['frame']
Rmat = np.column_stack([fr['right_unit'], fr['anterior_unit'], fr['superior_unit']])
center = np.array(fr['center_mm'])
assert np.allclose(Rmat.T @ Rmat, np.eye(3), atol=1e-10)
arches = [ct.prepare((d / (j + '.stl')).read_bytes(), meta['arches'][j], j, Rmat, center) for j in ['upper', 'lower']]
(U, ul, uc, uf, us) = arches[0]
(L, ll, lc, lf, ls) = arches[1]
ub = ct.boxes(U)
lb = ct.boxes(L)
origin = np.minimum(ub[:, :2].min(0), lb[:, :2].min(0)) - 0.8
hi = np.maximum(ub[:, 3:5].max(0), lb[:, 3:5].max(0))
shape = np.ceil((hi - origin) / 0.8).astype(np.int64) + 2
(ptr, ids) = ct.buckets(lb, origin, 0.8, shape)
(best, wi, xy, n, ne) = ct.extrema(U, L, ul, ll, uc, lc, ub, lb, origin, 0.8, shape, ptr, ids, True)
m = ct.make_map(U, L, ul, ll, uc, lc)
(pairs, features) = ct.summaries(m, best, wi, xy, U, L, uc, lc, uf, lf)
np.savez_compressed(d / 'contact_map.npz', **m, minimum_gap=best)
cmd = [sys.executable, str(x7 / 'code/diagnose.py'), '--upper', str(d / 'upper.stl'), '--lower', str(d / 'lower.stl'), '--out', str(out / 'X7_PREDICTION.json')]
with (out / 'X7.log').open('w') as log:
    p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if p.returncode:
        raise RuntimeError('X7 failed')
freeze(out / 'FROZEN_PREDICTIONS.json', dict(patient_id=_release_expand('@DENTAL_CASE_C@'), source=src, files=[artifact(d / 'labels+landmarks.json'), artifact(d / 'contact_map.npz'), artifact(out / 'X7_PREDICTION.json')], report_observed_in_this_run=False, physical_measurement='NOT_RUN', prior_exposure='X7 test and X21 upstream panel already evaluated; descriptive replay'))
parser_script = R / 'code/report_compare.py'
subprocess.run([sys.executable, str(parser_script), str(out), zpath], check=True)
v = np.array([center, center + [1, 2, 3], center + [-7, 5, 11]])
rt = float(np.max(np.abs((v - center) @ Rmat @ Rmat.T + center - v)))
cached = load(x21 / _release_expand('raw/cases/@DENTAL_SURFACE_ID_C@.json'))
dump(out / 'X21_REPLAY.json', dict(pairs=pairs, features=features, stats=dict(upper=us, lower=ls), cached_keys=list(cached), map=artifact(d / 'contact_map.npz'), physical_contact='UNKNOWN'))
dump(out / 'FRAME.json', dict(patient_id=_release_expand('@DENTAL_CASE_C@'), native_frame='original registered IOS pair, no anatomical calibration', unit='mm', local_to_native_rotation=Rmat.tolist(), native_origin_xyz_mm=center.tolist(), equations=['local=(native-center)@R', 'native=local@R.T+center'], roundtrip_error_mm=rt, determinant=float(np.linalg.det(Rmat)), X7_frame=load(out / 'X7_PREDICTION.json')['frame'], X7_scope='scalar categories in own published axis-selected frame; coordinates not joined without explicit transform', cross_modal_transform=None))
dump(out / 'IOS_RESULTS.json', dict(patient_id=_release_expand('@DENTAL_CASE_C@'), claim_type='capability', source=src, modalities=modalities, segmentation=artifact(d / 'labels+landmarks.json'), contact_pair_count=len(pairs), features=features, coordinate_roundtrip_mm=rt, source_pose='SIMULTANEOUS', FDI_accuracy='UNKNOWN on Bite2Text', runtime_s=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, data_dir=str(d)))
print('IOS done', len(pairs), flush=True)
