from dental_release.paths import expand as _release_expand
from common import *
import numpy as np
import importlib.util
lock = read(ROOT / 'SOURCE_LOCK.json')
work = Path(lock['X77_work'])
scene = np.load(work / 'scene.npz')
task = np.load(work / 'TASK_ARRAYS.npz')
points = np.c_[task['xy'], task['preparation_z']] + scene['base_mm']
DATA.mkdir(parents=True, exist_ok=True)
np.savez_compressed(DATA / 'X77_preparation_points.npz', points_CBCT_xyz_mm=points, lower_IOS_to_CBCT=scene['lower_IOS_to_CBCT'], local_to_cbct_rotation=scene['local_to_cbct_rotation'], base_mm=scene['base_mm'])
rows = read(X69 / 'raw/zenodo_zip_directory.json')
selected = [r for r in rows if r['name'].startswith('Demo/Demo_1/') and r['bytes']]
base = Path(read(X77 / 'PATIENT_SELECTION.json')['source_base'])
groups = {}
for r in selected:
    k = '/'.join(r['name'].split('/')[:-1])
    v = groups.setdefault(k, dict(members=0, bytes=0, local_members=0))
    v['members'] += 1
    v['bytes'] += r['bytes']
    v['local_members'] += int((base.parent.parent / r['name']).is_file())
decisions = [dict(patient_id=lock['same_patient'], query=q, status='UNKNOWN', nominal_value=None, calibrated_value=None, resolution='PER_POINT', registration_error_mm=None, tissue_error_mm=None, reason=reason) for (q, reason) in [('pulp_distance', 'No same-patient pulp observation or independent tissue/IOS error'), ('bone_gray_profile', 'No local same-patient source image or verified implant site'), ('canal_distance_and_annotation_spread', 'No same-patient canal masks; other patients are inadmissible')]]
spec = importlib.util.spec_from_file_location('x78_existing_guards', X77 / 'code/guards.py')
guards = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guards)
faults = []
for tissue in ('pulp', 'bone', 'canal'):
    for (kind, record) in [('missing', None), ('foreign_patient', dict(patient_id='P1', tissue=tissue, frame='CBCT_native', path=_release_expand('@DENTAL_WORK_ROOT@/X73-canal-wall-full/P1_points.npz')))]:
        try:
            guards.admit_tissue(record, tissue, lock['same_patient'])
            rejected = False
            reason = ''
        except ValueError as e:
            rejected = True
            reason = str(e)
        faults.append(dict(tissue=tissue, injected=kind, rejected=rejected, reason=reason))
try:
    guards.registration_bound(dict(kind='surface_residual_p95', epsilon_mm=0.1))
    rejected = False
except ValueError:
    rejected = True
faults.append(dict(injected='surface_residual_as_registration_truth', rejected=rejected))
assert all((f['rejected'] for f in faults))
r = dict(claim_type='capability', status='SAME_PATIENT_TISSUE_NOT_ESTABLISHED', patient_id=lock['same_patient'], preparation_point_count=len(points), point_file=dict(path=str(DATA / 'X77_preparation_points.npz'), sha256=sha(DATA / 'X77_preparation_points.npz')), source_directory_groups=groups, decisions=decisions, missing_tissues=3, requested_tissues=3, rejection_fraction=1.0, rejection_reason='Missing local matched semantic observation', external_referent=read(ROOT / 'PREREG_R1_SOURCE_RECOVERY.json')['external_referent'], modelled=['X77 virtual preparation; metric scale conditional on STL'], proven=['Coordinate transport of stored X77 nodal preparation'], calibrated=[], unknown=['Pulpa/bone/canal', 'IOS independent TRE', 'metric scale', 'real osteotomy'], full_patient_gate=False)
dump(ROOT / 'raw/R1_RESULTS.json', r)
dump(ROOT / 'raw/CONTROLS_R1.json', faults)
dump(ROOT / 'HANDOFF_R1.md.json', dict(outcome=r['status'], next_operation='R2 changes representation: test full native tissue consumer on separate externally annotated P1; Demo1 goal remains open'))
(ROOT / 'HANDOFF_R1.md').write_text('# R1 handoff\n\nDemo1 has no local pulp/bone/canal observations in the inspected original source roots. Three tissue questions remain UNKNOWN. X77 nodal preparation is transported into the original CBCT frame and retained. Archive directory shows nonlocal 511 image slices, jaw predictions and jaw STL, not a local raw calibrated volume. Next: actual all-tissue operator on explicitly separate P1 validation frame; no foreign anatomy assigned to Demo1.\n')
state('R1_DECIDED', '3/3 Demo1 tissue questions UNKNOWN', 'Run R2 on independent source-labelled P1 validation frame')
print(r['status'], len(points), 'query nodes')
