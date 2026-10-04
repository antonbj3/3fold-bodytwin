from common import *
import sys, trimesh
sys.path.insert(0, str(V4 / 'payload/participant_code'))
from task_io import attach
from legacy.geometry import shell
from whole import three_mf

def run():
    records = []
    panel = read(ROOT / 'PANEL_SELECTION.json')
    for tag in ['volume_1.0', 'volume_1.1']:
        for case in panel['selected']:
            key = case['case_key']
            ts = [t for t in read(V4 / 'payload/public/tasks' / f'{key}.json') if t['family'] in panel['families']]
            with np.load(DATA / 'predictions' / tag / f'{key}.npz') as a:
                for (i, t) in enumerate(ts):
                    if str(a['status'][i]) != 'DESIGN':
                        continue
                    with np.load(V4 / 'payload/public' / t['geometry_file']) as s:
                        t = attach(t, dict(s))
                    z = a['outer_' + str(i)]
                    inner = a['inner_' + str(i)]
                    (v, f) = shell(t['xy'], z, inner, t['faces'])
                    mesh = trimesh.Trimesh(v, f, process=True)
                    if mesh.volume < 0:
                        mesh.invert()
                    dest = DATA / 'roof_exports' / tag / t['task_id']
                    dest.mkdir(parents=True, exist_ok=True)
                    mesh.export(dest / 'specimen.stl')
                    three_mf(mesh, dest / 'specimen.3mf')
                    records.append(dict(task_id=t['task_id'], participant=tag, frame=t['frame'], units='mm', watertight=bool(mesh.is_watertight), volume_mm3=float(mesh.volume), scope='Restricted occlusal shell specimen; not a complete axial crown or clinically approved restoration', files=inventory(dest), relative_directory=str(dest.relative_to(DATA))))
    dump(ROOT / 'raw/ROOF_EXPORTS.json', records)
    freeze(ROOT / 'FROZEN_ROOF_EXPORTS.json', dict(records=records, measurement='NOT_RUN', source_prediction_hashes={tag: sha(ROOT / f'FROZEN_PREDICTIONS_{tag}.json') for tag in ['volume_1.0', 'volume_1.1']}, physical_manufacturability='UNKNOWN', time_scale='HANDOVER', resolution='PER_POINT'))
    print('roof shell specimens', len(records), 'watertight', sum((r['watertight'] for r in records)), 'bytes', budget())
if __name__ == '__main__':
    run()
