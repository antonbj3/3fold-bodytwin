from common import *
from annotated_pairs import read_obj
import trimesh

def run():
    if not (ROOT / 'PREREG_COMPONENT_AUDIT.json').exists():
        freeze(ROOT / 'PREREG_COMPONENT_AUDIT.json', dict(claim_type='capability', capability='Verify that library component selection does not silently add measured tooth surface', obstacle='Trimesh split defaults to repair=True; negative dropped-area fraction seen in first low-quality cohort', changed_operation='Compare default split against explicit repair=False on all36 annotated teeth, before deciding whether calibration surfaces were altered', consumer='Annotated variability envelope provenance', metrics=dict(added_faces_allowed=0, area_difference_tolerance_mm2=1e-10), decision='If synthesized faces exist, flag calibration sensitivity and require separate preserved strict recomputation. Do not overwrite frozen envelope.', strongest_equally_informed_control='Original mesh plus strict face subset', falsifier='Any added face changes measurement-surface provenance', full_cost='36 small tooth component splits, no new datasets/fit; inherited read cost unknown', external_referent=dict(kind='published_dataset', locator='https://osf.io/xctdy/', compared_quantity='original annotated triangle counts and areas', refutes_us=True), resolution='PER_SURFACE_REGION'))
    rows = []
    for src in read(ROOT / 'raw/ANNOTATED_PAIRS.json')['sources']:
        (v, f) = read_obj(src['path'])
        labels = np.array(read(src['labels_path'])['labels'])
        for fdi in [11, 14, 16, 21, 24, 26]:
            tri = v[f[np.all(labels[f] == fdi, axis=1)]]
            m = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=True)
            normal = max(m.split(only_watertight=False), key=lambda a: a.area)
            strict = max(m.split(only_watertight=False, repair=False), key=lambda a: a.area)
            rows.append(dict(case_key=src['case_key'], fdi=fdi, default_faces=len(normal.faces), strict_faces=len(strict.faces), area_difference_mm2=normal.area - strict.area, pass_no_added_faces=len(normal.faces) == len(strict.faces) and abs(normal.area - strict.area) < 1e-10))
    dump(ROOT / 'raw/COMPONENT_AUDIT.json', dict(rows=rows, passed=sum((x['pass_no_added_faces'] for x in rows)), total=len(rows)))
    return rows
if __name__ == '__main__':
    run()
