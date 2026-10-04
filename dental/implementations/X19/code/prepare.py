"""Read selected LOCAL meshes only. Store bounded geometry statistics, no mesh copies."""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
import sys, json, time, hashlib, zipfile
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1]
X11 = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts/LANE_X11_TOOTH_SEG'))
sys.path.insert(0, str(X11 / 'code'))
from morphology import mesh_obj, mesh_stl, frame

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def write(n, x):
    (R / n).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def geometry(v, lab, fr, landmarks=None):
    assert len(v) == len(lab)
    z = np.array(fr['occlusal_unit'])
    c = np.array(fr['center_mm'])
    teeth = {}
    for k in np.unique(lab):
        k = int(k)
        if k == 0:
            continue
        p = v[lab == k].astype(float)
        h = p @ z
        if len(p) < 20:
            continue
        cent = p.mean(0)
        tip = p[h >= np.quantile(h, 0.97)].mean(0)
        radial = cent - c
        radial -= z * np.dot(z, radial)
        norm = np.linalg.norm(radial)
        if norm < 1e-08:
            continue
        buccal = radial / norm
        if landmarks and str(k) in landmarks:
            buccal = np.array(landmarks[str(k)]['buccal_unit'])
            buccal /= np.linalg.norm(buccal)
        teeth[str(k)] = {'fdi': k, 'tooth_type': 'incisor' if k % 10 <= 2 else 'canine' if k % 10 == 3 else 'premolar' if k % 10 <= 5 else 'molar', 'centroid_mm': cent.tolist(), 'contact_mm': tip.tolist(), 'crown_height_mm': float(np.quantile(h, 0.95) - np.quantile(h, 0.05)), 'buccal_unit': buccal.tolist(), 'labelled_vertices': len(p), 'moment_origin': 'crown centroid; CoR UNKNOWN'}
    return teeth

def main():
    t = time.perf_counter()
    arches = []
    errors = []
    manifest = []
    visual = {}
    m = json.load(open(R / 'inputs/TEETH3DS_SELECTION.json'))
    for row in m['records']:
        (v, f) = mesh_obj(row['obj'])
        ann = json.load(open(row['annotation']))
        lab = np.array(ann['labels'], dtype=np.int32)
        (_, fr) = frame(v, f, row['jaw'], center_mode='mean', fast=True)
        key = 'T' + hashlib.sha256(row['patient'].encode()).hexdigest()[:8]
        teeth = geometry(v, lab, fr)
        arches.append({'case': key, 'dataset': 'Teeth3DS', 'jaw': row['jaw'], 'label_status': 'source expert FDI labels', 'frame': fr, 'teeth': teeth, 'vertices': len(v), 'faces': len(f)})
        manifest.extend(({'path': row[n], 'sha256': sha(row[n]), 'bytes': Path(row[n]).stat().st_size} for n in ['obj', 'annotation']))
        if len(visual) < 2:
            ix = np.linspace(0, len(v) - 1, min(6000, len(v)), dtype=int)
            visual[key + '_' + row['jaw']] = {'points': v[ix].tolist(), 'labels': lab[ix].tolist(), 'frame': fr}
    m = json.load(open(R / 'inputs/BITE2TEXT_SELECTION.json'))
    with zipfile.ZipFile(m['zip']) as zf:
        for row in m['pairs']:
            predpath = Path(m['labels_root']) / row['case_id'] / 'labels+landmarks.json'
            pred = json.load(open(predpath))
            key = 'B' + hashlib.sha256(row['case_id'].encode()).hexdigest()[:8]
            manifest.append({'path': str(predpath), 'sha256': sha(predpath), 'bytes': predpath.stat().st_size})
            for jaw in ['upper', 'lower']:
                data = zf.read(row[jaw])
                p = pred['arches'][jaw]
                assert hashlib.sha256(data).hexdigest() == p['input_sha256']
                (v, f) = mesh_stl(data, fast=True)
                labp = Path(p['labels']['path'])
                assert sha(labp) == p['labels']['sha256']
                lab = np.load(labp)['labels']
                teeth = geometry(v, lab, p['frame'], p['landmarks']['teeth'])
                arches.append({'case': key, 'dataset': 'Bite2Text', 'jaw': jaw, 'label_status': 'X11 predicted target FDI; correctness UNKNOWN', 'frame': p['frame'], 'teeth': teeth, 'vertices': len(v), 'faces': len(f)})
                manifest.extend([{'path': str(labp), 'sha256': sha(labp), 'bytes': labp.stat().st_size}, {'zip': m['zip'], 'member': row[jaw], 'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}])
    cost = {'wall_seconds': time.perf_counter() - t, 'arches': len(arches), 'mesh_reads': len(arches), 'mesh_copies': 0, 'labelled_teeth': sum((len(a['teeth']) for a in arches)), 'peak_rss_unit': 'KiB'}
    import resource
    cost['peak_rss_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    write('inputs/ARCH_GEOMETRY.json', {'schema': 'x19-geometry-v1', 'units': 'mm', 'arches': arches, 'errors': errors, 'source_manifest': manifest})
    write('inputs/VISUAL_GEOMETRY.json', visual)
    write('raw/PREP_COST.json', cost)
    print(cost)
if __name__ == '__main__':
    main()
