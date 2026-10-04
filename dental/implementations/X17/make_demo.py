"""Export the measured field, a figure and printable lumen models. Not a MAD."""
import json, hashlib, struct, time
from pathlib import Path
import numpy as np
from scipy import ndimage
from skimage.measure import marching_cubes
from run_r1 import ROOT, DATA, write, state, freeze
EXPORT = ROOT / 'export'
EXPORT.mkdir(exist_ok=True)

def stl(path, v, f):
    tri = v[f]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-20)
    dtype = np.dtype([('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attr', '<u2')])
    arr = np.zeros(len(f), dtype=dtype)
    arr['normal'] = n
    arr['vertices'] = tri
    with path.open('wb') as fp:
        fp.write(b'X17 research airway lumen; units mm; no clinical device'.ljust(80, b' '))
        fp.write(struct.pack('<I', len(f)))
        fp.write(arr.tobytes())
    vol = abs(np.einsum('ij,ij->i', tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)
    return float(vol)

def warp(v_zyx, mask, spacing, delta, ap_sign):
    v = v_zyx.copy()
    z = np.clip(np.rint(v[:, 0] / spacing[0]).astype(int), 0, len(mask) - 1)
    x = np.clip(np.rint(v[:, 2] / spacing[2]).astype(int), 0, mask.shape[2] - 1)
    for i in range(len(mask)):
        ids = np.where(z == i)[0]
        for col in np.unique(x[ids]):
            ii = ids[x[ids] == col]
            ys = np.where(mask[i, :, col])[0]
            if not len(ys):
                continue
            ymin = (ys[0] - 0.5) * spacing[1]
            ymax = (ys[-1] + 0.5) * spacing[1]
            fraction = (v[ii, 1] - ymin) / (ymax - ymin) if ap_sign > 0 else (ymax - v[ii, 1]) / (ymax - ymin)
            v[ii, 1] += ap_sign * delta[i] * np.clip(fraction, 0, 1)
    return v

def main():
    start = time.monotonic()
    r = [json.loads((ROOT / f'raw/round_R{i}.json').read_text()) for i in range(1, 5)]
    cases = json.loads((ROOT / 'raw/cases_R1.json').read_text())
    profiles = json.loads((ROOT / 'raw/profiles.json').read_text())
    witness = json.loads((ROOT / 'raw/witnesses_R3.json').read_text())
    case = cases[0]
    name = case['case']
    d = np.load(case['crop_file']['path'])
    mask = d['mask']
    sp = d['spacing']
    offset = d['offset']
    p = profiles[name]
    lm = case['landmarks']
    pad = np.pad(mask, 1)
    phi = (ndimage.distance_transform_edt(~pad, sampling=sp) - ndimage.distance_transform_edt(pad, sampling=sp)).astype(np.float32)
    sdf = DATA / (name + '_airway_sdf.npz')
    np.savez_compressed(sdf, phi_mm=phi, spacing_zyx_mm=sp, padded_offset_zyx=offset - 1, sign_convention=np.array('negative inside'))
    (vv, ff, _, _) = marching_cubes(phi, level=0, spacing=sp)
    vv -= sp
    models = []
    basevol = stl(EXPORT / 'airway_baseline_mm.stl', vv[:, ::-1], ff[:, ::-1])
    models.append({'name': 'airway_baseline_mm.stl', 'mesh_volume_mm3': basevol, 'profile_volume_mm3': case['volume_mm3'], 'status': 'measured label -> interpolated SDF -> mesh; annotation extent'})
    w = witness[name + '_nonresponders']
    eligible = np.array(p['eligible'], bool)
    a = np.array(p['area_mm2'])
    width = np.array(p['lateral_support_mm'])
    fullz = np.array(p['z_mm'])
    for key in ['low', 'high']:
        delta = np.zeros(len(mask))
        increments = np.array(w['area_' + key + '_mm2']) - np.array(w['area0_mm2'])
        delta[eligible] = increments / width[eligible]
        moved = warp(vv, mask, sp, delta, lm['ap_sign'])
        filename = f'airway_equal_volume_{key}_mm.stl'
        vol = stl(EXPORT / filename, moved[:, ::-1], ff[:, ::-1])
        target = case['volume_mm3'] + w['budget_area_sum_mm2'] * sp[0]
        models.append({'name': filename, 'mesh_volume_mm3': vol, 'exact_profile_volume_mm3': target, 'relative_mesh_to_profile_volume_error': abs(vol - target) / target, 'status': 'conditional geometry witness, not observed anatomy; raster/mesh interpolation prevents exact matching volume in STL'})
    write('export/mesh_manifest.json', {'models': models, 'coordinate_system': 'xyz = reversed scanner array zyx; units mm; not LPS/RAS', 'not_a_splint': 'These printable models show the lumen/volume ambiguity. No clinical MAD design exported.', 'sdf': {'path': str(sdf), 'bytes': sdf.stat().st_size, 'sha256': hashlib.sha256(sdf.read_bytes()).hexdigest()}})
    import subprocess, sys, os
    subprocess.run([sys.executable, str(ROOT / 'plot_demo.py')], env=dict(os.environ, PYTHONNOUSERSITE='1'), check=True)
    results = {'lane': 'X17-sleep-apnea', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'overall_outcome': 'MEASURED_GEOMETRY_CAPABILITY; INDIVIDUAL_RESPONSE_UNKNOWN; NO_10X_CLAIM', 'rounds': r, 'external_referent': r[0]['external_referent'], 'prediction_freezes': [f'FROZEN_PREDICTIONS_R{i}.json' for i in range(1, 5)], 'clinical_CSAmin': 'UNKNOWN: fixed annotation trim and array plane not clinical palatal/epiglottis ROI', 'large_arrays': [c['crop_file'] for c in cases if 'crop_file' in c] + [{'path': str(sdf), 'sha256': hashlib.sha256(sdf.read_bytes()).hexdigest(), 'bytes': sdf.stat().st_size}], 'new_information': 'Actual expert pharynx profiles measures; separate post AP /lateral diameters in external reference .', 'new_operation': 'Scenario pose-to-profile map; equal-volume extremal profiles and conditional acquisition set; conventional controls match mathematics.', 'empirical_validation': 'One published cohort, aggregate source-aware reconstruction / group transfer. No paired TF2 post scan. No clinical prediction validated.', 'data_license': {'ToothFairy2': 'CC-BY-SA 4.0 as verified in local dataset.json', 'Shi_2023_article': 'CC-BY 4.0 as in local XML license'}, 'cost': {'demo_export_seconds': time.monotonic() - start, 'thread_limit': 4, 'gpu': False, 'own_data_bytes': sum((x.stat().st_size for x in DATA.glob('*') if x.is_file())), 'preparation_seconds': 'UNMEASURED', 'physical_measurements_performed': 0}, 'figure': ['figure.png', 'figure.pdf'], 'exports': models, 'graph_feedback': 'Coverage proposal only; working rank unavailable (JSONDecodeError). No source/generated graph edited.'}
    write('results.json', results)
    fs = {f'R{i}': json.loads((ROOT / f'FROZEN_PREDICTIONS_R{i}.json').read_text()) for i in range(1, 5)}
    write('FROZEN_PREDICTIONS.json', {'rounds': fs, 'prospective_laboratory_prediction': 'NOT_PERFORMED; baseline geometry and post-pose acquisition required before physical validation'})
    manifest = []
    for f in sorted(ROOT.rglob('*')):
        if f.is_file() and '__pycache__' not in str(f) and (f.name not in ['ARTIFACT_MANIFEST.json', 'CURRENT_WORK_STATE.json']):
            manifest.append({'path': str(f.relative_to(ROOT)), 'bytes': f.stat().st_size, 'sha256': hashlib.sha256(f.read_bytes()).hexdigest()})
    write('ARTIFACT_MANIFEST.json', manifest)
    state('Four rounds and lab demo saved', 'Individual response UNKNOWN; population transfer FAIL; affine method TIE', 'Independent review; next construction: paired, registered airway response at known incremental jaw pose, with repeated baseline')
    print('Saved figure.png, figure.pdf, results.json, three STL lumen models, and measured airway SDF', flush=True)
if __name__ == '__main__':
    main()
