"""Prerequisite, source identity and complete cohort; no preparation queries."""
from dental_release.paths import expand as _release_expand
import copy, hashlib, json, os, resource, time, zipfile
from pathlib import Path
import numpy as np
from archive_io import index, nifti
from tf2_io import ZIP, ROOT as TFROOT, read_mha_bytes
from verify_local_artifacts import verify
ROOT = Path(__file__).resolve().parents[1]
CANONICAL_DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X95-ct-pulp-prep-port'))
DATA = Path(os.environ.get('X95_DATA_DIR', str(CANONICAL_DATA)))

def source_path(path):
    p = Path(path)
    return DATA / p.relative_to(CANONICAL_DATA) if p.is_relative_to(CANONICAL_DATA) else p

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def write(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2) + '\n')

def main():
    started = time.monotonic()
    package = ROOT / 'sources/x12_bound'
    original_sha = sha(package / 'raw/R8_teeth.csv')
    good = verify(package)
    csv = (package / 'raw/R8_teeth.csv').read_bytes()
    lines = csv.decode().splitlines()
    fields = lines[0].split(',')
    k = fields.index('horn_boundary_min_mm')
    row = lines[1].split(',')
    old = row[k]
    row[k] = '99.0'
    lines[1] = ','.join(row)
    try:
        (package / 'raw/R8_teeth.csv').write_text('\n'.join(lines) + '\n')
        bad = verify(package)
    finally:
        (package / 'raw/R8_teeth.csv').write_bytes(csv)
    guard = {'valid_pass': not good, '99mm_corruption_rejected': str(package / 'raw/R8_teeth.csv') in bad, 'original_distance_mm': float(old), 'injected_distance_mm': 99.0, 'restoration_hash_exact': sha(package / 'raw/R8_teeth.csv') == original_sha, 'reviewer_guard_sha256': sha(ROOT / 'code/verify_local_artifacts.py'), 'source_job_status': 'INTE_REDO_UNCHANGED'}
    write(ROOT / 'raw/PREREQUISITE_GUARD.json', guard)
    assert all((guard[x] for x in ['valid_pass', '99mm_corruption_rejected', 'restoration_hash_exact'])), guard
    pairs = json.loads((package / 'raw/R6_pairs.json').read_text())
    pair = next((r for r in pairs if r['case'] == 'P48'))
    cal = next((r for r in json.loads((package / 'raw/R8_calibration.json').read_text()) if r['case'] == 'P48'))
    assert cal['accepted'] and cal['calibration']['selected_map'] == 'left_right_swap'
    records = {r['name']: r for r in index()}
    (image, ihash) = nifti(records['Pulpy3D/P48/data.nii.gz'])
    (pulpimg, phash) = nifti(records['Pulpy3D/P48/gt_instance.nii.gz'])
    with zipfile.ZipFile(ZIP) as z:
        ib = z.read(f'{TFROOT}/imagesTr/ToothFairy2P_048_0000.mha')
        lb = z.read(f'{TFROOT}/labelsTr/ToothFairy2P_048.mha')
    (ct, spacing, ih) = read_mha_bytes(ib)
    (tooth, lspacing, lh) = read_mha_bytes(lb)
    (perm, flip) = pair['image_transforms'][0]

    def pose(a):
        return a.transpose(perm)[tuple((slice(None, None, -1) if f else slice(None) for f in flip))]
    exact = np.array_equal(pose(np.asanyarray(image.dataobj)), ct)
    rawpulp = pose(np.asanyarray(pulpimg.dataobj))
    lut = np.arange(65536, dtype=np.uint16)
    for fdi in list(range(31, 39)) + list(range(41, 49)):
        lut[fdi] = fdi + 10 if fdi < 40 else fdi - 10
    pulp = lut[rawpulp]
    positive = rawpulp > 0
    swapped_containment = float(np.mean(tooth[positive] == pulp[positive]))
    wrong_containment = float(np.mean(tooth[positive] == rawpulp[positive]))
    identity = {'pair': pair, 'calibration': cal, 'full_array_equal': exact, 'n_voxels': ct.size, 'source_hash_match': ihash == pair['pulpy_CT_sha256'] and phash == pair['pulpy_pulp_sha256'] and (hashlib.sha256(ib).hexdigest() == pair['tf2_CT_sha256']) and (hashlib.sha256(lb).hexdigest() == pair['TF2_tooth_sha256']), 'grid_exact': spacing == lspacing and all((ih.get(k) == lh.get(k) for k in ['DimSize', 'TransformMatrix', 'Offset'])), 'physical_scale_mm': list(spacing), 'scale_valid': np.allclose(spacing, [0.3] * 3, rtol=0, atol=1e-12), '1mm_scale_rejected': not np.allclose([1.0] * 3, spacing, rtol=0, atol=1e-12), 'mapped_containment': swapped_containment, 'wrong_map_containment': wrong_containment, 'wrong_map_rejected': wrong_containment < 0.95, 'pulp_source_labels': np.unique(rawpulp).tolist(), 'tooth_source_labels': np.unique(tooth).tolist(), 'source_header': lh, 'wall_s': time.monotonic() - started, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024}
    if (ROOT / 'FROZEN_COHORT.json').exists():
        frozen_identity = json.loads((ROOT / 'raw/IDENTITY.json').read_text())
        assert identity['pair'] == frozen_identity['pair'] and identity['source_header'] == frozen_identity['source_header']
        write(ROOT / 'raw/IDENTITY_REPLAY.json', identity)
    else:
        write(ROOT / 'raw/IDENTITY.json', identity)
    assert exact and identity['source_hash_match'] and identity['grid_exact'] and identity['scale_valid']
    assert identity['1mm_scale_rejected'] and identity['wrong_map_rejected'] and (swapped_containment == 1)
    fdis = [int(v) for v in np.unique(tooth) if 11 <= v <= 48 and 1 <= int(v) % 10 <= 8]
    cohort = []
    for fdi in fdis:
        coords = np.argwhere(tooth == fdi)
        pc = np.argwhere(pulp == fdi)
        lo = np.maximum(coords.min(0) - 2, 0)
        hi = np.minimum(coords.max(0) + 3, tooth.shape)
        sl = tuple((slice(int(a), int(b)) for (a, b) in zip(lo, hi)))
        t = tooth[sl] == fdi
        p = pulp[sl] == fdi
        outside = int(len(pc) - np.sum(p))
        truncated = bool(np.any(coords.min(0) == 0) or np.any(coords.max(0) == np.array(tooth.shape) - 1))
        path = DATA / f'P48_FDI{fdi}_source.npz'
        np.savez_compressed(path, tooth=t, pulp=p, spacing=np.array(spacing), origin_zyx=lo)
        cohort.append({'case': 'P48', 'canonical_fdi': fdi, 'source_pulpy_fdi': fdi + 10 if fdi < 40 else fdi - 10, 'source_path': str(CANONICAL_DATA / path.name), 'sha256': sha(path), 'tooth_voxels': int(t.sum()), 'pulp_voxels': int(len(pc)), 'pulp_outside_roi_voxels': outside, 'domain_truncated': truncated, 'status': 'AVAILABLE' if len(pc) >= 30 and outside == 0 and (not truncated) else 'UNKNOWN', 'reason': 'MISSING_PULP' if not len(pc) else 'DOMAIN_TRUNCATED' if truncated else 'PULP_OUTSIDE_ROI' if outside else None})
    payload = {'frozen_utc': __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(), 'case': 'P48', 'tf2_case': 'ToothFairy2P_048', 'map': 'left_right_swap', 'frame': 'native TF2 global ZYX centres; mm physical frame in IDENTITY.json', 'cohort': cohort, 'identity_sha256': sha(ROOT / 'raw/IDENTITY.json'), 'scope': 'Every available permanent source tooth label; missing same-case pulp becomes UNKNOWN, never borrowed', 'numeric_predictions': 'NOT_YET_FROZEN; no preparation measured', 'measurement_type': 'published annotation inventories, not physical laboratory measurement'}
    dest = ROOT / 'FROZEN_COHORT.json'
    if dest.exists():
        frozen = json.loads(dest.read_text())
        assert frozen['cohort'] == cohort, 'Frozen cohort changed'
    else:
        write(dest, payload)
        (ROOT / 'FROZEN_COHORT.sha256').write_text(sha(dest) + '\n')
    write(ROOT / 'CURRENT_WORK_STATE.json', {'lane': 'X95-ct-pulp-prep-port', 'milestone': 'IDENTITY_AND_COHORT_FROZEN', 'latest_gate': 'PASS_DIGITAL_PREREQUISITE', 'next_operation': 'freeze regional solid and witness operators before questions', 'source_X12_status': 'INTE_REDO_UNCHANGED', 'cohort_teeth': len(cohort), 'available_pulp_teeth': sum((r['status'] == 'AVAILABLE' for r in cohort))})
    print(json.dumps({'guard': guard, 'identity_exact': exact, 'cohort': cohort, 'cost_s': time.monotonic() - started}))
if __name__ == '__main__':
    main()
