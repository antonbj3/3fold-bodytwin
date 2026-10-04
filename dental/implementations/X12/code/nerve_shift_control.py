from dental_release.paths import expand as _release_expand
import json, sys, time, resource, zipfile
import numpy as np
from scipy.signal import fftconvolve
from pathlib import Path
from archive_io import index, nifti
from pair_atlas import transformations, transform, ROOT
sys.path.insert(0, _release_expand('@DENTAL_INPUT_ROOT@/workspace/cells/geometry'))
from tf2_io import ZIP, ROOT as TFROOT, read_mha_bytes
start = time.time()
lookup = {r['name']: r for r in index()}
(im, _) = nifti(lookup['Pulpy3D/P2/gt_ian.nii.gz'])
source = np.asanyarray(im.dataobj) > 0
with zipfile.ZipFile(ZIP) as z:
    (lab, _, _) = read_mha_bytes(z.read(TFROOT + '/labelsTr/ToothFairy2P_002.mha'))
target = (lab == 3) | (lab == 4)
pad = np.pad(target, 20)
n1 = int(source.sum())
n2 = int(target.sum())
poses = []
for (perm, flip) in transformations(source.shape, target.shape):
    a = transform(source, perm, flip)
    coords = np.argwhere(a)
    lo = coords.min(0)
    hi = coords.max(0) + 1
    sl = tuple((slice(int(l), int(h)) for (l, h) in zip(lo, hi)))
    shape = a[sl].shape
    small = a[sl].astype(np.float32)
    region = pad[tuple((slice(int(l), int(h) + 40) for (l, h) in zip(lo, hi)))].astype(np.float32)
    overlap = fftconvolve(region, small[::-1, ::-1, ::-1], mode='valid')
    best = int(np.rint(float(overlap.max())))
    at = np.array(np.unravel_index(overlap.argmax(), overlap.shape)) - 20
    valid = np.all(lo + at >= 0) and np.all(hi + at <= np.array(target.shape))
    dice = 2 * best / (n1 + n2)
    poses.append({'permutation': perm, 'flip': flip, 'best_shift_zyx': at.tolist(), 'intersection_voxels': best, 'Dice': dice, 'preserves_all_source_voxels': bool(valid), 'gate': bool(valid and dice >= 0.95)})
    print(poses[-1], flush=True)
result = {'case': 'P2', 'poses': poses, 'any_gate_pass': any((p['gate'] for p in poses)), 'outcome': 'FAIL' if not any((p['gate'] for p in poses)) else 'PASS', 'method': 'Integer correlation over all41^3boundedtranslations for each8shape-compatibleaxis-signposes; float32FFTintersection rounded to integer; not a physical accuracy certificate', 'wall_s': time.time() - start, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024}
(ROOT / 'raw/R3_P2_SHIFT_CONTROL.json').write_text(json.dumps(result, indent=2) + '\n')
print(result['outcome'], result['wall_s'], result['peak_RSS_MiB'])
