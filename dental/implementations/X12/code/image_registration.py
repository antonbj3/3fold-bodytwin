"""Image-only cavity contrast registration; whole-tooth labels are held out."""
from dental_release.paths import expand as _release_expand
import itertools, json, sys, time, zipfile
from pathlib import Path
import numpy as np
from archive_io import nifti
from pair_atlas import ROOT, transformations, transform_coords
sys.path.insert(0, _release_expand('@DENTAL_INPUT_ROOT@/workspace/cells/geometry'))
from tf2_io import ZIP, ROOT as TFROOT, read_mha_bytes

def pose_scores(image, c, shifts, radius=3):
    out = []
    shape = np.array(image.shape)
    for start in range(0, len(shifts), 128):
        s = shifts[start:start + 128]
        points = c[None, :, :] + s[:, None, :]
        valid = np.all((points >= radius) & (points < shape - radius), axis=(1, 2))
        points = points[valid]
        sv = s[valid]
        if not len(points):
            continue
        center = image[tuple(points.transpose(2, 0, 1))]
        collar = np.zeros_like(center, dtype=float)
        for axis in range(3):
            for sign in [-1, 1]:
                q = points.copy()
                q[:, :, axis] += sign * radius
                collar += image[tuple(q.transpose(2, 0, 1))] / 6
        score = np.mean(collar - center, axis=1)
        out.extend(((float(v), tuple((int(x) for x in shift))) for (v, shift) in zip(score, sv)))
    return sorted(out, reverse=True)

def register(image, pulp):
    samples = []
    for fdi in list(range(31, 39)) + list(range(41, 49)):
        c = np.argwhere(pulp == fdi)
        if len(c):
            samples.append(c[np.linspace(0, len(c) - 1, min(64, len(c)), dtype=int)])
    if not samples:
        return None
    source = np.concatenate(samples)
    best = []
    coarse = np.array(list(itertools.product(range(-20, 21, 4), repeat=3)), dtype=int)
    for (perm, flip) in transformations(pulp.shape, image.shape):
        c = transform_coords(source, pulp.shape, perm, flip)
        scores = pose_scores(image, c, coarse)
        if not scores:
            continue
        b = np.array(scores[0][1])
        fine = np.array([q for q in itertools.product(*(range(max(-20, int(v) - 4), min(20, int(v) + 4) + 1) for v in b))], dtype=int)
        scores = pose_scores(image, c, fine)
        best.extend(({'score': v, 'shift': s, 'perm': perm, 'flip': flip} for (v, s) in scores[:2]))
    best.sort(key=lambda x: x['score'], reverse=True)
    return {'best': best[:2], 'margin': best[0]['score'] - best[1]['score'], 'n_points': len(source)} if best else None

def main():
    lookup = {r['name']: r for r in json.load(open(ROOT / 'raw/local_records.json'))}
    out = []
    t0 = time.time()
    with zipfile.ZipFile(ZIP) as z:
        for pid in ['P1', 'P2', 'P4', 'P5', 'P9']:
            (pimg, _) = nifti(lookup['Pulpy3D/' + pid + '/gt_instance.nii.gz'])
            p = np.asanyarray(pimg.dataobj)
            case = 'ToothFairy2P_' + pid[1:].zfill(3)
            raw = z.read(TFROOT + '/imagesTr/' + case + '_0000.mha')
            (image, sp, h) = read_mha_bytes(raw)
            fitted = register(image, p)
            frozen = {'case': pid, 'fitted_before_tooth_validation': True, 'registration': fitted}
            out.append(frozen)
            (ROOT / 'raw/R5_image_registration_fits.json').write_text(json.dumps(out, indent=2) + '\n')
            (tooth, _, _) = read_mha_bytes(z.read(TFROOT + '/labelsTr/' + case + '.mha'))
            if fitted:
                b = fitted['best'][0]
                c = np.argwhere(p > 0)
                v = p[tuple(c.T)]
                mapped = transform_coords(c, p.shape, b['perm'], b['flip']) + np.array(b['shift'])
                valid = np.all((mapped >= 0) & (mapped < np.array(image.shape)), axis=1)
                inside = np.zeros(len(c), dtype=bool)
                inside[valid] = tooth[tuple(mapped[valid].T)] == v[valid]
                frozen['overall_containment'] = float(np.mean(inside))
                frozen['per_FDI'] = {str(fdi): float(np.mean(inside[v == fdi])) for fdi in np.unique(v)}
                frozen['gate'] = bool(fitted['best'][0]['score'] > 0 and fitted['margin'] >= 0.5 and (frozen['overall_containment'] >= 0.95))
            print(json.dumps(frozen), flush=True)
            (ROOT / 'raw/R5_image_registration_validation.json').write_text(json.dumps(out, indent=2) + '\n')
    print('wall_s', time.time() - t0)
if __name__ == '__main__':
    main()
