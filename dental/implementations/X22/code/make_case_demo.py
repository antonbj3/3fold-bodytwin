import argparse
import numpy as np
from common import *

def main():
    p = argparse.ArgumentParser(description='Generate known voxel lesion truth and local DRRs from a segmented research tooth')
    p.add_argument('--input', required=True, help='NPZ labels 0outside 1enamel 2dentin 3pulp, spacing_mm and origin_mm; or mask with explicit --shell-mm')
    p.add_argument('--shell-mm', type=float, default=None)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    q = np.load(a.input, allow_pickle=False)
    spacing = np.asarray(q['spacing_mm']).ravel()
    if spacing.size == 1:
        spacing = np.repeat(spacing, 3)
    if spacing.size != 3 or not np.allclose(spacing, spacing[0]) or np.any(spacing <= 0):
        raise ValueError('Positive isotropic spacing required; resample externally and track error')
    h = float(spacing[0])
    origin = np.asarray(q['origin_mm']) if 'origin_mm' in q else np.zeros(3)
    if a.shell_mm is not None:
        if not 0.15 <= a.shell_mm < 3:
            raise ValueError('Shell closure must be 0.15..3mm')
        mask = q['mask'] if 'mask' in q else q['labels'] > 0
        (labels, _) = tissue(mask, h, a.shell_mm)
        provenance = 'EXPLICIT_SYNTHETIC_TISSUE_CLOSURE'
    else:
        labels = q['labels']
        if not np.all(np.isin(labels, [0, 1, 2, 3])):
            raise ValueError('Unrecognized tissue class')
        provenance = 'CALLER_SUPPLIED_TISSUE_SEGMENTATION_UNVERIFIED'
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    mu = coefficient(labels)
    rng = np.random.default_rng(22203)
    records = []
    for site in ['proximal', 'occlusal']:
        (_, center) = choose_site(labels, site, h, origin, 1.2)
        for name in CLASSES[1:]:
            loss = 0.4 * lesion(labels, h, site, center, name)
            truth = verify_truth(labels, loss, h, site)
            if truth['actual_class'] != name:
                raise ValueError(f'Finite-resolution placement failed: {name} -> {truth}')
            for angle in [-10, 0, 10]:
                tau = project(mu - loss * (mu - 0.02683), h, angle)
                expected = transmission(tau, h)
                fn = out / f'{site}_{name}_{angle}deg.npz'
                np.savez_compressed(fn, image=noisy_image(expected, 20000, rng).astype(np.float32), expected=expected, mineral_loss=loss, labels=labels, spacing_mm=spacing, origin_mm=origin)
                records.append({'path': str(fn), 'sha256': sha(fn), 'site': site, 'depth_class': name, 'angle_deg': angle, **truth})
    write(out / 'manifest.json', {'input_path': a.input, 'input_sha256': sha(a.input), 'tissue_provenance': provenance, 'physics': 'uncalibrated effective40keV; no clinical realism claim', 'records': records})
    print('Saved', len(records), 'DRRs with voxel truth to', out)
if __name__ == '__main__':
    main()
