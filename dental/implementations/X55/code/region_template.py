"""Generate a face-region map in the CAD's explicit +Z crown frame.

CAD region masks or an intaglio reference mesh are preferred. With no intaglio
reference, geometry-based labels are a draft that must be visually checked.
"""
import argparse
import numpy as np
from metrology import load_mesh, closest, sha, write_json

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--margin-z-mm', required=True, type=float)
    p.add_argument('--margin-band-mm', default=0.45, type=float)
    p.add_argument('--intaglio-reference')
    p.add_argument('--intaglio-tolerance-mm', default=0.25, type=float)
    p.add_argument('--units', required=True, choices=['mm'])
    a = p.parse_args()
    m = load_mesh(a.design, a.units)
    x = m.triangles_center
    n = m.face_normals
    z = x[:, 2] - a.margin_z_mm
    labels = np.full(len(x), 'axial', dtype='<U8')
    if a.intaglio_reference:
        ref = load_mesh(a.intaglio_reference, a.units)
        (_, dist, _) = closest(ref, x)
        inner = dist <= a.intaglio_tolerance_mm
        source = 'CAD inner/preparation reference proximity; verify tolerance and masks'
    else:
        radial = x[:, :2] - m.centroid[:2]
        inner = (np.sum(radial * n[:, :2], axis=1) < 0) | (n[:, 2] < -0.35) & (z > a.margin_band_mm)
        source = 'DRAFT geometry heuristic; manual verification required'
    labels[inner] = 'intaglio'
    labels[~inner & (n[:, 2] > 0.35) & (z > 2)] = 'occlusal'
    labels[z <= a.margin_band_mm] = 'marginal'
    write_json(a.output, dict(design_sha256=sha(a.design), face_regions=labels.tolist(), frame='CAD +Z crown axis; margin plane z supplied by user', definition=source, region_validation_status='USER_VERIFICATION_REQUIRED', margin_z_mm=a.margin_z_mm, margin_band_mm=a.margin_band_mm, counts={k: int((labels == k).sum()) for k in np.unique(labels)}))
    print(a.output)
if __name__ == '__main__':
    main()
