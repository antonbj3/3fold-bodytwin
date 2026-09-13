"""Exact shared faces of a voxel label grid; no manifold guarantee or repair.

Label 0 is exterior. Coordinates are voxel corners, with origin zero and pitch
one. Each interface appears once. Its normal points from back to front. This
module performs no artifact writes and does not initialize a compute device.
"""
import numpy as np


def label_interfaces(labels):
    labels = np.asarray(labels)
    if labels.ndim != 3 or labels.dtype != np.uint8 or min(labels.shape) < 1:
        raise ValueError('expected a nonempty three-dimensional uint8 label array')
    padded = np.pad(labels, 1)
    quads, front, back = [], [], []
    for axis in range(3):
        left = [slice(1, -1)] * 3
        right = left.copy()
        left[axis], right[axis] = slice(0, -1), slice(1, None)
        a, b = padded[tuple(left)], padded[tuple(right)]
        mask = a != b
        base = np.argwhere(mask)
        u = np.eye(3, dtype=np.int64)[(axis + 1) % 3]
        v = np.eye(3, dtype=np.int64)[(axis + 2) % 3]
        quads.append(base[:, None, :] + np.array([u * 0, u, u + v, v])[None, :, :])
        back.append(a[mask])
        front.append(b[mask])
    corners = np.concatenate(quads)
    vertices, inverse = np.unique(corners.reshape(-1, 3), axis=0, return_inverse=True)
    if len(vertices) > np.iinfo(np.int32).max:
        raise ValueError('vertex indices exceed int32 capacity')
    quad_indices = inverse.reshape(-1, 4)
    faces = quad_indices[:, [0, 1, 2, 0, 2, 3]].reshape(-1, 3)
    return dict(vertices=vertices.astype(np.float64), faces=faces.astype(np.int32),
                front=np.repeat(np.concatenate(front), 2).astype(np.int32),
                back=np.repeat(np.concatenate(back), 2).astype(np.int32))


if __name__ == '__main__':
    from pathlib import Path
    import runpy
    import sys
    probe = Path(__file__).resolve().parents[3] / 'probes/optics/label_interface_mesh_probe.py'
    sys.path.insert(0, str(probe.parent))
    runpy.run_path(str(probe), run_name='__main__')
