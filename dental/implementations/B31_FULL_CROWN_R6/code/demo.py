from score6 import *
from construct_f import repair_outer
import subprocess

def run():
    st = time.perf_counter()
    locks = read(R6 / 'DEPENDENCIES.json')
    for row in locks['files']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Frozen dependency changed ' + row['path'])
    (a, b) = first()
    old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
    (v, f, info) = repair_outer(old, 'DEMO_' + now().replace(':', '').replace('.', ''))
    prior = npz(D6 / 'D_exterior.npz')
    same_shape = v.shape == prior['vertices'].shape and f.shape == prior['faces'].shape
    error = float(np.max(np.linalg.norm(v - prior['vertices'], axis=1))) if same_shape else None
    faces_equal = bool(same_shape and np.array_equal(f, prior['faces']))
    si = intersect(v, f, 'DEMO_repaired_exterior')
    assert error == 0 and faces_equal and (si['status'] == 'PASS')
    rawpath = D6 / ('DEMO_outer_' + str(int(time.time())) + '.npz')
    np.savez_compressed(rawpath, vertices=v, faces=f)
    out = dict(claim_type='capability', geometry_coordinate_max_error_mm=error, faces_identical=faces_equal, self_intersections=si['count'], frozen_reconstruction_path=rawpath, sha256=sha(rawpath), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, scope='Fresh reconstruction of priority-molar exterior; whole crown remains subject to all failed gates', passed=True)
    save(R6 / 'raw/DEMO_REPLAY.json', out)
    print(clean(out))
if __name__ == '__main__':
    run()
