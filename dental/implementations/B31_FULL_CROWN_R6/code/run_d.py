from check_intersections import *

def run():
    (a, b) = first()
    old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
    (v, f) = compact(old['vertices'], old['faces'][old['roles'] == 0])
    intersect(v, f, 'B_exterior_compact')
    p = subprocess.run([str(D6 / 'repair'), str(D6 / 'B_exterior_compact.mesh'), str(D6 / 'D_exterior.mesh')], capture_output=True, text=True)
    save(R6 / 'raw/D_REPAIR_RETRY.json', dict(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr))
    if p.returncode:
        raise RuntimeError(p.stderr)
    with (D6 / 'D_exterior.mesh').open() as h:
        (nv, nf) = map(int, h.readline().split())
        v = np.loadtxt(h, max_rows=nv)
        f = np.loadtxt(h, dtype=int)
    dest = D6 / 'D_exterior.npz'
    np.savez_compressed(dest, vertices=v, faces=f)
    freeze(R6 / 'FROZEN_PREDICTIONS_D.json', dict(rows=[dict(key=a['key'], mesh_path=dest, mesh_sha256=sha(dest))], prereg_sha256=sha(R6 / 'PREREG_D.json'), code_sha256=sha(Path(__file__))))
    x = intersect(v, f, 'D_exterior_check')
    save(R6 / 'raw/D_INTERSECTIONS.json', x)
    print('D', x.get('count'), p.stdout, flush=True)
if __name__ == '__main__':
    run()
