from construct_f import *

def repair_outer(old, key):
    (v, f) = compact(old['vertices'], old['faces'][old['roles'] == 0])
    intersect(v, f, key + '_outer_before')
    p = subprocess.run([str(D6 / 'repair'), str(D6 / (key + '_outer_before.mesh')), str(D6 / (key + '_outer_after.mesh'))], capture_output=True, text=True, timeout=30)
    if p.returncode:
        raise ValueError('repair aborted ' + p.stderr[:250])
    with (D6 / (key + '_outer_after.mesh')).open() as h:
        (n, nf) = map(int, h.readline().split())
        v = np.loadtxt(h, max_rows=n)
        f = np.loadtxt(h, dtype=int)
    return (v, f, dict(returncode=p.returncode, backend_output=p.stdout, stderr=p.stderr))
