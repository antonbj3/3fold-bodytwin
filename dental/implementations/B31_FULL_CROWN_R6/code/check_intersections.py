from local import *
import subprocess

def intersect(v, f, key):
    dest = D6 / (key + '.mesh')
    with dest.open('w') as h:
        h.write(f'{len(v)} {len(f)}\n')
        np.savetxt(h, v, fmt='%.17g')
        np.savetxt(h, f, fmt='%d')
    p = subprocess.run([str(D6 / 'intersections'), str(dest)], capture_output=True, text=True)
    if p.returncode:
        return dict(status='UNKNOWN', error=p.stderr, exit=p.returncode)
    x = json.loads(p.stdout)
    x.update(status='PASS' if not x['count'] else 'FAIL', backend='CGAL5.4 EPIK exact predicates on stored binary64 vertices', mesh_sha256=sha(dest))
    return x

def run():
    (a, b) = first()
    m = npz(a['mesh_path'])
    old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
    c = npz(D / 'C_IMPLICIT_CACHE' / (a['key'] + '.npz'))
    out = {}
    for (name, v, f) in [('source_cap', m['vertices'], m['faces'][m['roles'] == 1]), ('source_closed', c['prep_vertices'], c['prep_faces']), ('B_exterior', old['vertices'], old['faces'][old['roles'] == 0]), ('B_intaglio', old['vertices'], old['faces'][old['roles'] == 1]), ('B_shell', old['vertices'], old['faces'])]:
        x = intersect(v, f, name)
        out[name] = x
        save(R6 / 'raw/SOURCE_INTERSECTIONS.json', out)
        print(name, x.get('count'), x['status'], flush=True)
if __name__ == '__main__':
    run()
