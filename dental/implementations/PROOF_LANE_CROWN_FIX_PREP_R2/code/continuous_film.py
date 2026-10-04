from common import *
from construction import exact_mesh

def one(row):
    m = load(row['mesh_path'])
    stem = Path(row['mesh_path']).stem
    (v, f) = (m['vertices'], m['faces'])
    active = v[f[m['roles'] == 1]]
    ip = D / (stem + '_FILM_active.mesh')
    meshwrite(ip, *triangles_mesh(active))
    pv = m.get('prep_cap_vertices', m['prep_vertices'])
    pf = m.get('prep_cap_faces', m['prep_faces'])
    pp = D / (stem + '_FILM_prep.mesh')
    meshwrite(pp, pv, pf)
    lower = exact_mesh(pp, ip, 0.04, stem + '_FILM_lower', True)
    z = subprocess.run([str(D / 'exact_gap_upper'), str(pp), str(ip)], capture_output=True, text=True, timeout=240)
    if z.returncode:
        raise ValueError('UPPER_GAP_EXIT ' + str(z.returncode) + ' ' + z.stderr[:200])
    upper = json.loads(z.stdout)
    upper.update(support_path=str(pp), support_sha256=sha(pp), query_path=str(ip), query_sha256=sha(ip))
    path = R / 'raw' / (stem + '_FILM_UPPER.json')
    dump(path, upper)
    out = dict(key=row['key'], material=row['material'], round=row['round'], continuous_pass=bool(lower['all_pass'] and upper['all_pass']), interval_mm=[0.04, 0.06], lower=lower, upper=dict(path=str(path), sha256=sha(path), **upper), scope='Every point on emitted active inner facets has nearest-distance to complete stored preparation cap/core within40–60um. Shoulder ramp and physical seated film separate.', resolution='PER_SURFACE_REGION')
    dump(R / 'raw' / (stem + '_FILM.json'), out)
    return out

def run(tag):
    rows = []
    for r in read(R / f'FROZEN_PREDICTIONS_{tag}.json')['rows']:
        if r['status'] != 'GENERATED':
            continue
        try:
            x = one(r)
        except Exception as e:
            x = dict(key=r['key'], material=r['material'], round=tag, continuous_pass=False, error=repr(e))
            __import__('traceback').print_exc()
        rows.append(x)
        dump(R / f'RESULTS_FILM_{tag}.json', dict(claim_type='capability', rows=rows))
        print(r['key'], r['material'], x['continuous_pass'], x.get('error'), flush=True)
if __name__ == '__main__':
    run(sys.argv[1])
