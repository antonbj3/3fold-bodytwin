from dental_release.paths import expand as _release_expand
from fc_common import *
import trimesh, resource
from scipy.spatial import cKDTree
sys.path.insert(0, _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V5/deps'))
from homolog_plugin import Plugin
from prescan_solid import wall_bound, norm_interval
F = functional()

def sufficiency():
    xy = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [3, 0], [3, 1], [2, 1]], float)
    f = np.array([[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]])
    za = np.array([0] * 4 + [0.25] * 4)
    zb = np.array([0.25] * 4 + [0] * 4)
    zref = np.zeros(8)
    opp = np.array([0.0625] * 4 + [0.3125] * 4)
    ref = np.c_[xy, zref][f]
    sm = F.make_mesh(ref)
    scores = []
    fn = []
    for z in [za, zb]:
        tr = np.c_[xy, z][f]
        err = max(np.quantile(F.distance(sm, F.sample(tr)), 0.95), np.quantile(F.distance(F.make_mesh(tr), F.sample(ref)), 0.95))
        scores.append(float(err))
        fn.append(F.functional_map(xy, f, opp - z))
    return dict(resolution='PER_SURFACE_REGION', summary='whole-surface bidirectional p95 mm', states_p95_mm=scores, summary_identity_error_mm=abs(scores[0] - scores[1]), both_pass_035=all((s <= 0.35 for s in scores)), function=fn, penetration_difference_mm2=fn[1]['penetration_area_mm2'] - fn[0]['penetration_area_mm2'], contact_area_difference_mm2=fn[1]['contact_area_mm2'] - fn[0]['contact_area_mm2'], minimal_sufficient_extension='spatial signed gap field on supported faces; scalar p95/contact area not sufficient', external_referent=dict(kind='our_own_fixture', locator='code/validate_demo.py:sufficiency', compared_quantity='two exact equal source p95 values with distinct spatial contact/interference', refutes_us=True))

def run():
    st = time.perf_counter()
    checks = []
    pairrows = []

    def ck(name, value, detail=None):
        checks.append(dict(name=name, pass_gate=bool(value), detail=detail))
    total = 0
    for tag in ['R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
        fr = read(ROOT / f'FROZEN_PREDICTIONS_{tag}.json')
        base = DATA / (tag + '_predictions')
        for (rel, meta) in fr['files'].items():
            ck('frozen_output_' + tag + '/' + rel, sha(base / rel) == meta['sha256'])
            total += 1
    for (rel, meta) in read(ROOT / 'FROZEN_PUBLIC_INPUTS.json')['files'].items():
        ck('public_' + rel, sha(DATA / 'public' / rel) == meta['sha256'])
    try:
        Plugin().generate(dict(public={}, method='homolog_rigid', kind='shell', hidden_reference=np.zeros(3)))
        rejected = False
    except ValueError:
        rejected = True
    ck('hidden_reference_argument_rejected', rejected)
    from contracts import validate_external
    t = dict(task_id='x', frame='site_x')
    ck('NaN_native_geometry_rejected', not validate_external(t, dict(t, units='mm', points=np.array([[np.nan, 0, 0]]))))
    ck('false_zero_support_contact_rejected', not F.validate_function(dict(status='GEOMETRIC_ONLY', support_mm2=0, contact_regions=2, contact_area_mm2=1, penetration_area_mm2=0)))
    rec = read(DATA / 'public/RECORDS.json')[0]
    a = npz(V4 / 'payload/whole_private' / rec['key'] / 'reference.npz')
    p = npz(V4 / 'payload/whole_inputs' / rec['key'] / 'preparation.npz')
    tri = a['source_triangles']
    tri = tri[tri.mean(1)[:, 2] >= p['margin_z']]
    pts = F.sample(tri)
    sm = F.make_mesh(tri)
    zero = float(np.quantile(F.distance(sm, pts), 0.95))
    shift = float(np.quantile(F.distance(sm, pts + [0, 0, 2.0]), 0.95))
    ck('exact_source_accepted', zero < 1e-07, zero)
    ck('shifted_source_rejected', shift > 0.35, shift)
    s = sufficiency()
    ck('summary_identity_exact', s['summary_identity_error_mm'] == 0)
    ck('identical_p95_distinct_function', s['penetration_difference_mm2'] != 0 and s['contact_area_difference_mm2'] != 0)
    for tag in ['R6', 'R7']:
        for rec in read(DATA / (tag + '_predictions') / 'RECORDS.json'):
            if rec['status'] != 'EXPORTED':
                continue
            base = DATA / (tag + '_predictions') / rec['participant'] / rec['key']
            q = npz(base / 'mesh.npz')
            cert = read(base / 'CERTIFICATE.json')
            (v, f, roles) = (q['vertices'], q['faces'], q['face_roles'])
            c = np.array(cert['axis_top'])
            r = cert['cavity_radius_mm']
            ext = v[f[roles == 0]]
            inn = v[f[roles == 1]]
            (vv, ff) = trimesh.remesh.subdivide_to_size(v, f[roles == 0], max_edge=0.25, max_iter=10)
            (bound, _) = wall_bound(vv[ff], c, r)
            bound -= 1e-08
            (bad, _) = wall_bound(vv[ff], c, r + 0.75)
            bad -= 1e-08
            ck('wall_full_face_' + tag + '_' + rec['key'], bound >= 0.5, bound)
            ck('inflated_cavity_rejected_' + tag + '_' + rec['key'], bad < 0.5, bad)
            mesh = trimesh.Trimesh(v, f, process=False)
            ck('removed_face_rejected_' + tag + '_' + rec['key'], not trimesh.Trimesh(v, f[:-1], process=False).is_watertight)
            nv = -np.cross(inn[:, 1] - inn[:, 0], inn[:, 2] - inn[:, 0])
            length = np.linalg.norm(nv, axis=1)
            normal = nv / length[:, None]
            nz = normal[:, 2]
            b = np.einsum('ij,ij->i', normal, inn[:, 0] - c)
            rdie = r - 0.08
            minb = float(b.min())
            film = 0.08 / r * b
            ck('all_cavity_halfspaces_downward_' + tag + '_' + rec['key'], np.min(nz) >= -1e-10, float(nz.min()))
            ck('analytic_die_fits_all_halfspaces_' + tag + '_' + rec['key'], rdie < minb, float(minb - rdie))
            ck('common_support_film_' + tag + '_' + rec['key'], float(film.min()) >= cert['normal_film_interval_mm'][0] - 1e-10 and float(film.max()) <= 0.08 + 1e-10, [float(film.min()), float(film.max())])
            ck('bad_die_radius_rejected_' + tag + '_' + rec['key'], not r + 0.2 < minb)
            dp = inn.reshape(-1, 3) - c
            dp[:, 2] = np.maximum(dp[:, 2], 0)
            normhi = norm_interval(dp)[1]
            ck('polyhedral_cavity_inside_analytic_' + tag + '_' + rec['key'], normhi.max() <= r + 1e-10, float(normhi.max() - r))
            p = npz(DATA / 'public' / (rec['key'] + '.npz'))
            world = v @ p['source_R'].T + p['source_base']
            stl = trimesh.load(base / 'crown.stl', process=True)
            delta = float(max(cKDTree(stl.vertices).query(world, workers=1)[0].max(), cKDTree(world).query(stl.vertices, workers=1)[0].max()))
            ck('STL_roundtrip_' + tag + '_' + rec['key'], delta <= 2e-05, delta)
            ck('STL_closed_' + tag + '_' + rec['key'], stl.is_watertight and stl.is_winding_consistent and (stl.volume > 0))
            pairrows.append(dict(round=tag, key=rec['key'], resolution='PER_POINT', wall_lower_mm=bound, wrong_radius_wall_lower_mm=bad, film_support_interval_mm=[float(film.min()), float(film.max())], stl_roundtrip_max_mm=delta, analytic_die_clearance_to_faceted_cavity_mm=minb - rdie, ideal_withdrawal='all nonbasal cavity outward normals n_z>=0; every downward translation remains inside their halfspaces', proof_scope='finite polyhedral/analytic geometry; physical field/film/holder UNKNOWN'))
    out = dict(checks=checks, count=len(checks), passed=sum((c['pass_gate'] for c in checks)), all_pass=all((c['pass_gate'] for c in checks)), pairs=pairrows, sufficiency=s, files_checked=total, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'raw/VALIDATION.json', out)
    print('validation', out['passed'], '/', out['count'], flush=True)
    return out
if __name__ == '__main__':
    run()
