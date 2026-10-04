from core import *
from evaluate import contracts, mf
import zipfile, xml.etree.ElementTree as ET

def run():
    tests = []

    def add(name, ok, detail):
        tests.append(dict(name=name, passed=bool(ok), detail=detail))
    good = next((x for x in read(R / ('RESULTS_R2_PILOT.json' if (R / 'RESULTS_R2_PILOT.json').exists() else 'RESULTS_R2.json'))['rows'] if x['status'] == 'GENERATED'))
    m = load(good['mesh_path'])
    (qv, qf) = (m['prep_vertices'], m['prep_faces'])
    (v, f) = (m['vertices'], m['faces'])
    base = float(m['base'])
    M = m['M']
    basal = np.max(abs(qv[qf] @ M[2] - base), axis=1) < 1e-07
    c = good['certificates']['Q_dilated_inside_native']
    cert = read(c['path'])
    queries = np.loadtxt(cert['query_path'], skiprows=1)
    for (name, qq) in [('radius_100mm', np.array([[*queries[0, :3], 100.0]])), ('outside_center', np.array([[*queries[0, :3] + 100, 0.01]])), ('wall_10mm', np.array([[*queries[0, :3], queries[0, 3] + 10.0]]))]:
        rr = certify(cert['support_path'], qq, 'FAULT_' + name)
        add(name, not rr['all_pass'], rr)
    nn = normals_exact(qv, qf[~basal], -M[2])
    add('reversed_insertion', not nn['pass_exact'], nn)
    (iv, iff) = triangles_mesh(v[f[m['roles'] == 1]])
    gg = Distance(qv, qf).query(sample(v[f[m['roles'] == 1]], 128))[0]
    add('film_75um', np.max(abs(gg - 0.075)) > 0.01, dict(max_error_mm=float(abs(gg - 0.075).max())))
    broken = trimesh.Trimesh(v, f[1:], process=False)
    add('deleted_facet', not broken.is_watertight, dict(watertight=bool(broken.is_watertight)))
    invalid = False
    try:
        mf.write(D / 'FAULT_regions.3mf', v, f, np.full(len(f), 99), {})
    except ValueError:
        invalid = True
    add('unknown_face_region', invalid, 'Actual3MF writer rejects out-of-contract face role')
    src = good['exports']['three_mf']
    wrong = D / 'FAULT_wrong_units.3mf'
    with zipfile.ZipFile(src) as z, zipfile.ZipFile(wrong, 'w') as w:
        for n in z.namelist():
            data = z.read(n)
            if n.endswith('.model'):
                data = data.replace(b'unit="millimeter"', b'unit="meter"')
            w.writestr(n, data)
    rejected = False
    try:
        mf.readback(wrong)
    except ValueError:
        rejected = True
    add('wrong_export_unit', rejected, 'Actual3MF reader called')
    ex = v[f[m['roles'] == 0]]
    shape = shape_score(ex + 100, m['target'])
    add('translated_shape', shape['p95_mm'] > read(R / 'PREREG_A.json')['metrics']['shape_p95_mm'][good['family']], shape)
    port = dict(case_key=good['key'], geometry_sha256=good['mesh_sha256'], frame='tooth_local', unit='mm', quantity='surface', resolution='PER_POINT', observation_state='CAD', time_scale='SIMULTANEOUS', epistemic='PROVEN_UNDER_EXPLICIT_ASSUMPTIONS')
    bad = dict(port, case_key='different')
    jj = contracts.join(port, bad)
    add('cross_specimen_port', jj['status'] == 'FAIL', jj)
    pred = dict(key=good['key'], design_mesh_sha256=good['mesh_sha256'], quantity='force', unit='N', frame='local', measurement_state='loaded', prediction_interval=None)
    obs = dict(pred, value=100, absolute_error_bound=1)
    obs.pop('prediction_interval')
    cc = contracts.compare_observation(pred, obs)
    add('unprovenanced_observation', cc['status'] == 'UNKNOWN_MISSING_MEASUREMENT_PROVENANCE', cc)
    radii = [np.array([1, 0.875, 0.75, 0.625]), np.array([1, 0.6875, 0.9375, 0.625])]
    summ = [np.array([3.0, r[0], r[-1], r.mean()]) for r in radii]
    profiles = []
    for r in radii:
        corners = np.array([[1, 1], [-1, 1], [-1, -1], [1, -1]])
        vv = np.concatenate([np.c_[corners * x, np.full(4, j)] for (j, x) in enumerate(r)])
        ff = []
        for j in range(3):
            for i in range(4):
                k = (i + 1) % 4
                ff.extend([[4 * j + i, 4 * j + k, 4 * (j + 1) + k], [4 * j + i, 4 * (j + 1) + k, 4 * (j + 1) + i]])
        ff = np.array(ff)
        mm = trimesh.Trimesh(vv, ff, process=False)
        nn = mm.face_normals
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            lp = linprog([0, 0, 0, -1], A_ub=np.c_[-nn, np.ones(len(nn))], b_ub=np.zeros(len(nn)), bounds=[(-1, 1)] * 3 + [(None, None)], method='highs', options={'threads': 1})
        profiles.append(dict(radius=r, positive_rebound_mm=float(np.maximum(np.diff(r), 0).max()), lp_margin=float(lp.x[3]), axial_exact=normals_exact(vv, ff, [0, 0, 1])))
    sf = dict(summary_names=['height_mm', 'basal_radius_mm', 'top_radius_mm', 'mean_radius_mm'], summaries=summ, identity_error=float(abs(summ[0] - summ[1]).max()), array_equal=bool(np.array_equal(*summ)), downstream_difference_mm=profiles[1]['positive_rebound_mm'] - profiles[0]['positive_rebound_mm'], profiles=profiles, minimum_extension='Maximum positive radius rebound for this fixed radial family. General3D requires oriented facet normals and incidence.', scope='Own exact discriminating fixture; not anatomical external facit')
    add('summary_does_not_determine_insertion', sf['array_equal'] and sf['identity_error'] == 0 and (sf['downstream_difference_mm'] == 0.25) and profiles[0]['axial_exact']['pass_exact'] and (not profiles[1]['axial_exact']['pass_exact']), sf)
    from wall_certificate import membership
    cr = read(good['certificates']['cavity_wall_inside_outer']['path'])
    qq = np.loadtxt(cr['query_path'], skiprows=1)
    badtri = m['vertices'][m['faces'][m['roles'] == 1]][:2] + 100
    cover = membership(badtri, qq, good['material_contract']['wall_mm'])
    add('emitted_triangle_outside_wall_cover', not cover['all_covered'], cover)
    for p in R.glob('PREREG*.json'):
        add('prereg_hash_' + p.stem, sha(p) == p.with_suffix('.sha256').read_text().split()[0], str(p))
    out = dict(tests=tests, passed=sum((x['passed'] for x in tests)), count=len(tests), all_pass=all((x['passed'] for x in tests)), sufficiency=sf)
    dump(R / 'raw/CONTROLS.json', out)
    print('CONTROLS', out['passed'], '/', out['count'])
    return out
if __name__ == '__main__':
    run()
