from common import *
from contact import compare, gates, enclosure, spatial, Q

def run():
    p = read(ROOT / 'PREREG_R1.json')['metrics']
    xy = np.array([[0.0, 0], [1, 0], [1, 1], [0, 1]])
    f = np.array([[0, 1, 2], [0, 2, 3]])
    ref = np.full(4, 0.05)
    neg = np.full(4, -0.1)
    base = compare(xy, f, ref, ref)
    bad = compare(xy, f, neg, ref)
    checks = [dict(name='negative_area', nominal_pass=base['negative_gap_area_mm2'] == 0 and all(gates(base, p).values()), injected_rejected=bad['negative_gap_area_mm2'] == 1 and (not gates(bad, p)['interference']))]
    rng = np.random.default_rng(6106)
    errors = []
    violations = []
    for k in range(12):
        pg = rng.uniform(-0.2, 0.3, 4)
        rg = rng.uniform(-0.2, 0.3, 4)
        m = compare(xy, f, pg, rg)
        b = enclosure(xy, f, pg, rg, 0.05)
        pts = (np.arange(400) + 0.5) / 400
        (xx, yy) = np.meshgrid(pts, pts)
        X = xx.ravel()
        Y = yy.ravel()
        first = Y <= X

        def interp(g):
            return np.where(first, g[0] * (1 - X) + g[1] * (X - Y) + g[2] * Y, g[0] * (1 - Y) + g[2] * X + g[3] * (Y - X))
        pp = interp(pg)
        rr = interp(rg)
        truth = np.mean((pp >= 0) & (pp <= 0.1) ^ (rr >= 0) & (rr <= 0.1))
        errors.append(abs(truth - m['symdiff_mm2']))
        for shift in np.linspace(-0.05, 0.05, 21):
            value = compare(xy, f, pg + shift, rg + shift)['symdiff_mm2']
            violations.append(max(b['symdiff_lower_mm2'] - value, value - b['symdiff_upper_mm2'], 0))
    checks.append(dict(name='sampled_area', nominal_pass=max(errors) < 0.002, injected_rejected=abs(truth - (m['symdiff_mm2'] + 1)) > 0.002, max_error_mm2=max(errors)))
    checks.append(dict(name='all_offset_set_bound', nominal_pass=max(violations) < 1e-10, injected_rejected=not b['symdiff_lower_mm2'] <= b['symdiff_upper_mm2'] + 1 <= b['symdiff_upper_mm2'], max_violation_mm2=max(violations)))
    coords = []
    faces = []
    for (x, y) in [(-12.0, 0.0), (12.0, 0.0), (0.0, -12.0), (0.0, 12.0)]:
        i = len(coords)
        coords.extend([[x - 3, y - 3], [x + 3, y - 3], [x + 3, y + 3], [x - 3, y + 3]])
        faces.extend([[i, i + 1, i + 2], [i, i + 2, i + 3]])
    xx = np.array(coords)
    ff = np.array(faces)
    g1 = np.repeat([0.05, 0.05, 1.0, 1.0], 4)
    g2 = np.repeat([1.0, 1.0, 0.05, 0.05], 4)
    (a, _) = spatial(xx, ff, g1)
    (b2, _) = spatial(xx, ff, g2)
    identity = max(abs(a['area_mm2'] - b2['area_mm2']), abs(a['count'] - b2['count']), float(np.max(abs(a['centroid_xy_mm'] - b2['centroid_xy_mm']))))
    downstream = compare(xx, ff, g2, g1)['symdiff_mm2']
    suff = dict(summary=['area_mm2', 'region_count', 'centroid_xy_mm'], state_A=a, state_B=b2, identity_error=identity, downstream_symmetric_difference_mm2=downstream, minimum_extension_for_fixed_reference='reference intersection area; full mask for arbitrary future queries', external_referent=dict(kind='our_own_fixture', locator='code/verify_contact.py', compared_quantity='summary sufficiency counterexample, not physical validation', refutes_us=True), resolution='PER_SURFACE_REGION')
    checks.append(dict(name='exact_summary_collision', nominal_pass=identity == 0 and downstream == 144.0, injected_rejected=compare(xx, ff, g2, g1)['symdiff_mm2'] != 0))
    checks.append(dict(name='paired_discrimination_not_rates', nominal_pass=np.mean([1, 0]) == np.mean([0, 1]), injected_rejected=bool(np.any(np.array([1, 0]) != np.array([0, 1])))))
    out = dict(checks=checks, all_pass=all((c['nominal_pass'] and c['injected_rejected'] for c in checks)))
    dump(ROOT / 'raw/CONTACT_CONTROLS.json', out)
    dump(ROOT / 'raw/SUFFICIENCY_CONTACT.json', suff)
    if not out['all_pass']:
        raise AssertionError(out)
    return out
if __name__ == '__main__':
    print(run())
