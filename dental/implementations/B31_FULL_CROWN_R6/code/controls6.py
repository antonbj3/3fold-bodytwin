from score6 import *
from construct_g import retained_faces

def run():
    st = time.perf_counter()
    rows = []
    folder = D6 / 'controls'
    folder.mkdir(exist_ok=True)

    def add(name, good, bad, **kw):
        rows.append(dict(name=name, good_pass=bool(good), injected_rejected=not bool(bad), pass_control=bool(good) and (not bool(bad)), **kw))
    base = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]]])
    q = base.mean(1)
    gap_a = np.array([0.03125, 0.03125])
    gap_b = np.array([0.015625, 0.046875])
    a = base.copy()
    b = base.copy()
    a[:, :, 2] += gap_a[:, None]
    b[:, :, 2] += gap_b[:, None]
    da = closest(base, a.mean(1))[1]
    db = closest(base, b.mean(1))[1]
    suff = dict(mean_a_mm=float(da.mean()), mean_b_mm=float(db.mean()), identity_error_mm=float(abs(da.mean() - db.mean())), minimum_a_mm=float(da.min()), minimum_b_mm=float(db.min()), downstream_difference_mm=float(da.min() - db.min()), gate_a=x49.classify_interval([float(da.min())] * 2, 0.025), gate_b=x49.classify_interval([float(db.min())] * 2, 0.025), minimum_extension='Minimum gap distinguishes this pair. General seating/flow still needs a spatial field and constitutive data.', resolution=dict(mean='PER_TOOTH', minimum='PER_SURFACE_REGION'), external_referent=read(R6 / 'PREREG_CONTROLS.json')['external_referent'])
    assert suff['identity_error_mm'] == 0
    save(R6 / 'raw/SUFFICIENCY.json', suff)
    add('gap', da.min() >= 0.025, db.min() >= 0.025)
    wgood = base.copy()
    wbad = base.copy()
    wgood[:, :, 2] += 0.625
    wbad[:, :, 2] += 0.375
    add('wall', closest(base, wgood.mean(1))[1].min() >= 0.5, closest(base, wbad.mean(1))[1].min() >= 0.5)
    (a0, b0) = first()
    public = npz(b0['public_path'])
    target = npz(b0['private_path'])['target']
    dgood = fast_nearest(target, sample(target, 128))[1].max()
    moved = sample(target, 128) + np.array([100.0, 0.0, 0.0])
    dbad = fast_nearest(target, moved)[1].min()
    add('shape', dgood <= 0.437, dbad <= 0.437, good_mm=dgood, bad_mm=dbad)
    cv = public['margin_curve']
    cvbad = cv.copy()
    cvbad[0, 0] += 0.125
    add('margin_identity', np.array_equal(cv, cv.copy()), np.array_equal(cv, cvbad))
    tv = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    tf = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]])
    tri = tv[tf]
    valid = folder / 'valid.stl'
    ascii_stl(valid, tri)
    dat = x49.G.load(valid, 1.0)
    h = x49.mesh_health({'crown': dat})
    good = h['status'] == 'PASS'
    wrong = folder / 'normal.stl'
    lines = valid.read_text().splitlines()
    i = next((i for (i, l) in enumerate(lines) if l.startswith('facet normal')))
    xyz = np.array(list(map(float, lines[i].split()[2:])))
    lines[i] = 'facet normal ' + ' '.join(map(str, -xyz))
    wrong.write_text('\n'.join(lines) + '\n')
    add('STL_stored_normal', good, x49.mesh_health({'crown': x49.G.load(wrong, 1.0)})['status'] == 'PASS')
    path = folder / 'open.stl'
    ascii_stl(path, tri[:-1])
    add('watertight', good, x49.mesh_health({'crown': x49.G.load(path, 1.0)})['status'] == 'PASS')
    path = folder / 'winding.stl'
    badtri = tri.copy()
    badtri[0] = badtri[0, ::-1]
    ascii_stl(path, badtri)
    add('winding', good, x49.mesh_health({'crown': x49.G.load(path, 1.0)})['status'] == 'PASS')
    add('units', x49.scale_rule(dat, 'mm')['status'] == 'PASS', x49.scale_rule(x49.G.load(valid, 1000.0), 'mm')['status'] == 'PASS')
    crossv = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.25, 0.25, -1.0], [0.25, 0.25, 1.0], [0.75, 0.25, 0.0]])
    cf = np.array([[0, 1, 2], [3, 4, 5]])
    si_good = intersect(tv, tf, 'control_si_valid')
    si_bad = intersect(crossv, cf, 'control_si_cross')
    add('self_intersection', si_good['status'] == 'PASS', si_bad['status'] == 'PASS', bad_pairs=si_bad.get('count'))
    add('source_retention', retained_faces(tv, tf, tv, tf).all(), retained_faces(tv, tf[:-1], tv, tf).all())
    pp = R6 / 'FROZEN_PREDICTIONS_F_FIRST.json'
    origsha = sha(pp)
    copy = folder / 'tampered.json'
    copy.write_bytes(pp.read_bytes() + b' ')
    add('hash_binding', sha(pp) == origsha, sha(copy) == origsha)
    out = dict(controls=rows, count=len(rows), all_pass=all((r['pass_control'] for r in rows)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, sufficiency_file='raw/SUFFICIENCY.json')
    save(R6 / 'raw/CONTROLS.json', out)
    print(clean(out))
if __name__ == '__main__':
    run()
