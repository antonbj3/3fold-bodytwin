from oriented_crown import *
import oriented_crown as oc
from score_joint import intersections, shape_metrics
from geometry import curve_error

def repaired_close(rec):
    S = load(rec['private_path'])['target']
    (vv, ff) = meshread(D / (rec['key'] + '_closed.mesh'))
    A = vv[ff[len(S):]].copy()
    dest = D / (rec['key'] + '_G_closed.mesh')
    proc = subprocess.run([str(D.parent / 'PROOF_LANE_FULL_CROWN_R6/repair'), str(D / (rec['key'] + '_closed.mesh')), str(dest)], capture_output=True, text=True)
    if proc.returncode:
        raise ValueError('WHOLE_SUPPORT_REPAIR_RUNTIME ' + str(proc.returncode) + ' ' + proc.stderr[:300])
    (v, f) = meshread(dest)
    mm = trimesh.Trimesh(v, f, process=False)
    if mm.volume < 0:
        f = f[:, ::-1]
    si = intersections(v, f, rec['key'] + '_G_support')
    source = {tuple(sorted((tuple(x) for x in t))) for t in S}
    closure = {tuple(sorted((tuple(x) for x in t))) for t in A}
    native = []
    for t in v[f]:
        k = tuple(sorted((tuple(x) for x in t)))
        native.append(True if k in source else False if k in closure else None)
    unknown = np.array([x is None for x in native])
    p = v[f[unknown]].mean(1)
    n_new = fast_nearest(S, p)[1] < fast_nearest(A, p)[1] if len(p) else []
    native = np.array([False if x is None else x for x in native])
    native[unknown] = n_new
    if not si.get('status') == 'PASS' or not mm.is_watertight:
        raise ValueError('REPAIRED_SUPPORT_STILL_INVALID ' + json.dumps(si))
    repaired = v[f[native]].copy()
    error = shape_metrics(repaired, S)
    outer = trimesh.Trimesh(v, f[native], process=False)
    ll = loops(outer)
    p = load(rec['public_path'])
    margin = curve_error(v[ll[0]], p['margin_curve']) if len(ll) == 1 else None
    pr = read(R / 'PREREG_C.json')['metrics']
    info = dict(intersection_count=si.get('count'), closed=bool(mm.is_watertight), native_retained_count=sum((tuple(sorted((tuple(x) for x in t))) in source for t in repaired)), source_count=len(S), volume_mm3=float(abs(mm.volume)), original_shape_error=error, original_margin_error_mm=margin, changed_facets=int(unknown.sum()), native_repair_allowed=True)
    dump(R / 'raw' / ('G_REPAIR_' + rec['key'] + '.json'), info)
    if error['p95_mm'] > pr['shape_p95_mm'][rec['family']] or margin is None or margin > 0.025:
        raise ValueError('SOURCE_REPAIR_EXCEEDS_ORIGINAL_GATES ' + json.dumps(info))
    return (v, f, native, info)

def run():
    p = R / 'PREREG_G.json'
    assert not p.exists()
    x = read(R / 'PREREG_D.json')
    x.pop('frozen_utc')
    x.update(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), changed_operation='Permit local mesh repair of the3 refused digital supports, then rerun identical D search and Boolean.100% source-facet identity is replaced by the original R4 p95 shape and25um margin error gates; report every changed facet and actual shape/margin discrepancy. Use existing full-support CGAL repair with3 internal iterations.', obstacle='F proves2 native meshes self-intersect, impossible under exact retention; third has closure-only intersections. Geometry repair within original tolerances is a different admissible contract.', selection='The3 original D refusals, no replacement teeth', strongest_equally_informed_control='D/E exact-retention construction on same original source surfaces')
    dump(p, x)
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')
    dump(R / 'DECOMPOSITION_G.json', dict(idea='Exact mesh retention is not necessary for preserving measured shape within its original acceptance tolerance', operation=x['changed_operation'], leaves=[dict(status='DERIVED_UNDER_ASSUMPTIONS', statement='Exact obstruction is removed only by explicitly allowing a different surface', stopping_argument='Retain all original metrics; no claim of true anatomical correction without observations'), dict(status='CONSTITUTIVE_CLOSURE', statement='Local CGAL repair', stopping_argument='Numerical repair, not new anatomical measurement')]))
    idx = {r['key']: r for r in read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']}
    bad = [r for r in read(R / 'RESULTS_D.json')['rows'] if r['status'] == 'REJECTED']
    rows = []
    oc.close = repaired_close
    for r in bad:
        try:
            out = oc.generate_d(idx[r['key']])
            out['round'] = 'G'
            oldpath = Path(out['mesh_path'])
            newpath = D / (r['key'] + '_G.npz')
            oldpath.replace(newpath)
            out.update(mesh_path=str(newpath), mesh_sha256=sha(newpath))
            out['reference_record'] = idx[r['key']]
        except Exception as e:
            out = dict(key=r['key'], family=r['family'], status='REJECTED', reason=repr(e))
            print(out, flush=True)
        rows.append(out)
        dump(R / 'raw/G_GENERATION.json', rows)
    p = R / 'FROZEN_PREDICTIONS_G.json'
    dump(p, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), rows=rows, prereg_sha256=sha(R / 'PREREG_G.json'), code_sha256=sha(Path(__file__)), scope='Native repair explicitly allowed under original shape/margin tolerances'))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')
if __name__ == '__main__':
    run()
