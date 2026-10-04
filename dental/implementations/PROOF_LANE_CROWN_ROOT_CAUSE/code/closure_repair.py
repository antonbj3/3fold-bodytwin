from oriented_crown import *
from score_joint import intersections
import oriented_crown as oc

def repaired_close(rec):
    source = load(rec['private_path'])['target']
    path = D / (rec['key'] + '_closed.mesh')
    dest = D / (rec['key'] + '_E_closed.mesh')
    proc = subprocess.run([str(D.parent / 'PROOF_LANE_FULL_CROWN_R6/repair_protected'), str(path), str(dest), str(len(source))], capture_output=True, text=True)
    if proc.returncode:
        raise ValueError('PROTECTED_REPAIR_RUNTIME ' + str(proc.returncode) + ' ' + proc.stderr[:200])
    (v, f) = meshread(dest)
    mm = trimesh.Trimesh(v, f, process=False)
    if mm.volume < 0:
        f = f[:, ::-1]
    source_set = {tuple(sorted((tuple(x) for x in t))) for t in source}
    native = np.array([tuple(sorted((tuple(x) for x in t))) in source_set for t in v[f]])
    si = intersections(v, f, rec['key'] + '_E_support')
    info = dict(border_count=1, intersection_count=si.get('count'), original_facets=len(source), final_facets=len(f), closed=bool(mm.is_watertight), native_retained_count=int(native.sum()), source_count=len(source), watertight=bool(mm.is_watertight), volume_mm3=float(abs(mm.volume)), repair=json.loads(proc.stdout), status=si['status'])
    dump(R / 'raw' / ('E_REPAIR_' + rec['key'] + '.json'), info)
    if si['status'] != 'PASS' or not mm.is_watertight or (not mm.is_winding_consistent) or (native.sum() != len(source)):
        raise ValueError('SOURCE_PRESERVING_REPAIR_FAILED ' + json.dumps(info))
    return (v, f, native, info)

def run():
    p = R / 'PREREG_E.json'
    assert not p.exists()
    base = read(R / 'PREREG_D.json')
    base.pop('frozen_utc')
    base.update(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), changed_operation='One protected CGAL patch-repair invocation,7 internal iterations, for the3 D closure failures only. All native faces must remain exactly, all geometry gates and D search unchanged.', obstacle='Three computational closing patches self-intersect. Their clinical identity is still unknown even if repaired.', selection='Exactly the3 D closure refusals; retain D passes and failures without replacement', strongest_equally_informed_control='Original failed CGAL triangulation of the same native boundaries')
    dump(p, base)
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')
    dump(R / 'DECOMPOSITION_E.json', dict(idea='Repair artificial closure without changing observed anatomy', operation=base['changed_operation'], leaves=[dict(status='DERIVED_UNDER_ASSUMPTIONS', statement='Exact native-face set retention is necessary for same-data comparison', stopping_argument='It does not give clinical validity to the artificial patch'), dict(status='CONSTITUTIVE_CLOSURE', statement='Artificial cap and selected patch repair', stopping_argument='Real complete preparation/finish line is absent')]))
    idx = {r['key']: r for r in read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']}
    bad = [r for r in read(R / 'RESULTS_D.json')['rows'] if r['status'] == 'REJECTED']
    rows = []
    oc.close = repaired_close
    for r in bad:
        try:
            out = oc.generate_d(idx[r['key']])
            out['round'] = 'E'
        except Exception as e:
            out = dict(key=r['key'], family=r['family'], status='REJECTED', reason=repr(e))
            print(out, flush=True)
        rows.append(out)
        dump(R / 'raw/E_GENERATION.json', rows)
    p = R / 'FROZEN_PREDICTIONS_E.json'
    dump(p, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), rows=rows, prereg_sha256=sha(R / 'PREREG_E.json'), code_sha256=sha(Path(__file__)), scope='No native geometry change; computational closure only'))
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')
if __name__ == '__main__':
    run()
