from fc_common import *
import argparse, resource
F = functional()

def run(tag):
    root = DATA / (tag + '_predictions')
    fr = read(ROOT / ('FROZEN_PREDICTIONS_' + tag + '.json'))
    records = read(root / 'RECORDS.json')
    out = []
    st = time.perf_counter()
    for rec in records:
        r = dict(rec)
        if rec['status'] == 'EXPORTED':
            rel = rec['participant'] + '/' + rec['key'] + '/mesh.npz'
            assert sha(root / rel) == fr['files'][rel]['sha256']
            m = npz(root / rel)
            r.update(F.score_mesh(rec, m['vertices'], m['faces'], m['face_roles']))
            r['kind'] = rec.get('kind', 'shell')
            if r['reconstruction_p95_mm'] <= 0.35:
                a = npz(V4 / 'payload/whole_private' / rec['key'] / 'reference.npz')
                p = npz(V4 / 'payload/whole_inputs' / rec['key'] / 'preparation.npz')
                tri = a['source_triangles']
                tri = tri[tri.mean(1)[:, 2] >= p['margin_z']]
                ext = m['vertices'][m['faces'][m['face_roles'] == 0]]
                err = max(np.quantile(F.distance(F.make_mesh(tri), F.sample(ext, 4096)), 0.95), np.quantile(F.distance(F.make_mesh(ext), F.sample(tri, 4096)), 0.95))
                r['refined_p95_mm'] = float(err)
            print(tag, rec['key'], rec['participant'], round(r['reconstruction_p95_mm'], 4), 'closed', r['digital_closed_shell'], flush=True)
        out.append(r)
        dump(ROOT / 'raw' / f'{tag}_SCORE_CHECKPOINT.json', out)
    summary = {}
    for method in sorted({r['participant'] for r in out}):
        summary[method] = {}
        for fam in ['molar_crown', 'premolar_crown', 'anterior_crown']:
            rr = [r for r in out if r['participant'] == method and r['family'] == fam]
            sc = [r for r in rr if r.get('reconstruction_p95_mm') is not None]
            summary[method][fam] = dict(requested=len(rr), scored=len(sc), passes=sum((r['reconstruction_p95_mm'] <= 0.35 and r['digital_closed_shell'] for r in sc)), median_p95_mm=float(np.median([r['reconstruction_p95_mm'] for r in sc])) if sc else None)
    dump(ROOT / 'rounds' / f'{tag}.json', dict(round=tag, claim_type='capability', external_referent=REFERENT, rows=out, summary=summary, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, dropout_fraction=sum((r['status'] != 'SCORED' for r in out)) / len(out)))
    state(tag + '_DECIDED', 'All frozen outputs scored in unchanged v5', 'Inspect regional source of failure and choose next changed construction')
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('tag')
    a = p.parse_args()
    run(a.tag)
