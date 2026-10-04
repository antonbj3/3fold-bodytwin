from cadlib import *
from margin import *
import resource, scipy
from scipy.optimize import minimize_scalar, linprog
sys.path.insert(0, str(BASE / 'LANE_X52_PREP_ERROR_LIBRARY/code'))
x52 = module('readonly_x52_generate', BASE / 'LANE_X52_PREP_ERROR_LIBRARY/code/generate.py')
cone = module('readonly_insertion_cone', BASE / 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py')

def predict():
    st = time.perf_counter()
    rows = []
    inputs = []
    refs = {}
    tasks = sorted((DATA.parent / 'PROOF_LANE_GENCAD_V2/public').glob('*_crown_normal.json'))
    tasks = [p for p in tasks if 'molar_crown' in p.name or 'premolar_crown' in p.name]
    for p in tasks:
        key = p.stem
        die = DATA.parent / 'PROOF_LANE_GENCAD_V2/lab_pairs' / key / 'preparation_die.stl'
        if not die.exists():
            rows.append(dict(key=key, cohort='V2', status='UNAVAILABLE_SOURCE'))
            continue
        m = trimesh.load_mesh(die, process=True)
        d = detect(m)
        rec = dict(key=key, cohort='V2', **d, source=str(die), source_sha256=sha(die))
        b = m.vertices[np.isclose(m.vertices[:, 2], m.vertices[:, 2].min(), atol=1e-08)]
        b = b[np.linalg.norm(b[:, :2], axis=1) > 1e-06]
        b = b[np.argsort(np.arctan2(b[:, 1], b[:, 0]))]
        refs[key] = b
        if d['status'] == 'PROPOSED':
            roi = np.all(m.vertices[m.faces, 2] > d['points'][:, 2].min() + 1e-06, axis=1)
            roi = m.face_normals[:, 2] > -1e-06
            rec['insertion'] = cone.classify(m.face_normals[roi])
            rec['roi_scope'] = 'Normals with n_z>-1e-6 for known scan occlusal frame; a closed die basal support is excluded, not an unrestricted whole-surface proof'
        rows.append(rec)
        inputs.append(dict(path=die, sha256=sha(die)))
        print('V2', key, d['status'], flush=True)
    cases = read(BASE / 'LANE_X52_PREP_ERROR_LIBRARY/FROZEN_PREDICTIONS.json')['cases']
    for design in ['D1', 'D2', 'D3']:
        (g, a, tooth) = x52.base(design)
        for c in [c for c in cases if c['source_design'] == design]:
            t0 = time.perf_counter()
            key = c['case_id']
            (P, I, C, M, note) = x52.fields(g, a, tooth, c['fault'], c['severity'])
            m = x52.mesh(P, g, 2)
            dest = DATA / 'R1_preps' / (key + '.npz')
            dest.parent.mkdir(exist_ok=True)
            np.savez_compressed(dest, vertices=m.vertices, faces=m.faces)
            d = detect(m)
            val = m.vertices[:, 2] - x52.sample(np.broadcast_to(M, g['shape']), g, m.vertices)
            ref = scalar_contour(m, val)
            if ref:
                refs[key] = ref[0]['points']
            rows.append(dict(key=key, cohort='X52', fault=c['fault'], severity=c['severity'], source_design=design, **d, source=dest, source_sha256=sha(dest), generation_note=note, source_pitch_mm=0.12, mesh_step_mm=0.24, reference_available=bool(ref), seconds=time.perf_counter() - t0))
            print('X52', key, d['status'], flush=True)
            dump(ROOT / 'raw/R1_PARTIAL.json', rows)
        budget()
    np.savez_compressed(DATA / 'R1_references.npz', **refs)
    dump(ROOT / 'raw/R1_PREDICTIONS.json', rows)
    freeze(ROOT / 'FROZEN_PREDICTIONS_R1.json', dict(claim_type='capability', rows_sha256=sha(ROOT / 'raw/R1_PREDICTIONS.json'), prereg_sha256=sha(ROOT / 'PREREG_R1.json'), references_path=DATA / 'R1_references.npz', references_sha256=sha(DATA / 'R1_references.npz'), code={str(p): sha(p) for p in (ROOT / 'code').glob('*.py')}, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, reference_scored=False))
    state('R1_PREDICTIONS_FROZEN', len(rows), 'Score private rule-based margins')

def score():
    rows = read(ROOT / 'raw/R1_PREDICTIONS.json')
    refs = np.load(DATA / 'R1_references.npz')
    out = []
    for r in rows:
        r = dict(r)
        if r['status'] == 'PROPOSED' and r['key'] in refs:
            r['margin_error'] = errors(np.array(r['points']), refs[r['key']])
            r['margin_gate'] = 'PASS' if r['margin_error']['max_um'] <= 25 else 'FAIL'
        else:
            r['margin_gate'] = 'UNKNOWN'
        out.append(r)
    summary = {c: dict(requested=sum((r['cohort'] == c for r in out)), proposed=sum((r['cohort'] == c and r['status'] == 'PROPOSED' for r in out)), pass_count=sum((r['cohort'] == c and r['margin_gate'] == 'PASS' for r in out)), fail_count=sum((r['cohort'] == c and r['margin_gate'] == 'FAIL' for r in out)), unknown_count=sum((r['cohort'] == c and r['margin_gate'] == 'UNKNOWN' for r in out))) for c in ['V2', 'X52']}
    dump(ROOT / 'raw/R1_SCORED.json', out)
    dump(ROOT / 'RESULTS_R1.json', dict(claim_type='capability', summary=summary, external_referent=read(ROOT / 'PREREG_R1.json')['external_referent'], rows_path='raw/R1_SCORED.json', physical_margin_reference='NONE'))
    state('R1_COMPLETE', summary, 'Change feature selection/representation, preserve all failed margins')
    (ROOT / 'HANDOFF_R1.md').write_text('R1 closed high-curvature cycles: ' + json.dumps(summary) + '. The basal support can win over a finish line. This is synthetic-reference geometry; no annotated clinical preparation. Next construction must remove support-surface semantic ambiguity rather than changing the 25um gate.\n')
    print(summary)
if __name__ == '__main__':
    if not (ROOT / 'FROZEN_PREDICTIONS_R1.json').exists():
        predict()
    score()
