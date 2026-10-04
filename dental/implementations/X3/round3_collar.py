"""R3 information intervention: observed natural cervical collar, new patients.
Different input task from full tooth deletion. No withheld upper tooth surface
is passed to a method; the observed1mm collar is deliberately available.
"""
from crownbench import *
from ring_method import edge_components
import copy, resource, collections

def freeze():
    if (P / 'PREREG_R3.json').exists():
        return
    parent = json.load(open(P / 'RESULTS_R2.json'))
    root = pathlib.Path(json.load(open(P / 'PREREG_R1.json'))['data']['root'])
    original = json.load(open(P / 'SPLIT.json'))
    paths = {p.stem: p for p in root.glob('data_part_*/[ul]*/*/*.json')}
    both = {k.rsplit('_', 1)[0] for k in paths if k.endswith('_lower') and k.rsplit('_', 1)[0] + '_upper' in paths}
    pool = set()
    for name in ['testing_lower.txt', 'testing_upper.txt']:
        pool |= {s.rsplit('_', 1)[0] for s in (root / 'Teeth3DS_train_test_split' / name).read_text().splitlines()}
    excluded = set(sum(original['patients'].values(), []))
    cases = sorted(pool & both - excluded, key=lambda c: hashlib.sha256(('X3-crown-v1/' + c).encode()).hexdigest())[:4]
    records = []
    for case in cases:
        for jaw in ['upper', 'lower']:
            labels = paths[case + '_' + jaw]
            m = labels.with_suffix('.obj')
            records.append({'case': case, 'jaw': jaw, 'labels': str(labels), 'mesh': str(m), 'labels_sha256': sha(labels), 'mesh_sha256': sha(m)})
    put(P / 'SPLIT_R3.json', {'frozen_utc': now(), 'cases': cases, 'scans': records, 'disjoint_from_R1_R2_train': not set(cases) & excluded})
    pre = {'round': 'R3_OBSERVED_CERVICAL_COLLAR', 'frozen_utc': now(), 'capability': 'Quantify whether one locally observed cervical tooth band makes missing-upper-crown placement self-calibrating on a distinct natural IOS cohort. This is a partial crown-completion task, not the full-deletion benchmark.', 'parent': {'results_sha256': sha(P / 'RESULTS_R2.json'), 'decision': parent['decision']}, 'obstacle': 'A gingival mask rim is an indirect positional anchor. Asymmetric cervical tooth shape/axis remains unknown without local tooth-surface observations.', 'changed_operation': 'Retain only target vertices at or below median natural cervical boundary height plus1.0mm along visible arch normal. Remove all target vertices above this and incident faces. Methods may read observed collar; hidden upper surface is evaluator-only. Fit reflected homolog collar to observed collar by robust two-sided six-parameter ICP, then optionally SDF roundtrip.', 'consumer': 'A researcher comparing information needed for crown completion, before a prepared-abutment/registered-antagonist study.', 'selection': 'Four NEW official test-union patients, both annotated arches; same SHA256 rank with all R1/R2/train patients excluded. Frozen scans/hash SPLIT_R3.json.32 scheduled teeth; same FDI2,3,5,6.', 'observations': 'Only actual surviving natural tooth surface within1mm cervical band; synthetic cutting operation on independently acquired IOS. Not a clinical preparation.', 'methods': ['mirror', 'collar_mesh', 'collar_sdf'], 'fit': 'At most512 collar/template surface area samples. Source patch<=source median cervical boundary height+1mm. Translation+-3mm,axial rotation+-0.4rad,xy scales0.70..1.30. Bidirectional nearest-point residual /0.30mm; regularizer0.35 translations,0.35/0.20 angle/scales;soft_l1,max80eval. Visible lower band only; no fit to hidden points.', 'metrics': json.load(open(P / 'PREREG_R1.json'))['metrics'], 'gates': {'data': '>=3 of4 patients,>=3 tooth classes and>=80% of32 scheduled outputs', 'primary': 'patient-balanced mean occlusal_p95<=0.90*mirror; upper95% paired bootstrap CI<0,1000resamples seed6103', 'absolute': 'mean occlusal_p95<=0.50mm', 'proximal': 'mean local proximal p95 error<=mirror+0.05mm', 'strong_control': 'collar_mesh gets identical observations and fit; SDF vs mesh tolerance0.05mm', 'numerical': 'original train SDF roundtrip<=0.30mm', 'decision': 'All anatomical gates required for positive. Any evaluable failed gate NEGATIVE. Observed collar is new information, not an algorithmic novelty claim. Small independent cohort limits population inference.'}, 'falsifiers': ['collar anchoring improves lower surface but not upper held-out region', 'contralateral shape mismatch survives a good cervical alignment', 'geometry fit cannot cover80% of scheduled cohort', 'SDF error defeats direct mesh control'], 'external_referent': {'kind': 'published_dataset', 'locator': 'https://arxiv.org/abs/2210.06094 ; local hashes inSPLIT_R3.json', 'compared_quantity': 'Withheld upper natural crown surface and original neighbor proximity on NEW published IOS scans; observed lower1mm band is an explicit input intervention', 'refutes_us': True}, 'full_cost': {'preparation': '8source hash/imports and32 collar masks', 'fit': 'collar extraction,robust ICP and SDF per eligible tooth, charged to both standalone methods', 'discovery': 'prior R1/R2 failures retained', 'validation': 'all predicted surfaces frozen before hidden reference reads; exact distance kernel certified;patient bootstrap; borrowed metric-injection correctness checked R1/R2 and all thresholds reject injected values', 'queries': 'all32 scheduled outputs incl failure', 'fallback': 'no silent fallback; registered-antagonist/preparation/clinical quality UNKNOWN'}, 'resources': {'threads': 1, 'estimated_ram_gb': 1.5, 'lane_disk_limit_bytes': 3000000000}, 'split_sha256': sha(P / 'SPLIT_R3.json'), 'source_code_sha256': sha(pathlib.Path(__file__))}
    put(P / 'PREREG_R3.json', pre)
    put(P / 'DECOMPOSITION_R3.json', {'chain': 'Observed natural cervical surface -> target-local pose/width information -> template registration -> hidden upper crown prediction -> real withheld-surface evaluation', 'leaves': [{'name': 'observed_collar', 'status': 'EXTERNALLY_MEASURED', 'equation': 'S_obs={x inS_target :n_arch dot x<=median(n_arch dot cervical_boundary)+1mm}', 'stop': 'Acquired IOS points transformed by synthetic cut. This information does not exist in full tooth deletion or necessarily on a real prepared tooth.'}, {'name': 'registration', 'status': 'CONSTITUTIVE_CLOSURE', 'equation': 'min_p sum rho(||T_p(S_source_collar)-NN(S_obs)||/.30)^2 + reverse + prior', 'stop': 'Contralateral lower band is not guaranteed homologous; the target upper morphology remains unobserved. Fitted cervical accuracy does not establish cusp accuracy.'}, {'name': 'SDF', 'status': 'DERIVED_UNDER_ASSUMPTIONS', 'equation': 'negative-inside0.20mm signed EDT then zero mesh, identical fit to mesh control', 'stop': 'Numerical error remains; no material or microscopic margin claim.'}, {'name': 'clinical_preparation_antagonist', 'status': 'UNKNOWN', 'equation': 'Real abutment-margin/bite transform ports missing', 'stop': 'Cannot transfer intact-collar success to clinical crowns without those observations.'}]})
    (P / 'HANDOFF_R2.md').write_text('# R2 closed\n\n' + json.dumps(parent['decision'], indent=2) + '\n\nNext construction isR3 observed natural cervical collar on4new patients. This changes the physical information/input contract. PREREG_R3.json/SPLIT_R3.json are frozen before import or measurements. Full-deletion results are not overwritten or pooled with partial-crown inputs.\n')
    state('R3_PREREG_FROZEN', parent['decision'], 'measure observed collar and generate all32 predictions before opening hidden upper surfaces')

def collar_fit(mv, mf, obs, of, F):
    c = mv.mean(0)
    lv = (mv - c) @ F
    (_, bi) = boundary(mv, mf)
    cut = float(np.median(lv[bi, 2]) + 1.0)
    sf = mf[np.all(lv[mf, 2] <= cut, axis=1)]
    if len(sf) < 30 or len(of) < 30:
        raise ValueError('INSUFFICIENT_COLLAR_SURFACE')
    source = sampled(mv, sf, 512, 6103)
    target = sampled(obs, of, 512, 6103)
    sl = (source - c) @ F
    tl = (target - c) @ F
    tt = cKDTree(tl)

    def warp(p, z):
        out = z.copy()
        out[:, 0] *= p[4]
        out[:, 1] *= p[5]
        co = np.cos(p[3])
        si = np.sin(p[3])
        xy = out[:, :2].copy()
        out[:, 0] = co * xy[:, 0] - si * xy[:, 1]
        out[:, 1] = si * xy[:, 0] + co * xy[:, 1]
        return out + p[:3]

    def residual(p):
        b = warp(p, sl)
        near = tl[tt.query(b, workers=1)[1]]
        back = b[cKDTree(b).query(tl, workers=1)[1]]
        return np.r_[(b - near).ravel() / 0.3, (tl - back).ravel() / 0.3, 0.35 * p[:3], 0.35 * p[3] / 0.2, 0.35 * (p[4:] - 1) / 0.2]
    p0 = np.r_[np.clip(tl.mean(0) - sl.mean(0), -2.5, 2.5), 0, 1, 1]
    fit = least_squares(residual, p0, bounds=([-3, -3, -3, -0.4, 0.7, 0.7], [3, 3, 3, 0.4, 1.3, 1.3]), loss='soft_l1', f_scale=1, max_nfev=80, diff_step=0.001, ftol=1e-06, xtol=1e-06, gtol=1e-06)
    return (warp(fit.x, lv) @ F.T + c, mf, {'params': fit.x.tolist(), 'cost': float(fit.cost), 'nfev': int(fit.nfev), 'success': bool(fit.success), 'source_collar_faces': len(sf), 'observed_collar_faces': len(of), 'input': 'observed lower collar only'})

def run():
    start = time.perf_counter()
    freeze()
    records = []
    pubrows = []
    hidden = []
    for rec in json.load(open(P / 'SPLIT_R3.json'))['scans']:
        assert sha(rec['mesh']) == rec['mesh_sha256']
        assert sha(rec['labels']) == rec['labels_sha256']
        (v, f) = load_obj(rec['mesh'])
        labels = np.array(json.load(open(rec['labels']))['labels'], dtype=np.int16)
        teeth = tooth_dict(v, f, labels)
        q = 1 if rec['jaw'] == 'upper' else 3
        for i in TARGETS:
            k = q * 10 + i
            key = rec['case'] + '_' + rec['jaw'] + '_' + str(k)
            row = {'key': key, 'case': rec['case'], 'jaw': rec['jaw'], 'fdi': k, 'methods': {}, 'eligible': False}
            if any((z not in teeth or len(teeth[z][1]) < 200 for z in [k, swap(k), k - 1, k + 1])):
                row['failure'] = 'MISSING_TOOTH_OR_NEIGHBOR'
                records.append(row)
                continue
            try:
                (vv, ff, ll) = mask(v, f, labels, k)
                visible = tooth_dict(vv, ff, ll)
                n = axis(visible)
                (rn, rc, pl) = reflection(visible, q)
                F = local_frame(visible, k, n)
                (tv, tf) = teeth[k]
                (_, bi) = boundary(tv, tf)
                cut = float(np.median(tv[bi] @ n) + 1.0)
                keep = (labels != k) | (v @ n <= cut) & (labels == k)
                inds = np.flatnonzero(keep)
                lut = np.full(len(v), -1, dtype=np.int32)
                lut[inds] = np.arange(len(inds))
                fv = f[np.all(keep[f], axis=1)]
                cv = v[inds]
                cf = lut[fv]
                cl = labels[inds]
                op = DATA / 'public' / 'R3' / (key + '.npz')
                op.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(op, vertices=cv, faces=cf, labels=cl, target_fdi=k, observed_collar_max_height_mm=cut, arch_normal=n)
                (obs, of) = select(cv, cf, cl, k)
                assert np.all(obs @ n <= cut + 1e-10)
                row.update(eligible=True, context=str(op), context_sha256=sha(op), observed_collar_vertices=len(obs), cut_height_mm=cut, frame=F.tolist())
                hp = DATA / 'hidden' / 'R3' / (key + '.npz')
                write_surface(hp, tv, tf)
                hidden.append({'key': key, 'reference': str(hp), 'reference_sha256': sha(hp), 'context': str(op)})
                (src, sf) = visible[swap(k)]
                mirror = reflect(src, rn, rc)
                sf = sf[:, ::-1]
                center = mirror.mean(0)
                ts = time.perf_counter()
                (new, nf, fit) = collar_fit(mirror, sf, obs, of, F)
                fit['seconds'] = time.perf_counter() - ts
                ts = time.perf_counter()
                (sv, ss, phi, sd) = sdf_roundtrip(new, nf, F, center)
                sd['seconds'] = time.perf_counter() - ts
                for (name, (pv, pf, meta)) in {'mirror': (mirror, sf, {}), 'collar_mesh': (new, nf, fit), 'collar_sdf': (sv, ss, dict(sd, shared_fit=fit))}.items():
                    path = DATA / 'predictions' / 'R3' / name / (key + '.npz')
                    write_surface(path, pv, pf)
                    row['methods'][name] = {'path': str(path), 'sha256': sha(path), 'details': meta}
            except Exception as e:
                row['failure'] = type(e).__name__ + ': ' + str(e)
            records.append(row)
        print('R3_GENERATE', rec['case'], rec['jaw'], 'seconds', round(time.perf_counter() - start, 1), flush=True)
    put(DATA / 'public' / 'R3_INDEX.json', [{k: v for (k, v) in r.items() if k != 'methods'} for r in records])
    put(DATA / 'hidden' / 'R3_INDEX.json', hidden)
    idx = DATA / 'predictions' / 'R3' / 'INDEX.json'
    put(idx, records)
    fr = {'round': 'R3', 'frozen_utc': now(), 'prereg_sha256': sha(P / 'PREREG_R3.json'), 'split_sha256': sha(P / 'SPLIT_R3.json'), 'index_sha256': sha(idx), 'predictions_index': str(idx), 'prediction_files': [{'key': r['key'], 'method': m, 'path': v['path'], 'sha256': v['sha256']} for r in records for (m, v) in r['methods'].items()], 'before_reference_evaluation': True, 'generation_seconds': time.perf_counter() - start, 'source_code_sha256': sha(pathlib.Path(__file__)), 'peakrss_kb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    put(P / 'FROZEN_PREDICTIONS_R3.json', fr)
    state('R3_PREDICTIONS_FROZEN', {'scheduled': 32, 'completed': sum((bool(r['methods']) for r in records))}, 'evaluate hidden natural upper surfaces; no fit to hidden geometry')
    import crownbench as C
    from exact_distance import exact_distances
    C.distances = exact_distances
    from evaluate import metric, balanced, TYPES, controls
    lookup = {r['key']: r for r in hidden}
    rows = []
    evstart = time.perf_counter()
    probe = None
    with open(P / 'RAW_METRICS_R3.jsonl', 'w', buffering=1) as out:
        for r in records:
            row = {k: r[k] for k in ['key', 'case', 'jaw', 'fdi', 'eligible']}
            row['tooth_type'] = TYPES[r['fdi'] % 10]
            row['metrics'] = {}
            if not r['methods']:
                row['failure'] = r.get('failure')
                rows.append(row)
                out.write(json.dumps(row) + '\n')
                continue
            h = lookup[r['key']]
            assert sha(h['reference']) == h['reference_sha256']
            ref = read_surface(h['reference'])
            with np.load(h['context']) as d:
                cv = d['vertices']
                cf = d['faces']
                cl = d['labels']
            k = r['fdi']
            neigh = [select(cv, cf, cl, k - 1), select(cv, cf, cl, k + 1)]
            F = np.array(r['frame'])
            n = F[:, 2]
            for (m, det) in r['methods'].items():
                assert sha(det['path']) == det['sha256']
                row['metrics'][m] = metric(read_surface(det['path']), ref, neigh, n, r['key'])
            if probe is None:
                probe = (ref, neigh, n, F, r['key'])
            rows.append(row)
            out.write(json.dumps(row) + '\n')
            if len(rows) % 4 == 0:
                print('R3_EVALUATE', len(rows), 'seconds', round(time.perf_counter() - evstart, 1), flush=True)
    methods = ['mirror', 'collar_mesh', 'collar_sdf']
    tables = []
    for typ in ['all', 'incisor', 'canine', 'premolar', 'molar']:
        sub = rows if typ == 'all' else [r for r in rows if r['tooth_type'] == typ]
        for m in methods:
            values = balanced(sub, m, 'occlusal_p95_mm')
            if values:
                tables.append({'tooth_type': typ, 'method': m, 'teeth': sum((m in r['metrics'] for r in sub)), 'patients': len(values), **{q: float(np.mean(list(balanced(sub, m, q).values()))) for q in ('surface_rms_mm', 'occlusal_p95_mm', 'proximal_patch_p95_mm', 'contact_location_mm')}})
    a = balanced(rows, 'collar_sdf', 'occlusal_p95_mm')
    b = balanced(rows, 'mirror', 'occlusal_p95_mm')
    cases = sorted(set(a) & set(b))
    av = np.array([a[c] for c in cases])
    bv = np.array([b[c] for c in cases])
    diff = av - bv
    ii = np.random.default_rng(6103).integers(0, len(cases), (1000, len(cases)))
    ci = np.percentile(diff[ii].mean(1), [2.5, 97.5])
    proxA = balanced(rows, 'collar_sdf', 'proximal_patch_p95_mm')
    proxB = balanced(rows, 'mirror', 'proximal_patch_p95_mm')
    pd = float(np.mean([proxA[c] - proxB[c] for c in cases]))
    control = balanced(rows, 'collar_mesh', 'occlusal_p95_mm')
    cd = float(np.mean([a[c] - control[c] for c in cases]))
    completed = sum(('collar_sdf' in r['metrics'] for r in rows))
    gc = controls(*probe) if probe else None
    gates = {'data': len(cases) >= 3 and completed / 32 >= 0.8, 'primary_relative_and_bootstrap': av.mean() <= 0.9 * bv.mean() and ci[1] < 0, 'absolute_0.50mm': av.mean() <= 0.5, 'proximal_noninferiority_0.05mm': pd <= 0.05, 'numerical': json.load(open(P / 'TRAIN_FIT.json'))['roundtrip']['passes'], 'metric_controls': bool(gc and gc['all_pass'])}
    decision = {'decision': 'POSITIVE_PARTIAL_CROWN_COMPLETION' if all(gates.values()) else 'NEGATIVE', 'gates': {k: 'PASS' if v else 'FAIL' for (k, v) in gates.items()}, 'completed': completed, 'scheduled': 32, 'test_patients': len(cases), 'mean_occlusal_p95_mm': float(av.mean()), 'mirror_mean_occlusal_p95_mm': float(bv.mean()), 'ratio_to_mirror': float(av.mean() / bv.mean()), 'paired_patient_bootstrap95ci_mm': ci.tolist(), 'paired_proximal_patch_difference_mm': pd, 'candidate_minus_equal_information_mesh_mm': cd, 'algorithmic_comparison': 'TIE_WITHIN_0.05mm' if abs(cd) <= 0.05 else 'NO_NOVELTY_CLAIM'}
    put(P / 'METRIC_CONTROLS_R3.json', gc)
    put(P / 'RESULTS_R3.json', {'round': 'R3', 'evaluated_utc': now(), 'external_referent': json.load(open(P / 'PREREG_R3.json'))['external_referent'], 'decision': decision, 'tables': tables, 'frozen_predictions_sha256': sha(P / 'FROZEN_PREDICTIONS_R3.json'), 'raw_metrics_sha256': sha(P / 'RAW_METRICS_R3.jsonl'), 'prereg_sha256': sha(P / 'PREREG_R3.json'), 'generation_seconds': fr['generation_seconds'], 'evaluation_seconds': time.perf_counter() - evstart, 'peakrss_kb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'failures': [{'key': r['key'], 'failure': r.get('failure')} for r in records if r.get('failure')], 'scope': 'Different input from full deletion: natural cervical collar observed. No prepared abutment,registered antagonist or clinical quality. Four new patients: limited population inference.', 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
    state('R3_COMPLETE', decision, 'bind evidence; next true preparation-margin and labeled registered bite inputs')
    print(json.dumps(decision), flush=True)
if __name__ == '__main__':
    run()
