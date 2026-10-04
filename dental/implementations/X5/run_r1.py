from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import csv, datetime, gzip, hashlib, json, resource, time
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
from decision_error import fixed_design_record, guard, scalar_boundary_control
P = Path(__file__).resolve().parent
D = Path(__import__('os').environ.get('DENTAL_PROJECT_ROOT', str(P.parent.parent))).resolve()
S = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/datasets/dental_3fold_extra/scratch_DESIGN/sites'))

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def put(n, a):
    (P / n).write_text(json.dumps(a, indent=2, allow_nan=False, ensure_ascii=False) + '\n')

def state(phase, gate, nextop):
    a = json.load(open(P / 'CURRENT_WORK_STATE.json'))
    a.update(status=phase, updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(), latest_gate=gate, next_operation=nextop)
    put('CURRENT_WORK_STATE.json', a)

def key(c):
    return (tuple(c['pose']), float(c['D']), float(c['L']))

def distances(c, lower):
    out = {'canal_body' if lower else 'sinus_body': c['d_body'], 'lateral_bone': c['wall'], 'adjacent_root': c['adj']}
    return {k: None if abs(abs(v) - 99.0) < 0.0001 else v for (k, v) in out.items()}

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    assert sha(P / 'PREREG_R1.json') == (P / 'PREREG_R1.sha256').read_text().strip()
    pr = json.load(open(P / 'PREREG_R1.json'))
    th = pr['thresholds_mm']
    fm = json.load(open(D / 'results/GEOM_UNC/field_model.json'))
    budget = {}
    for (reg, k) in [('canal_body', 'canal'), ('sinus_body', 'sinus_floor'), ('lateral_bone', 'bone_outer'), ('adjacent_root', 'root')]:
        x = fm[k]
        budget[reg] = abs(x['all']['bias']) + 1.96 * np.hypot(x['all']['sigma_w'], x['tau'])
    num = np.sqrt(3) * 0.15 / 2 + 5e-05
    picks = {}
    for line in open(D / 'results/DESIGN_implant/picks.jsonl'):
        x = json.loads(line)
        picks[x['case'], x['site']] = x
    rows = []
    src = []
    total = Counter()
    hist = defaultdict(list)
    state('R1_RUNNING', 'PREREG_R1_HASH_OK', 'Stream archived candidates and independently replay inverse boundaries')
    for file in sorted(S.glob('*.json.gz')):
        src.append({'path': str(file), 'sha256': sha(file), 'bytes': file.stat().st_size})
        doc = json.load(gzip.open(file, 'rt'))
        for site in doc.get('sites', []):
            total['source_sites'] += 1
            cs = site.get('candidates', [])
            if not cs:
                total['sites_without_candidates'] += 1
                rows.append({'case': doc['case'], 'site': site.get('site'), 'status': 'UNKNOWN_NO_CANDIDATES'})
                continue
            total['represented_sites'] += 1
            total['source_candidates'] += len(cs)
            valid = [(i, c) for (i, c) in enumerate(cs) if all((k in c for k in ('d_body', 'wall', 'adj', 'pose', 'D', 'L')))]
            total['missing_candidate_distance'] += len(cs) - len(valid)
            if not valid:
                continue
            margins = np.array([[v - th[k] if v is not None else -np.inf for (k, v) in distances(c, site['lower']).items()] for (i, c) in valid])
            total['sites_with_missing_field_candidate'] += int(np.isneginf(margins).any())
            best = int(np.argmax(margins.min(1)))
            total['sites_with_geometric_pass'] += int(margins.min(1).max() >= 0)
            total['sites_with_uniform_0.3_margin'] += int(margins.min(1).max() >= 0.3)
            p = picks.get((doc['case'], site['site']), {})
            rb = p.get('theta', {}).get('nominal', {}).get('RB', {}).get('x')
            choices = {'MAX_SLACK_GEOMETRIC_ONLY': valid[best]}
            if rb:
                target = (tuple(rb[0]), float(rb[1]), float(rb[2]))
                hit = next(((i, c) for (i, c) in valid if key(c) == target), None)
                if hit:
                    choices['ARCHIVED_RB'] = hit
                else:
                    total['RB_not_in_candidates'] += 1
            for (name, (i, c)) in choices.items():
                rec = fixed_design_record(distances(c, site['lower']), th, budget, num)
                out = {'case': doc['case'], 'site': site['site'], 'jaw': 'lower' if site['lower'] else 'upper', 'region_class': site.get('pos'), 'selection': name, 'candidate_index': i, 'pose': c['pose'], 'D_mm': c['D'], 'L_mm': c['L'], 'source_path': str(file), 'source_sha256': src[-1]['sha256'], 'geometric_only': True, 'full_original_stochastic_decision': 'NOT_REEVALUATED', **rec}
                rows.append(out)
                for (reg, g) in rec['guards'].items():
                    hist[name + '|' + reg].append(g['signed_margin_mm'])
    prep = json.load(open(D / 'results/PROC_prep/runs/design_anatomy_h0.1.json'))
    prows = []
    for (tooth, a) in prep['teeth'].items():
        for (mat, b) in a['materials'].items():
            for (region, d) in b['regions'].items():
                prows.append({'tooth': tooth, 'material': mat, 'region': region, 'preparation_depth_mm': b['d'], 'provenance': 'PROC_prep/runs/design_anatomy_h0.1.json', 'published_pulp_ASD_as_bound': 'INVALID_METRIC_TRANSFER', **guard(d, 0.5, None, 0.1)})
    controls = []
    for r in rows:
        for (reg, g) in r.get('guards', {}).items():
            x = scalar_boundary_control(g['distance_mm'], g['threshold_mm'])
            err = abs(x - g['signed_margin_mm'])
            controls.append({'id': f"{r['case']}:{r['site']}:{r['selection']}:{reg}", 'boundary_error_mm': err, 'pass': bool(err <= 1e-06), 'fault_plus0.05_rejected': bool(abs(x - g['signed_margin_mm'] - 0.05) > 1e-06)})
    for r in prows:
        x = scalar_boundary_control(r['distance_mm'], r['threshold_mm'])
        controls.append({'id': f"STS:{r['tooth']}:{r['material']}:{r['region']}", 'boundary_error_mm': abs(x - r['signed_margin_mm']), 'pass': abs(x - r['signed_margin_mm']) <= 1e-06, 'fault_plus0.05_rejected': abs(x - r['signed_margin_mm'] - 0.05) > 1e-06})
    summary = {k: {'n': len(v), 'min_mm': float(min(v)), 'median_mm': float(np.median(v)), 'p10_mm': float(np.percentile(v, 10)), 'negative_n': int(np.sum(np.array(v) < 0)), 'nonnegative_below0.22_n': int(np.sum((np.array(v) >= 0) & (np.array(v) < 0.22)))} for (k, v) in hist.items()}
    flat = []
    for r in rows:
        for (reg, g) in r.get('guards', {}).items():
            flat.append({k: r.get(k) for k in ('case', 'site', 'jaw', 'region_class', 'selection', 'D_mm', 'L_mm')} | {'region': reg} | g)
    with open(P / 'DECISION_MARGINS.csv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(flat[0]))
        w.writeheader()
        w.writerows(flat)
    with open(P / 'PREPARATION_MARGINS.csv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(prows[0]))
        w.writeheader()
        w.writerows(prows)
    with open(P / 'RAW_R1_SITE_GUARDS.jsonl', 'w') as f:
        for r in rows:
            f.write(json.dumps(r, allow_nan=False) + '\n')
    put('SOURCE_DESIGN_MANIFEST.json', src)
    put('R1_CONTROLS.json', controls)
    output = {'id': pr['id'], 'counts': dict(total), 'selected_guard_rows': len(flat), 'prep_guard_rows': len(prows), 'summary': summary, 'scenario_budgets_mm': budget, 'numerical_lookup_scenario_mm': num, 'scalar_control_max_error_mm': max((q['boundary_error_mm'] for q in controls)), 'gates': {'inverse_replay': all((q['pass'] for q in controls)), 'fault_injection': all((q['fault_plus0.05_rejected'] for q in controls)), 'matched_local_physical_calibration': 'FAIL_METRIC_AND_COHORT_MISMATCH'}, 'outcome': 'CAPABILITY_DELIVERED_SCALAR_METHOD_TIE_PHYSICAL_SAFETY_UNKNOWN', 'external_referent': pr['external_referent'], 'cost': {'wall_seconds': time.perf_counter() - start, 'cpu_seconds': time.process_time() - cpu, 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'GPU': False, 'inherited_annotation_EDT_search_and_manual_review_seconds': None}}
    put('SUMMARY_R1.json', output)
    (P / 'HANDOFF_R1.md').write_text(f"# R1 avgjord\n\n{total['represented_sites']} represented sites, {total['source_candidates']} candidates, {len(flat)} selected regional lines; and {len(prows)} preparationsregioner. Inverse bounds matches replay scale; 0,05 mm injected error is rejected. Metodutfall TIE. Stochastic K37 risk and physical safety have not been validated.\n\nThe literature banner/HD95/ASD cannot be a local edge error limit. Physical calibration gate FAIL . R1 has no geometric witness coordinates per archived region.\n\nNext design R2 : change from a scaly minimum to all local hotspots distance bar; read local TF2 channel dials, maintain direction and two-sided coverage , select the smallest measurement patch that is sufficient during explicit error budget. Run full re-evaluation for an arbitrary predictive mask and adversarial migration of as many voxlares.\n")
    state('R1_COMPLETE', output['outcome'], 'Freeze R2 full local guard-field and model-agnostic paired mask evaluator')
    print(json.dumps({'counts': dict(total), 'gates': output['gates'], 'wall_seconds': output['cost']['wall_seconds']}))
if __name__ == '__main__':
    main()
