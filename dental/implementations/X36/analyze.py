from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, csv, hashlib, re, time, datetime, resource, xml.etree.ElementTree as ET
import numpy as np
from meta_engine import design, sampling, fit, egger
P = Path(__file__).resolve().parent
XML = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))

def read(name):
    return json.loads((P / name).read_text())

def dump(name, obj):
    (P / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def text(e):
    return ' '.join(' '.join(e.itertext()).split())

def write_csv(name, rows):
    keys = list(dict.fromkeys((k for r in rows for k in r)))
    with (P / name).open('w') as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(rows)

def grid(tb):
    cells = {}
    rows = []
    for (i, tr) in enumerate(tb.findall('.//tr')):
        col = 0
        for e in tr:
            while (i, col) in cells:
                col += 1
            for a in range(int(e.get('rowspan', '1'))):
                for b in range(int(e.get('colspan', '1'))):
                    cells[i + a, col + b] = text(e)
            col += int(e.get('colspan', '1'))
        rows.append([cells[i, c] for c in range(max((j for (k, j) in cells if k == i)) + 1)])
    return rows

def num(s):
    return float(re.search('\\d+(?:\\.\\d+)?', s).group())

def pairs(s):
    return [(float(a.replace(',', '')), float(b.replace(',', ''))) for (a, b) in re.findall('([\\d,]+(?:\\.\\d+)?)\\s*±\\s*([\\d,]+(?:\\.\\d+)?)', s)]

def validate_unit(row):
    expected = {'cement': 'um', 'crown': 'N'}[row['dataset']]
    if row['unit'] != expected:
        raise ValueError(f"Unit mismatch: {row['row_id']}, {row['unit']} vs {expected}")
    return True

def prepare():
    rows = []
    audit = []
    sources = {}
    gap = list(csv.DictReader((P / 'inputs/cement_measurements.csv').open()))
    srcs = {s['study']: s for s in read('inputs/cement_source_manifest.json')['sources']}
    for r in gap:
        s = r['study']
        path = Path(r['source_path'])
        if s not in sources:
            root = ET.parse(path).getroot()
            sources[s] = dict(study=s, source=str(path), sha256=sha(path), permissions=text(root.find('.//permissions')))
            sources[s]['tables'] = {t.get('id'): grid(t) for t in root.findall('.//table-wrap')}
        assert sources[s]['sha256'] == r['source_sha256'], s
        tb = sources[s]['tables'][r['source_table']]
        (i, j) = (int(r['source_row']), int(r['source_column']))
        cell = tb[i][j]
        (mu, sd) = (float(r['measured_mean_um']), float(r['reported_sd_um']))
        pp = pairs(cell)
        if pp:
            (em, es) = pp[-1]
        else:
            em = num(cell)
            if s == 'PMC10721348':
                z = next((a for a in tb[i + 1:] if a[0] == 'SD' or 'standard deviation' in a[0].lower()))
                es = num(z[1])
            else:
                es = num(tb[i][j + 1])
        ok = abs(em - mu) <= 0.011 and abs(es - sd) <= 0.011
        assert ok, (r['row_id'], em, es, mu, sd)
        reason = None
        if r['region'] not in ['marginal', 'axial', 'occlusal']:
            reason = 'different_region_not_primary'
        elif r.get('sd_basis') == 'measurement_locations':
            reason = 'location_SD_not_specimen_SD'
        elif s == 'PMC10582242' and r['region'] == 'occlusal' and (r['arm'] == 'layer_100'):
            reason = 'inconsistent_mean_SD_minmax'
        cad = r['cad'].lower()
        family = '3Shape' if '3shape' in cad else 'Exocad' if 'exocad' in cad else 'CEREC' if 'cerec' in cad else 'Other'
        n = int(r['n_specimens'])
        sp = float(r['internal_spacer_um'])
        rows.append(dict(**r, dataset='cement', outcome=mu, variance=sd ** 2 / n, variance_no_n=sd ** 2, cluster=r['source_family'], CAD_family=family, cemented=int(r['state'] == 'cemented'), CAD_um=sp, eligible=reason is None, exclusion_reason=reason or '', unit='um', resolution_level='PER_SURFACE_REGION', time_scale='SIMULTANEOUS', sampling_precision_status='CONDITIONAL_SD_BASIS_NOT_FULLY_VERIFIED', locator=f"https://doi.org/{r['doi']}; {s} #{r['source_table']}, expanded row {i}, col {j}", original_cell_verified=ok))
        audit.append(dict(row_id=r['row_id'], mean_error=abs(em - mu), SD_error=abs(es - sd), original_sha_pass=True))
    for r in rows:
        group = [a['CAD_um'] for a in rows if a['study'] == r['study'] and a['eligible']]
        center = float(np.mean(group)) if group else r['CAD_um']
        r['CAD_within_per10um'] = (r['CAD_um'] - center) / 10
        r['CAD_between_per10um'] = (center - 50) / 10
    crown = []
    for r in read('inputs/crown_groups.json'):
        s = r['study']
        path = P.parent / _release_expand('X1B') / r['source']
        assert sha(path) == r['source_sha256'], s
        primary = XML / (s + '.xml')
        if primary.exists():
            root = ET.parse(primary).getroot()
            perm = text(root.find('.//permissions'))
            allpairs = pairs(' '.join((text(t) for t in root.findall('.//table-wrap'))))
            if s == 'PMC10560138':
                paragraph = next((p for p in root.findall('.//body//p') if text(p).startswith('Descriptive statistics for fracture')))
                allpairs = pairs(text(paragraph))
        else:
            root = None
            perm = 'UNKNOWN: inspect source permission before redistribution'
            allpairs = pairs(path.read_text())
        (mu, sd) = (r['mean_N'], r['sd_N'])
        origpair = any((abs(a - mu) < 0.011 and abs(b - sd) < 0.011 for (a, b) in allpairs))
        if s == 'PMC4764450':
            origpair = any((abs(a * 1000 - mu) < 0.011 and abs(b * 1000 - sd) < 0.011 for (a, b) in allpairs))
        if s == 'PMC7274823':
            content = text(root) if root is not None else path.read_text()
            pp = [(float(a.replace(',', '')), float(b)) for (a, b) in re.findall('([\\d,]+(?:\\.\\d+)?)\\s*\\(([\\d.]+)\\)', content)]
            origpair = any((abs(a - mu) < 0.011 and abs(b - sd) < 0.011 for (a, b) in pp))
        if s not in sources:
            sources[s] = dict(study=s, source=str(primary if primary.exists() else path), sha256=sha(primary if primary.exists() else path), permissions=perm)
        mat = r['material']
        if mat == 'UNKNOWN' and s == 'PMC7274823':
            mat = 'zirconia_unspecified'
        reason = 'different_endpoint_critical_splitting' if s == 'PMC4764450' else ''
        if not origpair:
            reason = 'inherited_numeric_pair_not_in_original_table'
        n = r['n']
        v = sd ** 2 / (n * mu ** 2)
        if s == 'PMC10478297':
            mat = 'zirconia_reported5Y_product_conflict'
        crown.append(dict(**r, dataset='crown', outcome=float(np.log(mu)), variance=v, variance_no_n=v * n, cluster=s, arm=r['row_id'], material_class=mat, angle30=int(r['angle_deg'] == 30), log_thickness=float(np.log(r['thickness_mm'])) if r['thickness_mm'] is not None else None, eligible=not reason, thickness_eligible=not reason and r['thickness_mm'] is not None, endpoint_exclusion=reason, resolution_level='PER_TOOTH', time_scale='SIMULTANEOUS', local_primary_path=sources[s]['source'], local_primary_sha256=sources[s]['sha256'], license=perm, original_pair_verified=origpair))
        audit.append(dict(row_id=r['row_id'], original_pair_found=origpair, original_sha_pass=True, excluded_from_R1=not origpair))
    for r in rows + crown:
        validate_unit(r)
    dump('SOURCE_AUDIT.json', dict(rows=len(audit), all_pass=all((a.get('original_pair_found', True) for a in audit)), checks=audit, note='Single-producer replay of original table cells; no independent extraction review; crown check confirms pair presence in primary tables, locator placement still needs second reviewer. Nonmatching rows excluded, never corrected in R1.'))
    for s in sources.values():
        s.pop('tables', None)
    dump('SOURCE_PERMISSIONS.json', list(sources.values()))
    write_csv('EXTRACTION_CEMENT.csv', rows)
    write_csv('EXTRACTION_CROWN.csv', crown)
    dump('NORMALIZED_DATA.json', dict(cement=rows, crown=crown))
    return (rows, crown)
GBASE = [('CAD_within_per10um', 'numeric', None), ('CAD_between_per10um', 'numeric', None), ('region', 'categorical', 'marginal')]
CBASE = [('log_thickness', 'numeric', None), ('material_class', 'categorical', '3Y'), ('angle30', 'numeric', None)]

def model(rows, terms, rho=0, no_n=False, profile=True, control=False):
    (X, names) = design(rows, terms)
    out = fit([r['outcome'] for r in rows], X, sampling(rows, rho, no_n), [r['cluster'] for r in rows], names, profile, control)
    out['row_ids'] = [r['row_id'] for r in rows]
    out['terms'] = terms
    out['rho_sampling'] = rho
    out['no_n_division'] = no_n
    out['resolution_level'] = 'POPULATION'
    out['physical_quantity_resolution'] = rows[0]['resolution_level']
    return out

def predict_loso(rows, terms, kind):
    clusters = sorted(set((r['cluster'] for r in rows)))
    records = []
    for held in clusters:
        tr = [r for r in rows if r['cluster'] != held]
        te = [r for r in rows if r['cluster'] == held]
        for r in te:
            reason = []
            for (key, cat, _) in terms:
                if cat == 'categorical' and r[key] not in {a[key] for a in tr}:
                    reason.append('unseen_category:' + key + '=' + r[key])
            if reason:
                records.append(dict(row_id=r['row_id'], held_cluster=held, status='UNKNOWN_UNSEEN_CATEGORY', reason=reason))
                continue
        known = [r for r in te if not any((a['row_id'] == r['row_id'] for a in records))]
        if not known:
            continue
        (X, names) = design(tr, terms)
        f = fit([r['outcome'] for r in tr], X, sampling(tr), [r['cluster'] for r in tr], names, False)
        if 'beta' not in f:
            records.extend((dict(row_id=r['row_id'], held_cluster=held, status=f['status']) for r in known))
            continue
        b = np.array(f['beta'])
        C = np.array(f['coefficient_covariance'])
        for r in known:
            x = np.array([1.0] + [r[k] if '=' not in k else float(r[k.split('=')[0]] == k.split('=')[1]) for k in names[1:]])
            pred = float(x @ b)
            var = float(x @ C @ x + f['tau2_study'] + f['omega2_cell'] + r['variance'])
            if kind == 'cement':
                base = r['CAD_um']
            else:
                bystudy = [np.mean([a['outcome'] - 2 * a['log_thickness'] for a in tr if a['cluster'] == s]) for s in sorted({a['cluster'] for a in tr})]
                base = float(np.mean(bystudy) + 2 * r['log_thickness'])
            (lo, hi) = (pred - 1.95996398454 * np.sqrt(var), pred + 1.95996398454 * np.sqrt(var))
            records.append(dict(row_id=r['row_id'], held_cluster=held, status='SCORED', observed=r['outcome'], prediction=pred, practice=base, pi95=[float(lo), float(hi)], covered=bool(lo <= r['outcome'] <= hi), doi=r['doi'], external_referent_kind='independent_measurement', resolution_level=r['resolution_level']))
    scored = [r for r in records if r['status'] == 'SCORED']
    summary = dict(eligible_rows=len(rows), scored_rows=len(scored), unknown_rows=len(rows) - len(scored), limitation='Conditional normal intervals and sparse source-family training; retrospective holdout, not independent extraction or prospective laboratory validation')
    if scored:
        sgroups = sorted({r['held_cluster'] for r in scored})
        mse = np.mean([np.mean([(r['prediction'] - r['observed']) ** 2 for r in scored if r['held_cluster'] == s]) for s in sgroups])
        bmse = np.mean([np.mean([(r['practice'] - r['observed']) ** 2 for r in scored if r['held_cluster'] == s]) for s in sgroups])
        coverage = np.mean([np.mean([r['covered'] for r in scored if r['held_cluster'] == s]) for s in sgroups])
        summary.update(study_balanced_RMSE=float(np.sqrt(mse)), practice_study_balanced_RMSE=float(np.sqrt(bmse)), RMSE_ratio=float(np.sqrt(mse / bmse)), study_balanced_coverage=float(coverage), gate_pass=bool(mse <= 0.8 ** 2 * bmse and coverage >= 0.9), scored_clusters=len(sgroups))
    return dict(summary=summary, rows=records)

def funnels(gap, crown):
    out = {}
    for (label, rows) in [('cement_marginal_discrepancy', [r for r in gap if r['eligible'] and r['region'] == 'marginal']), ('crown_log_force', [r for r in crown if r['eligible']])]:
        records = []
        for c in sorted({r['cluster'] for r in rows}):
            a = [r for r in rows if r['cluster'] == c]
            vals = np.array([r['outcome'] - (r['CAD_um'] if label.startswith('cement') else 0) for r in a])
            se = np.mean(np.sqrt([r['variance'] for r in a]))
            records.append(dict(cluster=c, value=float(np.mean(vals)), se_upper_correlation1=float(se), n_rows=len(a), dois=sorted({r['doi'] for r in a})))
        out[label] = dict(study_summaries=records, egger=egger([r['value'] for r in records], [r['se_upper_correlation1'] for r in records]), resolution_level='POPULATION', summary_method='Equal-weight study-family mean; SE perfect-positive-correlation bound, not empirically identified SE; descriptive only')
    return out

def flow(gap, crown):
    selected_c = set((r['study'] for r in crown))
    selected_g = set((r['study'] for r in gap))
    out = {}
    screen = []
    for (kind, fn, selected, rows, key) in [('cement', 'cement_corpus.jsonl', selected_g, gap, 'eligible'), ('crown', 'crown_corpus.jsonl', selected_c, crown, 'thickness_eligible')]:
        corpus = [json.loads(x) for x in (P / 'inputs' / fn).read_text().splitlines()]
        ids = [x['pmcid'] for x in corpus]
        unique = set(ids)
        included = set((r['study'] for r in rows if r[key]))
        additional = selected - unique
        out[kind] = dict(local_records=len(corpus), additional_local_primary_reports=len(additional), additional_report_ids=sorted(additional), total_identified_local_records=len(corpus) + len(additional), unique_records=len(unique) + len(additional), duplicate_records=len(corpus) - len(unique), inherited_extracted_reports=len(selected), not_extracted_or_eligibility_not_assessed=len(unique - selected), extracted_rows=len(rows), primary_rows=sum((r[key] for r in rows)), primary_reports=len(included), reports_not_primary=len(selected - included), rejected_row_fraction=float(1 - sum((r[key] for r in rows)) / len(rows)), original_database_search_total=None, fulltext_retrieval_failures=None, dual_screening=False, dual_extraction=False, scope='Reconstructed flow for selected local corpus secondary synthesis, NOT a completed systematic review. Unassessed reports are not exclusions.')
        for x in corpus:
            screen.append(dict(dataset=kind, pmcid=x['pmcid'], doi=x.get('doi'), title=x.get('title'), selection_status='INHERITED_EXTRACTION' if x['pmcid'] in selected else 'UNASSESSED_NOT_EXTRACTED', eligibility_status='PRIMARY_INCLUDED' if x['pmcid'] in included else 'UNASSESSED_OR_NOT_PRIMARY'))
        for s in sorted(additional):
            a = next((r for r in rows if r['study'] == s))
            screen.append(dict(dataset=kind, pmcid=s, doi=a['doi'], title='Additional inherited local primary report; see source', selection_status='INHERITED_EXTRACTION_ADDITIONAL_LOCAL_SOURCE', eligibility_status='PRIMARY_INCLUDED' if s in included else 'NOT_PRIMARY'))
    write_csv('SCREENING_LEDGER.csv', screen)
    return out

def main():
    start = time.perf_counter()
    (gap, crown) = prepare()
    prep = time.perf_counter() - start
    g = [r for r in gap if r['eligible']]
    c = [r for r in crown if r['thickness_eligible']]
    ca = [r for r in crown if r['eligible']]
    models = {}
    specs = [('cement_full', g, GBASE + [('material', 'categorical', 'zirconia'), ('CAD_family', 'categorical', '3Shape'), ('cemented', 'numeric', None)]), ('cement_reduced', g, GBASE), ('cement_material', g, GBASE + [('material', 'categorical', 'zirconia')]), ('cement_system', g, GBASE + [('CAD_family', 'categorical', '3Shape')]), ('crown_full', c, CBASE), ('crown_all_material_setup', ca, CBASE[1:])]
    for (name, rows, terms) in specs:
        t = time.perf_counter()
        models[name] = model(rows, terms, control=True)
        models[name]['seconds'] = time.perf_counter() - t
        dump('MODEL_' + name + '.json', models[name])
        print(name, models[name]['status'], models[name]['support'], flush=True)
    sensitivities = {}
    for rho in [0.5, 0.9]:
        sensitivities[f'cement_rho_{rho}'] = model(g, GBASE, rho=rho)
    sensitivities['cement_no_n_division'] = model(g, GBASE, no_n=True)
    sensitivities['cement_drop_possible_shared_paper'] = model([r for r in g if r['study'] != 'PMC10971874'], GBASE)
    dump('SENSITIVITY_R1.json', sensitivities)
    loso = {'cement': predict_loso(g, GBASE, 'cement'), 'crown': predict_loso(c, CBASE, 'crown')}
    dump('LOSO_R1.json', loso)
    funnel = funnels(gap, crown)
    dump('FUNNEL_DATA.json', funnel)
    fl = flow(gap, crown)
    dump('PRISMA_FLOW.json', fl)
    outcome = dict(claim_type='information_link', round='R1', review_state='PENDING_INDEPENDENT_REVIEW', verdict='TRANSFERABLE_FULL_ADJUSTED_INFERENCE_NOT_SUPPORTED', models=models, sensitivity=sensitivities, holdout=loso, funnel=funnel, flow=fl, external_referent=dict(kind='independent_measurement', locator='EXTRACTION_CEMENT.csv; EXTRACTION_CROWN.csv DOI/table/cell; LOSO_R1.json held primary DOI', compared_quantity='Published group-mean regional gap (um) and crown fracture force (N)', refutes_us=True), cost=dict(preparation_seconds=prep, elapsed_seconds=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, threads=4, prior_collection_and_extraction_cost='UNKNOWN', agent_reading_and_discovery_cost='UNKNOWN', physical_measurements=0, questions=0), numeric_levels={'gap': 'PER_SURFACE_REGION', 'fracture_force': 'PER_TOOTH', 'coefficients_and_heterogeneity': 'POPULATION', 'unmeasured_sampling_correlation': 'PHENOMENOLOGICAL'}, debts=['Specimen-level SD definition/covariance unresolved; sensitivity is not validation', 'Unassessed corpus records and lack of independent dual extraction preclude systematic review claim', 'No causal loading-angle, support, or full-region spacer effect', 'No new clinical risk or pointwise film prediction'])
    dump('RESULTS_R1.json', outcome)
    dump('CURRENT_WORK_STATE.json', dict(lane='X36-meta-regression', phase='R1_COMPLETE', latest_gate=outcome['verdict'], next_operation='Freeze within-study contrast construction R2; preserve R1', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    (P / 'HANDOFF_R1.md').write_text('R1 completed: provenance audited and random-study-intercept models fitted where full rank. Full inferential claim fails independent-study support gates. See RESULTS_R1.json for unchanged criteria and all numerical outcomes. Next: within-study endpoint contrasts that cancel study offsets; frozen separately before computing.\n')
    print('Holdout', json.dumps({k: v['summary'] for (k, v) in loso.items()}))
if __name__ == '__main__':
    main()
