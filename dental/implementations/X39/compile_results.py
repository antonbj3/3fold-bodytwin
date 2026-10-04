"""Compile existing evidence. No refit, inference search, or new bootstrap."""
from dental_release.paths import expand as _release_expand
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '4'
import copy
import csv
import datetime
import hashlib
import json
import math
import resource
import sys
import time
import zipfile
from collections import Counter
from pathlib import Path
P = Path(__file__).resolve().parent
R = P.parent
FIELDS = ['overbite', 'overjet', 'crossbite', 'midlines', 'spee', 'molar_right', 'molar_left', 'canine_right', 'canine_left']
NAMES = dict(zip(FIELDS, ['Overbite', 'Overjet', 'Crossbite', 'Midlines', 'Curve of Spee', 'Molar right', 'Molar left', 'Canine right', 'Canine left']))
TOL = 1e-12
Z = 1.959963984540054
CHECKS = []

def read(lane, name):
    return json.loads((R / lane / name).read_text())

def write(name, obj):
    (P / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def checked(name, valid, mutated):
    CHECKS.append(dict(name=name, passes=bool(valid), injected_wrong_value_rejected=not bool(mutated)))
    if not valid or mutated:
        raise AssertionError(name)

def parity(name, observed, expected):
    predicate = lambda value: abs(value - expected) <= TOL
    checked(name, predicate(observed), predicate(observed + 0.05))

def wilson(k, n):
    if not 0 <= k <= n or n < 1:
        raise ValueError('Invalid conflict denominator')
    p = k / n
    den = 1 + Z * Z / n
    centre = (p + Z * Z / (2 * n)) / den
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / den
    return [centre - half, centre + half]

def verify_wilson(k, n, interval):
    (a, b, c) = (n + Z * Z, -(2 * k + Z * Z), k * k / n)
    discriminant = b * b - 4 * a * c
    expected = [(-b - math.sqrt(discriminant)) / (2 * a), (-b + math.sqrt(discriminant)) / (2 * a)]
    pred = lambda pair: len(pair) == 2 and all((abs(u - v) <= TOL for (u, v) in zip(pair, expected)))
    checked(f'Wilson {k}/{n}', pred(interval), pred([interval[0] + 0.02, interval[1]]))

def table(name, title, columns, rows, note):
    with (P / 'tables' / (name + '.csv')).open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    text = f'### {title}\n\n| ' + ' | '.join((label for (key, label, fmt) in columns)) + ' |\n'
    text += '| ' + ' | '.join(('---' for _ in columns)) + ' |\n'
    for row in rows:
        text += '| ' + ' | '.join((fmt(row[key]) for (key, label, fmt) in columns)) + ' |\n'
    text += '\n' + note + '\n'
    (P / 'tables' / (name + '.md')).write_text(text)

def interval_text(x):
    return f'{x[0]:.3f}–{x[1]:.3f}'

def outcome_summary(rows, ref, pred):
    n = len(rows)
    hit = sum((r[ref] == r[pred] for r in rows))
    classes = sorted(set((r[ref] for r in rows)))
    ba = sum((sum((r[ref] == c and r[pred] == c for r in rows)) / sum((r[ref] == c for r in rows)) for c in classes)) / len(classes)
    p_o = hit / n
    (actual, predicted) = (Counter((r[ref] for r in rows)), Counter((r[pred] for r in rows)))
    p_e = sum((actual[c] * predicted[c] for c in set(actual) | set(predicted))) / (n * n)
    labels = sorted(set(actual) | set(predicted))
    matrix = [[sum((r[ref] == a and r[pred] == b for r in rows)) for b in labels] for a in labels]
    return dict(n=n, hits=hit, accuracy=p_o, balanced_accuracy=ba, kappa=(p_o - p_e) / (1 - p_e), confusion_labels=labels, confusion=matrix)

def main():
    start = time.perf_counter()
    write('CURRENT_WORK_STATE.json', dict(lane='X39-manuscript-bite', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='COMPILING_LOCKED_EVIDENCE', latest_gate='PREREG_REPORTING frozen', next_operation='Verify estimands; regenerate tables and figures'))
    contract = json.loads((P / 'PREREG_REPORTING.json').read_text())
    expected = (P / 'PREREG_REPORTING.json.sha256').read_text().strip()
    checked('Frozen reporting contract', sha(P / 'PREREG_REPORTING.json') == expected, '0' + expected[1:] == expected)
    original_sources = json.loads((P / 'INPUT_LOCK.json').read_text())['sources']
    additional_sources = json.loads((P / 'REPORTING_SOURCE_ADDENDUM.json').read_text())['sources']
    for source in original_sources + additional_sources:
        observed = sha(Path(source['path']))
        corrupted = ('1' if observed[0] != '1' else '2') + observed[1:]
        checked('Source ' + source['path'], observed == source['sha256'], corrupted == source['sha256'])
    manifest = read(_release_expand('X7'), 'raw/DATA_MANIFEST.json')
    partitions = [set(manifest[k]) for k in ['train', 'calibration', 'test']]
    disjoint = lambda sets: not any((a & b for (i, a) in enumerate(sets) for b in sets[i + 1:]))
    bad = copy.deepcopy(partitions)
    bad[2].add(next(iter(bad[0])))
    checked('Patient-case split disjointness', disjoint(partitions), disjoint(bad))
    stats = json.loads((R / 'LANE_X35_PACKAGE_STATISTICS/RESULTS_R1_V2.json').read_text())
    intervals = {(m['demo'], m['metric']): m for m in stats['metrics']}
    baseline = read(_release_expand('X7'), 'raw/RESULTS_R1.json')
    predictions = {r['case_id']: r for r in read(_release_expand('X7'), 'raw/PREDICTIONS_R1.json') if r['split'] == 'test'}
    labels = read(_release_expand('X7'), 'raw/TEST_LABELS.json')
    (base_rows, confusion_rows) = ([], [])
    for field in FIELDS:
        rows = [dict(case_id=c, y=labels[c][0]['labels'][field], pred=predictions[c]['fields'][field]['point']) for c in manifest['test'] if labels[c] and labels[c][0]['labels'].get(field) is not None]
        summary = outcome_summary(rows, 'y', 'pred')
        source = baseline['metrics'][field]['point']
        for metric in ['n', 'accuracy', 'balanced_accuracy', 'kappa']:
            parity('X7 ' + field + ' ' + metric, summary[metric], source[metric])
        for metric in ['accuracy', 'kappa']:
            parity('X35 X7 ' + field + ' ' + metric, summary[metric], intervals['X7', field + ' ' + metric]['point'])
        base_rows.append(dict(field=NAMES[field], n=summary['n'], hits=summary['hits'], accuracy=summary['accuracy'], ci95=intervals['X7', field + ' accuracy']['cluster_ci95'], balanced_accuracy=source['balanced_accuracy'], kappa=source['kappa'], kappa_ci95=intervals['X7', field + ' kappa']['cluster_ci95'], excluded_from_198=198 - summary['n'], excluded_fraction=(198 - summary['n']) / 198, resolution_level='PER_ARCH -> POPULATION', reference_version='X7 R1 original first report'))
        confusion_rows.append(dict(model='Arch-profile X7 R1', field=field, **summary))
    table('table2_arch_profile', 'Table 2. Arch-profile category agreement, original reference', [('field', 'Field', str), ('n', 'n', str), ('hits', 'Correct', str), ('accuracy', 'Agreement', lambda x: f'{x:.3f}'), ('ci95', '95% patient-case bootstrap CI', interval_text), ('balanced_accuracy', 'Balanced accuracy', lambda x: f'{x:.3f}'), ('kappa', 'κ', lambda x: f'{x:.3f}'), ('kappa_ci95', 'κ 95% CI', interval_text)], base_rows, 'Reference: category extracted from the first archive-sorted report, provided that report contains a readable category. Crossbite uses the historical parser, superseded in Table 3B. Fixed-model descriptive intervals are copied from X35 (10,000 whole-patient-case bootstrap draws). Outcomes are PER_ARCH; rates are POPULATION summaries conditional on this internal sample. Field exclusions are missing/unreadable first-report categories, not failed scans.')
    opposition = read(_release_expand('X21'), 'raw/RESULTS_R2.json')
    comparison = read(_release_expand('X21'), 'raw/COMPARISON_ROWS_R2.json')
    main_rows = []
    for field in ['overbite', 'overjet', 'molar_right', 'molar_left']:
        rr = [r for r in comparison if r['field'] == field]
        summary = outcome_summary(rr, 'reference_first', 'prediction')
        ctrl = outcome_summary(rr, 'reference_first', 'x7_prediction')
        for metric in ['n', 'accuracy', 'balanced_accuracy']:
            parity('X21 ' + field + ' ' + metric, summary[metric], opposition['fields'][field]['first_report'][metric])
        ci = intervals['X21', field + ' first-report accuracy']
        contrast = intervals['X21', field + ' paired accuracy gain vs X7']
        delta = summary['accuracy'] - ctrl['accuracy']
        parity('X35 X21 ' + field + ' mean', summary['accuracy'], ci['point'])
        parity('X35 X21 ' + field + ' paired contrast', delta, contrast['point'])
        main_rows.append(dict(field=NAMES[field], n=summary['n'], hits=summary['hits'], accuracy=summary['accuracy'], ci95=ci['cluster_ci95'], comparator_hits=ctrl['hits'], comparator_accuracy=ctrl['accuracy'], delta=delta, delta_ci95=contrast['cluster_ci95'], balanced_accuracy=summary['balanced_accuracy'], majority_accuracy=opposition['fields'][field]['practice_majority']['accuracy'], excluded_from_198=198 - summary['n'], excluded_fraction=(198 - summary['n']) / 198, resolution_level='PER_ARCH -> POPULATION', reference_version='X21 R2 first readable report'))
        confusion_rows.append(dict(model='Contact plus opposition X21 R2', field=field, **summary))
    table('table3a_opposition', 'Table 3A. Contact plus opposition, shared internal test cases', [('field', 'Field', str), ('n', 'n', str), ('hits', 'Correct', str), ('accuracy', 'Agreement', lambda x: f'{x:.3f}'), ('ci95', '95% patient-case bootstrap CI', interval_text), ('comparator_accuracy', 'X7 in same subset', lambda x: f'{x:.3f}'), ('delta', 'Paired difference', lambda x: f'{x:+.3f}'), ('delta_ci95', '95% paired bootstrap CI', interval_text), ('balanced_accuracy', 'Balanced accuracy', lambda x: f'{x:.3f}')], main_rows, 'Reference: first readable report per field, an explicit change from Table 2. X7 predictions are rescored in each identical subset. All models are fixed; no new training or bootstrap. Intervals copied from X35. These are descriptive contrasts after repeated testing of this cohort, not confirmatory superiority tests. The 0.078 overbite difference uses 192 shared cases, not 191 plus 192 independent patients.')
    corrected = read(_release_expand('X21'), 'raw/RESULTS_R5.json')
    rr = read(_release_expand('X21'), 'raw/COMPARISON_ROWS_R5.json')
    corrected_rows = []
    for (model, key, source_key) in [('Corrected-label model', 'prediction', 'corrected_label_model'), ('Frozen opposition model', 'frozen_R2', 'frozen_R2_model'), ('Frozen arch-profile model', 'frozen_X7', 'frozen_X7_model')]:
        rows = [dict(y=r['reference_first'], pred=r[key]['point'] if key == 'prediction' else r[key]) for r in rr]
        summary = outcome_summary(rows, 'y', 'pred')
        for metric in ['n', 'accuracy', 'balanced_accuracy']:
            parity('Corrected crossbite ' + model + ' ' + metric, summary[metric], corrected['scores'][source_key][metric])
        ci = wilson(summary['hits'], summary['n'])
        verify_wilson(summary['hits'], summary['n'], ci)
        corrected_rows.append(dict(model=model, n=summary['n'], hits=summary['hits'], accuracy=summary['accuracy'], ci95=ci, balanced_accuracy=summary['balanced_accuracy'], resolution_level='PER_ARCH -> POPULATION', reference_version='X21 R5 corrected global category'))
        confusion_rows.append(dict(model=model, field='crossbite_corrected', **summary))
    table('table3b_corrected_crossbite', 'Table 3B. Crossbite against the corrected textual reference', [('model', 'Model', str), ('n', 'n', str), ('hits', 'Correct', str), ('accuracy', 'Agreement', lambda x: f'{x:.3f}'), ('ci95', '95% Wilson CI', interval_text), ('balanced_accuracy', 'Balanced accuracy', lambda x: f'{x:.3f}')], corrected_rows, 'No refit of the two frozen comparators. R5 corrects global presence/absence; R6 separately preserves anterior/lateral scope. These category outcomes do not validate named teeth. R5 retains its failed global mixed-scope injection gate. One event per patient-case; Wilson intervals condition on the extracted reference and fixed predictions.')
    conflicts = []
    for field in ['overbite', 'overjet', 'molar_right', 'molar_left']:
        source = opposition['fields'][field]
        conflicts.append(dict(field=NAMES[field], population='Internal test, first-readable version', conflicts=source['discordant_report_cases'], n=source['multi_report_cases'], source='X21 R2'))
    q = corrected['test_report_disagreement_corrected']
    conflicts.append(dict(field='Crossbite, corrected', population='Internal test, corrected reference', conflicts=q['n_discordant_cases'], n=q['n_multi_report_cases'], source='X21 R5'))
    all_disagreement = read(_release_expand('X7'), 'raw/REPORT_DISAGREEMENT.json')
    for field in FIELDS:
        if field == 'crossbite':
            q = corrected['all_cohort_disagreement_corrected']
            n = q['n_multi_report_cases']
            k = q['n_discordant_cases']
            version = 'X21 R5 corrected'
        else:
            q = all_disagreement[field]
            n = q['n_patients']
            k = sum((len(set(p['all_values'])) > 1 for p in q['pairs']))
            version = 'X7 R2'
            parity('Whole-cohort conflict ' + field, k / n, q['all_report_discordance'])
        conflicts.append(dict(field=NAMES[field] + (', corrected' if field == 'crossbite' else ''), population='All eligible cases with multiple readable reports', conflicts=k, n=n, source=version))
    for row in conflicts:
        row['fraction'] = row['conflicts'] / row['n']
        row['ci95'] = wilson(row['conflicts'], row['n'])
        verify_wilson(row['conflicts'], row['n'], row['ci95'])
        row['resolution_level'] = 'PER_ARCH -> POPULATION'
        parity('Conflict denominator ' + row['field'] + ' ' + row['population'], row['fraction'], row['conflicts'] / row['n'])
        pred = lambda n: n == row['n']
        checked('Reject 989 denominator ' + row['field'] + ' ' + row['population'], pred(row['n']), pred(989))
    table('table4_report_conflicts', 'Table 4. Report conflict conditional on multiple readable reports', [('field', 'Field', str), ('population', 'Population', str), ('conflicts', 'Conflicts', str), ('n', 'n with ≥2 reports', str), ('fraction', 'Fraction', lambda x: f'{x:.3f}'), ('ci95', '95% Wilson CI', interval_text)], conflicts, 'Conflict means at least two distinct extracted categories among all readable reports for one patient-case. Reader IDs and repeat-reading design are unavailable; these are report conflicts, not identified interrater error. Whole-cohort and internal-test subsets are dependent and are not pooled. These observations are an empirical ambiguity context, not a universal 7–16% accuracy floor.')
    geometry = read(_release_expand('X21'), 'results.json')['population_geometry']
    flows = [dict(stage='Advertised dataset', included=1000, excluded='UNKNOWN', denominator=1000, reason='Public description; local archive differs by six cases, reason unestablished', resolution_level='POPULATION'), dict(stage='Local paired scan entries', included=994, excluded=6, denominator=1000, reason='Six advertised cases not in local paired inventory; do not label clinical exclusions', resolution_level='PER_ARCH'), dict(stage='Cases with English IOS reports', included=989, excluded=5, denominator=994, reason='No English IOS report for five local paired cases', resolution_level='PER_ARCH'), dict(stage='Training / calibration / internal test', included='593 / 198 / 198', excluded=0, denominator=989, reason='Case-level SHA256 partition; all reports retained within case', resolution_level='PER_ARCH'), dict(stage='Tooth-labelled maps', included=993, excluded=1, denominator=994, reason='Both source STL members empty for one report-ineligible case', resolution_level='PER_TOOTH'), dict(stage='Projected overlap >0.03 mm', included=877, excluded='N/A', denominator=993, reason='Geometry caution; retained rather than discarded', resolution_level='PER_ARCH'), dict(stage='Assigned upper / lower triangles', included='49.3% / 58.3% median', excluded='50.7% / 41.7% median', denominator='per-case triangle count', reason='Faces require ≥2 identical nonzero vertex labels; median unassigned fractions, no physical area interpretation', resolution_level='PER_POINT')]
    table('table1_flow', 'Table 1. Data accounting and geometric coverage', [('stage', 'Stage', str), ('included', 'Included / retained', str), ('excluded', 'Excluded / unassigned', str), ('denominator', 'Denominator', str), ('reason', 'Reason or scope', str)], flows, 'Different branches have different denominators. The five report-ineligible cases and the empty-geometry case overlap; these losses must not be added. Category exclusions appear in Tables 2–3 and their CSV columns. No exclusion count is substituted for a clinical eligibility assessment.')
    write('raw/CONFUSION_MATRICES.json', confusion_rows)
    with (P / 'tables/test_category_distributions.csv').open('w') as handle:
        writer = csv.writer(handle)
        writer.writerow(['model', 'field', 'reference_category', 'count', 'n', 'resolution_level'])
        for q in confusion_rows:
            for (i, a) in enumerate(q['confusion_labels']):
                writer.writerow([q['model'], q['field'], a, sum(q['confusion'][i]), q['n'], 'PER_ARCH -> POPULATION'])
    with (P / 'tables/confusion_matrices.csv').open('w') as handle:
        writer = csv.writer(handle)
        writer.writerow(['model', 'field', 'reference_category', 'predicted_category', 'count', 'resolution_level'])
        for q in confusion_rows:
            for (i, a) in enumerate(q['confusion_labels']):
                for (j, b) in enumerate(q['confusion_labels']):
                    writer.writerow([q['model'], q['field'], a, b, q['confusion'][i][j], 'PER_ARCH'])
    members = []
    chosen = sorted(manifest['test'], key=lambda c: hashlib.sha256(('X39-source:' + c).encode()).hexdigest())[:2] + [_release_expand('@DENTAL_CASE_ID@')]
    with zipfile.ZipFile(manifest['zip']) as archive:
        for case in chosen:
            for report in labels[case]:
                observed = hashlib.sha256(archive.read(report['member'])).hexdigest()
                checked('External report ' + report['member'], observed == report['sha256'], observed == ('0' if observed[0] != '0' else '1') + observed[1:])
                members.append(dict(member=report['member'], sha256=observed, case_id=case))
        correction = read('PROOF_LANE_XREVIEW_BATCH2', 'evidence/X21_ellipsis_correction.json')
        blob = archive.read(correction['member'])
        checked('Independent ellipsis source identity', hashlib.sha256(blob).hexdigest() == correction['source_sha256'], hashlib.sha256(blob + b'wrong').hexdigest() == correction['source_sha256'])
        assert b'as is overbite' in blob.lower()
        members.append(dict(member=correction['member'], sha256=correction['source_sha256'], known_error='Missing increased-overbite label, training case; frozen test numbers unchanged'))
    write('raw/EXTERNAL_REFERENCE_LOCATORS.json', dict(archive=manifest['zip'], members=members, no_text_exported=True, reference_quantity='Report categories, not numeric anatomical truth'))
    metadata = json.loads((P / 'sources/DOI_VERIFICATION.json').read_text())
    for ref in metadata['records']:
        if ref.get('identity_verified'):
            pred = lambda value: value.lower() == ref['doi'].lower()
            checked('DOI ' + ref['name'], pred(ref['doi']), pred('10.0000/wrong'))
            assert sha(P / ref['metadata_file']) == ref['sha256']
    write('raw/REPORTING_NUMBERS.json', dict(arch_profile=base_rows, opposition=main_rows, corrected_crossbite=corrected_rows, report_conflicts=conflicts, flows=flows))
    make_figures(main_rows, conflicts, base_rows)
    r2 = read(_release_expand('X7'), 'raw/RESULTS_R2.json')
    r4 = read(_release_expand('X7'), 'raw/RESULTS_R4.json')
    joint = read(_release_expand('X7'), 'raw/JOINT_DIAGNOSTIC.json')
    out = dict(lane='X39-manuscript-bite', claim_type=['information_link', 'capability'], review_state='PENDING_INDEPENDENT_REVIEW', scope='Retrospective manuscript compilation; no refit/new bootstrap/measurement', arch_profile=base_rows, opposition=main_rows, corrected_crossbite=corrected_rows, report_conflicts=conflicts, geometry=geometry, sets=dict(fieldwise=r2['macro'], joint_coverage=r4['joint_patient_coverage'], joint_mean_size=r4['mean_field_set_size'], joint_singleton=r4['singleton_fraction'], joint_gate=r4['primary_gate'], fieldwise_simultaneous_coverage=joint['all_available_fields_all_reports_covered']), original_gates=dict(X7_R1=baseline['primary_gate'], X7_R2=r2['primary_gate'], X7_R4=r4['primary_gate'], X21_R2=opposition['primary_gate'], X21_R5=corrected['primary_gate'], X21_R6=read(_release_expand('X21'), 'raw/RESULTS_R6.json')['primary_gate']), publication_ready=False, submission_gaps=['Author identities, funding, conflicts and authorship contribution statements', 'Dataset access terms; CC BY-NC-SA is brief-supplied, exact Bite2Text license unverified', 'Secondary-use ethics/exemption determination', 'Independent clinician label adjudication and new external/untouched test cohort', 'Target-domain FDI/landmark and matched contact reference'], full_cost=dict(new_fits=0, new_bootstrap_draws=0, compiler_seconds=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, preparation_discovery='UNKNOWN unmetered document work', upstream_cost='Not charged as zero; see locked source costs', measurement='NOT_RUN', fallback='UNKNOWN lab acquisition cost'), controls=dict(n=len(CHECKS), all_pass=all((q['passes'] and q['injected_wrong_value_rejected'] for q in CHECKS))), resolution_levels=['PER_POINT', 'PER_TOOTH', 'PER_ARCH', 'POPULATION'], time_scale='SIMULTANEOUS', claim_scope='Category agreement and mathematical projected mesh; no physical/clinical diagnosis claim', large_arrays=[], graph_target=None, proposed_graph_target='DENT-VAL-BITE2TEXT-REPORT-AGREEMENT')
    out['phenomenological_debts'] = [dict(quantity='Inherited rigid vertical pose perturbation', values_mm=[0.05, 0.1], resolution_level='PHENOMENOLOGICAL', status='DECLARED_SCENARIO_NOT_MEASUREMENT', replacement_measurement='Repeated matched bite registrations and independent reference-scanner/landmark errors for the same arches', source='X21 PREREG_R2 representation.uncertainty', physical_error_inference=False)]
    write('VALIDATION.json', dict(checks=CHECKS, pass_=out['controls']['all_pass'], scope='Reporting/source/statistical controls; no clinical validation'))
    write('results.json', out)
    write('CURRENT_WORK_STATE.json', dict(lane='X39-manuscript-bite', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='EVIDENCE_TABLES_VERIFIED', latest_gate=f'{len(CHECKS)} checks passed and all injected wrong values rejected', next_operation='Complete manuscript, guideline checklists and publication gap handoff'))
    print(f"Compiled existing evidence: {len(CHECKS)} checks PASS; no model refits or bootstrap. {out['full_cost']['compiler_seconds']:.2f} s.")

def make_figures(opposition, conflicts, baseline):
    vendor = R / 'LANE_X7_BITE2TEXT/code/vendor_plotting'
    if vendor.is_dir():
        sys.path.insert(0, str(vendor))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    (fig, axs) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1.05, 1]})
    ax = axs[0]
    for (i, row) in enumerate(opposition):
        (lo, hi) = row['ci95']
        ax.errorbar(row['accuracy'], i, xerr=[[row['accuracy'] - lo], [hi - row['accuracy']]], fmt='o', color='#16697a', capsize=3)
        ax.plot(row['comparator_accuracy'], i, 's', color='#ac513b', markersize=5)
    ax.set_yticks(range(len(opposition)), [r['field'] for r in opposition])
    ax.invert_yaxis()
    ax.set_xlim(0.35, 1)
    ax.set_xlabel('First-readable-report agreement')
    ax.set_title('A  Fixed-model internal comparison')
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], marker='o', color='#16697a', label='Contact + opposition; 95% bootstrap'), Line2D([0], [0], marker='s', linestyle='', color='#ac513b', label='Arch-profile in identical cases')], loc='lower right', fontsize=8)
    ax = axs[1]
    rr = conflicts[:5]
    for (i, row) in enumerate(rr):
        (lo, hi) = row['ci95']
        ax.errorbar(row['fraction'], i, xerr=[[row['fraction'] - lo], [hi - row['fraction']]], fmt='o', color='#55456f', capsize=3)
        ax.text(0.49, i, f"{row['conflicts']}/{row['n']}", va='center', ha='right', fontsize=9)
    ax.set_yticks(range(len(rr)), [r['field'] for r in rr])
    ax.invert_yaxis()
    ax.set_xlim(0, 0.5)
    ax.set_xlabel('Any conflict among readable reports')
    ax.set_title('B  Report ambiguity; 95% Wilson')
    for ax in axs:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='x', alpha=0.2)
    fig.text(0.5, 0.01, 'Same previously analysed case cohort. Report conflicts are not an identified clinician error floor.', ha='center', fontsize=9)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    save(fig, 'figure2_agreement_and_conflict')
    (fig, ax) = plt.subplots(figsize=(9, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    boxes = [(5, 9, 'Advertised release: 1,000 cases\nLocal paired entries: 994 (6 absent; cause unknown)'), (2.4, 6.9, 'Report branch: 989 cases\n5/994 without English IOS reports (0.5%)'), (2.4, 4.7, 'Case partition: 593 train / 198 calibration / 198 test\nAll reports from a case remain in one partition'), (2.4, 2.4, 'Field-specific test counts: 178–195 (arch-profile)\nMissing/unreadable first-report labels: 1.5–10.1%\nFirst-readable opposition subsets: 190–193*'), (7.6, 6.9, 'Geometry branch: 993 pairs\n1/994 with empty STL sources (0.1%)'), (7.6, 4.7, '38,280 projected-overlap tooth pairs\nMedian unassigned faces: upper 50.7%, lower 41.7%'), (7.6, 2.4, '877/993 projected overlaps >0.03 mm (88.3%)\n440 refined bounds / 12 cases: width <0.005 mm\nPhysical contact/force and target landmarks: UNKNOWN')]
    for (x, y, text) in boxes:
        ax.text(x, y, text, ha='center', va='center', fontsize=9, bbox=dict(boxstyle='round,pad=.6', facecolor='#f0f4f6', edgecolor='#66828d'))
    for (xy, xytext) in [((2.4, 7.55), (5, 8.35)), ((7.6, 7.55), (5, 8.35)), ((2.4, 5.4), (2.4, 6.25)), ((2.4, 3.15), (2.4, 4.05)), ((7.6, 5.4), (7.6, 6.25)), ((7.6, 3.15), (7.6, 4.05))]:
        ax.annotate('', xy=xy, xytext=xytext, arrowprops=dict(arrowstyle='->', color='#66828d'))
    ax.text(5, 0.55, '*Corrected crossbite: 184 cases. Branch losses overlap; they are not summed.\nNo photographs or new clinical measurements analysed.', ha='center', fontsize=9)
    fig.tight_layout()
    save(fig, 'figure1_case_flow')

def save(fig, name):
    import matplotlib.pyplot as plt
    for suffix in ['png', 'pdf', 'svg']:
        fig.savefig(P / 'figures' / (name + '.' + suffix), dpi=220, bbox_inches='tight')
    plt.close(fig)
if __name__ == '__main__':
    main()
