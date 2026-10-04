"""Reader-facing statistics, one table per demo, machine table, and exportable figure."""
import collections
import csv
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parent

def fmt(x):
    if x is None:
        return 'UNKNOWN'
    return f'{x:.5g}' if isinstance(x, (int, float)) else str(x)

def interval(x):
    return '[' + ', '.join((fmt(y) for y in x)) + ']' if x is not None else 'UNKNOWN'

def locator(f):
    ext = f.get('external_referent', {})
    return ext.get('locator', f.get('external_locator', 'Locator in demo external_referents / SOURCE_MANIFEST.'))

def figure(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lookup = {(m['demo'], m['metric']): m for m in result['metrics']}
    selections = [('X3', 'R3 upper-surface p95 mean', 1, 'mm'), ('GENCAD_V2', 'test digital PASS constraint_optimizer', 100, '% tasks'), ('X15', 'R4 survival of R3-admissible queries', 100, '% queries'), ('X18', 'median informed/practice area-error ratio', 1, 'ratio'), ('X22', 'test AUC', 1, 'AUC'), ('GENCAD_V3', 'molar_crown digital PASS constraint_optimizer', 100, '% tasks'), ('X31', 'fully_guided apex-scenario decision change', 100, '% site decisions')]
    selections = [(d, l if (d, l) in lookup else next((m['metric'] for m in result['metrics'] if m['demo'] == d and m.get('primary'))), scale, u) for (d, l, scale, u) in selections]
    (fig, axs) = plt.subplots(2, 4, figsize=(16, 7))
    for (ax, (d, label, scale, u)) in zip(axs.flat, selections):
        m = lookup[d, label]
        for (y, key, color) in [(1, 'naive_row_ci95', '#b46a25'), (0, 'cluster_ci95', '#21658b')]:
            bounds = m[key]
            p = m.get('naive_row_point', m['point']) if y == 1 else m['point']
            ax.hlines(y, bounds[0] * scale, bounds[1] * scale, color=color, lw=3)
            ax.plot(p * scale, y, 'o', color=color)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(['Whole cases', 'Row IID'])
        ax.set_ylim(-0.5, 1.5)
        ax.grid(axis='x', alpha=0.2)
        ax.set_xlabel(u)
        name = d + ' (archived)' if d == 'GENCAD_V2' else d
        ax.set_title(f"{name}: {m['n_records']} rows / {m['k_clusters']} cases", fontsize=10)
    ax = axs.flat[-1]
    c = next((c for c in result['round2']['paired_case_contrasts'] if 'collar_mesh' in c['metric']))
    ax.axvspan(-0.05, 0.05, color='#9dc4aa', alpha=0.45, label='Original equivalence margin')
    ax.hlines(0, *c['case_t90'], color='#21658b', lw=3)
    ax.plot(c['point_mm'], 0, 'o', color='#21658b')
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([0])
    ax.set_yticklabels(['Case t90%'])
    ax.set_title('X3: equivalence not established (K=4)', fontsize=10)
    ax.set_xlabel('SDF minus mesh error (mm)')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(axis='x', alpha=0.2)
    fig.suptitle('Frozen predictions: conditional intervals, not clinical validation', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    dest = ROOT / 'figures'
    dest.mkdir(exist_ok=True)
    fig.savefig(dest / 'cluster_intervals.png', dpi=180)
    fig.savefig(dest / 'cluster_intervals.pdf')
    plt.close(fig)

def build(result, checklists):
    rows = result['metrics']
    demos = result['demos']
    r2 = result['round2']
    nlimited = [d['id'] for d in demos if d['support'].startswith(('POPULATION_INSUFFICIENT', 'TRANSPORT_INSUFFICIENT'))]
    unknown = [d['id'] for d in demos if d['support'].startswith('UNKNOWN')]
    result['sample_size_limited_demos'] = nlimited
    result['unidentified_independence_or_variance_demos'] = unknown
    lines = ['# Demo package statistics: independent cases and open evidence limits', '', 'All **35 retained demos in the R3 snapshot** have a table below: **34 active demos and archived GenCAD v2**. X33/v3 were added during analysis; the previous 33 reference records retain their semantics. The original R1 scope of 33 cases and detected hash drift are preserved. Recalculation retains outcome definitions and failed gates. References are hash-bound local measurements/annotations with external locators; statistical references are separate algebra, binomial tails and rank sums. No new physical measurement was made.', '', '| Decided question | Result | Independent unit |', '|---|---|---|', '| X3: mean tooth surface p95 error | 0.7303 mm; row-IID95% [0.5770-0.9210], case bootstrap [0.4973-1.0902] | 32 teeth / 4 scan cases |', '| GenCAD v2 (archived): digital PASS, constraint_optimizer | 81.94%; case bootstrap95% [75.35-88.19] | 288 test tasks / 8 test cases; 576 tasks across all splits |', '| GenCAD v2 (archived): all assigned designs PASS in one case | 0/8; exact two-sided95% [0-36.94%] under IID cases | Same 8 cases, different estimand from task fraction |', '| GenCAD v3 (active): actual test size | 588 patient/cases per family; 2352 dependent level tasks, 27 recalculated family/method fractions | Within Bite2Text; original already uses case bootstrap |', '| X3: mesh equivalence within +/-0.05 mm | Delta=0.0060 mm; case-t90% [-0.05545-0.06749], TOST p=0.0954: not established | 4 paired cases, Gaussian assumption |', '| X25: no 100 N failures throughout specimen horizon | 0/3; one-sided95% upper risk 63.16% under IID specimens | 3 full survivors; batch identity UNKNOWN |', '', '**Insufficient evidence for broad population/transfer claims:** ' + ', '.join(nlimited) + '. Reasons are given per demo; low n is not a universal rule that large effects cannot be detected. **Unidentified specimen variance/clustering:** ' + ', '.join(unknown) + '. A valid negative result, export identity or constructed counterexample does not require an invented population size.', '', '![Conditional intervals and X3 equivalence test](figures/cluster_intervals.png)', '', '## Reading the numbers', '', '**Before** explicitly recalculates endpoint rows as IID; it is not automatically the original reported interval. **After** resamples entire case/study blocks with replacement 10000 times (frozen seed 6103501). All paired methods, report values and repeated queries follow the same block. Identified original intervals are also shown separately. Original case bootstrap is retained as correctly scoped practice. Small Monte Carlo differences when K=N are not a biological cluster effect.', '', 'Clustering must preserve the estimand. Equal case weighting uses the mean of case means. R4 row IID retains each row frozen inverse count within its original block and normalizes sampled weight mass; before/after therefore share the equal-study point estimate. Earlier unweighted row sensitivity for X13/X1B is preserved in RESULTS_R3_COMPLETE.json and does not count as a pure cluster effect. X16 retains original equal weights for seven family folds while resampling their **five** shared source studies. Medians/median ratios are recalculated over all rows in sampled blocks.', '', '**K** counts observed blocks, not verified independent patients. TF2 volumes and STS prefixes are case proxies; cross-dataset identity and repeat visits are unresolved. **n_eff(rho=1)** = (sum m)^2/sum m^2 is concentration sensitivity if all rows within a block are perfectly dependent, not an estimated ICC or universal effective patient n. Equal case weighting reports K. UNKNOWN lineage remains UNKNOWN.', '', '**MDE80/alpha0.05**: two-sided noncentral t for hypothetical independent Gaussian case means, reported in their SD units. When raw SD exists for the same mean estimand, a native-unit planning scenario is shown. Native MDE for RMSE, ratios, medians, kappa and AUC is UNKNOWN; the d multiplier is a general planning reference, not power for these metrics. It is not retrospective observed power or proof that an observed effect exists/is absent; gates are not adjusted afterward. K<10 is flagged as pilot because empirical bootstrap resamples only a few observed anatomies/studies.', '', 'Percentile intervals are marginal 95%, conditional on delivered predictions, selection and annotations. They exclude refitting, model selection, shared cross-validation fits, label bias, systematic scanner error and physical transport. Shared data/fits make demos and folds dependent; package count does not count independent validations. No package-wide confirmatory test series is declared.', '', 'R2 exact **all-design-PASS per case** intervals prevent false zero uncertainty when empirical bootstrap gives [0,0]. This changes to a declared case event and is not an interval for task fraction. X3 exact sign flip has 16 possible sign patterns: smallest two-sided p is 0.125; this assumes symmetric independent errors, not randomized treatment. t/TOST gives another answer under stronger Gaussian assumptions.', '', 'Every numeric row states input/output resolution. The information edge uses the finest common level; downstream aggregates are separate. This reanalysis is SIMULTANEOUS for the same frozen outcomes. It creates no new biological HANDOVER link. Each uncalibrated population/batch closure is recorded as a debt in results.json.', '', '## Reporting sources and DOI', '', 'CLAIM2024: [Tejani et al.,10.1148/ryai.240300](https://pubs.rsna.org/doi/pdf/10.1148/ryai.240300), table on printed page 3, 44 items. Applies to medical imaging AI; geometry tools use explicitly adapted questions. The guideline is not a quality grade.', '', 'TRIPOD+AI: [Collins et al.,10.1136/bmj-2023-078378](https://pmc.ncbi.nlm.nih.gov/articles/PMC11019967/#tbl2), Table 2, 27 main items/52 subitems. Clinical prediction models; material/CAD forecasts fall outside formal scope.', '', 'STARD2015: [Bossuyt et al.,10.1136/bmj.h5527](https://doi.org/10.1136/bmj.h5527); [official checklist](https://www.equator-network.org/wp-content/uploads/2015/03/STARD-2015-checklist.pdf), 30 main items/34 subitems. 10.1136/bmjopen-2016-012799 is explanation/elaboration, not the 2015 statement. No categorical/synthetic demo is relabeled as real diagnostic validation.', '', 'CRIS: [Krithikadatta et al.,10.4103/0972-0707.136338](https://pmc.ncbi.nlm.nih.gov/articles/PMC4127685/). The 2014 article is a concept proposal calling for development and validation; it has no finalized numbered checklist. Here 8 explicitly marked **LOCAL** thematic prompts are used, not fabricated official CRIS certification.', '', 'Statistical methods: [Field & Welsh2007, cluster bootstrap](https://doi.org/10.1111/j.1467-9868.2007.00593.x); [Clopper & Pearson1934, binomial bounds](https://doi.org/10.1093/biomet/26.4.404); [Schuirmann1987, TOST](https://doi.org/10.1007/BF01068419). Constitutive assumptions and limits of uncertainty are in DECOMPOSITION.json and PREREG_R1/R2.', '', 'Matrices cover every item as review material for demo artifacts. PARTIAL, MISSING, NOT_LOCATED and NOT_APPLICABLE have reasons/locators; no automated manuscript certification or quality score is assigned. Original study blinding/ethics cannot be attributed to this demo without evidence.', '']
    csvrows = []
    for d in demos:
        id = d['id']
        f = d['external_reference']
        a = d['attrition']
        cl = checklists[id]
        applicable = ', '.join((x['guideline'] for x in cl['checklists'])) or 'N/A (scoped geometry/identity)'
        lines.extend([f'## {id}', '', d['title'], '', '**Paketstatus:** ' + d['package_status'] + '.', '', '**Capability and evidence limit:** ' + d['support'], '', '**Oberoende enhet:** ' + d['independent_unit'] + '.', '**Externt facit:** ' + str(locator(f)) + '. Facitets storhet/scope: ' + str(f.get('external_referent', {}).get('compared_quantity', f.get('scope', 'se ursprungsdemo'))) + '.', '', f"""**Urval/attrition:** {fmt(a['retained'])}/{fmt(a['candidate'])} retained, {fmt(a['excluded'])} removed ({(fmt(100 * a['excluded_fraction']) + '%' if a['excluded_fraction'] is not None else 'UNKNOWN')}). {a['reason']}""", '', f"""**Checklist:** [{applicable}](CHECKLISTS/{id}.md): {cl['status']}. Formal compliance is not established.""", ''])
        if f.get('primary_value') is not None:
            lines.extend([f"Paketets facitpekare `{f.get('primary_pointer')}` = **{fmt(f['primary_value'])} {f.get('units', '')}**. " + ('This is a software test count; empirical design fractions below are separate.' if id == 'GENCAD_V2' else ''), ''])
        if d.get('original_intervals_note'):
            lines.extend([d['original_intervals_note'], ''])
        lines.extend(['| Metric (unit) | Point | Original95% | Before: row-IID95% | After: block95% | Rows / K / n_eff(rho=1) | MDE d / native | Resolution | Checklist status |', '|---|---:|---|---|---|---|---|---|---|'])
        for m in d['metrics']:
            na = m.get('naive_row_ci95')
            bc = m.get('cluster_ci95')
            orig = m.get('original_reported_ci95')
            after = interval(bc) if bc is not None else 'N/A' if m['inference'] == 'NOT_APPLICABLE' else 'UNKNOWN'
            before = interval(na) if na is not None else 'N/A' if m['inference'] == 'NOT_APPLICABLE' else 'UNKNOWN'
            md = fmt(m.get('mde_80pct_alpha05_cluster_SD')) + ' / ' + fmt(m.get('mde_native_approx'))
            resolution = m['input_resolution'] + ' → ' + m['resolution_level']
            lines.append(f"""| {m['metric']} ({m['unit']}) | {fmt(m['point'])} | {(interval(orig) if orig is not None else 'Unidentified / N/A')} | {before} | {after} | {fmt(m['n_records'])} / {fmt(m.get('k_clusters'))} / {fmt(m.get('kish_cluster_equivalent'))} | {md} | {resolution} | {('Incomplete' if cl['checklists'] else 'N/A scope')} |""")
            csvrows.append(dict(demo=id, metric=m['metric'], unit=m['unit'], point=m['point'], original95=json.dumps(orig), before95=json.dumps(na), after95=json.dumps(bc), n_records=m['n_records'], K=m.get('k_clusters'), n_eff_rho1=m.get('kish_cluster_equivalent'), MDE_d=m.get('mde_80pct_alpha05_cluster_SD'), MDE_native=m.get('mde_native_approx'), scope=m['inference'], checklist=cl['status'], support=d['support']))
        lines.append('')
        for m in d['metrics']:
            if m.get('primary') or m.get('note') or m.get('reported_specimen_t95_if_IID'):
                lines.append(f"- **{m['metric']}:** {m.get('note', '')}" + (' IID-provens t95-scenario ' + interval(m['reported_specimen_t95_if_IID']) + ' ' + m['unit'] + '; ingen identifierad batch-CI.' if m.get('reported_specimen_t95_if_IID') else ''))
        events = [e for e in r2['case_events'] if e['demo'] == id]
        if events:
            lines.extend(['', '| Separate case event: all tasks PASS | Accepted cases/K | Exact two-sided95% |', '|---|---:|---|'])
            for e in events:
                q = e['interval']
                lines.append(f"| {e['metric']} | {q['successes']}/{q['n']} | {interval(q['ci95'])} |")
            lines.extend(['', 'Event concerns digital complete PASS, including missing answers as non-PASS. Conditional on independent representative cases; it does not imply a clinical risk bound.'])
        if id == 'X3':
            lines.extend(['', '| Paired case test | Delta mm | t95% | Exact symmetric sign-flip p | TOST+/-0.05mm |', '|---|---:|---|---:|---|'])
            for c in r2['paired_case_contrasts']:
                lines.append(f"| {c['metric']} | {fmt(c['point_mm'])} | {interval(c['case_t95'])} | {fmt(c['exact_sign_flip_two_sided_p'])} | {(('PASS' if c['equivalence_pass'] else 'FAIL') + '; p=' + fmt(c['tost_p']) if 'equivalence_pass' in c else 'N/A')} |")
        lines.extend(['', 'Source: [' + id + ' README](' + d['path'] + '/README_DEMO.md:1); raw endpoint records in [ENDPOINT_RECORDS](raw/ENDPOINT_RECORDS.json); file hashes in [INPUT_LOCK](INPUT_LOCK.json).', ''])
    lines.extend(['## Shared evidence and next measurement', '', 'X7/X21 share Bite2Text cases; X8/X31/X32 share nerve geometry; X12/X30 and X31 reduction share tooth/pulp annotations; X5/X26/X32 reuse canal decisions; X13/X32 reuse margin group means. Participant linkage across STS/Teeth3DS/Bits2Bites in GenCAD is unverified. Case counts and CIs must therefore not be summed across demos.', '', 'Next construction is a batch/patient register with case -> patient/donor -> batch/center -> independent split, followed by predictions frozen before new cases and physical measurement in the same specimen region. X25 separate 0/29 plan for risk <=10% is an externally anchored IID specimen scenario; shared batches may make even 29 specimens insufficient. Power scenarios are planning, not executed experiments.', '', f"""Control outcome: {result['validation']['status']}, {result['validation']['n_controls']} controls that reject injected errors; see [VALIDATION.json](VALIDATION.json). An equally informed algebraic control should match a statistical reanalysis. No algorithmic improvement is claimed.""", '', 'Source graphs were unchanged. Graph reading view was STALE_INPUT; no permitted existing packet could be bound. GRAPH_FEEDBACK.json records scoped missing coverage with PENDING_INDEPENDENT_REVIEW.'])
    (ROOT / 'STATISTICS.md').write_text('\n'.join(lines) + '\n')
    with (ROOT / 'STATISTICS_TABLE.csv').open('w') as out:
        writer = csv.DictWriter(out, fieldnames=list(csvrows[0]))
        writer.writeheader()
        writer.writerows(csvrows)
    figure(result)
    return result
