from dental_release.paths import expand as _release_expand
import csv, json, hashlib, time, resource, sys, datetime
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_r1 import main as r1
from analyze_r2 import main as r2
from analyze_r3 import main as r3
ROOT = Path(__file__).resolve().parent
PKG = Path(_release_expand('@DENTAL_INPUT_ROOT@/artifacts/DEMO48_PACKAGE'))

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def load(p):
    return json.loads(p.read_text())

def fmt(x):
    return f'{x:.2f}'

def verify_inputs():
    for d in load(ROOT / 'raw/PACKAGE_SOURCE_MANIFEST.json'):
        snapshot = ROOT / 'raw/package_snapshot' / d['demo'] / Path(d['path']).name
        assert sha(snapshot) == d['sha256'], 'Frozen package input drift: ' + str(snapshot)
    for d in load(ROOT / 'SOURCE_LOCK.json'):
        assert sha(ROOT / d['path']) == d['sha256'], 'Source drift: ' + d['path']
    return True
coverage = {'framework_fracture': 'X1b : elastic response and fracture hypothesis; clinical frequency UNKNOWN', 'ceramic_fracture_loss': 'X1b only frame/dedicated crack origin; facade and long-term UNKNOWN', 'chipping': 'Missing: façade/interface /contact injury/aging', 'secondary_caries': 'X13 / X14 / X10 measures or models gap ; carie model missing', 'caries_loss': 'Missing: gap is not a caries forecast', 'endodontic': 'X12 pulpanity; vitality transfer missing', 'abutment_tooth_fracture': 'Missing: X1b tooth -/stop break rejects its crown mechanism', 'retention_loss': 'Missing: perfectly bonded cement is closure, no debonding/deduction model', 'aesthetic': 'Missing: clinical aesthetic outcome', 'discoloration': 'Missing: leakage / aging and aesthetic end point', 'soft_tissue': 'Missing: peri-implantation tissue /inflammation', 'bone_loss': 'Missing: implant biology; X15 is another orthodontic geometry', 'abutment_fracture': 'Missing: distance crime model', 'screw_fracture': 'Missing: screw preload / fatigue', 'screw_loosening': 'Missing: preload / relaxation', 'ceramic_fracture_or_chipping': 'X1b cannot distinguish this pooled end point; chipping missing', 'restoration_fracture_loss': 'X1b crown fracture is laboratory predecessor; implant support / lifetime missing', 'any_complication': 'No validated overall clinical outcome model'}

def make_tables(a, b, rules):
    comps = {x['demo']: x for x in a['comparisons']}
    records = []
    lines = ['Demo | Comparison rule | Changes per 100 (95 % bootstrap) | Riktning | Denominator and limitation', '--- | --- | --- | --- | ---']
    for rule in rules:
        d = rule['demo']
        x = comps.get(d)
        s = 'UNKNOWN'
        direction = 'UNKNOWN'
        n = None
        low = None
        hi = None
        if x:
            (low, hi) = x['intervals']['unit_rate_ci95'][0]
            n = x['n']
            s = f"{fmt(x['changed_per100'])} ({fmt(low)}–{fmt(hi)})"
            direction = 'mer restriktivt' if x['more_restrictive'] else 'no trades'
            if d == 'X5':
                s += '; uncalibrated scenario'
            if d == 'X13':
                s += '; group drug proxy, not patient case'
        record = dict(rule, decisions_per100=x['changed_per100'] if x else None, ci95_low=low, ci95_high=hi, n_units=n, geometric_direction=direction, clinical_direction='UNKNOWN', rule_locator='sources/iti2018.pdf' if d in ['X5', 'X8'] else f'raw/package_snapshot/{d}/README_DEMO.md; primary origin status in rule_status', reference_output=f'raw/package_snapshot/{d}/results.json')
        records.append(record)
        lines.append(f"{d} | {rule['practice_rule']} | {s} | {direction}; clinically safer/unnecessary conservative UNKNOWN | {('n=' + str(n) + '; ' if n else '')}{rule['reason']}")
    (ROOT / 'PRACTICE_CHANGES.md').write_text('\n'.join(lines) + '\n')
    with (ROOT / 'raw/practice_changes.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader()
        w.writerows(records)
    freq = list(csv.DictReader((ROOT / 'raw/review_frequencies.csv').open()))
    lines = ['The five-year risk is expressed per 100 crowns/abutment in the respective published model ( POPULATION). Each line has its own denominator. No summarizing to a unique breakdown.', '', 'Source/population | Error Mode and quantity compared | Frekvens/100 (95 % KI), n | Vad demona modellerar', '--- | --- | --- | ---']
    names = {'tooth_supported_metal_ceramic': 'Dentally supported metal ceramics', 'tooth_supported_leucite_or_lithium_disilicate': 'Tooth-supported leucitis/ LS2 pool', 'implant_supported_metal_ceramic': 'Implant supported metal ceramics', 'implant_supported_veneered_zirconia': 'Implant-supported facaded zirconia'}
    for x in freq:
        ci = 'CI not reported' if not x['ci95_low'] else f"{float(x['ci95_low']):g}–{float(x['ci95_high']):g}"
        lines.append(f"{names[x['population']]} ({x['source']}) | {x['source_quantity']} | {float(x['rate_5year_per100']):g} ({ci}), n={x['n_crowns_or_abutments']} | {coverage[x['failure_mode']]}")
    lines += ['', 'Sailer: [primary full text, Table 7, page 620 ](https://archive-ouverte.unige.ch/unige:79013), DOI 10.1016/j.dental.2015.02.011 . Pjetursson: [primary full text, Table 5, page 211 ](https://archive-ouverte.unige.ch/unige:123715), DOI 10.1111/clr.13306 . Local PDF , table images and cell extracts are available under sources/ and raw/SOURCE_CELL_CHECKS.json.', '', "The Zirkonia group of the Sailer original is not included: the 2016 is available, but the full text of the correction is UNVERIFIED here. Leucit/ LS2 is a pool and must not be described as pure LS2. Pjetursson's zirconia is façade; no new monolithic zirconia studies met the follow-up requirements. Chipping and fracture loss endpoints are different and must not be merged. Biological dental complications are not defined as corresponding implant complications."]
    (ROOT / 'FAILURE_MODES.md').write_text('\n'.join(lines) + '\n')
    return records

def figure(a, b):
    x = next((x for x in a['comparisons'] if x['demo'] == 'X8'))
    (lo, hi) = x['intervals']['unit_rate_ci95'][0]
    (f, axs) = plt.subplots(1, 2, figsize=(11, 4.7), gridspec_kw={'width_ratios': [1, 1.6]})
    axs[0].bar(['Per virtual site', 'Patients with\nany change'], [x['changed_per100'], x['patients_any_change_per100']], color=['#297a94', '#608b91'])
    (p_lo, p_hi) = x['intervals']['any_patient_change_ci95']
    axs[0].errorbar([0, 1], [x['changed_per100'], x['patients_any_change_per100']], yerr=[[x['changed_per100'] - lo, x['patients_any_change_per100'] - p_lo], [hi - x['changed_per100'], p_hi - x['patients_any_change_per100']]], fmt='none', ecolor='black', capsize=4)
    axs[0].set_ylabel('Changes per 100 eligible units')
    axs[0].set_ylim(0, 13)
    axs[0].set_title('X8: same 2 mm rule, whole body\n2581 virtual sites / 433 patients')
    axs[0].text(0.5, -0.25, 'Annotation geometry; injury benefit UNKNOWN', ha='center', transform=axs[0].transAxes, fontsize=9)
    causes = list(csv.DictReader((ROOT / 'raw/failure_causes.csv').open()))
    labels = ['Retention loss', 'Tooth loss', 'Tooth fracture', 'Changed restoration', 'Caries', 'Periapical disease', 'Bone loss', 'Major chipping', 'Mobility', 'Vitality loss', 'Aesthetics', 'Framework fracture', 'Unknown tooth loss', 'Poor fitting']
    vals = [100 * int(v['n']) / 230 for v in causes]
    colors = ['#b9bec4'] * len(vals)
    colors[11] = '#297a94'
    colors[13] = '#bd9440'
    axs[1].barh(range(len(vals)), vals, color=colors)
    axs[1].set_yticks(range(len(vals)), labels, fontsize=8)
    axs[1].invert_yaxis()
    axs[1].set_xlabel('Share of 230 attributed failures (%)')
    axs[1].set_title('External historical crown cohort\nX1b family match at most 1.74%')
    f.suptitle('Decision disagreement and failure-family scope are different quantities', fontsize=12)
    f.tight_layout()
    f.savefig(ROOT / 'figures/practice_and_failure_scope.png', dpi=170)
    f.savefig(ROOT / 'figures/practice_and_failure_scope.pdf')
    plt.close(f)

def main():
    t = time.perf_counter()
    verify_inputs()
    import extract_reviews
    original = extract_reviews.S[0]
    bad = list(original)
    bad[3] = 999.0
    extract_reviews.S[0] = tuple(bad)
    try:
        extract_reviews.extract()
        raise RuntimeError('Injected rate escaped source-cell check')
    except AssertionError:
        pass
    finally:
        extract_reviews.S[0] = original
    r1()
    r2()
    r3()
    c = load(ROOT / 'rounds/RESULTS_R3.json')
    a = load(ROOT / 'rounds/RESULTS_R1.json')
    b = load(ROOT / 'rounds/RESULTS_R2.json')
    rules = load(ROOT / 'DEMO_RULES.json')
    records = make_tables(a, b, rules)
    figure(a, b)
    comparisons = {x['demo']: x for x in a['comparisons']}
    x = comparisons['X8']
    five = comparisons['X5']
    gap = comparisons['X13']
    text = f"""# What decisions are changed, and which crown failures reach the model?\n\nThe new link makes demo decisions and clinical scope reviewable. Nervdemon Change {x['changed']}/{x['n']} virtual location decisions when the entire implant body is measured against the same annotated channel. At the same time, the external crown cohort shows that the spine fracture is {b['X1b_family_matched_failures']}/230 of the causes of failure. It distinguishes a testable fracture attempt from a general sustainability forecast.\n\nRun `./run_all.sh` from this directory. Offline, a frozen snapshot of the package, existing system-Python with NumPy /Matplotlib; maximum four threads, no GPU . The command controls source hash, recalculates the decisions/cluster intervals, reads the external cause table, tests error values and creates tables, figure and results.json . The run this comparative ability; the original demos heavy physics/geometric work is not run on new.\n\n| Question | Actually result | Evidence and resolution |\n|---|---|---|\n| 2 mm channel rule: what location decisions are changed? | {fmt(x['changed_per100'])}/100; 95 % patientbootstrap {fmt(x['intervals']['unit_rate_ci95'][0][0])}–{fmt(x['intervals']['unit_rate_ci95'][0][1])}; all more restrictive | PER_TOOTH, virtuell plan against TF2-etikett |\n| How many patients will have a slightly changed location decision? | {x['patients_with_any_change']}/{x['patients']}: {fmt(x['patients_any_change_per100'])}/100; {fmt(x['intervals']['any_patient_change_ci95'][0])}–{fmt(x['intervals']['any_patient_change_ci95'][1])} | PER_ARCH Aggregation of site decisions; no observed operational changes |\n| Scenario with an error budget, X5 | {five['changed']}/{five['n']}: {fmt(five['changed_per100'])}/100 ({fmt(five['intervals']['unit_rate_ci95'][0][0])}–{fmt(five['intervals']['unit_rate_ci95'][0][1])}); mer restriktivt | PER_TOOTH ; budget is uncalibrated, not measured new anatomy |\n| Marginal-CAD-spacer as spaltproxy, X13 | 0/22 gruppmedel byter 120 µm-klass | PER_SURFACE_REGION ; five source families, not 22 patients |\n| Which family of wrecks can X1b most hit? | 4/230 = {fmt(b['X1b_family_share_percent'])} % stommefraktur | POPULATION, extern disjunkt orsakstabell; mekanismmatchning may vara 0 |\n| How much clinical prediction is covered? | UNKNOWN; 0 validerade kliniska slutpunktsmodeller i paketet | POPULATION ; family matching is no measure of prediction |\n\nAll fifteen demos rules, wastes and unknown decisions are in [PRACTICE_CHANGES.md ](PRACTICE_CHANGES.md). Twelve lacks the prerequisites for a calculated comparison; this does not mean zero amended decisions. X5 and X13 have only defined scenario/proxy comparisons. Thus, there is no demo with a proven number of safer or unnecessarily conservative clinical decisions per 100 patients.\n\n[FAILURE_MODES.md ](FAILURE_MODES.md) has 40 verified frequency cells from two systematic overviews, with support type, material, five-year window, denominator and KI . Read the source's definition: Sailer Table 7 has separate chipping/carcase/tand fractures; Pjetursson Table 5 includes a combined ceramic fracture or chipping endpoint. They are not a common unique-haul distribution. [Sailer 2015 ](https://archive-ouverte.unige.ch/unige:79013), [Pjetursson 2018 ](https://archive-ouverte.unige.ch/unige:123715).\n\n![Decision change and error mode coverage](figures/practice_and_failure_scope.png)\n\nWhat does not holds : [ ITI : s documented 2 mm Convention](https://pure.uva.nl/ws/files/32522661/Wismeijer_et_al_2018_Clinical_Oral_Implants_Research.pdf) is not a nerve damage curve. No case-by-case clinical damage indicates that the prey is safer. 120 µm is a historical margin convention whose primary text from 1971 has not been verified here; it is used as an explicit benchmark, not cariesfacit. Internal gap , scannertrueness, maximum margin and regional funds are distinct quantities . X13 : s bootstrap provides 0 – 0 in this limited selection; it does not prove a zero frequency outside the selection.\n\nThe causal table comes from [an external historical cohort, Table 5 ]( https://pmc.ncbi.nlm.nih.gov/articles/PMC9546353/#eos12871-tbl-0005 : 1037 crowns/401 patients, 230 failures in 149 patients, mean follow-up 134,8 months. The materials and follow-up do not match a modern monolithic zirconium crown or implant prothetics. Patient ID by cause is missing: statistical patient interval is UNKNOWN. Boundary [0 ;{fmt(b['X1b_family_share_percent'])}] % only describes how many of these breakdowns can be in X1b : s exact mechanism when crack origin/material is unknown . It's not a 95 %-KI. Poor-passing family adds a maximum of one cause (all passerer together maximum{fmt(b['X1b_plus_fit_precursor_family_upper_percent'])} %); gap must not be automatically counted as retention, decay or inflammation.\n\nGone: X8 rejects 125/2706 original location lines (4,62%) with incomplete previous geometry and deduplicates 7743 wizard copies; X5 uses 2036 channel controls of 2342 original sites and rejects second control regions/selection; X13 rejects 126 / 148 measurement cells ( 85,14 %) for wrong region or missing explicit margin setting. Overview Table 40 / 80 rejected material cells: 30 outside recognized groups, 10 Sailer-Zirconia cells with unsolved full text correction from 2016 . The table is a targeted source selection, no new systematic overview.\n\nThe check with the same information gives the same threshold. Wrong number, wrong device/quantity , missing locator, wrong cause label and overlapping complications as unique breakdowns are rejected. This is an information link with separate practice comparison, no algorithm victory. Frozen replay expectations are retrospective; no new physical predictions or measurements have been created. Status PENDING_INDEPENDENT_REVIEW.\n\nData/licenses: only existing local demoartefacts are used. ToothFairy2 CC BY-SA4.0; Bits2Bites CC BY - NC - SA (exact version according to source); Mandible defects CC BY4.0 ; STS-Tooth3D CC BY4.0. Teeth3DS the origin conditions of the package indicate CC BY-NC-ND4.0 with older conflicting data ; see [the OWNERSHIP of the package ](../DEMO48_PACKAGE/OWNERSHIP.md). Bite2Text and Pulpy3D extensions are UNKNOWN. These second datasets are read here only through small existing result/README , no new mesh/CT extracted. X13 : s Article Licences are specified per source in its README. Sailer- PDF is copyright protected, Pjetursson CC BY-NC4.0, the external cause of the CC BY-NC-ND4.0 ; the documents are primary source data for private research. No common new distribution licence is adopted.\n\nThe next crucial design is a patient-bound, independently reviewed accident-cause and crack-origin list as well as measured retention after aging in the same kronbatch. It would test the largest open prosthetic family in this cohort: retention loss. See [HANDOFF.md ](HANDOFF.md) for exact continuation and graphation blocked binding .\n"""
    text += f"\nThe third construction reread all 2706 original location keys: they come from 442 represented cases, and all 125 non-response lacks L_plan . The loss is therefore not a random loss of measurement value. The original cohort decision change frequency is UNKNOWN when the plan is undefined. Only if a binary decision is hypothetically assigned to each failure will the logical complement limits be obtained {c['site_completion_bound_per100'][0]:.2f}–{c['site_completion_bound_per100'][1]:.2f} per 100 original sites and {c['patient_any_change_completion_bound_per100'][0]:.2f}–{c['patient_any_change_completion_bound_per100'][1]:.2f} per 100 represented cases with any change (POPULATION, PER_TOOTH and PER_ARCH respectively). These are not confidence intervals or surgical outcomes. Paketets480 declared scans and 479 actual archive case is a different denominator than 442 case with placebars. Se rounds/RESULTS_R3.json.\n"
    (ROOT / 'README_DEMO.md').write_text(text)
    (ROOT / 'README.md').write_text('Run `./run_all.sh` . Klinikerunderlag: [README_DEMO.md](README_DEMO.md). All 15 demon: [PRACTICE_CHANGES.md ](PRACTICE_CHANGES.md). Felmoder: [FAILURE_MODES.md](FAILURE_MODES.md).\n')
    cost = {'run_wall_seconds': time.perf_counter() - t, 'process_cpu_seconds': resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime, 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads_max': 4, 'GPU': False, 'preparation_external_source_reading_and_literature_discovery_seconds': None, 'fit': 'none', 'full_discovery_and_manual_validation_labor': 'UNKNOWN', 'inherited_demo_preparation_fit_cost': 'See original demo costs; not measured again', 'queries': '15 frozen outputs, 3 raw decision sources, original K3 keys, 2 review tables, 1 external cause table', 'fallback': 'Explicit UNKNOWN, no dataset or physical acquisition'}
    out = {'lane': 'X32-practice-changes', 'claim_type': 'information_link', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'capability': 'Paired practice-rule decision counts with clustered intervals; externally grounded clinical failure-family scope', 'clinical_benefit_per100': 'UNKNOWN', 'primary_result': x, 'case_level_comparisons': a['comparisons'], 'all_demos': records, 'failure_coverage': b, 'incompleteness_analysis': c, 'external_referent': b['external_referent'], 'external_referents': [{'kind': 'published_dataset', 'locator': 'https://toothfairy2.grand-challenge.org/dataset/', 'compared_quantity': 'Annotated canal geometry vs virtual implant clearance, not injury', 'refutes_us': True}, {'kind': 'external_review', 'locator': 'https://doi.org/10.1016/j.dental.2015.02.011', 'compared_quantity': '5-year complication risks Table7, own denominators', 'refutes_us': True}, {'kind': 'external_review', 'locator': 'https://doi.org/10.1111/clr.13306', 'compared_quantity': '5-year implant-supported single-crown complication risks Table5', 'refutes_us': True}], 'controls': {'same_information': a['same_information_control'], 'R1_faults': a['fault_injections'], 'R2_faults': b['fault_injections'], 'R3_faults': c['fault_injections'], 'source_cell_check_count': 40, 'injected_source_rate_rejected': True}, 'cost': cost, 'input_bindings': a['input_bindings'], 'source_lock_sha256': sha(ROOT / 'SOURCE_LOCK.json'), 'prereg_hashes': {p.name: sha(p) for p in sorted(ROOT.glob('PREREG_*.json'))}, 'dropout': {'demos_without_usable_paired_practice_comparison': 12, 'all_demos': 15, 'fraction': 0.8, 'X8_original_sites': 2706, 'X8_incomplete_sites': 125, 'X8_incomplete_fraction': 125 / 2706, 'X13_source_cells': 148, 'X13_rejected_cells': 126, 'X13_rejected_fraction': 126 / 148, 'review_cells': b['review_exclusions']}, 'large_arrays_over_50MB': [], 'granularity_policy': 'Edges use finest shared level; consumers aggregate', 'physical_clinical_claims': 'No new clinical endpoint validation or clinical recommendation', 'time_scale': 'SIMULTANEOUS decision comparison; clinical HANDOVER uncalibrated'}
    (ROOT / 'results.json').write_text(json.dumps(out, indent=2) + '\n')
    (ROOT / 'RESULTS.md').write_text(text)
    feedback = ROOT / 'GRAPH_FEEDBACK.json'
    if feedback.exists():
        fb = load(feedback)
        fb['sha256'] = sha(ROOT / 'results.json')
        fb['review_state'] = 'PENDING_INDEPENDENT_REVIEW'
        feedback.write_text(json.dumps(fb, indent=2) + '\n')
    replay3 = ROOT / 'FROZEN_PREDICTIONS_R3.json'
    if replay3.exists():
        claims = load(replay3)['expected_counts']
        actual = {'original_site_keys': c['represented_original_sites'], 'represented_case_keys': c['represented_original_patient_keys'], 'missing_L_plan': c['unmatched_rows_with_no_L_plan'], 'patients_with_unmatched': c['patients_with_unmatched_sites'], 'additional_potentially_changed': c['additional_potentially_changed_patients']}
        assert actual == claims, 'R3 frozen replay counts drift'
    replay = {'kind': 'RETROSPECTIVE_REPRODUCTION_EXPECTATIONS_NOT_PHYSICAL_PREDICTION', 'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'claims': {'X8_changed': 37, 'X8_n': 2581, 'X5_scenario_changed': 224, 'X5_n': 2036, 'X13_group_changes': 0, 'X13_n': 22, 'failure_denominator': 230, 'framework_failure_count': 4, 'review_cells': 40}, 'source_lock_sha256': sha(ROOT / 'SOURCE_LOCK.json'), 'future_physical_measurement': 'NONE; separate preregistration and prediction freeze required'}
    frozen = ROOT / 'FROZEN_PREDICTIONS.json'
    if frozen.exists():
        assert load(frozen)['claims'] == replay['claims'], 'Frozen replay claims drift'
    else:
        frozen.write_text(json.dumps(replay, indent=2) + '\n')
    assert x['changed'] == 37 and five['changed'] == 224 and (gap['changed'] == 0) and (b['X1b_family_matched_failures'] == 4)
    print('X32 COMPLETE: paired counts, source tables, falsifying probes and figure verified.')
if __name__ == '__main__':
    main()
