"""Bind later exact refinement and scope-aware source observations without changing old gates."""
import json
from contact import P, sha, write, state

def run():
    r = json.load(open(P / 'results.json'))
    extra = {k: json.load(open(P / 'raw' / f'RESULTS_{k}.json')) for k in ['R4', 'R5', 'R6']}
    r['rounds'].update(extra)
    r['report_observation_correction'] = '121 crossbite-negation source labels corrected, original X7/X21 labels/results retained. R6 separates anterior/lateral/unspecified statements. R5 failed blanket-region control preserved.'
    r['tight_gap_certificate'] = dict(n_cases=12, n_pairs=440, maximum_width_mm=extra['R4']['maximum_width_mm'], full_cohort_remaining='981 cases retain broad conservative R3 enclosures', domain='Rounded projected mesh, not patient truth')
    r['external_referent'] = extra['R5']['external_referent']
    r['costs']['exact_refinement_R4'] = extra['R4']['cost']
    r['costs']['report_repair_R5'] = dict(extraction=extra['R5']['source_repairs'], fit=json.load(open(P / 'raw/FIT_COST_R5.json')), evaluation_s=extra['R5']['evaluation_wall_s'])
    r['costs']['region_scope_R6_s'] = extra['R6']['wall_s']
    r['raw_negation_correction_control'] = '125 region-specific literal negations reject injected lateral presence'
    write(P / 'results.json', r)
    if (P / 'EXTERNAL_REFERENTS.json').exists():
        references = json.load(open(P / 'EXTERNAL_REFERENTS.json'))
        r['external_referents_supplement'] = references
        for k in ['R3', 'R4']:
            ref = r['rounds'][k]['external_referent']
            ref['frozen_locator'] = ref['locator']
            ref['locator'] = references['continuous_projection']['locator']
            ref['reference_metadata_appended_after_experiment'] = True
        write(P / 'results.json', r)
    a = extra['R4']
    b = extra['R5']
    c = extra['R6']
    orig = b['test_report_disagreement_original']
    cor = b['test_report_disagreement_corrected']
    all0 = b['all_cohort_disagreement_original']
    all1 = b['all_cohort_disagreement_corrected']
    lead = f"""# Dental gap certificates and actual findings of clinic reports\n\nRun `./run_all.sh` . [Contact and Report Figure](figures/contact_report_demo.png), [range width of certificates](figures/certificate_width.png), [tandparsmap as CSV ](raw/TOOTH_PAIR_MAP.csv).\n\n**993 scan pairs provide 38 280 tooth pairs with continuous projected gap intervals and source witnesses.** On the default 12 case panel, all 440 ranges are narrower than **0,005 mm**; complete rational tooth pairs and injected errors pass. The second 981 cases retain wider conservative certificates. Geometry physical , anatomical and clinical validity is still UNKNOWN .\n\nReading the external clinic facit revealed **121 incorrectly positive cross-beast labels** in X7 : s negative extraction. After a limited correction, the text conflicts in the test group drop from **{orig['n_discordant_cases']}/{orig['n_multi_report_cases']} ({orig['discordance']:.1%})** to **{cor['n_discordant_cases']}/{cor['n_multi_report_cases']} ({cor['discordance']:.1%})**. Among all {all1['n_multi_report_cases']} cases with multiple readable cross-bit reports are falling {all0['n_discordant_cases']} konflikter to {all1['n_discordant_cases']} ({all1['discordance']:.1%}). This applies to reporting categories with preserved range, not to a writer-identified inter-assessing noise floor. Author ID missing.\n\nFive of the 125 literal  A global positive label is then compatible with lateral absence. R5 : s global fault injection failed on these five ; R6 changes observation representation to separate regions and rejects errors polarity within each negative lateral statement. Total retained **{c['mixed_anterior_positive_lateral_negative_reports']}** mixed anterior/lateral reports. Original gate responses remain.\n\n| Corrected crossbite reference, 184 test cases | Balanced match | Exact match |\n|---|---:|---:|\n| New frozen contact/opposition model R5 | {b['scores']['corrected_label_model']['balanced_accuracy']:.3f} | {b['scores']['corrected_label_model']['accuracy']:.3f} |\n| Already frozen R2 , without retraining | {b['scores']['frozen_R2_model']['balanced_accuracy']:.3f} | {b['scores']['frozen_R2_model']['accuracy']:.3f} |\n| Already frozen X7 , without retraining | {b['scores']['frozen_X7_model']['balanced_accuracy']:.3f} | {b['scores']['frozen_X7_model']['accuracy']:.3f} |\n\nR5 : s category gates: BA , report coverage **{b['all_report_set_coverage']:.3f}**and size**{b['mean_set_size']:.2f}**passes; the preserved global negation control still makes the sanctuary**FAIL **. R6 : s region control and R4 : s narrow mathematical certificates **PASS**. No new algorithm superior is claimed; X7 provides a highly balanced match for corrected cross bites.\n\nWhere geometry is informative: R2 : s overbitt matches exactly the first readable report in **83,3%** and overjet in **74,1%**. In the same previously used test cohort, the training geometry model uncertainty combines with actual report conflict for overbit (AUC0,781, descriptive bootstrap 95% 0,642 – 0,904 ) and overjet ( 0,685,0,544 – 0,818 ). This supports an association with geometric ambiguity, **not a causal explanation**. M1 /Angle and per tooth cross bite holds weaker: the simple order rule gives 30,2% only specificity to the original explicit absence group. The Scopefelelet explains a large part of the previously calculated crossbeast conflict; remaining 35 conflicts throughout the cohort have not yet identified a geometric cause.\n\nThe limitation of a physical contact map is greater than the numerical error: **877/993** source poses have projected overlap> 0,03 mm. Gap cannot therefore be interpreted as real penetration, contact area or force. No independent contact image, registration uncertainty, goal - FDI or named point measurement available. The frozen next observation panel and reader are available in FROZEN_MEASUREMENT_PANEL.json and MATCHED_MEASUREMENT_SPEC.md ; the measurements are NOT_RUN .\n\nRaw source review: raw/REPORT_CORRECTIONS_R5.json, raw/MIXED_SCOPE_GEOMETRY_R6.json and raw/PRODUCER_SOURCE_REVIEW_R5.json. All sourceklasuler is anchored in exactly ZIP member and SHA256. The review of 12 statements is made by the producer, not independent review.\n\n"""
    old = (P / 'README_DEMO.md').read_text()
    old = old.replace('# Per tooth contact geometry and bite reports', '# Initial frozen R1 — R3 trials (before negative correction)')
    old = old.replace('Klinikerrapporternas oenighet may vara', 'The original R1/R2 cross-bets values below use the older negative extraction and may not be called actual cliniceroenity. Klinikerrapporternas oenighet may vara')
    text = lead + old
    for file in ['README_DEMO.md', 'RESULTS.md']:
        (P / file).write_text(text)
    hand = (P / 'HANDOFF.md').read_text()
    hand = 'X21 final handover: six designs executed, PENDING_INDEPENDENT_REVIEW .\n\n' + lead + hand
    hand += '\nThe next design: all-case exact container definition only if the consumer requires tighter geometrybounds; for physical /clinical claims, the frozen blinded FDI /point/loaded-bite port should be filled with independent data first. The report reader must continue to wear regional scope and unspecified scopeUNKNOWN.\n'
    (P / 'HANDOFF.md').write_text(hand)
    feedback = json.load(open(P / 'GRAPH_FEEDBACK.json'))
    feedback.update(sha256=sha(P / 'results.json'), outcome='993 conditional maps;12-case/440-pair5um projected certificates PASS; whole-bite report accuracy FAIL;121 negation corrections; regional source controls PASS; physical contact/target anatomy UNKNOWN')
    write(P / 'GRAPH_FEEDBACK.json', feedback)
    plot = json.load(open(P / 'raw/PLOT_DATA.json'))
    plot['r4'] = a
    plot['r5'] = b
    plot['r6'] = c
    write(P / 'raw/PLOT_DATA.json', plot)
    write(P / 'raw/CODE_SNAPSHOT.json', dict(files=[dict(path=str(p.relative_to(P)), sha256=sha(p)) for p in sorted((P / 'code').glob('*.py'))]))
    state('FINAL_VERIFICATION', 'R4/R6 scopedPASS; R1/R2/R3/R5 failed gates preserved; physical UNKNOWN', 'One-command reproduction and source-linked figure checks')
if __name__ == '__main__':
    run()
