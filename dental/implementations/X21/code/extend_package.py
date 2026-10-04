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
    lead = f"# Dental gap certificates and actual findings of clinic reports\n\nRun `./run_all.sh`. [Contact and Report Figure](figures/contact_report_demo.png), [intygens intervallbredd](figures/certificate_width.png), [tandparskarta som CSV](raw/TOOTH_PAIR_MAP.csv).\n\n**993 skanningspar ger 38 280 tooth pairs with continuous projected gap intervals and source witnesses.** On the default 12- The case panel is all. 440 range narrower than **0,005 mm**; complete rational tooth pairs checks and injected errors pass. The other 981 The cases retain broader conservative certificates. Geometry's physical, anatomical and clinical validity remains UNKNOWN.\n\nReading of the external clinics revealed **121 felaktigt positiva korsbettsetiketter** i X7:s negation extraction. After delimited correction, the text conflicts in the test group decline from **{orig['n_discordant_cases']}/{orig['n_multi_report_cases']} ({orig['discordance']:.1%})** till **{cor['n_discordant_cases']}/{cor['n_multi_report_cases']} ({cor['discordance']:.1%})**. Bland alla {all1['n_multi_report_cases']} cases with multiple readable crossbit reports dropping {all0['n_discordant_cases']} konflikter till {all1['n_discordant_cases']} ({all1['discordance']:.1%}). It concerns reporting categories with preserved scope, not a writer-identified inter-assessing noise floor.ID saknas.\n\nFem av de 125 bokstavliga ”no lateral”-reports also describe anterior presence. A global positive label is therefore compatible with lateral absence. R5the global error injection fell on these five; R6 changes observation representation to separate regions and rejects the wrong polarity within each negative lateral statement. **{c['mixed_anterior_positive_lateral_negative_reports']}** blandade anterior/lateral-rapporter. Original gate responses remain.\n\n| Korrigerat korsbettsfacit, 184 testfall | Balanced match | Exact match |\n|---|---:|---:|\n| Ny fryst kontakt/oppositionsmodell R5 | {b['scores']['corrected_label_model']['balanced_accuracy']:.3f} | {b['scores']['corrected_label_model']['accuracy']:.3f} |\n| Redan fryst R2, without re-training | {b['scores']['frozen_R2_model']['balanced_accuracy']:.3f} | {b['scores']['frozen_R2_model']['accuracy']:.3f} |\n| Redan fryst X7, without re-training | {b['scores']['frozen_X7_model']['balanced_accuracy']:.3f} | {b['scores']['frozen_X7_model']['accuracy']:.3f} |\n\nR5:s kategorigrindar: BA, report amount coverage **{b['all_report_set_coverage']:.3f}** and size **{b['mean_set_size']:.2f}** passing; the preserved global negation control still makes the sanctuary run **FAIL**. R6Regional control and R4:s narrow mathematical certificates **PASS**No new algorithm superior is claimed; X7 provides a highly balanced match for corrected cross bites.\n\nWhere the geometry is informative: R2:s overbett matches exactly with the first readable report in **83,3%** and overjet in **74,1%**. In the same previously used test cohort, the uncertainty of the trained geometry model coincides with actual report conflict for overbite. (AUC0,781, beskrivande bootstrap95%0,642–0,904) and overjet (0,685,0,544–0,818)This supports an association with geometric ambiguity, **not a causal explanation**. M1/Angle and toothy crossbite keeps weaker: the simple rule of order provides only30,2% specificity of the original explicit absenteeism group. The Scope error explains a large part of the previously calculated crossbeast conflict; remaining35 conflicts throughout the cohort have not yet been identified as a geometric cause.\n\nThe limitation of a physical contact map is greater than the numerical error: **877/993** source poses have projected overlaps>0,03mm. Gap therefore cannot be interpreted as real penetration, contact area or force. No independent contact image, registration uncertainty, target-FDI or named point measurement available. The frozen next observation panel and reader are available in FROZEN_MEASUREMENT_PANEL.json and MATCHED_MEASUREMENT_SPEC.mdThe measurements are: NOT_RUN.\n\nRaw source review: raw/REPORT_CORRECTIONS_R5.json, raw/MIXED_SCOPE_GEOMETRY_R6.json and raw/PRODUCER_SOURCE_REVIEW_R5.json. All source crystals are anchored in exactly ZIP-Member and SHA256. Granskningen av12 statements are made by the producer, not independent review.\n\n"
    old = (P / 'README_DEMO.md').read_text()
    old = old.replace("# Per tooth contact geometry and bite reports", "# Initial frozen R1 — R3 trials (before negative correction)")
    old = old.replace('Klinikerrapporternas oenighet kan vara', "The original R1/R2 cross-bets values below use the older negative extraction and may not be called actual cliniceroenity. Klinikerrapporternas oenighet may vara")
    text = lead + old
    for file in ['README_DEMO.md', 'RESULTS.md']:
        (P / file).write_text(text)
    hand = (P / 'HANDOFF.md').read_text()
    hand = "X21 final handover: six designs executed, PENDING_INDEPENDENT_REVIEW .\n\n" + lead + hand
    hand += "\nThe next design: all-case exact container definition only if the consumer requires tighter geometrybounds; for physical /clinical claims, the frozen blinded FDI /point/loaded-bite port should be filled with independent data first. The report reader must continue to wear regional scope and unspecified scopeUNKNOWN.\n"
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
