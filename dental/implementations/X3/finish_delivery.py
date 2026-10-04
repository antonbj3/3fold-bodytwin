from crownbench import *
import csv

def finish():
    r3 = json.load(open(P / 'RESULTS_R3.json'))
    result = json.load(open(P / 'results.json'))
    result['rounds']['R3'] = r3
    result['input_regimes'] = {'R1_R2': 'complete tooth deletion', 'R3': 'observed natural cervical1mm band; distinct partial-crown completion task, not pooled with full deletion'}
    result['external_referent_R3'] = r3['external_referent']
    result['large_data_bytes'] = disk_guard()
    result['numerical_kernel'] = {'K1': 'FAILED1e-10mm equality against Rtree; near-tie face-normal switch; failure retained', 'K2': 'PASS exact equality with exhaustive projection on3072 queries; same tolerance', 'verification_sha256': sha(P / 'NUMERICAL_KERNEL_VERIFICATION_K2.json'), 'aborted_Rtree_cost': json.load(open(P / 'ABORTED_RTREE_EVALUATION.json'))}
    put(P / 'results.json', result)
    with open(P / 'TABLES_R3.csv', 'w', newline='') as f:
        fields = ['tooth_type', 'method', 'teeth', 'patients', 'surface_rms_mm', 'occlusal_p95_mm', 'proximal_patch_p95_mm', 'contact_location_mm']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(r3['tables'])
    with open(P / 'PER_TOOTH_R3.csv', 'w', newline='') as f:
        fields = ['case', 'jaw', 'fdi', 'tooth_type', 'method', 'surface_rms_mm', 'occlusal_p95_mm', 'proximal_patch_p95_mm', 'contact_location_mm']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in map(json.loads, open(P / 'RAW_METRICS_R3.jsonl')):
            for (m, v) in r['metrics'].items():
                w.writerow(dict(case=r['case'], jaw=r['jaw'], fdi=r['fdi'], tooth_type=r['tooth_type'], method=m, **{k: v[k] for k in fields[5:]}))
    text = (P / 'README_DEMO.md').read_text()
    text += "\n## R3 : what an observed cervical tooth stripe adds\n\n"
    d = r3['decision']
    text += f"This is a separate task with more information: the bottom1mm above the natural cervical edge may be visible. Four new persons are separated from all training and R1/R2. {d['completed']}/32 teeth completed. Candidate's upper surface-p95 {d['mean_occlusal_p95_mm']:.3f} mot spegling {d['mirror_mean_occlusal_p95_mm']:.3f}mm, kvot {d['ratio_to_mirror']:.3f}; parad95%-KI [{d['paired_patient_bootstrap95ci_mm'][0]:.3f},{d['paired_patient_bootstrap95ci_mm'][1]:.3f}]mm. Utfall: **{d['decision']}**The results are not mixed with completely removed teeth.\n\n"
    text += "| Tooth type | Method | n | Surface RMS | Upper surface p95 | Neighbour region p95 |\n|---|---|---:|---:|---:|---:|\n"
    for t in r3['tables']:
        if t['tooth_type'] != 'all':
            text += f"| {t['tooth_type']} | {t['method']} | {t['teeth']} | {t['surface_rms_mm']:.3f} | {t['occlusal_p95_mm']:.3f} | {t['proximal_patch_p95_mm']:.3f} |\n"
    text += "\nA natural cervical strip is not a prepared stump. Four individuals may provide a disparate information test but do not carry a broad clinical generalisation. Minutes/split and outputhashar are located in PREREG_R3.json, SPLIT_R3.json and FROZEN_PREDICTIONS_R3.json. The same triangle and error injection controls are used.\n\n"
    text += "Registered antagonist arcs are now available local in Bits2Bites, according to the separate execution LANE_X2_OCCLUSION_B2B. They lack the FDI-tandeletics that this benchmark needs and are second individuals than Teeth3DS. Therefore, no other person's antagonist has been mounted on a Teeth3DS arc. The next concrete data port is an independent dental labeling of such a registered bite; then the same facit operator can measure the antagonist near geometry. Bits2Bites was not used in our numerical trials; no such data was copied.\n"
    if '--dataset-root /data/Teeth3DS' not in text:
        text = text.replace("Reference is published measurement geometry", "On another machine: `./run_all.sh --replay --dataset-root /data/Teeth3DS --data-dir /data/x3_outputs --output-dir /work/x3_demo` . The report directory must be new; each source file must match frozen SHA256 before numericals. The protocol operators, patients, measurements and thresholds are not changed.\n\n" + "Reference is published measurement geometry")
    (P / 'README_DEMO.md').write_text(text)
    with open(P / 'RESULTS.md', 'a') as f:
        f.write('\n## R3 — observed natural collar, separate input regime\n\n' + json.dumps(d, indent=2) + '\n\nSource tooth surfaces are published measurements from4new patients; lower1mm observed explicitly, upper anatomy withheld from fit. FROZEN_PREDICTIONS_R3.json precedes hidden-reference metrics. No clinical preparation/antagonist or new algorithmic superiority is claimed.\n')
    rows = []
    for (rn, rr) in result['rounds'].items():
        rows.append({'round': rn, 'source': 'R21/restorative and direction14', 'capability': json.load(open(P / f'PREREG_{rn}.json'))['capability'], 'changed_operation': json.load(open(P / f'PREREG_{rn}.json'))['changed_operation'], 'outcome': rr['decision'], 'evidence': f'RESULTS_{rn}.json', 'sha256': sha(P / f'RESULTS_{rn}.json')})
    rows.append({'round': 'NUMERICAL_K1', 'outcome': 'FAILED', 'evidence': 'NUMERICAL_KERNEL_FAILED_V1.json', 'next': 'exhaustive-projection comparator at unchanged1e-10mm tolerance'})
    rows.append({'round': 'NUMERICAL_K2', 'outcome': 'PASS', 'evidence': 'NUMERICAL_KERNEL_VERIFICATION_K2.json', 'scope': 'kernel exactness, not crown scientific validation'})
    put(P / 'ATTEMPTS.json', rows)
    feedback = json.load(open(P / 'GRAPH_FEEDBACK.json'))
    feedback.update(sha256=sha(P / 'results.json'), outcome={r: v['decision'] for (r, v) in result['rounds'].items()}, negative_result=any((v['decision']['decision'].startswith('NEGATIVE') for v in result['rounds'].values())))
    put(P / 'GRAPH_FEEDBACK.json', feedback)
    hand = '# X3 final handoff — PENDING_INDEPENDENT_REVIEW\n\nRun `./run_all.sh`; README_DEMO.md is the reader entry. `--replay` reruns the exact protocols with isolated data/output directories. External original IOS surfaces are hash-bound; generated fixture controls are separately named.\n\n'
    for (rn, rr) in result['rounds'].items():
        hand += rn + '\n' + json.dumps(rr['decision'], indent=2) + '\n\n'
    hand += 'Full-deletion R1/R2 and observed-collar R3 are different input regimes and must never be pooled. Generation failures and missing labels stay in scheduled denominators. Raw per-tooth metrics, frozen outputs and patient splits are retained. No thresholds changed after evaluation.\n\nNumerical issue: Rtree near-tie distance switch failed1e-10mm verification on one training query by1.391e-10mm. Exact centroid-radius search agrees with exhaustive projection on3072 queries; failed K1 and aborted slow-evaluator wall cost652.258s remain. This is an implementation repair, not dental algorithm novelty.\n\nMain open prerequisite: real prepared-margin and registered antagonist geometry with FDI labels for the same individual. Existing Bits2Bites registered pairs are a useful source, but X2 identifies missing tooth labels/support. Do not attach another individual to Teeth3DS. Next construction: independently tooth-label a frozen Bits2Bites subset, retain its supplied bite transform, then run leave-one-tooth-out plus antagonist distance-map error and lower/upper contact-patch sensitivity. For clinical crown completion, acquire one lab scan triplet of intact reference / prepared abutment / finished restoration, freeze outer-form/margin predictions before final measurement, and distinguish scanner/pose/CAM error. This measurement is NOT_RUN.\n\nGraph coverage proposal: DENT-VAL-NATURAL-CROWN-LEAVE-ONE-TOOTH-OUT under DENT-DESIGN-PROBLEM and DENT-VAL-BASELINE-COMPARISON. Existing material/fracture crown packet is unsuitable and has unlinked goal authority; no source/native/generated graph mutation, no fabricated dispatch receipt. GRAPH_FEEDBACK.json remains PENDING_INDEPENDENT_REVIEW.\n\nData root ' + str(DATA) + '; max3GB. Current bytes ' + str(disk_guard()) + '. One CPU thread, no GPU, observed peaks below4GB; no heavy run required. Licensing: author repository CC BY-NC-ND4.0, conflicting older local statements retained; no data publication.\n'
    (P / 'HANDOFF.md').write_text(hand)
    fr = json.load(open(P / 'FROZEN_PREDICTIONS_R3.json'))
    assert sha(P / 'PREREG_R3.json') == fr['prereg_sha256']
    assert sha(P / 'SPLIT_R3.json') == fr['split_sha256']
    assert fr['frozen_utc'] < r3['evaluated_utc']
    assert sha(fr['predictions_index']) == fr['index_sha256']
    assert sha(P / 'FROZEN_PREDICTIONS_R3.json') == r3['frozen_predictions_sha256']
    assert sha(P / 'RAW_METRICS_R3.jsonl') == r3['raw_metrics_sha256']
    for v in fr['prediction_files']:
        assert sha(v['path']) == v['sha256']
    orig = json.load(open(P / 'SPLIT.json'))
    new = json.load(open(P / 'SPLIT_R3.json'))
    assert not set(new['cases']) & set(sum(orig['patients'].values(), []))
    put(P / 'VERIFICATION_R3.json', {'verified_utc': now(), 'new_patients_disjoint': True, 'output_hashes_valid': True, 'frozen_before_hidden_reference_evaluation': True, 'different_input_contract_explicit': True, 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
    state('DELIVERY_COMPLETE_PENDING_REVIEW', {r: v['decision'] for (r, v) in result['rounds'].items()}, 'independent raw-data/code review; labeled registered bite or lab preparation triplet required for next scientific capability', active_operation='All three numerical constructions closed; delivery and independent review pending')
    print('R3_BOUND_AND_HANDOFF_WRITTEN', flush=True)
if __name__ == '__main__':
    finish()
