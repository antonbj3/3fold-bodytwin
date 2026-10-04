from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r2 = json.load(open(ROOT / 'rounds/R2.json'))
r3 = json.load(open(ROOT / 'rounds/R3.json'))
c = [r for r in r2['contrasts'] if r['control'] == 'control_learned' and r['family'] != 'all']
(fig, ax) = plt.subplots(1, 2, figsize=(11, 5))
for (a, k, label) in [(ax[0], 'rmse_mm', 'Native-height ΔRMSE (mm)'), (ax[1], 'contact_error_mm2', 'Static proximity-area Δerror (mm²)')]:
    y = np.arange(len(c))
    d = [x['paired'][k] for x in c]
    v = np.array([x['mean_difference'] if x else np.nan for x in d])
    lo = np.array([x['case_bootstrap_95'][0] if x else np.nan for x in d])
    hi = np.array([x['case_bootstrap_95'][1] if x else np.nan for x in d])
    a.errorbar(v, y, xerr=[v - lo, hi - v], fmt='o', color='#1b9e77')
    a.axvline(0, color='black', linewidth=0.6)
    a.set_yticks(y, [x['family'] for x in c] if a is ax[0] else [])
    a.set_xlabel(label)
    a.set_title('Kernel minus equal-data RF; negative improves')
fig.suptitle('DEVELOPMENT ONLY:104 cases; no heldout query')
fig.tight_layout()
fig.savefig(ROOT / 'figures/development.png', dpi=160)
plt.close(fig)
dump(ROOT / 'results.json', dict(claim_type='algorithm', status='FROZEN_TEST_PREDICTIONS_IN_PROGRESS', test_evaluation='NOT_RUN', hidden_queries=0, review_state='PENDING_INDEPENDENT_REVIEW', external_referent=json.load(open(ROOT / 'PREREG_R1.json'))['external_referent'], development={'R1_gate': 'FAIL', 'R2_gate': 'FAIL', 'R3_gate': 'FAIL', 'R2_strict_pass': 2900, 'R2_tasks': 3744, 'R2_pareto_vs_original_LP': True, 'R2_pareto_vs_equal_data_RF': False, 'R3_equal_data_contact_comparison': r3['contrasts'][-1]}, generation_failure={'path': 'raw/GENERATION_FAILED_01.log', 'public_file_hash_now_matches_release': True, 'cause': 'UNKNOWN', 'predictions_preserved': 279}, candidate_freeze_sha256=sha(ROOT / 'FROZEN_GENERATOR.json'), direct_conventional_control_freeze_sha256=sha(ROOT / 'FROZEN_CONTROL.json'), figure='figures/development.png', field_kernel='NOT_RUN'))
state('RESUMING_FROZEN_PREDICTIONS', {'all_development_parent_gates': 'FAIL', 'first_generation_completed': 279, 'failure_preserved': True}, 'Complete same candidate and frozen controls; freeze NPZ set; one trusted heldout query', hidden_queries=0)
(ROOT / 'RESULTS.md').write_text("Sourceed crown generator and separate contact field are built. They can be run with the same exact limitations as V3 and export nine development shells. Native IOS from published Bite2Text/Bits2Bites are reference observations; the test measurement is still NOT_RUN.\n\nR1:645/864 approved, the same as: LP; contact better, height errors worse. R2:2900/3744 approved, the same as: LP; height/contact better than original LP but the height gain against equally informed RF Insecure.3:separate contact field also helps RF; full the validation gate falls. Precise intervals in the rounds/ and the development figure figures/development.png.\n\nFrozen candidate and direct full-output kernel control is unchanged. A public file could temporarily not open after279 the case; the error and NPZThe prefix has been preserved, but the reboot is being replaced, and no hidden test question has been exhausted.\n\nFull milling access, approximate contacts, felling and measured mirror control are missing in the locked interface. PENDING_INDEPENDENT_REVIEW.\n")
(ROOT / 'HANDOFF.md').write_text('Checkpoint before first heldout measurement.\n\nAll3 development constructions retain their original failed parent gates. R3 primary is x42_contact_branch beta1, selected from development only. Frozen generator SHA' + sha(ROOT / 'FROZEN_GENERATOR.json') + '; extra equally informed full-output kernel control SHA' + sha(ROOT / 'FROZEN_CONTROL.json') + ". All208 dev donors refit;95/1872 donor family rows excluded. Test references have not been parsed by development/generation.\n\nOriginal isolated generation failed after279 persisted cases due a transient missing public scene; current SHA matches release; cause UNKNOWN. raw/GENERATION_FAILED_01.log, ISOLATION_FAILED_01.json, PARTIAL_PREDICTIONS_279.json preserve it. Resume unchanged participant on verified one-case snapshots is active; exact command/logs in COMMANDS.md. Finish pending717 case files, then generate additional full-output kernel control, freeze all12participant NPZs and use run_all.sh --first-measurement once. No query ledger yet. Do not tune after measurement or overwrite any freeze.\n\nNext construction after evaluation must change fixed preparation/full source-supported boundary in a new benchmark version. V3 roof pass is not a clinical restoration. Field port is documented; 9 CPU research shells exported, native GPU invocation NOT_RUN . Graph dispatch refused numerical_dispatch_not_allowed and missing best_baseline; no graph source/generated mutation. Missing generative-crown experiment coverage needs an explicit definition for review.\n")
print('development artifact and checkpoint saved')
