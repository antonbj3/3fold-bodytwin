import datetime, json, platform, resource, time
from common import *

def figure(net, stress, r2):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    (fig, axs) = plt.subplots(1, 3, figsize=(15, 5), gridspec_kw={'width_ratios': [1.6, 1, 1.4]})
    top = stress['all_variables_sorted'][:10]
    labels = [x['variable'].replace('_', ' ') for x in top][::-1]
    axs[0].barh(labels, [x['decision_count'] for x in top][::-1], color='#32647d')
    axs[0].set_xlabel('Distinct K consumers with OPEN/UNKNOWN edge')
    axs[0].set_title('Where evidence is binding')
    axs[0].text(0, -0.27, 'Direct curated dependencies; no effect-size ranking', transform=axs[0].transAxes, fontsize=8)
    counts = stress['status_counts']
    labs = ['TIGHT', 'OPEN', 'UNKNOWN']
    axs[1].bar(labs, [counts[k] for k in labs], color=['#3b8266', '#d69737', '#8b9298'])
    axs[1].set_title('78 scoped relations')
    axs[1].set_ylabel('Relations')
    axs[1].text(0, -0.27, 'TIGHT: scoped external comparison,\nnot validated patient prediction', transform=axs[1].transAxes, fontsize=8)
    cases = [r2['examples'][i] for i in [1, 2, 3, 4, 5]]
    axs[2].barh(['M01', 'M02 only', 'M01 + M02', 'M01 + M02 + M03', 'M01 + M02 + M04'], [c['count'] for c in cases], color='#6b6696')
    axs[2].set_title('AND dependencies change the answer')
    axs[2].set_xlabel('Structurally answerable X71 questions', fontsize=9)
    axs[2].text(0, -0.27, 'No physical measurement performed;\nM03 leaves hydraulic-model blocker open', transform=axs[2].transAxes, fontsize=8)
    fig.suptitle('Dental constraint network — evidence scope and next measurements', fontsize=15)
    fig.tight_layout(rect=[0, 0.1, 1, 0.95])
    fig.savefig(HERE / 'CONSTRAINT_NET_DEMO.png', dpi=160)
    fig.savefig(HERE / 'CONSTRAINT_NET_DEMO.pdf')
    plt.close(fig)

def main():
    start = time.perf_counter()
    net = read(HERE / 'CONSTRAINT_NET_DENTAL.json')
    s = read(HERE / 'CONSTRAINT_STRESS_MAP_DENTAL.json')
    r2 = read(HERE / 'RESULTS_R2.json')
    v = read(HERE / 'VALIDATION.json')
    figure(net, s, r2)
    dump('COST_RUN.json', {'R2_seconds': r2['runtime_seconds'], 'demo_seconds': time.perf_counter() - start, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'threads_max': 1, 'GPU_seconds': 0, 'recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'scope': 'Latest successful replay only; preparation, failed attempts and reasoning not included'})
    binding = read(HERE / 'BINDINGS_RECEIPT.json')
    suff = read(HERE / 'SUFFICIENCY_R1.json')
    result = {'lane': 'PROOF_LANE-constraint-net', 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'capability': 'Query every K01–K52 chain through typed quantities, scoped quantified relations, exact evidence and named measurement blockers', 'outcome': 'RUNNABLE_SCOPED_CONSTRAINT_NET; PHYSICAL_TRANSFER_OPEN', 'external_referent': {'kind': 'external_review', 'locator': 'REVIEW_CENSUS.json; exact source paths and hashes in SOURCE_MANIFEST.json', 'compared_quantity': 'Accepted review scopes and original numerical observations/contrasts; each TIGHT edge additionally names its primary published measurement/dataset', 'refutes_us': True}, 'network': {'variables': len(net['variables']), 'edges': len(net['edges']), 'chains': len(net['chains']), 'status_counts': s['status_counts'], 'count_resolution': 'PHENOMENOLOGICAL', 'count_scope': 'Exact counts within curated contract; completeness of all dental relations is UNKNOWN', 'replacement_measurement': 'Independent curation audit of omitted reviewed relation candidates'}, 'evidence_checks': v, 'current_source_checks': 'VALIDATION_LIVE_VALUES.json', 'whole_live_hash_failure': 'VALIDATION_LIVE_R1_ATTEMPT1_FAILED.json', 'sufficiency': suff, 'round2': {'portfolios': r2['exhaustive']['portfolios'], 'mismatches': r2['exhaustive']['set_vs_bitmask_mismatches'], 'decision_contracts': 64, 'physical_measurements': 0, 'physical_certificates': 0}, 'stress_top': s['top_stress_points'], 'x71_comparison': {k: z for (k, z) in s['x71_comparison'].items() if k != 'mapped_measurements'}, 'dropout': net['dropout'], 'binding': binding, 'negative_results': ['Scalar unresolved-edge/consumer counts insufficient: exact identity error0, downstream difference1', 'Whole current WORKING_VIEW SHA differs from frozen snapshot; exact queried values and reviewed result hashes are checked separately, full-hash failure retained', 'No physical end-to-end K chain validated by this network', 'X71 frozen input only; no finished priority ordering available for agreement test', 'Optical/cure branch and full rheology are outside this core network'], 'protocol_deviations': ['Formal DECOMPOSITION.json was saved after structural R1/R2 calculations. It is retrospective, not claimed preregistered. PREREG_R1 and PREREG_R2 were frozen before their computations.'], 'full_cost': {'preparation': 'Targeted local curation and source snapshots; exact historical preparation wall time UNKNOWN', 'fit': 0, 'discovery': 'Manual quantity/scope binding; human/agent reasoning time and tokens UNKNOWN', 'validation': 'Schema; exact evidence; review hashes;13 corruptions;8 source contrasts;32768 dependency combinations', 'questions': '52 K-chain and X71 named measurement queries', 'fallback': 'Keep UNKNOWN, request matched measurement contract through handoff; no lab or heavy run', 'execution_receipt': 'COST_RUN.json', 'threads_max': 1, 'GPU_seconds': 0, 'lab_cost': 'NOT_PERFORMED', 'independent_review_cost': 'UNKNOWN'}, 'large_arrays_over_50MB': [], 'intermediate_byte_limit': 3000000000, 'artifacts': {name: sha(HERE / name) for name in ['CONSTRAINT_NET_DENTAL.json', 'CONSTRAINT_STRESS_MAP_DENTAL.json', 'CONSTRAINT_NET_DEMO.png', 'CONSTRAINT_NET_DENTAL.schema.json', 'SOURCE_MANIFEST.json', 'CURATION_CONTRACT.json']}, 'next_construction': 'Type optical/rheology missing ports; independent review of full edge curation; join X71 final priority under matched specimen/protocol observation contracts', 'clinical_recommendation': False}
    dump('results.json', result)
    print(json.dumps({'outcome': result['outcome'], 'checks_pass': v['pass'] and v['mutation_tests']['pass'] and r2['checks_pass'], 'network': result['network'], 'figure': 'CONSTRAINT_NET_DEMO.png'}))
if __name__ == '__main__':
    main()
