from fc_common import *
import collections

def run():
    replay = read(ROOT / 'raw/REGENERATION_RECEIPT.json')
    demo = read(ROOT / 'raw/DEMO_RUN.json')
    assert replay['all_pass'] and demo['regenerated']
    result = read(ROOT / 'results.json')
    borrowed = []
    wf = read(V4 / 'FROZEN_WHOLE_INPUTS.json')
    for (kind, folder) in [('public_files', 'whole_inputs'), ('private_files', 'whole_private')]:
        for (rel, info) in wf[kind].items():
            path = V4 / 'payload' / folder / rel
            ok = sha(path) == info['sha256']
            borrowed.append(dict(path=str(path), sha256=sha(path), pass_gate=ok))
            assert ok, 'borrowed input drift'
    for (path, expected) in read(ROOT / 'FROZEN_SCORER.json')['files'].items():
        assert sha(path) == expected, 'scorer drift'
    dump(ROOT / 'raw/BORROWED_INPUT_INTEGRITY.json', dict(all_pass=True, files=borrowed, source_scorer_unchanged=True))
    attempts = []
    for tag in [f'R{i}' for i in range(1, 8)]:
        pr = read(ROOT / f'PREREG_{tag}.json')
        rr = read(ROOT / 'rounds' / f'{tag}.json')
        attempts.append(dict(round=tag, claim_type='capability', parent='PROOF_LANE_GENCAD_V4 / PROOF_LANE_GENCAD_V5 whole-crown frontier', changed_operation=pr['changed_operation'], prereg_sha256=sha(ROOT / f'PREREG_{tag}.json'), results_file=f'rounds/{tag}.json', outcome='ANATOMICAL_GATE_FAILED', negative_result=True, next_operation=pr.get('full_cost', {}).get('fallback', result['next_construction'])))
    dump(ROOT / 'ATTEMPTS.json', attempts)
    docs = ['README_DEMO.md', 'LAB_PROTOCOL.md', 'PLUGIN_PROTOCOL.md', 'NEXT_INPUT_CONTRACT.json', 'LIMITATIONS_SUPPLEMENT.json', 'GRAPH_COVERAGE_PROPOSAL.json', 'LAB_EXPORT_ENCLOSURES.json', 'FROZEN_EXPORT_ENCLOSURES.json', 'FROZEN_LAB_EXPORTS.json', 'FROZEN_PREDICTIONS.json', 'RELEASE_CODE.json', 'FROZEN_SCORER.json']
    freeze(ROOT / 'REVIEW_SNAPSHOT.json', dict(claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', results=result, results_sha256=sha(ROOT / 'results.json'), full_replay=replay, demo_run=demo, borrowed_integrity_sha256=sha(ROOT / 'raw/BORROWED_INPUT_INTEGRITY.json'), supplement=read(ROOT / 'LIMITATIONS_SUPPLEMENT.json'), actual_export_enclosures=read(ROOT / 'LAB_EXPORT_ENCLOSURES.json'), documents={name: sha(ROOT / name) for name in docs}, supplemental_code_sha256={p.name: sha(p) for p in [ROOT / 'code/quantization_enclosure.py', ROOT / 'code/finalize.py']}, scientific_admission=False))
    rsha = sha(ROOT / 'REVIEW_SNAPSHOT.json')
    feedback = dict(target_id='DENT-VAL-BASELINE-COMPARISON', result_file='results/PROOF_LANE_FULL_CROWN/REVIEW_SNAPSHOT.json', sha256=rsha, review_state='PENDING_INDEPENDENT_REVIEW', claim_type='capability', outcome='ANATOMICAL_GATE_NOT_BROKEN; COMPLETE_CONDITIONAL_RESEARCH_PAIRS_AND_SOURCE_SUPPORT_DIAGNOSIS', measured_quantity='Whole-exterior bidirectional sampledp95; spatial contact/interference; source/virtual facet support; conditional all-exterior-facet wall bound; digital export quantization', units='mm; mm2; counts and fractions', uncertainty=result['uncertainty'], population_regime='Retrospective6 case clusters,3 inferred tooth types;R2/R3 hidden-upper-form track andR4–R7 full-prescan/co-designed-preparation track kept separate', preregistered_gate=[f'PREREG_R{i}.json' for i in range(1, 8)] + ['PREREG_R4B.json', 'PREREG_VALIDATION.json', 'PREREG_EXPORTS.json'], baseline='Actual unchangedv4/v5 failures, rigid and affine measured homolog controls, known-source height envelope oracle. No algorithm-superiority claim.', negative_result=True, dispatch_scope='review; numerical capability-specific definition missing', resolution=result['resolution_schema'], edge_contracts=result['edge_contracts'], external_referent=REFERENT, missing_coverage_proposal='GRAPH_COVERAGE_PROPOSAL.json and NEXT_INPUT_CONTRACT.json', scientific_admission=False)
    dump(ROOT / 'GRAPH_FEEDBACK.json', feedback)
    n = sum((r.get('status') != 'FAILED' for r in replay['rows']))
    failed = sum((r.get('status') == 'FAILED' for r in replay['rows']))
    with (ROOT / 'HANDOFF.md').open('a') as f:
        f.write(f"\nFull replay PASS: {n} regenerated mesh outputs exactly identical; {failed} R6/R7 construction failures reproduced plus fullR2/R3 status lists (including6 R3 abstentions). Total actual-generator replay{replay['seconds']:.3f}s; full demo{demo['seconds']:.3f}s. Peak demo parentRSS{demo['peak_rss_MiB']:.3f}MiB. Canonical immutable review snapshot SHA256 {rsha}. Do not bind future mutable results.json to this receipt; snapshot embeds current results.\n")
    state('ROUND_COMPLETE_PENDING_INDEPENDENT_REVIEW', 'All anatomical0.35mm gates failed; actual generators reproduced exactly; selected research3MF exports valid', 'Use NEXT_INPUT_CONTRACT: independently annotated source/margin and real pre/post prep, then patient-disjoint full-surface prior. Retain all failed gates and information-track separation.')
    print(json.dumps(dict(snapshot_sha256=rsha, identical_mesh_outputs=n, reproduced_R6_R7_failures=failed, full_demo_seconds=demo['seconds'])))
if __name__ == '__main__':
    run()
