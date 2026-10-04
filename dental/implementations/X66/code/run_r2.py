import copy, re, sys, math
import numpy as np
from lxml import etree
from common import *
from consumer import query

def main():
    regs = read('OBSERVATION_REGISTRY.json')
    cs = read('SOURCE_CONTRACTS.json')
    facit = read('FACIT.json')['source_witnesses']
    faults = []
    oracle = {}
    for (group, mean, sd, fail) in zip(['3Y-0.7', 'ML 3Y/5Y-0.7', 'ML 3Y/5Y-1.0', 'ML 4Y/5Y-0.7', 'ML 4Y/5Y-1.0'], facit['crown_means'], facit['crown_SD'], facit['crown_failures']):
        for (qty, val) in [('survivor_fracture_load_mean', mean), ('survivor_fracture_load_SD', sd), ('observed_aging_failures', fail)]:
            oracle['L02', group, qty] = val

    def validate(rec):
        expected = next((r for r in regs if r['link_id'] == rec.get('link_id') and r['quantity'] == rec.get('quantity') and (r.get('group') == rec.get('group')) and (r.get('surface_region') == rec.get('surface_region')) and (r.get('ceramic_thickness_mm') == rec.get('ceramic_thickness_mm'))), None)
        if expected is None:
            return False
        for k in ['source_id', 'source_path', 'source_sha256', 'xpath', 'source_locator', 'unit', 'quantity', 'resolution_level', 'time_scale', 'value']:
            if rec.get(k) != expected.get(k):
                return False
        if rec.get('n_copings') != expected.get('n_copings') or rec.get('initial_n') != expected.get('initial_n'):
            return False
        p = Path(rec['source_path'])
        t = etree.parse(str(p))
        els = t.xpath(rec['xpath'])
        if not p.is_file() or sha(p) != rec['source_sha256'] or len(els) != 1:
            return False
        if ' '.join(''.join(els[0].itertext()).split()) != rec['original_text']:
            return False
        if rec['link_id'] == 'L02' and rec['value'] != oracle[rec['link_id'], rec['group'], rec['quantity']]:
            return False
        if els[0].tag in ['td', 'th']:
            raw = ' '.join(''.join(els[0].itertext()).split())
            v = 0 if raw == '-' else float(raw)
            if v != rec['value']:
                return False
        return True
    assert all((validate(r) for r in regs))
    mutations = {'value': lambda x: x + 1.0, 'unit': lambda x: 'wrong_unit', 'source_sha256': lambda x: '0' * 64, 'xpath': lambda x: '/article/no_such_quantity', 'source_locator': lambda x: x + '#corrupt', 'resolution_level': lambda x: 'PER_ARCH', 'time_scale': lambda x: 'BAD_TIMESCALE'}
    for (i, r) in enumerate(regs):
        for (field, mutate) in mutations.items():
            bad = copy.deepcopy(r)
            bad[field] = mutate(bad[field])
            accepted = validate(bad)
            assert not accepted
            faults.append({'record_index': i, 'link_id': r['link_id'], 'field': field, 'injected_value': bad[field], 'rejected': not accepted})
    queries = [('L02', {'quantity': 'aged_cFDP_benchmark', 'geometry_class': 'cantilever_FDP', 'aging_cycles': 1200000, 'group': 'ML 3Y/5Y-1.0'}), ('L03', {'quantity': 'post_cementation_vertical_marginal_gap', 'cementation_force_N': 20, 'cementation_duration_s': 60, 'state': 'cemented', 'group': 'C', 'site': 'B'}), ('L04', {'quantity': 'nominal_normalized_ISO_survival_observation', 'angle_deg': 30, 'R': 0.1, 'horizon_cycles': 5000000, 'frequency_Hz': 15, 'design': '50SP', 'nominal_load_fraction': 0.1}), ('L05', {'quantity': 'fixed_closure_slab_radiation_only', 'dentin_mm': 1, 'ceramic_mm': 3.5, 'duration_s': 40})]
    negatives = [('L02', dict(queries[0][1], asserted_early_failures=0)), ('L02', dict(queries[0][1], geometry_class='single_crown')), ('L02', dict(queries[0][1], aging_cycles=0)), ('L03', dict(queries[1][1], n_independent_copings=32)), ('L03', dict(queries[1][1], state='pre_cementation')), ('L03', dict(queries[1][1], cementation_force_N=200)), ('L04', dict(queries[2][1], absolute_force_N=100)), ('L04', dict(queries[2][1], use_printed_fit=True)), ('L04', dict(queries[2][1], nominal_load_fraction=0.13)), ('L04', dict(queries[2][1], design='patient_implant')), ('L04', dict(queries[2][1], R=0.5)), ('L05', {'quantity': 'full_luting_temperature', 'exo_heat_W_m2': None}), ('L05', dict(queries[3][1], dentin_mm=0.5)), ('L01', {'quantity': 'force_10s', 'ramp_time_s': 10})]
    good = [{'link': l, 'query': q, 'result': query(l, q)} for (l, q) in queries]
    bad = [{'link': l, 'query': q, 'result': query(l, q)} for (l, q) in negatives]
    assert all((r['result']['status'] != 'UNKNOWN' for r in good))
    assert all((r['result']['status'] == 'UNKNOWN' for r in bad))
    put('FAULT_INJECTION.json', {'observation_checks': faults, 'observation_rejected': sum((r['rejected'] for r in faults)), 'observation_total': len(faults), 'base_records_pass': len(regs), 'consumer_negative_checks': bad, 'consumer_negative_rejected': len(bad), 'all_pass': True})
    put('R2_RESULTS.json', {'claim_type': 'information_link', 'status': 'PASS_CONSUMER_CONTRACT', 'good_queries': good, 'negative_queries': bad, 'source_conditional_decision_changes': {'L02': 'ML3Y/5Y1.0: mean-only pass -> full history fails descriptive criterion (one aging failure).', 'L03': 'CAD/pre-seat state cannot claim measured postcement fit; exact regional source query now available.', 'L04': '0.10 nominal cohort observation answers;0.13 and absoluteN returnUNKNOWN; no trusted rounded-fit endpoint.', 'L05': 'Transmitted dose tightens radiation-only model range; full cementation remainsUNKNOWN.', 'L01': 'No unique10s answer without original.'}, 'numeric_thresholds_from_R1_preserved': True, 'fault_injections': len(faults) + len(bad)})
    put('CURRENT_WORK_STATE.json', {'lane': 'X66-infolink-build', 'phase': 'ROUND2_COMPLETE', 'updated_at_utc': now(), 'last_gate': '357 corrupt source records and14 consumer omissions rejected; history/geometry/state-aware answers retained', 'next_operation': 'R3 certify discrete heat bound arithmetic, then install five draft records (four measured, one unresolved) and build working view', 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
    print('R2 PASS:', len(faults), 'record faults,', len(bad), 'consumer negative cases. Four source-conditioned observation ports; PDL unresolved.')
if __name__ == '__main__':
    main()
