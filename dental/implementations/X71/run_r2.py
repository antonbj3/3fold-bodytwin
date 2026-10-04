"""Conservative portfolios and exact matched enumeration over frozen assumptions."""
import itertools, json, resource, time
from fractions import Fraction
from planning import *
start = time.perf_counter()
verify_lock()
reg = read('PREREG_R2.json')
if reg['input_lock_sha256'] != sha(HERE / 'INPUT_LOCK.json'):
    raise ValueError('PREREG_INPUT_DRIFT')
ds = read('inputs/DECISIONS.json')
ms = read('inputs/MEASUREMENTS.json')
ass = read('inputs/LAB_ASSUMPTIONS.json')
allp = list(portfolios(ds, ms))

def ref_ready(ids):
    ports = {s for m in ms if m['id'] in ids for s in m['supplies']}
    return sorted((d['id'] for d in ds if d.get('count_in_value', True) and d['requires'] and all((s in ports for s in d['requires'])) and (len(d['other_blockers']) == 0)))
control_errors = []
for p in allp:
    if sorted(p['evaluation_ready']) != ref_ready(p['ids']):
        control_errors.append(p['ids'])

def rankkey(p):
    return (-len(p['evaluation_ready']), p['new_specimens'], p['operator_hours'][1], p['ids'])

def feasible(p, cap, idx, manufacture=False):
    op = p['operator_hours'][idx] + (ass['shared_manufacture_operator_hours'][idx] if manufacture and p['ids'] else 0)
    machine = p['instrument_hours'][idx] + (ass['shared_manufacture_machine_hours'][idx] if manufacture and p['ids'] else 0)
    return op <= cap and machine <= cap

def select(cap, idx, manufacture=False):
    available = [p for p in allp if feasible(p, cap, idx, manufacture)]
    best = sorted(available, key=rankkey)[0]
    ties = [p for p in available if len(p['evaluation_ready']) == len(best['evaluation_ready'])]
    return dict(selected=best, max_count_ties=len(ties), ties=sorted(ties, key=rankkey), feasible_portfolios=len(available))
robust = select(40, 1)
optimistic = select(40, 0)
scratch = select(40, 1, True)
x40 = next((p for p in allp if p['ids'] == ['M01', 'M02', 'M04']))

def acceptance(p):
    return feasible(p, 40, 1) and len(p['evaluation_ready']) == len(ref_ready(p['ids']))
nominal = acceptance(robust['selected'])
budget_mutant = dict(robust['selected'])
budget_mutant['operator_hours'] = [999, 999]
budget_rejected = not acceptance(budget_mutant)
assert not admissible(['M02'], ms)
mutant = dict(robust['selected'])
mutant['evaluation_ready'] = list(mutant['evaluation_ready']) + ['D05']
dependency_rejected = not acceptance(mutant)
assert 'D20' not in evaluate(['M04'], ds, ms)['evaluation_ready']
rank = []
for m in ms:
    ids = {m['id']}
    while True:
        add = {x for mm in ms if mm['id'] in ids for x in mm['prerequisites']}
        if add <= ids:
            break
        ids |= add
    p = evaluate(ids, ds, ms)
    standalone = evaluate([m['id']], ds, ms)
    p['measurement'] = m['id']
    p['name'] = m['name']
    p['standalone_prereq_valid'] = admissible([m['id']], ms)
    p['individual_ready'] = standalone['evaluation_ready'] if p['standalone_prereq_valid'] else []
    p['value_per_operator_hour_interval'] = [len(p['evaluation_ready']) / p['operator_hours'][1], len(p['evaluation_ready']) / p['operator_hours'][0]]
    p['value_per_operator_hour_exact_enclosure'] = [str(Fraction(len(p['evaluation_ready']), p['operator_hours'][1])), str(Fraction(len(p['evaluation_ready']), p['operator_hours'][0]))]
    p['value_scope'] = 'Potential evaluation of measured specimen/held port; actual decisions0 before outcomes'
    rank.append(p)
rank.sort(key=lambda p: (-p['value_per_operator_hour_interval'][0], p['new_specimens'], p['ids']))
capacity = [dict(capacity_h=h, robust=select(h, 1)['selected'], optimistic=select(h, 0)['selected']) for h in range(20, 61)]
cost_order = []
for (a, b) in itertools.combinations(rank, 2):
    (la, ha) = map(Fraction, a['value_per_operator_hour_exact_enclosure'])
    (lb, hb) = map(Fraction, b['value_per_operator_hour_exact_enclosure'])
    cost_order.append(dict(a=a['measurement'], b=b['measurement'], a_dominates_all_costs=la > hb, b_dominates_all_costs=lb > ha, status='ORDER_DEPENDS_ON_COST' if not (la > hb or lb > ha) else 'STRICT_BOX_ORDER'))
out = dict(round='R2', claim_type='information_link', outcome='CONDITIONAL_ACQUISITION_PLAN_WITH_UNMEASURED_COSTS', robust_40h=robust, optimistic_40h=optimistic, from_scratch_40h=scratch, conventional_X40_pilot_prefix=x40, concordant_with_X40=robust['selected']['ids'] == x40['ids'], same_information_control=dict(portfolios_checked=len(allp), readiness_mismatches=control_errors, nominal_accepts=nominal, budget_live_mutation_rejected=budget_rejected, missing_prerequisite_rejected=not admissible(['M02'], ms), injected_false_ready_rejected=dependency_rejected, incompatible_foursupport_credit_rejected='D20' not in evaluate(['M04'], ds, ms)['evaluation_ready']), individual_packages=rank, cost_order_enclosures=cost_order, capacity_scenarios=capacity, physical_measurements_performed=0, actually_new_decisions=0, numerical_enclosure='Exact finite Boolean enumeration, integer endpoint sums and rational count/cost ratios over stated additive interval closure; no physical cost bound', physical_enclosure='MISSING; calibration/physical model truth and threshold crossing require matched laboratory observations', full_cost=dict(runtime_s=time.perf_counter() - start, max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, physical_prepare='NOT_RUN; ready specimens UNKNOWN', fit='NONE', discovery_labor='UNKNOWN', validation='controls only, not physical', money='UNKNOWN', fallback='new paired/held batches outside present ready-specimen plan'))
write('raw/RESULTS_R2.json', out)
write('HANDOFF_R2.json', dict(outcome=out['outcome'], next_operation='R3 destructive specimen-state and independent observation intake; test full72 fracture versus retention allocation', selected=robust['selected'], largest_obstacle='Costs/access/readiness unmeasured; a measurement set still does not record irreversible specimen state'))
if control_errors or not all([nominal, budget_rejected, dependency_rejected]):
    raise SystemExit('CONTROL_FAILURE')
print('robust', robust['selected']['ids'], 'ready', len(robust['selected']['evaluation_ready']), 'operator_h', robust['selected']['operator_hours'])
print('optimistic', optimistic['selected']['ids'], 'ready', len(optimistic['selected']['evaluation_ready']))
print('full chains', len(robust['selected']['complete_chains']), 'control mismatches', len(control_errors))
checkpoint('R2_COMPLETE', 'ALL_FINITE_CONTROL_PORTFOLIOS_MATCH', 'R3: ordered specimen ownership and prospective record intake')
