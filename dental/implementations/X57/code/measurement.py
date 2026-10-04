"""Decision-sufficiency on declared logical predicates, not lab costs."""
import itertools
PERIO_FIELDS = ['case_definition_confirmed', 'stage_iv_complexity', 'periodontal_tooth_loss_ge5', 'prior_stage_iv']
CARIES_FIELDS = ['ESE_deep_geometry_confirmed', 'pulp_state_compatible', 'irreversible_symptoms_absent', 'actual_pulp_exposure_absent']

def make_worlds(chain):
    fields = PERIO_FIELDS if chain == 'BD10' else CARIES_FIELDS
    worlds = []
    for vals in itertools.product([False, True], repeat=len(fields)):
        if chain == 'BD10':
            outcome = 'NOT_CONFIRMED_PERIODONTITIS_CASE' if not vals[0] else 'IV' if any(vals[1:]) else 'III'
        else:
            outcome = 'ELIGIBLE_RESEARCH_BRANCH' if all(vals) else 'OUTSIDE_RESEARCH_BRANCH'
        worlds.append({'state': dict(zip(fields, vals)), 'outcome': outcome})
    return (fields, worlds)

def collision(worlds, fields):
    seen = {}
    for w in worlds:
        key = tuple((w['state'][f] for f in fields))
        if key in seen and seen[key]['outcome'] != w['outcome']:
            return [seen[key], w]
        seen[key] = w
    return None

def minimum_separators(worlds, fields):
    for n in range(len(fields) + 1):
        result = [list(ss) for ss in itertools.combinations(fields, n) if collision(worlds, ss) is None]
        if result:
            return result
    raise ValueError('World identity does not determine outcome')

def minimum_answer_certificates(worlds, fields, state):
    target = next((w['outcome'] for w in worlds if w['state'] == state))
    for n in range(len(fields) + 1):
        answers = []
        for subset in itertools.combinations(fields, n):
            same = [w for w in worlds if all((w['state'][f] == state[f] for f in subset))]
            if {w['outcome'] for w in same} == {target}:
                answers.append({f: state[f] for f in subset})
        if answers:
            return answers
    raise ValueError('No answer certificate')

def independent_separator_check(worlds, fields):
    unresolved = []
    for (i, w) in enumerate(worlds):
        for v in worlds[i + 1:]:
            if w['outcome'] != v['outcome'] and (not any((w['state'][f] != v['state'][f] for f in fields))):
                unresolved.append([w, v])
    return unresolved

def answer(chain, measured):
    (fields, worlds) = make_worlds(chain)
    if any((k not in fields for k in measured)):
        raise ValueError('Unknown observation predicate')
    if any((type(x) is not bool for x in measured.values())):
        raise ValueError('Predicate must be literal Boolean')
    feasible = [w for w in worlds if all((w['state'][k] == v for (k, v) in measured.items()))]
    outcomes = sorted({w['outcome'] for w in feasible})
    missing = [f for f in fields if f not in measured]
    return {'possible_outcomes': outcomes, 'status': 'DETERMINATE_IN_DECLARED_MODEL' if len(outcomes) == 1 else 'UNKNOWN_NEEDS_MEASUREMENT', 'minimum_additional_predicate_sets': minimum_separators(feasible, missing), 'measurement_count_is_not_physical_cost': True, 'physical_measurement_error': 'UNKNOWN', 'clinical_validation': 'NOT_PERFORMED'}
