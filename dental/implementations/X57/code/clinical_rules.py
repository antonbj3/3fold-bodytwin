"""Conditional research predicates, with exact images on interval cells.

These functions do not diagnose disease or recommend treatment. Source-error
bounds are unknown. An interval is a stated observation scenario, never a CI.
"""
import itertools
STAGES = ['I', 'II', 'III', 'IV']

def validate_observation(o):
    q = o['quantity']
    expected = 'mm' if q in ['PD', 'CAL'] else 'root-length fraction semantic scenario' if q == 'RBL' else 'furcation grade' if q == 'FURCATION' else 'category'
    if o['unit'] != expected:
        raise ValueError('Observation has wrong unit for ' + q)
    if q == 'FURCATION' and (o['interval'] is None or any((v not in [1, 2, 3] for v in o['interval']))):
        raise ValueError('Invalid furcation grade')
    if o.get('timescale') != 'SIMULTANEOUS':
        raise ValueError('Wrong time contract')

def severity_rbl(value):
    if not 0 <= value <= 1:
        raise ValueError('RBL must be fraction of root length in [0,1]')
    return 1 if value < 0.15 else 2 if value <= 0.33 else 3

def severity_cal(value):
    if value < 0 or value > 20:
        raise ValueError('CAL mm outside declared domain')
    return 1 if value < 3 else 2 if value < 5 else 3

def severity_pd(value):
    if not 0 <= value <= 20:
        raise ValueError('PD mm outside declared domain')
    return 3 if value >= 6 else 1

def interval_image(interval, function, breaks):
    (lo, hi) = interval
    if lo > hi:
        raise ValueError('Reversed observation interval')
    cuts = sorted(set([lo, hi] + [b for b in breaks if lo <= b <= hi]))
    samples = cuts + [(a + b) / 2 for (a, b) in zip(cuts, cuts[1:])]
    return sorted(set((function(v) for v in samples)))

def periodontal(observations):
    minima = []
    supports = []
    for o in observations:
        if o.get('accepted') is not True or o.get('chain') != 'BD10':
            continue
        validate_observation(o)
        q = o['quantity']
        if q == 'PD':
            values = interval_image(o['interval'], severity_pd, [6.0])
        elif q == 'CAL':
            values = interval_image(o['interval'], severity_cal, [3.0, 5.0])
        elif q == 'RBL':
            values = interval_image(o['interval'], severity_rbl, [0.15, 0.33])
        elif q == 'FURCATION':
            values = [3] if o['interval'][0] >= 2 else [1]
        elif q == 'BONE_LOSS_REPORTED':
            values = [1]
        else:
            continue
        minima.append(min(values))
        supports.append(o['observation_id'])
    lower = max(minima, default=1)
    return {'conditional_stages': STAGES[lower - 1:], 'conditional_on': ['periodontitis case definition and etiology established', 'journal observations valid at relevant visit'], 'advanced_threshold_definitely_reported': lower >= 3, 'complete_stage': 'UNKNOWN', 'clinical_decision': 'MEASURE_CASE_DEFINITION_AND_STAGE_MODIFIERS', 'missing': ['nonadjacent site CAL with CEJ', 'exclude endodontic/fracture/third-molar causes', 'tooth-loss attribution', 'stage-IV complexity and prior treatment'], 'supports': supports, 'resolution': 'PER_ARCH', 'input_resolution': 'PER_SURFACE_REGION where stated, else PER_TOOTH or PER_ARCH', 'timescale': 'SIMULTANEOUS', 'uncertainty': 'Logical image conditional on text; physical/clinical measurement bounds UNKNOWN'}

def lesion_class(depth_fraction, barrier_present):
    if not 0 <= depth_fraction <= 1:
        raise ValueError('Lesion depth fraction outside [0,1]')
    if depth_fraction == 1 and (not barrier_present):
        return 'EXTREMELY_DEEP'
    if depth_fraction >= 0.75 and barrier_present:
        return 'DEEP'
    return 'OTHER_OR_INCONSISTENT'

def branch_eligible(deep_class, vital, irreversible_symptoms, actual_exposure):
    return deep_class == 'DEEP' and vital and (not irreversible_symptoms) and (not actual_exposure)

def caries(observations):
    relations = [o for o in observations if o.get('accepted') is True and o.get('chain') == 'BD08']
    lesion_sets = []
    for o in relations:
        validate_observation(o)
        possible = ['OTHER_OR_INCONSISTENT', 'DEEP', 'EXTREMELY_DEEP']
        if o['quantity'] == 'PULP_PENETRATION_REPORTED':
            possible = ['EXTREMELY_DEEP']
        lesion_sets.append({'observation_id': o['observation_id'], 'teeth': o['teeth'], 'conditional_lesion_classes': possible, 'resolution': o['resolution'], 'conditional_on': 'literal journal caries-front observation correct', 'actual_pulp_exposure': 'UNKNOWN', 'pulp_vitality': 'UNKNOWN'})
    return {'lesions': lesion_sets, 'pulp_penetration_definitely_reported': any((o['quantity'] == 'PULP_PENETRATION_REPORTED' for o in relations)), 'selective_excavation_eligibility': 'UNKNOWN', 'clinical_decision': 'MEASURE_LESION_BARRIER_AND_PULP_STATE', 'missing': ['lesion depth versus full dentin thickness', 'visible hard/firm barrier', 'tooth-specific pulp state and symptom duration', 'actual exposure and restorative feasibility'], 'resolution': 'PER_TOOTH', 'timescale': 'SIMULTANEOUS', 'cold_response_is_vitality': False, 'uncertainty': 'Narrative class is not anatomical, histological or therapeutic certification'}

def control_periodontal(observations):
    """Separately coded enumeration on source table, including unknown upgrades."""
    compatible = set(range(1, 5))
    for o in observations:
        if not o.get('accepted') or o.get('chain') != 'BD10':
            continue
        (low, high) = o['interval'] if o['interval'] is not None else (0, 0)
        q = o['quantity']
        if q == 'PD':
            bounds = [(0, 6, 1), (6, 20, 3)]
            required = 3 if low >= bounds[1][0] else 1
        elif q == 'CAL':
            required = next((s for (threshold, s) in [(5, 3), (3, 2), (0, 1)] if low >= threshold), 1)
        elif q == 'RBL':
            required = 1
            if low >= 0.15:
                required = 2
            if low > 0.33:
                required = 3
        elif q == 'FURCATION':
            required = 3 if low in [2, 3] else 1
        elif q == 'BONE_LOSS_REPORTED':
            required = 1
        else:
            continue
        compatible &= set(range(required, 5))
    return [STAGES[i - 1] for i in sorted(compatible)]

def control_caries(observations):
    output = []
    for o in observations:
        if not o.get('accepted') or o.get('chain') != 'BD08':
            continue
        classes = ['EXTREMELY_DEEP'] if o['quantity'] == 'PULP_PENETRATION_REPORTED' else ['OTHER_OR_INCONSISTENT', 'DEEP', 'EXTREMELY_DEEP']
        outcomes = sorted(set((c == 'DEEP' and v and (not s) and (not e) for (c, v, s, e) in itertools.product(classes, [False, True], [False, True], [False, True]))))
        output.append({'observation_id': o['observation_id'], 'classes': classes, 'eligibility_image': outcomes})
    return output
