"""Composition rules. Missing independent calibration is never filled from another case."""
import math

def join(producer, consumer):
    fields = ('case_key', 'geometry_sha256', 'frame', 'unit', 'quantity', 'resolution', 'observation_state', 'time_scale')
    if any((not producer.get(k) or not consumer.get(k) for k in fields)):
        return {'status': 'UNKNOWN', 'reason': 'Incomplete dimensioned port'}
    mismatch = [k for k in fields if producer[k] != consumer[k]]
    if mismatch:
        return {'status': 'FAIL', 'reason': 'Incompatible port', 'mismatches': mismatch}
    if producer.get('epistemic') not in ('EXTERNALLY_CALIBRATED', 'EXTERNALLY_MEASURED', 'PROVEN_UNDER_EXPLICIT_ASSUMPTIONS'):
        return {'status': 'UNKNOWN', 'reason': 'Required calibration/proof is absent'}
    return {'status': 'PASS', 'reason': 'Declared ports match; empirical scope still belongs to producer'}

def field_query(crown, graph=None, model=None):
    if graph is None:
        return dict(status='UNKNOWN', classification='UNKNOWN', force_interval_N=None, reason='No complete same-case, same-edited-mesh opposing-tooth gap graph; do not infer NEVER from missing edges')
    for k in ('case_key', 'geometry_sha256', 'frame'):
        if graph.get(k) != crown.get(k):
            raise ValueError('Field graph binding mismatch: ' + k)
    if not graph.get('complete_declared_pair_graph'):
        return dict(status='UNKNOWN', classification='UNKNOWN', force_interval_N=None, reason='Candidate graph coverage missing')
    answer = model.classify(graph['edges'], [crown['fdi']])
    return dict(status='CONDITIONAL_MODEL', **answer)

def compare_observation(frozen, obs):
    required = ('key', 'design_mesh_sha256', 'quantity', 'unit', 'frame', 'measurement_state')
    for k in required:
        if obs.get(k) != frozen.get(k):
            return dict(status='REJECTED_BINDING', field=k)
    if not obs.get('raw_measurement_sha256') or not obs.get('calibration_locator'):
        return dict(status='UNKNOWN_MISSING_MEASUREMENT_PROVENANCE')
    try:
        value = float(obs['value'])
        error = float(obs['absolute_error_bound'])
        if not math.isfinite(value) or not math.isfinite(error) or error < 0:
            raise ValueError()
    except (KeyError, ValueError, TypeError):
        return dict(status='REJECTED_MEASUREMENT_VALUE')
    pred = frozen.get('prediction_interval')
    if pred is None:
        return dict(status='UNKNOWN_NO_PHYSICAL_PREDICTION', observation_interval=[value - error, value + error])
    (lo, hi) = map(float, pred)
    if value - error > hi or value + error < lo:
        return dict(status='ALARM_DISJOINT_INTERVALS', observed=[value - error, value + error], predicted=[lo, hi])
    return dict(status='CONSISTENT_AT_DECLARED_ERRORS', observed=[value - error, value + error], predicted=[lo, hi], scientific_admission=False)
