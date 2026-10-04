"""Bind each decision's support to the patient/frame of its hashed theta registry."""

def validate_identity(row, registry):
    support = row.get('support', {})
    for key in ('patient_id', 'frame_id'):
        if not registry.get(key) or support.get(key) != registry[key]:
            raise ValueError('DECISION_IDENTITY_MISMATCH:' + key)
    if not row.get('id', '').startswith(registry['patient_id'] + ':'):
        raise ValueError('DECISION_ID_PATIENT_MISMATCH')
    return True

def force_balance(fe, expected_force):
    """Do not silently substitute zero for a missing solver result."""
    import math
    for answer in fe.values():
        if answer.get('status') == 'UNKNOWN_NO_CONTACT':
            continue
        if answer.get('status') != 'SIMULATED':
            return False
        for key in ('force_balance_N', 'total_load_N'):
            if key not in answer or not math.isfinite(answer[key]):
                return False
        if abs(answer['force_balance_N']) >= 1e-06:
            return False
        if abs(answer['total_load_N'] - expected_force) >= 1e-06:
            return False
    return True
