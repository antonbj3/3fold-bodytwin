import numpy as np

def observation_gate(pred, observed, tol):
    return bool(np.isfinite([pred, observed, tol]).all() and tol >= 0 and (abs(pred - observed) <= tol))

def reconstruction_tiers(error_mm):
    return {str(t): observation_gate(0, error_mm, t) and error_mm >= 0 for t in [0.35, 0.5, 0.75, 1.0, 1.5]}

def validate_external(task, result):
    return bool(result.get('task_id') == task.get('task_id') and result.get('frame') == task.get('frame') and (result.get('units') == 'mm') and isinstance(result.get('points'), np.ndarray) and (result['points'].ndim == 2) and (result['points'].shape[1] == 3) and (len(result['points']) > 0) and np.isfinite(result['points']).all() and (np.max(np.abs(result['points'])) <= 100))
