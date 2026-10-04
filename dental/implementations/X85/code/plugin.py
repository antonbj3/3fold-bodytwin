"""GenCAD dictionary -> dictionary plugin. Does not modify the benchmark."""
from common import *
from regional import inverse
from geometry import basis, roof_checks, ball_access
import hashlib

def generate(task):
    cal = task.get('regional_calibration')
    if not cal:
        return dict(task_id=task['task_id'], status='ABSTAIN', reason='UNKNOWN: same-specimen regional actuation calibration absent')
    if cal.get('case') != task['case_key'] or cal.get('frame') != task['frame'] or cal.get('geometry_sha256') != task['geometry_sha256']:
        raise ValueError('Calibration geometry/frame identity mismatch')
    B = np.asarray(task['actuation_basis'])
    expected = task['basis_sha256']
    z = np.asarray(task['unadjusted_z'])
    inner = np.asarray(task['inner_z'])
    allow = z - inner - task['requirements']['wall_mm']
    computed_basis = hashlib.sha256(B.astype('<f8').tobytes()).hexdigest()
    computed_geometry = hashlib.sha256(np.c_[task['xy'], z, inner].astype('<f8').tobytes() + np.asarray(task['faces']).astype('<i8').tobytes()).hexdigest()
    if cal.get('basis_sha256') != expected or expected != computed_basis:
        raise ValueError('Actuation basis identity mismatch')
    if computed_geometry != cal['geometry_sha256']:
        raise ValueError('Actual reference geometry identity mismatch')
    if task.get('height_operation') == 'signed_generation':
        from round7 import inverse_signed
        u = np.array(cal.get('final_coefficient_error_mm', np.zeros(B.shape[1])))
        AB = np.asarray(task['A'] @ B)
        rhs = task['obstacle_b'] - task['A'] @ z - abs(AB) @ u
        ans = inverse_signed(cal, B, task['weights'], allow - B @ u, extra_A=AB, extra_b=rhs)
    else:
        ans = inverse(cal, B, task['weights'], allow)
    if ans['status'] not in ['CONDITIONAL_TARGET_DESIGN', 'CONDITIONAL_SIGNED_TARGET']:
        return dict(task_id=task['task_id'], status='ABSTAIN', reason=ans['status'], diagnostics=ans)
    final = z + B @ np.asarray(ans['coefficients_mm'])
    checks = roof_checks(task, final, inner)
    tool = ball_access(task['xy'], task['faces'], final)
    return dict(task_id=task['task_id'], status='DESIGN', units='mm', frame=task['frame'], outer_vertices=np.c_[task['xy'], final].tolist(), inner_vertices=np.c_[task['xy'], inner].tolist(), faces=np.asarray(task['faces']).tolist(), diagnostics=dict(inverse=ans, geometry=checks, milling=tool, manufacturing_status='CONDITIONAL_ROOF_ONLY' if checks['wall_pass'] and checks['antagonist_pass'] and checks['intaglio_milling_pass'] and (tool['status'] == 'CONDITIONAL_IDEAL_TOOL_PASS') else 'UNKNOWN_OR_FAIL', physical_force='UNKNOWN'))
