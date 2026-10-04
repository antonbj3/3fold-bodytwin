from dental_release.paths import expand as _release_expand
import json, hashlib, datetime
from pathlib import Path
H = Path(__file__).resolve().parents[1]
p = json.loads((H / 'PREREG_R1.json').read_text())
ids = sorted([int(x.name) for x in Path(_release_expand('@DENTAL_WORK_ROOT@/X11/targets/Bits2Bites')).iterdir() if x.is_dir()], key=lambda x: hashlib.sha256(('X18:' + str(x)).encode()).hexdigest())
p.update(round='R2', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), parents=['R1 contact-recovery FAIL, continuous nonpenetration PASS'], test_cases=ids[20:32], changed_operation='Replace vertex clipping+global translation by exact local polygon-vertex inequalities; solve a morphology-regularized convex roof obstacle with optional cusp lift', lift_candidates_mm=[0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0], morphology_regularizer=0.1, dual_maxiter=2000, dual_gtol=1e-09, selection=dict(rule='First feasible lift candidate with at least1 patch, at most4 patches, contact area0.2..5mm2. If none, select first nonempty contact map, else zero lift. Contact-recovery decisions unchanged.', patch_min=1, patch_max=4, area_min_mm2=0.2, area_max_mm2=5.0, status='CONSTITUTIVE_CLOSURE, not a published clinical norm'), strongest_equally_informed_control='L1 minimum morphology-change linear program over same full triangle-intersection constraints, donor and selected lift. No algorithm-superiority claim; actual costs and scores retained.', decision='Same contact-recovery and nonpenetration thresholds as R1 on12 new case IDs; target FDI validity and physical stress remain UNKNOWN', cost={**p['cost'], 'questions': 'QP optimization and LP equal-information check charged separately; no original contacts for selection'})
p['external_referent']['locator'] += ';12 new untouched target contact cases'
(H / 'PREREG_R2.json').write_text(json.dumps(p, indent=2) + '\n')
d = json.loads((H / 'DECOMPOSITION_R1.json').read_text())
d['changed_equation'] = 'For every projected triangle polygon vertex q: sum_i w_i(q) z_i <= z_U(q)-0.025mm. Minimize (z-z_prior)^T(I+0.1 L^T L)(z-z_prior)/2; z_prior=z_donor+lift.'
d['leaves'].append(dict(leaf='contact-target selection', status='CONSTITUTIVE_CLOSURE', relation='First antagonist-obstacle candidate satisfying prereg patch/area constraints', stopping_argument='No physical norm enforces this exact count/area; external original contacts can refute'))
(H / 'DECOMPOSITION_R2.json').write_text(json.dumps(d, indent=2) + '\n')
print(p['test_cases'])
