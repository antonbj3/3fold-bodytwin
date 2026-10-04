"""Use inherited continuous projected triangle operator without its CLI imports.

AST selects the original geometry functions verbatim. Normal-coordinate
triangles verify the affine patch gap only: no full-aligner collision claim.
"""
import ast, hashlib, json
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1]
SOURCE = R.parents[1] / 'results/LANE_NEXT_P_OCCLUSION_VALIDATION/code/continuous_gap.py'

def functions():
    tree = ast.parse(SOURCE.read_text())
    selected = ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ['cross', 'plane', 'inside', 'exact_batch']], type_ignores=[])
    scope = {'np': np}
    exec(compile(selected, str(SOURCE), 'exec'), scope)
    return scope

def projected_gap(gap, scope):
    tooth = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    shell = tooth.copy()
    shell[:, :, 2] += gap
    up = scope['plane'](shell)
    lp = scope['plane'](tooth)
    return scope['exact_batch'](shell, tooth, up, lp, np.array([[0, 0]]))[0]

def main():
    scope = functions()
    pr = json.loads((R / 'PREREG_R1.json').read_text())
    rows = json.loads((R / 'raw/PREDICTIONS_R1.json').read_text())['rows']
    checked = []
    for row in rows:
        for region in row['regional']['regional']:
            gap = region['gap_mm']
            g = projected_gap(gap, scope)
            checked.append(dict(case=row['case'], active_fdi=row['active_fdi'], contact_fdi=region['fdi'], region=region['region'], model_gap_mm=gap, inherited_operator_gap_mm=g, absolute_difference_mm=abs(g - gap), source_kind='our_own_fixture'))
    injected = projected_gap(-0.01, scope)
    result = dict(source=str(SOURCE), sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(), checks=checked, maximum_gap_difference_mm=max((r['absolute_difference_mm'] for r in checked)), corruption_gap_mm=injected, corruption_rejected=max(0, -injected) > pr['gates']['nonpenetration_mm'], scope='Continuous projected triangle verification in each patch normal frame; fixture-only, fixed affine gap; no anatomical 3D collision certificate')
    (R / 'raw/INHERITED_GAP_CHECK.json').write_text(json.dumps(result, indent=2) + '\n')
    if not result['corruption_rejected'] or result['maximum_gap_difference_mm'] > 1e-12:
        raise RuntimeError('Inherited gap check failed')
    print('Inherited projected-gap operator:', len(checked), 'patch checks, corruption rejected')
if __name__ == '__main__':
    main()
