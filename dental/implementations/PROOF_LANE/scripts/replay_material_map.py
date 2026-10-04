"""Execute the existing O replay with only its final output path redirected.

The scientific code is not duplicated. The AST IO change is narrow and checked;
the original source and its inputs are immutable. No predecessor state changes.
"""
import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[k] = '4'
import ast, sys, json, hashlib
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gencad_bench.io import dump, sha, DENTAL
SOURCE = DENTAL / 'results/LANE_NEXT_O_MATERIAL_FEASIBILITY/code/replay.py'
dest = ROOT / 'raw/material_map_replay.json'
tree = ast.parse(SOURCE.read_text())
changed = []

class Redirect(ast.NodeTransformer):

    def visit_Call(self, node):
        self.generic_visit(node)
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'write_text':
            old = ast.unparse(node.func.value)
            if old != "R / 'raw/independent_replay.json'":
                raise ValueError('Unexpected writer ' + old)
            node.func.value = ast.Call(func=ast.Name(id='Path', ctx=ast.Load()), args=[ast.Constant(str(dest))], keywords=[])
            changed.append(old)
        return node
tree = Redirect().visit(tree)
ast.fix_missing_locations(tree)
assert len(changed) == 1
ns = {'__file__': str(SOURCE), '__name__': 'gencad_inherited_material_replay'}
exec(compile(tree, str(SOURCE), 'exec'), ns)
x = json.loads(dest.read_text())
dump(ROOT / 'raw/material_map_reuse_receipt.json', dict(source_path=str(SOURCE), source_sha256=sha(SOURCE), redirected_writer=changed, scientific_operations_changed=False, input_json_sha256=sha(SOURCE.parents[1] / 'raw/r1_results.json'), result_sha256=sha(dest), witnesses_replayed=x['witnesses_replayed'], failures=x['failures'], constructors_replayed=x['constructors_replayed'], destructive_radius_controls_rejected=x['destructive_radius_controls_rejected'], arithmetic='Inherited floating distance bounds with 1e-7 mm checks; not newly formalized', scope='K/O digital necessary inclusion; raw tooth/pulp anatomy and clinical truth retain original limits'))
