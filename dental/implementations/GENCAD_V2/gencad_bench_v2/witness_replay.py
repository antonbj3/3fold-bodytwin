"""Replay inherited tissue/material impossibility witnesses; redirect the sole writer locally."""
import ast, sys, time, collections
from pathlib import Path
from .common import *

def run():
    source = DENTAL / 'results/LANE_NEXT_O_MATERIAL_FEASIBILITY/code/replay.py'
    input_path = source.parents[1] / 'raw/r1_results.json'
    dest = ROOT / 'raw/MATERIAL_WITNESS_REPLAY.json'
    freeze(ROOT / 'PREREG_WITNESS_REPLAY.json', dict(claim_type='capability', capability='Keep the prior LS2/tissue impossibility instrument executable alongside new IOS tasks', operation='Run existing independently written replay on unchanged voxel boxes; redirect only output', source_sha256=sha(source), input_sha256=sha(input_path), gate='All 224 inherited witnesses and two constructive controls recheck; destructive radius controls reject; no new clinical inference', external_referent=dict(kind='published_code', locator=str(source), compared_quantity='digital set-inclusion and cone/voxel necessary conditions, not new pulp measurements', refutes_us=True)))
    tree = ast.parse(source.read_text())
    changed = []

    class Redirect(ast.NodeTransformer):

        def visit_Call(self, node):
            self.generic_visit(node)
            if isinstance(node.func, ast.Attribute) and node.func.attr == 'write_text':
                old = ast.unparse(node.func.value)
                if old != "R / 'raw/independent_replay.json'":
                    raise ValueError('Unexpected source writer ' + old)
                node.func.value = ast.Call(func=ast.Name(id='Path', ctx=ast.Load()), args=[ast.Constant(str(dest))], keywords=[])
                changed.append(old)
            return node
    tree = Redirect().visit(tree)
    ast.fix_missing_locations(tree)
    if len(changed) != 1:
        raise ValueError('Source writer count changed')
    ns = {'__file__': str(source), '__name__': 'gencad_v2_readonly_material_instrument'}
    exec(compile(tree, str(source), 'exec'), ns)
    r = read(dest)
    receipt = dict(source_path=str(source), source_sha256=sha(source), input_sha256=sha(input_path), result_sha256=sha(dest), witnesses=r['witnesses_replayed'], constructors=r['constructors_replayed'], destructive_controls_rejected=r['destructive_radius_controls_rejected'], by_material=dict(collections.Counter((a['material'] for a in r['rows']))), failures=r['failures'], seconds=r['wall_s'], scope='inherited voxel proof instrument, not a patient join, not newly formalized or a new external measurement')
    if r['failures'] or r['destructive_radius_controls_rejected'] != r['constructors_replayed']:
        raise ValueError('material witness replay failed')
    dump(ROOT / 'raw/MATERIAL_REPLAY_RECEIPT.json', receipt)
    return receipt
if __name__ == '__main__':
    print(run())
