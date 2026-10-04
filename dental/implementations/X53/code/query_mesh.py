"""Single-command tool query on a local mesh, optionally including fixture STL."""
import argparse
from common import *
from collision import *
from run_r1 import lib
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mesh')
    ap.add_argument('--point', nargs=3, type=float, required=True)
    ap.add_argument('--normal', nargs=3, type=float, required=True)
    ap.add_argument('--axes', type=int, choices=[4, 5], default=5)
    ap.add_argument('--library', choices=['vhf', 'ceramill'], default='vhf')
    ap.add_argument('--scale', type=float, default=0.8)
    ap.add_argument('--fixture')
    ap.add_argument('--tool-card', help='Measured tool library JSON with same schema; replaces scenario library')
    ap.add_argument('--output', default='QUERY.json')
    args = ap.parse_args()
    m = trimesh.load_mesh(args.mesh, process=True)
    if not m.is_watertight or not m.is_winding_consistent or m.volume <= 0:
        raise ValueError('mesh must be a closed oriented positive solid')
    if args.fixture:
        m = trimesh.util.concatenate([m, trimesh.load_mesh(args.fixture, process=True)])
    n = np.array(args.normal)
    ln = np.linalg.norm(n)
    if ln <= 0:
        raise ValueError('zero normal')
    n /= ln
    p = np.array(args.point)
    s = Scene(m.vertices, m.faces)
    if args.tool_card:
        tools = read(args.tool_card)['tools']
        for t in tools:
            for k in ['diameter_mm', 'neck_reach_mm', 'shank_mm', 'gauge_mm', 'holder_diameter_mm', 'holder_length_mm']:
                t[k] *= args.scale
    else:
        tools = lib(args.library, args.axes, args.scale)
    rows = [dict(tool_id=t['id'], nominal_diameter_mm=t['diameter_mm'] / args.scale, **search(s, p, n, t, args.axes, offsets=(0, 0.05, 0.1, 0.2, 0.4))) for t in tools]
    found = [x for x in rows if x['status'] == 'FOUND' and x['offset_mm'] == 0]
    dump(args.output, dict(claim_type='capability', resolution='PER_POINT', mesh_sha256=sha(args.mesh), fixture_sha256=sha(args.fixture) if args.fixture else None, scale=args.scale, point=p, normal=n, rows=rows, smallest_found_mm=min((x['nominal_diameter_mm'] for x in found), default=None), status='CONDITIONAL_DECLARED_ASSEMBLY_AND_POSE_GRID', full_machine_status='UNKNOWN_WORKSPACE_AND_CONTINUOUS_POSE_COVERAGE'))
