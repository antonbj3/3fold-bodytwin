"""Three-unit geometry and fresh CalculiX solve; absolute brittle strength is absent.

Reuses X1B's actual C3D10 crown/die tie, normalized loads and surface stress operator.
The repeated tooth is a CAD fixture, never claimed to be the published bridge geometry.
"""
import os, sys, time, subprocess, re
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates, label
import trimesh
from common import R, read, dump, sha, state
sys.path.insert(0, str(R / 'vendor'))
import crown_design_fe as FE
from crown_fit_geometry import Grid, surface_mesh
FE.GMSH_SCRIPT = str(R / 'vendor/crown_design_gmsh.py')

def build(height, width, spacing):
    a = np.load(R / 'inputs/D1_model.npz')
    g = read('inputs/D1_grid.json')
    old = Grid(g['origin'], np.array(g['origin']) + g['h'] * (np.array(g['shape']) - 1), g['h'])
    grid = Grid([-18.8, -7.8, -3.1], [18.8, 7.8, 8.5], spacing)
    (X, Y, Z) = grid.mesh()
    P = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    shape = grid.shape

    def shifted(field, s):
        q = P.copy()
        q[:, 0] -= s
        return old.sample(a[field], q, cval=12.0).reshape(shape).astype('float32')
    c = np.minimum(shifted('crown', -11.0), shifted('crown', 11.0))
    pontic = np.maximum(shifted('outer', 0.0), g['z_m'] - Z)
    c = np.minimum(c, pontic)
    for s in [-5.5, 5.5]:
        Q = np.stack([np.abs(X - s) - 3.5, np.abs(Y) - width / 2, np.abs(Z - 3.4) - height / 2], axis=-1)
        box = np.linalg.norm(np.maximum(Q, 0.0), axis=-1) + np.minimum(np.max(Q, axis=-1), 0.0)
        c = np.minimum(c, box)
    die = np.minimum(shifted('die', -11.0), shifted('die', 11.0))
    ncc = int(label(c < 0)[1])
    ndie = int(label(die < 0)[1])
    if ncc != 1 or ndie != 2:
        raise ValueError(f'Geometry components bridge{ncc}/dies{ndie}, expected1/2')
    V = a['dV']
    T = a['dT']
    v = np.vstack([V + [-11.0, 0, 0], V + [11.0, 0, 0]])
    t = np.vstack([T, T + len(V)])
    return (c.astype('float32'), die.astype('float32'), grid, (v, t), g['z_m'])

def solve(height, width, spacing, tag):
    w = R / 'raw/fe' / tag
    w.mkdir(exist_ok=True, parents=True)
    start = time.monotonic()
    (c, d, grid, dieVT, z_m) = build(height, width, spacing)
    (Vs, Fs) = surface_mesh(c, grid)
    surf = trimesh.Trimesh(Vs, Fs, process=False)
    surf.merge_vertices(digits_vertex=5)
    surf.update_faces(surf.unique_faces())
    surf.update_faces(surf.nondegenerate_faces(height=1e-06))
    surf.remove_unreferenced_vertices()
    (Vs, Fs) = (surf.vertices, surf.faces)
    target_faces = 16000 if spacing >= 0.24 else 24000
    if len(surf.faces) > target_faces:
        surf = surf.simplify_quadric_decimation(face_count=target_faces, aggression=5)
        surf.remove_unreferenced_vertices()
        (Vs, Fs) = (surf.vertices, surf.faces)
    if not surf.is_watertight:
        raise ValueError('Bridge STL is not watertight')
    surf.export(w / 'bridge.stl')
    FE.tet_mesh(str(w / 'bridge.stl'), str(w / 'mesh.npz'), 0.65, 0.15)
    a = np.load(w / 'mesh.npz')
    m = FE.Model((a['V'], a['T']), dieVT, c, d, grid, z_m, z_m - 3.0, tie_tol=0.18)
    if not len(m.fix_nodes) or not len(m.Fc_tied) or (not len(m.Fd_master)):
        raise ValueError('Missing bridge end supports/tie')
    m._outer_pts = m._outer_pts[np.abs(m._outer_pts[:, 0]) < 4.0]
    cases = []
    diags = []
    for (name, direction) in [('axial', [0, 0, -1]), ('offaxis30', [0, 0.5, -np.sqrt(3) / 2])]:
        (loads, diag) = FE.load_case(m, dict(kind='fossa' if name == 'axial' else 'directed', R=2.5, d=direction, center=[0, 0, 5.0], a_fixed=0.75, sticking=True), 210000.0, 0.27)
        err = float(np.linalg.norm(np.array(diag['total_force_per_N']) - direction))
        if err > 1e-05:
            raise ValueError('Normalized force is not1N: ' + str(err))
        cases.append((name, loads))
        diags.append(diag)
    FE.write_deck(str(w / 'solve.inp'), m, (210000.0, 0.27), (200000.0, 0.3), cases)
    deck = (w / 'solve.inp').read_text().splitlines()
    start_nodes = deck.index('*NODE, NSET=NALL') + 1
    allV = np.vstack([m.Vc, m.Vd])
    deck[start_nodes:start_nodes + len(allV)] = [f'{i + 1},{x:.12f},{y:.12f},{z:.12f}' for (i, (x, y, z)) in enumerate(allV)]
    (w / 'solve.inp').write_text('\n'.join(deck) + '\n')
    manifest = read('INPUT_MANIFEST.json')
    solver = Path(manifest['solver'])
    locked = [x for x in manifest['items'] if x['local'] is None][0]
    if sha(solver) != locked['sha256']:
        raise ValueError('Solver changed')
    env = dict(os.environ, LD_LIBRARY_PATH=manifest['solver_library_path'])
    for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'CCX_NPROC_RESULTS', 'CCX_NPROC_EQUATION_SOLVER', 'CCX_NPROC_STIFFNESS']:
        env[k] = '4'
    cmd = ['/usr/bin/time', '-v', '-o', str(w / 'time.txt'), str(solver), '-i', 'solve']
    with (w / 'solve.log').open('w') as f:
        p = subprocess.run(cmd, cwd=w, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=240)
    if p.returncode or 'Job finished' not in (w / 'solve.log').read_text():
        raise ValueError('Solver failed: preserved ' + str(w / 'solve.log'))
    ga = FE.model_geometry_arrays(m)
    stresses = FE.read_frd_stress(str(w / 'solve.frd'), m.nc + len(m.Vd))
    out = []
    if len(stresses) != 2:
        raise ValueError('Expected2 load cases')
    np.savez_compressed(w / 'geometry.npz', cf=ga['cf'], nf=ga['nf'], af=ga['af'], tied=ga['tied'], Vc=m.Vc, Tc=m.Tc)
    for ((name, _), S, diag) in zip(cases, stresses, diags):
        (Sf, Se) = FE.face_elem_stress(m, S)
        s = FE.tangential_sigma1(Sf.astype(float), ga['nf'])
        cf = ga['cf']
        connector = np.abs(np.abs(cf[:, 0]) - 5.5) < 1.25
        contact = np.zeros(len(cf), bool)
        for p in diag['patches']:
            contact |= (np.linalg.norm(cf - p['p'], axis=1) < 3 * p['a_mm']) & (np.abs(cf[:, 0]) < 4.0)
        regions = {'connector_left': connector & (cf[:, 0] < 0), 'connector_right': connector & (cf[:, 0] > 0), 'pontic_contact': contact, 'pontic_other': (np.abs(cf[:, 0]) < 4.0) & ~contact, 'end_crown': (np.abs(cf[:, 0]) >= 4.0) & ~connector}
        regional = {}
        for (region, mask) in regions.items():
            if not mask.any():
                raise ValueError('Empty surface region ' + region)
            regional[region] = dict(max_MPa_per_N=float(s[mask].max()), p99_MPa_per_N=float(np.quantile(s[mask], 0.99)), area_mm2=float(ga['af'][mask].sum()), point_max_mm=cf[mask][np.argmax(s[mask])].tolist(), resolution='PER_SURFACE_REGION', mode_closure='tangential surface tensile stress; excludes contact crack mechanics')
        np.savez_compressed(w / (name + '_stress.npz'), Sf=Sf, s1=s)
        out.append(dict(name=name, load=diag, regions=regional, absolute_fracture_N=None, absolute_status='UNKNOWN_STRENGTH_AND_CONTACT_DAMAGE'))
    sdf_err = float(np.quantile(np.abs(grid.sample(c, Vs)), 0.99))
    txt = (w / 'time.txt').read_text()
    rss = int(re.search('Maximum resident set size \\(kbytes\\): (\\d+)', txt)[1]) * 1024
    result = dict(tag=tag, height_mm=height, width_mm=width, nominal_area_mm2=height * width, field_spacing_mm=spacing, nodes=m.nc + len(m.Vd), tet10=len(m.Tc) + len(m.Td), load_cases=out, surface_p99_error_mm=sdf_err, watertight=True, bridge_components=1, die_components=2, fix_nodes=len(m.fix_nodes), mesh_cleanup_merge_mm=1e-05, deck_coordinate_decimals=12, surface_target_faces=target_faces, missed_tie_node_warnings=(w / 'solve.log').read_text().count('no tied MPC'), source_geometry='X1B STS tooth D1 outer/shell/die fields; repeated tooth prototype, not original study', seconds=time.monotonic() - start, max_solver_rss_bytes=rss, hashes={p.name: sha(p) for p in [w / 'bridge.stl', w / 'mesh.npz', w / 'solve.inp', w / 'solve.frd', w / 'geometry.npz']})
    dump(f'raw/fe/{tag}/FE.json', result)
    print(tag, 'nodes', result['nodes'], 'RSS_MB', rss / 1000000.0, 'seconds', round(result['seconds'], 2), flush=True)
    return result

def run():
    state('R2_fresh_bridge_FE', 'running; absolute fracture capacity UNKNOWN', 'solve two9mm2 shapes at two field spacings; retain stress fields and mesh convergence gate')
    runs = []
    for (h, b) in [(4.0, 2.25), (3.0, 3.0)]:
        for spacing in [0.24, 0.18]:
            tag = f'h{h:g}_b{b:g}_s{spacing:g}'
            runs.append(solve(h, b, spacing, tag))
            state('R2_fresh_bridge_FE', tag + ' solved', 'continue remaining mesh/shape solve')
    checks = []
    for (h, b) in [(4.0, 2.25), (3.0, 3.0)]:
        (old, new) = [x for x in runs if x['height_mm'] == h]
        for (i, lc) in enumerate(['axial', 'offaxis30']):
            for region in old['load_cases'][i]['regions']:
                a = old['load_cases'][i]['regions'][region]['p99_MPa_per_N']
                z = new['load_cases'][i]['regions'][region]['p99_MPa_per_N']
                e = abs(a - z) / z
                checks.append(dict(height_mm=h, load_case=lc, region=region, relative_p99_difference=e, pass_gate=e <= 0.05))
    out = dict(runs=runs, refinement_checks=checks, refinement_gate_pass=all((c['pass_gate'] for c in checks)), note='Grid/surface/tetrahedral refinement together; perfect ties and fixed load patch remain closures. Regionalp99 diagnostic; does not certify singular maxima.')
    dump('raw/BRIDGE_FE.json', out)
    return out
if __name__ == '__main__':
    run()
