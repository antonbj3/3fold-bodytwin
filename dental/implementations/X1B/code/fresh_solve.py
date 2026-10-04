"""Actually rerun the vendored CPU solver on one unchanged anatomical deck."""
import os, sys, shutil, subprocess, time, re
from pathlib import Path
import numpy as np
from common import R, read, dump, sha
sys.path.insert(0, str(R / 'vendor'))
import crown_design_fe as FE
from crown_fit_geometry import Grid

def solve(work):
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    sid = 'L005_lo_M1_k6_thin_3Y'
    shutil.copy2(R / f'inputs/fe/{sid}.inp', work / 'solve.inp')
    env = dict(os.environ, LD_LIBRARY_PATH=str(R / 'vendor/solver/lib'))
    for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'CCX_NPROC_RESULTS', 'CCX_NPROC_EQUATION_SOLVER', 'CCX_NPROC_STIFFNESS']:
        env[k] = '4'
    start = time.monotonic()
    cmd = ['/usr/bin/time', '-v', '-o', str(work / 'time.txt'), str(R / 'vendor/solver/ccx'), '-i', 'solve']
    with (work / 'solve.log').open('w') as f:
        p = subprocess.run(cmd, cwd=work, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=120)
    if p.returncode or 'Job finished' not in (work / 'solve.log').read_text():
        raise ValueError('Fresh CPU solve failed; preserved solve.log')
    a = np.load(R / 'inputs/geometry/D1_model.npz')
    g = read('inputs/geometry/D1_grid.json')
    grid = Grid(g['origin'], np.array(g['origin']) + g['h'] * (np.array(g['shape']) - 1), g['h'])
    m = FE.Model((a['cV'], a['cT']), (a['dV'], a['dT']), a['crown'], a['die'], grid, g['z_m'], g['z_m'] - 3)
    stresses = FE.read_frd_stress(str(work / 'solve.frd'), m.nc + len(m.Vd))
    old = np.load(R / f'inputs/fe/{sid}.npz')
    errors = {}
    if len(stresses) != 2:
        raise ValueError('Expected two fresh load cases')
    for (name, S) in zip(['axial', 'offaxis30'], stresses):
        (Sf, Se) = FE.face_elem_stress(m, S)
        truth = old[f'Sf__{name}']
        errors[name] = float(np.linalg.norm(Sf - truth) / np.linalg.norm(truth))
    tm = (work / 'time.txt').read_text()
    rss = int(re.search('Maximum resident set size \\(kbytes\\): (\\d+)', tm)[1]) * 1024
    artifacts = [dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size) for p in [work / 'solve.inp', work / 'solve.frd', work / 'solve.log', work / 'time.txt']]
    result = dict(fresh_solve=True, solver='CalculiX2.23 MKL PARDISO vendored', deck_sha256=sha(work / 'solve.inp'), solver_sha256=sha(R / 'vendor/solver/ccx'), load_cases=2, relative_tensor_L2_errors=errors, gate_pass=all((e <= 0.05 for e in errors.values())), seconds=time.monotonic() - start, max_rss_bytes=rss, artifacts=artifacts, scope='Fresh solve of unchanged saved geometry/deck; no new CBCT segmentation or mesh generation; independent physical validation remains absent')
    dump('raw/FRESH_SOLVE.json', result)
    if not result['gate_pass']:
        raise ValueError('Fresh tensor response exceeds frozen numerical gate')
    return result
