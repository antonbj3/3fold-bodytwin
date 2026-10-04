"""Freshly solve a pinned complete bridge deck and compare stored surface tensors."""
import os, sys, subprocess, time, re, shutil
import numpy as np
from common import R, read, dump, sha
sys.path.insert(0, str(R / 'vendor'))
import crown_design_fe as FE

def run():
    tag = 'h4_b2.25_s0.24'
    base = R / 'raw/fe' / tag
    work = R / 'raw/replay'
    work.mkdir(exist_ok=True, parents=True)
    shutil.copyfile(base / 'solve.inp', work / 'solve.inp')
    m = read('INPUT_MANIFEST.json')
    env = dict(os.environ, LD_LIBRARY_PATH=m['solver_library_path'])
    for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'CCX_NPROC_RESULTS', 'CCX_NPROC_EQUATION_SOLVER', 'CCX_NPROC_STIFFNESS']:
        env[k] = '4'
    start = time.monotonic()
    with (work / 'solve.log').open('w') as f:
        p = subprocess.run(['/usr/bin/time', '-v', '-o', str(work / 'time.txt'), m['solver'], '-i', 'solve'], cwd=work, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=240)
    if p.returncode or 'Job finished' not in (work / 'solve.log').read_text():
        raise ValueError('Fresh bridge replay failed; logs retained')
    ga = np.load(base / 'geometry.npz')
    n = len(ga['Vc'])
    S = FE.read_frd_stress(str(work / 'solve.frd'), n)
    Fc = FE.boundary_faces10(ga['Tc'])
    errors = {}
    if len(S) != 2:
        raise ValueError('Expected2freshbridgecases')
    for (name, tensor) in zip(['axial', 'offaxis30'], S):
        sf = -1 / 9 * (tensor[Fc[:, 2]] + tensor[Fc[:, 3]] + tensor[Fc[:, 4]]) + 4 / 9 * (tensor[Fc[:, 5]] + tensor[Fc[:, 6]] + tensor[Fc[:, 7]])
        truth = np.load(base / (name + '_stress.npz'))['Sf']
        errors[name] = float(np.linalg.norm(sf - truth) / np.linalg.norm(truth))
    rss = int(re.search('Maximum resident set size \\(kbytes\\): (\\d+)', (work / 'time.txt').read_text())[1]) * 1024
    out = dict(tag=tag, fresh_solve=True, load_cases=2, seconds=time.monotonic() - start, max_solver_rss_bytes=rss, tensor_relative_L2_errors=errors, gate_pass=all((e <= 0.05 for e in errors.values())), deck_sha256=sha(work / 'solve.inp'), scope='Numerical replay only; parent mesh-convergence and physical calibration failures are unchanged.')
    dump('raw/FRESH_REPLAY.json', out)
    if not out['gate_pass']:
        raise ValueError('Fresh replay differs beyond frozen5% numerical tolerance')
    print('Fresh bridge replay', out['gate_pass'], 'seconds', round(out['seconds'], 2))
    return out
if __name__ == '__main__':
    run()
