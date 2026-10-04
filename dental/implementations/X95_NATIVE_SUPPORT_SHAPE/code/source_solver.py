from common import *
import sys, subprocess, shutil, re
sys.path.insert(0, str(X1 / 'vendor'))
import crown_design_fe as FE
from native_port import model

def run():
    verify_inputs()
    t0 = time.perf_counter()
    work = DATA / 'NATIVE_SOLVER'
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy2(X1 / 'inputs/fe/L005_lo_M1_k6_thin_3Y.inp', work / 'solve.inp')
    env = dict(os.environ, LD_LIBRARY_PATH=str(X1 / 'vendor/solver/lib'))
    for key in ['OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'CCX_NPROC_RESULTS', 'CCX_NPROC_EQUATION_SOLVER', 'CCX_NPROC_STIFFNESS']:
        env[key] = '4'
    cmd = ['/usr/bin/time', '-v', '-o', str(work / 'time.txt'), str(X1 / 'vendor/solver/ccx'), '-i', 'solve']
    with (work / 'solve.log').open('w') as f:
        proc = subprocess.run(cmd, cwd=work, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=900)
    log = (work / 'solve.log').read_text()
    out = {'command': cmd, 'returncode': proc.returncode, 'job_finished': 'Job finished' in log, 'cost': cost(t0), 'warning_lines': [l for l in log.splitlines() if any((w in l.lower() for w in ['warning', 'tie', 'slave', 'untied']))]}
    if proc.returncode == 0 and out['job_finished']:
        (a, g, m) = model()
        stresses = FE.read_frd_stress(str(work / 'solve.frd'), m.nc + len(m.Vd))
        old = np.load(X1 / 'inputs/fe/L005_lo_M1_k6_thin_3Y.npz')
        errors = {}
        source_forces = []
        for (name, S) in zip(['axial', 'offaxis30'], stresses):
            (Sf, Se) = FE.face_elem_stress(m, S)
            truth = old[f'Sf__{name}']
            errors[name] = float(np.linalg.norm(Sf - truth) / np.linalg.norm(truth))
        out['relative_tensor_L2_error'] = errors
        out['native_numeric_parity_pass'] = len(stresses) == 2 and all((e <= 0.05 for e in errors.values()))
        out['source_solver_peak_RSS_MiB'] = int(re.search('Maximum resident set size \\(kbytes\\): (\\d+)', (work / 'time.txt').read_text())[1]) / 1024
    else:
        out['native_numeric_parity_pass'] = False
    out['artifacts'] = [{'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in work.iterdir() if p.is_file()]
    dump(R / 'raw/NATIVE_SOLVER.json', out)
    print(json.dumps(clean(out), indent=2))

def cached_reference():
    verify_inputs()
    t0 = time.perf_counter()
    lock = read(SRC / 'SOURCE_NATIVE_REFERENCE_LOCK.json')
    if sha(lock['native_result_path']) != lock['native_result_sha256']:
        raise ValueError('Native reference result drift')
    out = read(lock['native_result_path'])
    if not lock['actual_source_solve_executed'] or not out['native_numeric_parity_pass']:
        raise ValueError('No successful actual native reference')
    work = DATA / 'NATIVE_SOLVER'
    work.mkdir(parents=True, exist_ok=True)
    artifacts = []
    for row in lock['artifacts']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Native reference artifact drift ' + row['path'])
        dest = work / Path(row['path']).name
        shutil.copy2(row['path'], dest)
        artifacts.append({'path': str(dest), 'sha256': sha(dest), 'bytes': dest.stat().st_size})
    out['actual_reference_cost'] = out['cost']
    out['cost'] = cost(t0)
    out['artifacts'] = artifacts
    out['replay_native_status'] = 'HASH_VERIFIED_CACHED_ACTUAL_REFERENCE_SHARED_RESOURCE_GATE_BLOCKED'
    out['fresh_native_solve_in_this_replay'] = False
    out['reference_result_sha256'] = lock['native_result_sha256']
    out['cache_scope'] = 'The actual fresh native solve happened in the original lane; this copied replay reuses that one reference, not independent evidence.'
    dump(R / 'raw/NATIVE_SOLVER.json', out)
    print(json.dumps(clean(out), indent=2))
if __name__ == '__main__':
    cached_reference() if '--cached-reference' in sys.argv else run()
