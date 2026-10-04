"""Executed adversarial checks; positive fixtures accompany refusal probes."""
import tempfile, subprocess, sys, time, os
from pathlib import Path
import numpy as np
from .common import *
from .sandbox import bundle, command
from .measurements import scalar_error

def run():
    start = time.perf_counter()
    log = ROOT / 'raw/controls_final.log'
    with log.open('w') as f:
        p = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT / 'tests'), '-v'], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
    if p.returncode:
        raise RuntimeError('control suite failed, ' + str(log))
    bundle()
    secret = DATA / 'private/sandbox_secret.txt'
    secret.parent.mkdir(parents=True, exist_ok=True)
    secret.write_text('private reference must not be read\n')
    with tempfile.TemporaryDirectory(dir=DATA) as td:
        inp = Path(td) / 'input'
        out = Path(td) / 'output'
        inp.mkdir()
        out.mkdir()
        (inp / 'public.txt').write_text('public')
        code = "import pathlib,socket,json; p=pathlib.Path; assert p('/inputs/public.txt').read_text()=='public'; assert not p(" + repr(str(secret)) + ").exists(); assert not p('/code/gencad_bench_v2/scorer.py').exists(); assert not p('/mnt').exists(); print('PUBLIC_READ_OK_PRIVATE_AND_SCORER_DENIED')"
        cmd = command(inp, out, ['-c', code])
        probe = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    if probe.returncode:
        raise RuntimeError('sandbox control failed ' + probe.stderr)
    tests = []
    for (name, correct, bad, truth, tol) in [('fracture_force_N', [1000.0], [1000000.0], [1000.0], 250.0), ('cement_um', [80.0], [80000.0], [80.0], 20.0), ('sinter_angle_deg', [0.5], [50.0], [0.5], 0.2)]:
        good = scalar_error(correct, truth)['rmse'] <= tol
        rejected = scalar_error(bad, truth)['rmse'] > tol
        if not good or not rejected:
            raise ValueError('measurement control does not discriminate')
        tests.append(dict(check=name, positive=good, injected_error_rejected=rejected))
    result = dict(unit_suite_pass=True, unit_test_count=9, unit_test_log=str(log), sandbox=dict(returncode=probe.returncode, stdout=probe.stdout, argv=cmd, private_read_rejected=True, scope='unprivileged plugin mount/network namespace; kernel escapes and privileged host outside threat model'), measurement_faults=tests, seconds=time.perf_counter() - start)
    dump(ROOT / 'raw/CONTROL_REPORT.json', result)
    return result
