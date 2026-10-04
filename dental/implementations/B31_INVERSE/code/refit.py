"""Fresh optimizer replay in a small isolated directory; frozen originals untouched."""
from common import *
import subprocess, time

def run():
    start = time.perf_counter()
    for r in read(ROOT / 'INPUT_LOCK.json'):
        if sha(r['path']) != r['sha256']:
            raise RuntimeError('Source drift ' + r['path'])
    total = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    if total > 2500000000:
        raise RuntimeError('Insufficient lane intermediate budget for replay')
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    work = DATA / 'refit' / stamp
    work.mkdir(parents=True)
    (work / 'data').mkdir()
    for name in ['PREREG_B.json', 'PREREG_B.sha256']:
        (work / name).write_bytes((ROOT / name).read_bytes())
    env = os.environ.copy()
    env.update(PROOF_LANE_INVERSE_RUN_ROOT=str(work), PROOF_LANE_INVERSE_RUN_DATA=str(work / 'data'))
    with (work / 'RUN.log').open('w') as log:
        subprocess.run([sys.executable, str(SOURCE_ROOT / 'code/inverse.py')], env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    observed = read(work / 'raw/B_FIT.json')
    reference = read(ROOT / 'raw/B_FIT.json')
    err = float(np.max(abs(np.array(observed['coefficients_mm']) - reference['coefficients_mm'])))
    a = np.load(DATA / 'B_FIELDS.npz')
    b = np.load(work / 'data/B_FIELDS.npz')
    field_err = max((float(np.max(abs(a[k] - b[k])) / max(np.max(abs(a[k])), 1.0)) for k in a.files))
    out = dict(work=work, coefficients_error_mm=err, fields_relative_error=field_err, pass_gate=err <= 1e-08 and field_err <= 1e-08, seconds=time.perf_counter() - start, physical_validation=False)
    save(ROOT / 'raw/REFIT_REPLAY.json', out)
    print(json.dumps(clean(out)))
    if not out['pass_gate']:
        raise RuntimeError('Optimizer replay mismatch')
if __name__ == '__main__':
    run()
