from fc_common import *
import subprocess
out = DATA / 'R3_predictions'
out.mkdir(exist_ok=True)
cmd = read(ROOT / 'FROZEN_GENERATOR_R2.json')['argv']
cmd = [str(out) if x == str(DATA / 'R2_predictions') else '/runner/preparation_plugin.py' if x == '/runner/homolog_plugin.py' else x for x in cmd]
freeze(ROOT / 'FROZEN_GENERATOR_R3.json', dict(code_sha256={n: sha(ROOT / 'code' / n) for n in ['preparation_plugin.py', 'homolog_plugin.py']}, public_sha256=sha(ROOT / 'FROZEN_PUBLIC_INPUTS.json'), prereg_sha256=sha(ROOT / 'PREREG_R3.json'), argv=cmd, withheld_reference_mounted=False))
st = time.perf_counter()
with (ROOT / 'raw/generate_R3.log').open('w') as f:
    r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
dump(ROOT / 'raw/PROCESS_R3.json', dict(argv=cmd, exit_code=r.returncode, seconds=time.perf_counter() - st))
if r.returncode:
    raise RuntimeError('R3 failed, preserve log')
freeze(ROOT / 'FROZEN_PREDICTIONS_R3.json', dict(files={str(p.relative_to(out)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in out.rglob('*') if p.is_file()}, generator_sha256=sha(ROOT / 'FROZEN_GENERATOR_R3.json'), reference_queries_before_freeze=0))
state('R3_PREDICTIONS_FROZEN', 'Public-only preparation-conditioned full surfaces produced', 'Score unchanged v5 per-tooth vector')
