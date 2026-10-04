from common import *
files = {str(p.relative_to(ROOT)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in [ROOT / 'PREREG_CONTROL_EXTENSION.json', ROOT / 'EVAL_CONFIG.json', ROOT / 'code/control_extension.py', ROOT / 'code/full_kernel_participant.py', ROOT / 'code/referee_v2.py']}
models = {p.name: {'sha256': sha(p), 'bytes': p.stat().st_size} for p in (DATA / 'FULL_KERNEL_CONTROL').glob('*.npz')}
dump(ROOT / 'FROZEN_CONTROL.json', dict(frozen_utc=now(), original_generator_sha256=sha(ROOT / 'FROZEN_GENERATOR.json'), candidate_modification=False, files=files, models=models, test_queries_at_freeze=0))
(ROOT / 'FROZEN_CONTROL.sha256').write_text(sha(ROOT / 'FROZEN_CONTROL.json') + '  FROZEN_CONTROL.json\n')
print('additional control frozen', sha(ROOT / 'FROZEN_CONTROL.json'))
