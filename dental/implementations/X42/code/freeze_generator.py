from common import *
from release_anchor import BENCHMARK_SHA256
import platform

def run():
    cfg = json.load(open(ROOT / 'FINAL_CONFIG.json'))
    files = {}
    for p in [ROOT / 'FINAL_CONFIG.json', ROOT / 'SPLITS.json', ROOT / 'PREREG_R1.json', ROOT / 'PREREG_R2.json', ROOT / 'PREREG_R3.json', *(ROOT / 'code').glob('*.py')]:
        files[str(p.relative_to(ROOT))] = {'sha256': sha(p), 'bytes': p.stat().st_size}
    models = {str(p.relative_to(DATA / 'FINAL_MODELS')): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in (DATA / 'FINAL_MODELS').glob('*/*.npz')}
    dump(ROOT / 'FROZEN_GENERATOR.json', dict(frozen_utc=now(), claim_type='algorithm', benchmark_sha256=BENCHMARK_SHA256, files=files, models=models, model_root=str(DATA / 'FINAL_MODELS'), selected_candidate=cfg['selected_candidate'], test_queries_at_freeze=0, model_selection_after_test='FORBIDDEN', fit_metadata_sha256=sha(ROOT / 'raw/FINAL_FIT.json'), scope='Frozen participant + predeclared referee aggregate extension; test reference bytes not read by development scripts'))
    (ROOT / 'FROZEN_GENERATOR.sha256').write_text(sha(ROOT / 'FROZEN_GENERATOR.json') + '  FROZEN_GENERATOR.json\n')
    state('GENERATOR_FROZEN', {'R3': 'no superiority vs equally informed control; full gate failed in development'}, 'Run isolated generation on public inputs only; freeze every output before one reference query', generator_sha256=sha(ROOT / 'FROZEN_GENERATOR.json'))
    print('generator frozen', sha(ROOT / 'FROZEN_GENERATOR.json'))
if __name__ == '__main__':
    run()
