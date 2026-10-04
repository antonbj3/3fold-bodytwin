from common import *
from release_anchor import BENCHMARK_SHA256
from scorer import decode
from integrity import predictions

def verify_generator():
    f = json.load(open(ROOT / 'FROZEN_GENERATOR.json'))
    if sha(ROOT / 'FROZEN_GENERATOR.json') != (ROOT / 'FROZEN_GENERATOR.sha256').read_text().split()[0]:
        raise ValueError('generator anchor mismatch')
    for (name, r) in f['files'].items():
        if sha(ROOT / name) != r['sha256']:
            raise ValueError('frozen generator changed: ' + name)
    for (name, r) in f['models'].items():
        if sha(DATA / 'FINAL_MODELS' / name) != r['sha256']:
            raise ValueError('frozen learned model changed: ' + name)
    if (ROOT / 'FROZEN_CONTROL.json').exists():
        c = json.load(open(ROOT / 'FROZEN_CONTROL.json'))
        if sha(ROOT / 'FROZEN_CONTROL.json') != (ROOT / 'FROZEN_CONTROL.sha256').read_text().split()[0]:
            raise ValueError('strong control anchor mismatch')
        for (name, r) in c['files'].items():
            if sha(ROOT / name) != r['sha256']:
                raise ValueError('strong control code/config changed: ' + name)
        for (name, r) in c['models'].items():
            if sha(DATA / 'FULL_KERNEL_CONTROL' / name) != r['sha256']:
                raise ValueError('full kernel model changed: ' + name)
    return f

def run():
    verify_generator()
    base = DATA / 'PREDICTIONS'
    cfg = json.load(open(ROOT / 'EVAL_CONFIG.json'))
    keys = [r['key'] for r in split()['metadata']]
    files = {}
    for name in cfg['participants']:
        actual = {p.stem for p in (base / name).glob('*.npz')}
        if actual != set(keys):
            raise ValueError('missing or extra cases ' + name)
        for key in keys:
            p = base / name / (key + '.npz')
            with np.load(p, allow_pickle=False) as a:
                decode(dict(a), tasks(key))
            files[name + '/' + p.name] = {'sha256': sha(p), 'bytes': p.stat().st_size}
    manifest = {'frozen_utc': now(), 'benchmark_sha256': BENCHMARK_SHA256, 'generator_sha256': sha(ROOT / 'FROZEN_GENERATOR.json'), 'control_extension_sha256': sha(ROOT / 'FROZEN_CONTROL.json'), 'participants': cfg['participants'], 'files': files}
    text = json.dumps(manifest, separators=(',', ':'), allow_nan=False) + '\n'
    p = base / 'FROZEN_PREDICTIONS.json'
    if p.exists():
        raise ValueError('refuse refreeze')
    if len(text.encode()) > 2000000:
        raise ValueError('prediction manifest exceeds official intake bound')
    p.write_text(text)
    receipt = sha(p)
    dump(ROOT / 'FROZEN_PREDICTIONS.json', dict(frozen_utc=manifest['frozen_utc'], manifest_path=str(p), manifest_sha256=receipt, generator_sha256=manifest['generator_sha256'], benchmark_sha256=BENCHMARK_SHA256, files=len(files), cases=len(keys), participants=cfg['participants'], reference_query_count_at_freeze=0, prediction_bytes=sum((r['bytes'] for r in files.values())), resolution='PER_POINT', physical_measurement='NOT_RUN'))
    predictions(base, receipt, keys, BENCHMARK_SHA256)
    state('PREDICTIONS_FROZEN', {'cases': len(keys), 'fields': len(files), 'test_queries': 0}, 'Run single lane-owned trusted referee query; models cannot change', prediction_sha256=receipt)
    print('predictions frozen', receipt, 'bytes', sum((r['bytes'] for r in files.values())))
if __name__ == '__main__':
    run()
