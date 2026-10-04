from common import *
for p in list(ROOT.glob('PREREG*.json')) + list(ROOT.glob('FROZEN_PREDICTIONS*.json')):
    assert (ROOT / (p.name + '.sha256')).read_text().strip() == sha(p), p
lock = ROOT / 'REPRODUCTION_INPUTS_LOCK.json'
assert sha(lock) == (ROOT / (lock.name + '.sha256')).read_text().strip(), lock
for (path, record) in read(lock.name)['inputs'].items():
    assert sha(path) == record['sha256'], path
for (path, record) in read('SOURCE_MANIFEST.json').items():
    assert sha(path) == record['sha256'], path
print('Frozen files unchanged')
