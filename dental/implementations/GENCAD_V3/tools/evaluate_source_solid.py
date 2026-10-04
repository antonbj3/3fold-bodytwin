import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *
from integrity import safe_bytes, load_npz, h
from source_solid_replay import run
r = int(sys.argv[1])
fn = 'FROZEN_CROWNS_R' + str(r) + '.json'
fr = read(ROOT / fn)
assert digest(fr['payload']) == fr['payload_sha256']
files = {'payload/' + k: v['sha256'] for (k, v) in fr['payload']['files'].items()}
files[fn] = sha(ROOT / fn)
if r == 7:
    prior = read(ROOT / 'FROZEN_CROWNS_R6.json')['payload']
    files.update({'payload/' + k: v['sha256'] for (k, v) in prior['files'].items()})
for name in ['PREREG_R' + str(r) + '.json', 'code/source_solid_replay.py', 'code/crown_track.py', 'code/preparation.py', 'code/integrity.py', 'code/util.py']:
    files[name] = sha(ROOT / name)
freeze(ROOT / ('rounds/R' + str(r) + '_VALIDATION_INPUTS.json'), files)

class Frozen:
    payload = PAYLOAD

    def bytes(self, name):
        (root, rel) = (PAYLOAD, name[8:]) if name.startswith('payload/') else (ROOT, name)
        b = safe_bytes(root, rel)
        if h(b) != files[name]:
            raise ValueError('frozen bytes differ')
        return b

    def json(self, name):
        return json.loads(self.bytes(name))

    def npz(self, name):
        return load_npz(self.bytes(name))
out = run(Frozen(), r)
dump(ROOT / ('rounds/R' + str(r) + '.json'), out)
print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
