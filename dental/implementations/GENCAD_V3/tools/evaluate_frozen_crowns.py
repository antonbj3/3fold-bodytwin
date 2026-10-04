import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *
from integrity import safe_bytes, load_npz, h
from crown_track import evaluate
files = {}
for fn in ['FROZEN_PREPARATIONS_R3.json', 'FROZEN_CROWNS.json']:
    fr = read(ROOT / fn)
    assert digest(fr['payload']) == fr['payload_sha256']
    files.update({'payload/' + k: v['sha256'] for (k, v) in fr['payload']['files'].items()})
    files[fn] = sha(ROOT / fn)
for name in ['PREREG_R4.json', 'code/crown_track.py', 'code/preparation.py', 'code/integrity.py', 'code/util.py']:
    files[name] = sha(ROOT / name)
freeze(ROOT / 'rounds/R4_VALIDATION_INPUTS.json', files)

class Frozen:
    payload = PAYLOAD

    def bytes(self, name):
        (root, rel) = (PAYLOAD, name[8:]) if name.startswith('payload/') else (ROOT, name)
        b = safe_bytes(root, rel)
        if h(b) != files[name]:
            raise ValueError('R4 frozen bytes differ')
        return b

    def json(self, name):
        return json.loads(self.bytes(name))

    def npz(self, name):
        return load_npz(self.bytes(name))
out = evaluate(Frozen())
dump(ROOT / 'rounds/R4.json', out)
print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
