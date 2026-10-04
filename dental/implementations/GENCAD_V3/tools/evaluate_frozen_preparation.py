"""Early R2 test, after immutable predictions; same evaluator as final demo."""
import sys, hashlib, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *
from integrity import safe_bytes, load_npz, h
from preparation import evaluate
fr = read(ROOT / 'FROZEN_PREPARATIONS.json')
assert digest(fr['payload']) == fr['payload_sha256']
assert sha(ROOT / 'PREREG_R2.json') == fr['payload']['prereg_sha256']
files = {'payload/' + k: v['sha256'] for (k, v) in fr['payload']['files'].items()}
for name in ['FROZEN_PREPARATIONS.json', 'PREREG_R2.json', 'code/preparation.py', 'code/integrity.py', 'code/util.py']:
    files[name] = sha(ROOT / name)
freeze(ROOT / 'rounds/R2_VALIDATION_INPUTS.json', files)

class Frozen:
    payload = PAYLOAD

    def bytes(self, name):
        if name.startswith('payload/'):
            b = safe_bytes(PAYLOAD, name[len('payload/'):])
        else:
            b = safe_bytes(ROOT, name)
        if h(b) != files[name]:
            raise ValueError('R2 frozen bytes differ')
        return b

    def json(self, name):
        return json.loads(self.bytes(name))

    def npz(self, name):
        return load_npz(self.bytes(name))
out = evaluate(Frozen())
dump(ROOT / 'rounds/R2.json', out)
print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
