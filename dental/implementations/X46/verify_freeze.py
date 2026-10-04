import hashlib
import json
from pathlib import Path
root = Path(__file__).resolve().parent
sha = lambda p: hashlib.sha256((root / p).read_bytes()).hexdigest()
one = json.loads((root / 'PREREG_RECIPES.json').read_text())
two = json.loads((root / 'R2_FREEZE.json').read_text())
for (file, expected) in [('review.py', one['code_sha256']), ('PREREG_R1.json', one['prereg_sha256']), ('SAMPLE.json', one['sample_sha256']), ('repairs.py', two['code_sha256']), ('PREREG_R2.json', two['prereg_sha256'])]:
    assert sha(file) == expected, f'Frozen input drift: {file}'
env = json.loads((root / 'EXECUTION_FREEZE.json').read_text())
assert sha('execute_review.py') == env['execution_wrapper_sha256']
print('Six frozen code/input/execution hashes verified')
