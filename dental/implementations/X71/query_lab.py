import argparse, json
from pathlib import Path
from planning import HERE, sha, read
from measurement_port import validate
p = argparse.ArgumentParser()
p.add_argument('--record', required=True)
p.add_argument('--output')
p.add_argument('--freeze', default=str(HERE / 'FROZEN_PREDICTIONS.json'))
a = p.parse_args()
record = json.loads(Path(a.record).read_text())
freeze = json.loads(Path(a.freeze).read_text())
out = validate(record, freeze, sha(Path(a.freeze)))
if a.output:
    Path(a.output).write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
