"""CLI for the corrected conditional DC observation port; no file mutation."""
import json, sys
from pathlib import Path
from lab_port import intake
from lab_port_guard import checked_intake
P = Path(__file__).resolve().parents[1]
if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python3 code/lab_port_reviewed.py observation.json')
    print(json.dumps(checked_intake(json.loads(Path(sys.argv[1]).read_text()), intake, json.loads((P / 'FROZEN_PREDICTIONS.json').read_text())), indent=2))
