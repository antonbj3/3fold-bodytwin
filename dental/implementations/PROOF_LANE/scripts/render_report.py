import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gencad_bench.reporting import report
report(json.loads((ROOT / 'results.json').read_text()))
