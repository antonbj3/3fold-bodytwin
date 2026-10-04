"""Expose already locked numerical evidence omitted from the reporting JSON."""
from pathlib import Path
import json
P = Path(__file__).resolve().parent
q = json.loads((P / 'results.json').read_text())
q['review_numeric_evidence'] = json.loads((P / 'REVIEW_NUMERIC_EVIDENCE.json').read_text())
(P / 'results.json').write_text(json.dumps(q, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
