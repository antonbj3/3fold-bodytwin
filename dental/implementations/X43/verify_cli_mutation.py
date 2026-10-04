"""A real CLI fault probe on this lane's own curated map, restored byte-exactly in finally."""
import json, subprocess, datetime, xml.etree.ElementTree as ET
from pathlib import Path
p = Path(__file__).resolve().parent
f = p / 'raw/CURATED_EVIDENCE.json'
original = f.read_bytes()
try:
    d = json.loads(original)
    d['chains'][46]['proposed_status'] = 'PRESENT'
    f.write_text(json.dumps(d, ensure_ascii=False, indent=2) + '\n')
    bad = subprocess.run(['./run_all.sh'], cwd=p, text=True, capture_output=True)
    (p / 'raw/CLI_INVALID_VALIDATION.json').write_bytes((p / 'VALIDATION.json').read_bytes())
    (p / 'raw/CLI_WRONG_PRESENT.stdout.txt').write_text(bad.stdout)
    (p / 'raw/CLI_WRONG_PRESENT.stderr.txt').write_text(bad.stderr)
    if bad.returncode == 0:
        raise RuntimeError('Wrong full-chain PRESENT was accepted')
finally:
    f.write_bytes(original)
good = subprocess.run(['./run_all.sh'], cwd=p, text=True, capture_output=True)
(p / 'raw/FINAL_RUN.stdout.json').write_text(good.stdout)
(p / 'raw/FINAL_RUN.stderr.txt').write_text(good.stderr)
ET.parse(p / 'status_counts.svg')
receipt = {'checked_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'wrong_K47_PRESENT_exit': bad.returncode, 'wrong_input_error': bad.stderr.strip(), 'restored_input_exit': good.returncode, 'svg_xml': 'PASS', 'restoration_byte_exact': f.read_bytes() == original}
(p / 'raw/CLI_MUTATION_RECEIPT.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')
print(json.dumps(receipt, ensure_ascii=False, indent=2))
if good.returncode:
    raise SystemExit(good.returncode)
