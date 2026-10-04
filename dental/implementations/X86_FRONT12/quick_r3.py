"""Execute existing R3 adversaries and table audit, without all72 rescoring."""
import copy
import importlib.util
import json
import sys
import tempfile
import time
from pathlib import Path

def main():
    root = Path(sys.argv[1])
    batch = Path(sys.argv[2])
    started = time.monotonic()
    sys.path.insert(0, str(root / 'code'))
    from validation import run as validate
    from scorer_controls import run as scorer_controls
    witness = validate()
    scorer = scorer_controls()
    spec = importlib.util.spec_from_file_location('x86b_original_audit', batch / 'audit_checks.py')
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    current = audit.audit('PROOF_LANE_FULL_CROWN_R3', root)
    failures = []
    original = json.loads((root / 'results.json').read_text())
    for (field, value) in [('median_p95_mm', 999.0), ('joint_pass', 1)]:
        corrupted = copy.deepcopy(original)
        table = next((t for t in corrupted['table'] if not t['round'].startswith('R2')))
        table[field] = value
        with tempfile.TemporaryDirectory(dir=root / 'raw') as tmp:
            Path(tmp, 'results.json').write_text(json.dumps(corrupted))
            try:
                audit.audit('PROOF_LANE_FULL_CROWN_R3', tmp)
            except AssertionError:
                failures.append({'field': field, 'value': value, 'rejected': True})
            else:
                raise AssertionError('Injected result passed original R3 audit: ' + field)
    receipt = {'demo_id': 'PROOF_LANE_FULL_CROWN_R3', 'status': 'PASS', 'pass': True, 'execution_scope': 'QUICK_EXISTING_CONTROLS: frozen72 records hashed; original spatial/scorer mutations executed; independent table audit rerun. No full72 rescoring or generator replay.', 'wall_seconds': time.monotonic() - started, 'table_audit': current, 'sufficiency': witness, 'scorer_controls': scorer, 'injected_table_faults': failures, 'physical_measurements': 0}
    (root / 'QUICK_X86B_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('R3 QUICK PASS: original controls, exact witness, table audit,2 injected values rejected', flush=True)
if __name__ == '__main__':
    main()
