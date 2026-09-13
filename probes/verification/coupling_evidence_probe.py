"""Replay synthetic coupling records against an explicit abstention contract."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from bodytwin.framework.coupling_evidence_v1 import CouplingEvidence, self_test


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    root = Path(__file__).resolve().parents[2]
    path = root / 'reports/coupling_schema_control.json'
    raw = path.read_bytes()
    report = json.loads(raw)
    contract = CouplingEvidence(('vitamin_d3_days', '25_OH_D_days', '1_25_OH2_D_hours'))
    def evaluate():
        return [dict(arm=row['arm'], leg=i, source_pass=leg['coupling'],
                     verdict=asdict(contract.assess(leg['comparisons'])))
                for row in report['rows'] for i,leg in enumerate(row['legs'])]
    first, second = evaluate(), evaluate()
    expected = {'absent': 'ABSTAIN', 'empty': 'ABSTAIN', 'partial': 'ABSTAIN',
                'complete_match': 'PASS', 'complete_mismatch': 'FAIL'}
    gates = self_test()
    gates['archived_schema_parity'] = len(first) == 10 and all(row['verdict']['state'] == expected[row['arm']] for row in first)
    gates['complete_verdict_repeat'] = first == second
    result = dict(input_sha256=hashlib.sha256(raw).hexdigest(), verdicts=first,
                  source_pass_abstentions=sum(row['source_pass'] and row['verdict']['state'] == 'ABSTAIN' for row in first),
                  gates=gates, overall_pass=all(gates.values()))
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(gates, sort_keys=True))
    return 0 if result['overall_pass'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
