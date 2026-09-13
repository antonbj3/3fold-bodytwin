"""Strict complete-chain verification beside the frozen measurement observer."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

REPORT = 'reports/wording_refresh/coupled_measurement.json'
PIN = 'd1abf2fd61870b818a15e59eb8e0fdbc9fb3e5aada78a75fcedcec9f5783cc5c'


def assess(record, expected):
    if type(record) is not dict or record.keys() != expected.keys():
        raise ValueError('Exact record schema required')
    gates = record['gates']
    if type(gates) is not dict or gates.keys() != expected['gates'].keys() or any(v is not True for v in gates.values()) or record['overall_pass'] is not True:
        raise ValueError('Literal successful measurement gates required')
    rows = record['rows']
    if type(rows) is not list or len(rows) != 8:
        raise ValueError('Two complete four-stage chains required')
    for row, reference in zip(rows, expected['rows']):
        if type(row) is not dict or row.keys() != reference.keys():
            raise ValueError('Exact stage schema required')
        if type(row['leg']) is not int or row['leg'] != reference['leg'] or row['stage'] != reference['stage']:
            raise ValueError('Exact ordered stage inventory required')
        values = row['gates']
        if type(values) is not dict or values.keys() != reference['gates'].keys() or any(v is not True for v in values.values()) or row['overall'] is not True:
            raise ValueError('Exact literal source gate inventory required')
        if type(row['bytes']) is not int or row['bytes'] <= 0 or type(row['exit_code']) is not int or row['exit_code'] != 0:
            raise ValueError('Completed full artifact required')
    if record != expected:
        raise ValueError('Complete frozen artifacts must match')
    return True


def controls(expected):
    mutations = [lambda r:r['rows'].pop(),
                 lambda r:r['rows'].__setitem__(1, copy.deepcopy(r['rows'][0])),
                 lambda r:r['rows'][0]['gates'].update(extra=True),
                 lambda r:r['rows'][0]['gates'].__setitem__(next(iter(r['rows'][0]['gates'])), 1),
                 lambda r:r['rows'][0].__setitem__('overall', False),
                 lambda r:r['rows'][0]['gates'].__setitem__(next(iter(r['rows'][0]['gates'])), False),
                 lambda r:r['rows'][0].__setitem__('bytes', 0),
                 lambda r:r['rows'][0].__setitem__('sha256', '0'*64)]
    refused = 0
    for mutate in mutations:
        record = copy.deepcopy(expected)
        mutate(record)
        try:
            assess(record, expected)
        except ValueError:
            refused += 1
    return refused == len(mutations)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    root = Path(__file__).resolve().parents[2]
    raw = (root/REPORT).read_bytes()
    if hashlib.sha256(raw).hexdigest() != PIN:
        raise ValueError('Frozen measurement changed')
    expected = json.loads(raw)
    with tempfile.TemporaryDirectory(prefix='coupled-contract-') as directory:
        output = Path(directory)/'measurement.json'
        code = subprocess.run([sys.executable, str(root/'probes/verification/coupled_activation_measurement.py'),
                               '--output', str(output)], capture_output=True, timeout=820).returncode
        fresh = output.read_bytes() if output.exists() else b''
    accepted = False
    try:
        accepted = assess(json.loads(fresh), expected)
    except (ValueError, TypeError, KeyError):
        pass
    gates = dict(chain_complete=code == 0, strict_contract=accepted,
                 complete_report_identity=fresh == raw, malformed_refusals=controls(expected))
    result = dict(measurement_sha256=hashlib.sha256(fresh).hexdigest(), gates=gates, overall_pass=all(gates.values()))
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(gates, sort_keys=True))
    return 0 if result['overall_pass'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
