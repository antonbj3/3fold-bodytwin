"""Measure absent coupling evidence without changing source acceptance."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

SOURCE = 'src/bodytwin/cells/endocrine/vitamin_d_activation.py'
PIN = '3f23dd823ea265b75e5a47151dd08c9b9266c5289481998f80018a4f7f92f295'


def finite(value):
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, dict):
        return all(finite(v) for v in value.values())
    if isinstance(value, list):
        return all(finite(v) for v in value)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[2] / SOURCE
    if hashlib.sha256(source.read_bytes()).hexdigest() != PIN:
        raise ValueError('Frozen source changed')
    rows, blobs, checks = [], [], []
    with tempfile.TemporaryDirectory(prefix='optional-coupling-') as directory:
        for leg in range(2):
            root = Path(directory) / str(leg)
            root.mkdir()
            env = dict(os.environ, BODYTWIN_OUT=str(root), CUDA_VISIBLE_DEVICES='',
                       OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
            try:
                code = subprocess.run([sys.executable, str(source)], cwd=root,
                                      env=env, capture_output=True, timeout=200).returncode
            except subprocess.TimeoutExpired:
                code = 124
            path = root / 'vitamin_d_activation/vitamin_d_activation_results.json'
            blob = path.read_bytes() if path.exists() else b''
            record = json.loads(blob) if blob else {}
            gates = record.get('gates', {})
            cascade = record.get('activation_cascade_kinetics_reused', {})
            comparisons = cascade.get('independent_refetch_cross_check', {})
            observed = [v.get('exact_match') for v in comparisons.values()]
            available = sum(v is not None for v in observed)
            emitted = gates.get('S4_independent_refetch_matches_sibling_cache')
            typed = bool(gates) and all(type(v) is bool for v in gates.values())
            consistent = typed and type(record.get('overall_pass_strict_all')) is bool and record['overall_pass_strict_all'] == all(gates.values())
            evidence = emitted is not True or (cascade.get('sibling_file_found') is True and available > 0)
            checks.append((code == 0 and bool(blob), consistent, finite(record), evidence))
            blobs.append(blob)
            rows.append(dict(bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(),
                             exit_code=code, gates=gates,
                             overall=record.get('overall_pass_strict_all'),
                             sibling_found=cascade.get('sibling_file_found'),
                             comparison_count=len(observed), available_comparisons=available,
                             null_self_tests=[k for k,v in record.get('self_tests', {}).items() if v is None],
                             emitted_coupling_gate=emitted))
    gates = dict(children_complete=all(c[0] for c in checks),
                 complete_repeat=bool(blobs[0]) and blobs[0] == blobs[1],
                 typed_overall=all(c[1] for c in checks),
                 finite=all(c[2] for c in checks),
                 coupling_pass_has_evidence=all(c[3] for c in checks))
    report = dict(source=SOURCE, source_sha256=PIN, rows=rows, gates=gates,
                  overall_pass=all(gates.values()))
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(gates, sort_keys=True))
    return 0 if report['overall_pass'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
