"""Probe frozen coupling acceptance using explicit synthetic sibling schemas."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from optional_coupling_measurement import SOURCE, PIN, finite

KEYS = ('vitamin_d3_days', '25_OH_D_days', '1_25_OH2_D_hours')
MATCH = dict(zip(KEYS, (60.0, 15.0, 15.0)))
ARMS = {
    'absent': None,
    'empty': {},
    'partial': {'vitamin_d_kinetics': {'half_lives': {KEYS[0]: 60.0}}},
    'complete_match': {'vitamin_d_kinetics': {'half_lives': MATCH}},
    'complete_mismatch': {'vitamin_d_kinetics': {'half_lives': dict(MATCH, vitamin_d3_days=61.0)}},
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[2] / SOURCE
    if hashlib.sha256(source.read_bytes()).hexdigest() != PIN:
        raise ValueError('Frozen source changed')
    rows = []
    complete = repeat = typed = clean = supported = True
    with tempfile.TemporaryDirectory(prefix='coupling-schema-') as directory:
        for arm, sibling in ARMS.items():
            blobs, legs = [], []
            for leg in range(2):
                root = Path(directory) / arm / str(leg)
                root.mkdir(parents=True)
                if sibling is not None:
                    path = root / 'calcium_pth_vitd/calcium_pth_vitd_results.json'
                    path.parent.mkdir()
                    path.write_text(json.dumps(sibling))
                env = dict(os.environ, BODYTWIN_OUT=str(root), CUDA_VISIBLE_DEVICES='',
                           OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
                try:
                    code = subprocess.run([sys.executable, str(source)], cwd=root, env=env,
                                          capture_output=True, timeout=200).returncode
                except subprocess.TimeoutExpired:
                    code = 124
                path = root / 'vitamin_d_activation/vitamin_d_activation_results.json'
                blob = path.read_bytes() if path.exists() else b''
                record = json.loads(blob) if blob else {}
                gates = record.get('gates', {})
                cascade = record.get('activation_cascade_kinetics_reused', {})
                comparisons = cascade.get('independent_refetch_cross_check', {})
                evidence = set(comparisons) == set(KEYS) and all(type(comparisons[k].get('exact_match')) is bool for k in KEYS)
                coupling = gates.get('S4_independent_refetch_matches_sibling_cache')
                complete &= code == 0 and bool(blob)
                clean &= finite(record)
                typed &= bool(gates) and all(type(v) is bool for v in gates.values()) and type(record.get('overall_pass_strict_all')) is bool and record['overall_pass_strict_all'] == all(gates.values())
                supported &= coupling is not True or evidence
                blobs.append(blob)
                legs.append(dict(exit_code=code, bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(),
                                 gate_count=len(gates), passed=sum(v is True for v in gates.values()),
                                 coupling=coupling, evidence_complete=evidence,
                                 comparisons={k:v.get('exact_match') for k,v in comparisons.items()},
                                 overall=record.get('overall_pass_strict_all')))
            repeat &= bool(blobs[0]) and blobs[0] == blobs[1]
            rows.append(dict(arm=arm, legs=legs))
    indexed = {row['arm']: row['legs'] for row in rows}
    controls = all(v['coupling'] is True for v in indexed['complete_match']) and all(v['coupling'] is False for v in indexed['complete_mismatch'])
    gates = dict(children_complete=complete, complete_repeat=repeat, typed_overall=typed,
                 finite=clean, coupling_pass_has_complete_evidence=supported, complete_controls=controls)
    report = dict(source=SOURCE, source_sha256=PIN, rows=rows, gates=gates, overall_pass=all(gates.values()))
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(gates, sort_keys=True))
    return 0 if report['overall_pass'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
