"""Measure complete frozen dependency-chain records before adapter design."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from optional_coupling_measurement import finite

STAGES = (
    ('cardiovascular/cardiac_output_geometric.py', '8c55ca5dabe0c9a85c36212071687c266188a5ee1f8f2e733ab9586eb7db372b', 'cardiac_output_geometric'),
    ('renal/renal_filtration.py', 'd0d451282fec3b58b696ace32e9c8e79540f5e417a2472731ef23cb1128629ab', 'renal_filtration'),
    ('endocrine/calcium_pth_vitd.py', 'bb91034d05bc10001593dd5e8c7f755601f3ec9d986368b50147c3150fcd9e90', 'calcium_pth_vitd'),
    ('endocrine/vitamin_d_activation.py', '3f23dd823ea265b75e5a47151dd08c9b9266c5289481998f80018a4f7f92f295', 'vitamin_d_activation'),
)
KEYS = {'vitamin_d3_days', '25_OH_D_days', '1_25_OH2_D_hours'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    repo = Path(__file__).resolve().parents[2]
    for name, pin, _ in STAGES:
        if hashlib.sha256((repo/'src/bodytwin/cells'/name).read_bytes()).hexdigest() != pin:
            raise ValueError('Frozen source changed')
    rows, chains = [], []
    completed = numeric = clean = immutable = evidenced = True
    with tempfile.TemporaryDirectory(prefix='activation-chain-') as directory:
        for leg in range(2):
            root = Path(directory)/str(leg)
            root.mkdir()
            env = dict(os.environ, BODYTWIN_OUT=str(root), CUDA_VISIBLE_DEVICES='',
                       OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
            artifacts = {}
            for name, pin, stage in STAGES:
                try:
                    code = subprocess.run([sys.executable, str(repo/'src/bodytwin/cells'/name)],
                                          cwd=root, env=env, capture_output=True, timeout=200).returncode
                except subprocess.TimeoutExpired:
                    code = 124
                immutable &= all(path.exists() and path.read_bytes() == blob for path,blob in artifacts.items())
                path = root/stage/(stage+'_results.json')
                blob = path.read_bytes() if path.exists() else b''
                artifacts[path] = blob
                record = json.loads(blob) if blob else {}
                gates = record.get('gates', {})
                completed &= code == 0 and bool(blob)
                numeric &= bool(gates) and all(type(v) is bool and v for v in gates.values())
                clean &= finite(record)
                rows.append(dict(leg=leg, stage=stage, source_sha256=pin, exit_code=code,
                                 bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(), gates=gates,
                                 overall=record.get('overall_pass_strict_all', record.get('overall_pass', gates.get('overall_pass')))))
                if stage == 'vitamin_d_activation':
                    comparisons = record.get('activation_cascade_kinetics_reused', {}).get('independent_refetch_cross_check', {})
                    evidenced &= set(comparisons) == KEYS and all(v.get('exact_match') is True for v in comparisons.values())
            chains.append(list(artifacts.values()))
    gates = dict(children_complete=completed, complete_repeat=chains[0] == chains[1] and all(chains[0]),
                 source_gates_pass=numeric, finite=clean, upstream_immutable=immutable,
                 complete_coupling=evidenced)
    result = dict(rows=rows, gates=gates, overall_pass=all(gates.values()))
    a.output.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps(gates, sort_keys=True))
    return 0 if result['overall_pass'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
