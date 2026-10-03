"""Frozen C5 release-ledger consumer. NumPy only; no native model or ODE import.

python3 -B query.py --ledger <ledger.npz> --rate .03 --rate .3 --output <new.json>
The ledger is valid only for its acquired upstream model, forcing and decay fate.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import time
import numpy as np

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ledger',type=Path,required=True)
    ap.add_argument('--rate',type=float,action='append',required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    start=time.process_time()
    source=Path(__file__).with_name('cell.py')
    fn=[n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='terminal_mac']
    assert len(fn)==1
    env={'np':np}
    exec(compile(ast.Module(body=fn,type_ignores=[]),str(source),'exec'),env)
    with np.load(args.ledger,allow_pickle=False) as d:
        t=d['time_min'];c=d['C5_nM']
    if t.ndim!=1 or c.shape!=t.shape or len(t)<2 or not np.all(np.isfinite(t)) or not np.all(np.isfinite(c)):
        raise ValueError('invalid chronological ledger')
    if np.any(c< -1e-8) or np.any(np.diff(c)>1e-8):raise ValueError('C5 positivity/release monotonicity failed')
    if any(k<=0 or not np.isfinite(k) for k in args.rate):raise ValueError('strictly positive finite rates required')
    rows=[dict(k_MAC_per_min=k,**env['terminal_mac'](t,c,k)) for k in args.rate]
    out=dict(review='PENDING_INDEPENDENT_REVIEW',rows=rows,full_ODE_calls=0,source_model_imports=0,
             nodes=len(t),ledger_sha256=hashlib.sha256(args.ledger.read_bytes()).hexdigest(),
             operator_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cpu_s=time.process_time()-start,
             scope='MAC-only queries; rigorous history acquisition/quadrature enclosure UNKNOWN')
    with args.output.open('x') as f:json.dump(out,f,indent=2,allow_nan=False)

if __name__=='__main__':main()
