"""Capture values behind seven frozen required-gate failures."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

NAME='endocrine/glucose_meal_dallaman2007.py'
PIN='1c96a5970119e2ae2bc59b26f547b950ca425c35f102d476717de8ced97b1500'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells-root',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    root=a.cells_root or Path(__file__).resolve().parents[2]/'src/bodytwin/cells'
    script=root/NAME
    if hashlib.sha256(script.read_bytes()).hexdigest()!=PIN:raise ValueError('Frozen source identity mismatch')
    legs=[];records=[]
    with tempfile.TemporaryDirectory(prefix='glucose-gates-') as directory:
        for leg in (0,1):
            out=Path(directory)/str(leg);out.mkdir()
            env=os.environ.copy();env.update(BODYTWIN_OUT=str(out),CUDA_VISIBLE_DEVICES='',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            child=subprocess.run([sys.executable,str(script.resolve())],cwd=out,env=env,capture_output=True,text=True,timeout=200)
            if child.returncode!=0:raise ValueError('Numerical child failed')
            raw=(out/script.stem/f'{script.stem}_results.json').read_bytes();records.append(raw);data=json.loads(raw)
            target=data['recorded_original_claim'];rows=[]
            for point in (40,50,60,75,100):
                value=data['scenario4_egp_suppression_pct'][str(point)];expected=target[f'egp{point}']
                rows.append(dict(case=f'suppression_{point}',observed=value,target=expected,relative_error=abs(value-expected)/abs(expected)))
            for name,arm,key,target_key in [('endpoint','a_bolus_alone','G_at_600','a_end600'),('crossing','b_bolus_plus_basal','cross_below_200_t','b_cross200')]:
                value=data['scenario6_t1d_pump'][arm][key];expected=target[target_key]
                rows.append(dict(case=name,observed=value,target=expected,relative_error=None if value is None else abs(value-expected)/abs(expected)))
            legs.append(dict(rows=rows,full_record_hash=hashlib.sha256(raw).hexdigest(),
                             required_pass=data['required_gates_overall_pass'],failed_count=sum(v is False for v in data['gates'].values())))
    gates=dict(full_repeat=legs[0]==legs[1] and records[0]==records[1],
               captured_failures=all(r['failed_count']==7 and all(v['relative_error'] is None or v['relative_error']>.05 for v in r['rows']) for r in legs),
               original_required=all(r['required_pass'] is True for r in legs))
    result=dict(legs=legs,gates=gates,source_hash=PIN,scope='Frozen software target mismatch, not a clinical interpretation or changed acceptance band.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
