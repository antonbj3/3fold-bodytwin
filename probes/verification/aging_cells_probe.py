"""Observe required gates and complete repeat bytes of two unchanged cells."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def assess(data):
    gates=data.get('gates')
    if not isinstance(gates,dict) or not gates or any(type(v) is not bool for v in gates.values()):
        raise ValueError('Nonempty boolean gate dictionary required')
    required=[v for k,v in gates.items() if k.startswith('required_')]
    if not required or not all(required) or data.get('required_gates_overall_pass') is not True:
        raise ValueError('Required numerical gates did not pass')
    return dict(required_count=len(required),bonus_count=len(gates)-len(required),
                bonus_failed=sum(not v for k,v in gates.items() if not k.startswith('required_')))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    root=Path(__file__).resolve().parents[2]
    names=('cellular_senescence','telomere_attrition')
    rows=[];records=[]
    with tempfile.TemporaryDirectory(prefix='bodytwin-cell-observer-') as directory:
        for leg in (0,1):
            output=Path(directory)/str(leg);output.mkdir()
            env=os.environ.copy();env.update(BODYTWIN_OUT=str(output),CUDA_VISIBLE_DEVICES='',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            leg_rows=[];leg_records=[]
            for name in names:
                script=root/f'src/bodytwin/cells/aging/{name}.py'
                proc=subprocess.run([sys.executable,str(script)],cwd=output,env=env,capture_output=True,text=True,timeout=60)
                if proc.returncode!=0:raise ValueError('Cell process failed')
                record=(output/name/f'{name}_results.json').read_bytes()
                data=json.loads(record);assessment=assess(data)
                leg_rows.append(dict(module=script.relative_to(root).as_posix(),source_hash=hashlib.sha256(script.read_bytes()).hexdigest(),
                                     output_hash=hashlib.sha256(record).hexdigest(),output_bytes=len(record),**assessment))
                leg_records.append(record)
            records.append(leg_records);rows.append(leg_rows)
    bad=[{}, {'gates':{},'required_gates_overall_pass':True},
         {'gates':{'required_a':False},'required_gates_overall_pass':True},
         {'gates':{'required_a':'true'},'required_gates_overall_pass':True},
         {'gates':{'required_a':True},'required_gates_overall_pass':False},
         {'gates':{'bonus_a':True},'required_gates_overall_pass':True}]
    refused=0
    for data in bad:
        try:assess(data)
        except ValueError:refused+=1
    # A bonus failure is recorded, not silently promoted into a required gate.
    bonus=assess({'gates':{'required_a':True,'bonus_a':False},'required_gates_overall_pass':True})
    gates=dict(full_repeat=rows[0]==rows[1] and records[0]==records[1],required_gates=True,
               refusal_controls=refused==6,bonus_preserved=bonus['bonus_failed']==1)
    result=dict(rows=rows,gates=gates,scope='Unchanged source-defined required gate and full-record replay; bonus/abstention semantics retained, no fresh literature or clinical review.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
