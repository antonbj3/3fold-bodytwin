"""Pinned replay of an explicit embedded-overall gate contract."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from bodytwin.framework.numerical_gate_contract_v1 import GateContract

SOURCE='src/bodytwin/cells/gastrointestinal/deglutition_swallowing.py'
PIN='43078fbaed1dd0842b1413ba3cea6b7d798bbf778f4abe38f7808321696378d0'
OVERALL='overall_pass_strict_all'


def inventory(raw):
    if hashlib.sha256(raw).hexdigest()!=PIN:raise ValueError('Frozen source mismatch')
    functions=[n for n in ast.parse(raw).body if isinstance(n,ast.FunctionDef) and n.name=='compute_gates']
    if len(functions)!=1:raise ValueError('Explicit gate function missing')
    names=[]
    for node in functions[0].body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Subscript) and isinstance(target.value,ast.Name) and target.value.id=='gates':
                    names.append(ast.literal_eval(target.slice))
    if len(names)!=19 or names[-1]!=OVERALL or len(set(names))!=19:raise ValueError('Gate inventory changed')
    return tuple(sorted(names[:-1]))


def assess(record,names):
    gates=record.get('gates')
    if not isinstance(gates,dict) or set(gates)!=set(names)|{OVERALL}:raise ValueError('Exact embedded gate inventory required')
    return GateContract(names,names,'overall').assess({'gates':{k:gates[k] for k in names},'overall':gates[OVERALL]})


def controls():
    good={'gates':{'a':True,OVERALL:True}}
    bad=[{}, {'gates':{'a':True}}, {'gates':{'a':1,OVERALL:True}},
         {'gates':{'a':True,OVERALL:1}}, {'gates':{'a':False,OVERALL:True}},
         {'gates':{'a':True,OVERALL:False}}, {'gates':{'a':True,'extra':True,OVERALL:True}}]
    refused=0
    for record in bad:
        try:assess(record,('a',))
        except ValueError:refused+=1
    return refused==7 and assess(good,('a',)).accepted and not assess({'gates':{'a':False,OVERALL:False}},('a',)).accepted


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    root=Path(__file__).resolve().parents[2];source=root/SOURCE;names=inventory(source.read_bytes());rows=[];records=[]
    with tempfile.TemporaryDirectory(prefix='embedded-gates-') as directory:
        for leg in (0,1):
            out=Path(directory)/str(leg);out.mkdir()
            env=os.environ.copy();env.update(BODYTWIN_OUT=str(out),CUDA_VISIBLE_DEVICES='',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            try:code=subprocess.run([sys.executable,str(source)],cwd=out,env=env,capture_output=True,timeout=200).returncode
            except subprocess.TimeoutExpired:code=124
            path=out/'deglutition_swallowing/deglutition_swallowing_results.json';raw=path.read_bytes() if path.exists() else b''
            record=json.loads(raw) if raw else {};accepted=False
            try:accepted=assess(record,names).accepted
            except ValueError:pass
            values=record.get('gates',{})
            rows.append(dict(source_hash=PIN,returncode=code,output_bytes=len(raw),output_hash=hashlib.sha256(raw).hexdigest(),
                             gate_values=[values.get(k) for k in names],gate_key_hash=hashlib.sha256(json.dumps(names).encode()).hexdigest(),
                             accepted=accepted,finite=record.get('nan_inf_free_check') is True))
            records.append(raw)
    gates=dict(full_repeat=rows[0]==rows[1] and records[0]==records[1] and bool(records[0]),required=all(r['accepted'] for r in rows),finite=all(r['finite'] for r in rows),completed=all(r['returncode']==0 for r in rows),controls=controls())
    report=dict(rows=rows,gates=gates,scope='Pinned numerical gate replay, excluding intentional string sentinels from the source finite check; no biological or clinical revalidation.')
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
