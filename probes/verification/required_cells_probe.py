"""Replay pinned numerical cells with explicit all-required gate contracts."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

CELLS={
    'endocrine/glucose_meal_dallaman2007.py': ('1c96a5970119e2ae2bc59b26f547b950ca425c35f102d476717de8ced97b1500',19),
    'haematology/iron_hepcidin_dynamic_4state.py': ('404179f718a47a8f9a14c83ff47ab9e0b0264b61b5b6f0ac792acc496ac8c2b7',7)}


def keys_from_source(source):
    for node in ast.walk(ast.parse(source)):
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='GATES' for t in node.targets) and isinstance(node.value,ast.Dict):
            return sorted(ast.literal_eval(k) for k in node.value.keys)
    raise ValueError('Explicit source gate dictionary missing')


def assess(data,keys):
    gates=data.get('gates')
    if not isinstance(gates,dict) or sorted(gates)!=keys or any(type(v) is not bool for v in gates.values()):
        raise ValueError('Gate inventory or type mismatch')
    if data.get('required_gates_overall_pass') is not True or not all(gates.values()):
        raise ValueError('Required numerical gate failed')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells-root',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    root=a.cells_root or Path(__file__).resolve().parents[2]/'src/bodytwin/cells'
    rows=[];records=[]
    with tempfile.TemporaryDirectory(prefix='required-cells-') as directory:
        for leg in (0,1):
            out=Path(directory)/str(leg);out.mkdir()
            env=os.environ.copy();env.update(BODYTWIN_OUT=str(out),CUDA_VISIBLE_DEVICES='',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            leg_rows=[];leg_records=[]
            for name,(pin,count) in CELLS.items():
                script=root/name;raw=script.read_bytes()
                if hashlib.sha256(raw).hexdigest()!=pin:raise ValueError('Frozen cell source changed')
                keys=keys_from_source(raw.decode())
                if len(keys)!=count:raise ValueError('Unexpected source gate count')
                try:
                    child=subprocess.run([sys.executable,str(script.resolve())],cwd=out,env=env,capture_output=True,text=True,timeout=200)
                    code=child.returncode
                except subprocess.TimeoutExpired:
                    code=124
                path=out/script.stem/f'{script.stem}_results.json'
                record=path.read_bytes() if path.is_file() else b''
                data=json.loads(record) if record else {}
                accepted=False
                try:assess(data,keys);accepted=True
                except ValueError:pass
                values=data.get('gates',{})
                leg_rows.append(dict(module='src/bodytwin/cells/'+name,source_hash=pin,returncode=code,
                                     output_hash=hashlib.sha256(record).hexdigest(),output_bytes=len(record),
                                     gate_values=[values.get(k) for k in keys],gate_key_hash=hashlib.sha256(json.dumps(keys).encode()).hexdigest(),
                                     required_pass=accepted))
                leg_records.append(record)
                print(json.dumps(dict(leg=leg,module=name,returncode=code,required_pass=accepted)),flush=True)
            rows.append(leg_rows);records.append(leg_records)
    refused=0
    for bad in ({},{'gates':{'a':'true'},'required_gates_overall_pass':True},
                {'gates':{'a':False},'required_gates_overall_pass':True},
                {'gates':{'a':True},'required_gates_overall_pass':False}):
        try:assess(bad,['a'])
        except ValueError:refused+=1
    gates=dict(full_repeat=rows[0]==rows[1] and records[0]==records[1] and all(records[0]),
               required=all(r['required_pass'] for leg in rows for r in leg),
               completed=all(r['returncode']==0 for leg in rows for r in leg),refusals=refused==4)
    report=dict(rows=rows,gates=gates,scope='Pinned source-defined numerical gate replay with complete output hashes; no fresh literature/clinical review.')
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
