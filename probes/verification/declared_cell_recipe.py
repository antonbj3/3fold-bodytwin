"""Replay explicitly declared all-required cell contracts without schema inference."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

def load_spec(path,root):
    data=json.loads(path.read_text())
    if not isinstance(data,dict) or set(data)!={'cells'} or not isinstance(data['cells'],list) or not data['cells']:
        raise ValueError('Nonempty explicit cell list required')
    cells={}
    fields={'path','source_hash','gate_count','overall_key','gate_variable','output_filename','gate_key_hash'}
    for item in data['cells']:
        if not isinstance(item,dict) or set(item)!=fields:raise ValueError('Unexpected cell declaration fields')
        name=item['path']
        if type(name) is not str or not name or name in cells:raise ValueError('Unique relative cell path required')
        script=(root/name).resolve()
        if Path(name).is_absolute() or not script.is_relative_to(root.resolve()) or not script.is_file():
            raise ValueError('Cell source escapes declared root or is missing')
        for key in ('source_hash','gate_key_hash'):
            value=item[key]
            if type(value) is not str or len(value)!=64 or any(c not in '0123456789abcdef' for c in value):
                raise ValueError('Explicit SHA256 required')
        if type(item['gate_count']) is not int or item['gate_count']<=0:raise ValueError('Positive integer gate count required')
        for key in ('overall_key','gate_variable','output_filename'):
            if type(item[key]) is not str or not item[key]:raise ValueError('Explicit nonempty field required')
        filename=item['output_filename']
        if Path(filename).name!=filename or filename in ('.','..'):raise ValueError('Output filename must be a basename')
        raw=script.read_bytes()
        keys=keys_from_source(raw.decode(),item['gate_variable'])
        if hashlib.sha256(raw).hexdigest()!=item['source_hash'] or len(keys)!=item['gate_count'] or hashlib.sha256(json.dumps(keys).encode()).hexdigest()!=item['gate_key_hash']:
            raise ValueError('Frozen source or gate-key binding mismatch')
        cells[name]=(item['source_hash'],item['gate_count'],item['overall_key'],item['gate_variable'],filename)
    return cells




def keys_from_source(source,variable):
    matches=[]
    for node in ast.walk(ast.parse(source)):
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id==variable for t in node.targets) and isinstance(node.value,ast.Dict):
            matches.append(sorted(ast.literal_eval(k) for k in node.value.keys))
    if len(matches)!=1 or any(type(k) is not str for k in matches[0]) or len(set(matches[0]))!=len(matches[0]):
        raise ValueError('One unambiguous literal gate dictionary required')
    return matches[0]


def assess(data,keys,overall='required_gates_overall_pass'):
    from bodytwin.framework.numerical_gate_contract_v1 import GateContract
    verdict=GateContract(tuple(keys),tuple(keys),overall).assess(data)
    if not verdict.accepted:raise ValueError('Required numerical gate failed')


def controls(spec,root):
    import copy
    data=json.loads(spec.read_text());out={}
    with tempfile.TemporaryDirectory(prefix='cell-spec-controls-') as directory:
        path=Path(directory)/'spec.json'
        for kind in ('empty','duplicate','source_pin','key_pin','count_type','output_path','source_path','missing','extra'):
            changed=copy.deepcopy(data);item=changed['cells'][0]
            if kind=='empty':changed['cells']=[]
            elif kind=='duplicate':changed['cells'].append(copy.deepcopy(item))
            elif kind=='source_pin':item['source_hash']='0'*64
            elif kind=='key_pin':item['gate_key_hash']='0'*64
            elif kind=='count_type':item['gate_count']=True
            elif kind=='output_path':item['output_filename']='../report.json'
            elif kind=='source_path':item['path']='../../outside.py'
            elif kind=='missing':del item['overall_key']
            else:item['unknown']=True
            path.write_text(json.dumps(changed))
            try:load_spec(path,root);out[kind]=False
            except ValueError:out[kind]=True
    try:keys_from_source("gates={'a': True}\ngates={'b': False}",'gates');out['ambiguous']=False
    except ValueError:out['ambiguous']=True
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells-root',type=Path)
    p.add_argument('--spec',type=Path,required=True)
    p.add_argument('--controls-only',action='store_true')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    root=a.cells_root or Path(__file__).resolve().parents[2]/'src/bodytwin/cells'
    cells=load_spec(a.spec,root)
    if a.controls_only:
        report=dict(gates=controls(a.spec,root))
        a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
        return 0 if all(report['gates'].values()) else 2
    rows=[];records=[]
    with tempfile.TemporaryDirectory(prefix='required-cells-') as directory:
        for leg in (0,1):
            out=Path(directory)/str(leg);out.mkdir()
            env=os.environ.copy();env.update(BODYTWIN_OUT=str(out),CUDA_VISIBLE_DEVICES='',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            leg_rows=[];leg_records=[]
            for name,(pin,count,overall,variable,filename) in cells.items():
                script=root/name;raw=script.read_bytes()
                if hashlib.sha256(raw).hexdigest()!=pin:raise ValueError('Frozen cell source changed')
                keys=keys_from_source(raw.decode(),variable)
                if len(keys)!=count:raise ValueError('Unexpected source gate count')
                try:
                    child=subprocess.run([sys.executable,str(script.resolve())],cwd=out,env=env,capture_output=True,text=True,timeout=200)
                    code=child.returncode
                except subprocess.TimeoutExpired:
                    code=124
                path=out/script.stem/filename
                record=path.read_bytes() if path.is_file() else b''
                data=json.loads(record) if record else {}
                accepted=False
                try:assess(data,keys,overall);accepted=True
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
