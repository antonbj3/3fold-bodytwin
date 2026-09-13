"""Measure strict nested-gate acceptance separately from complete record repeat."""
import argparse,ast,hashlib,json,os,subprocess,sys,tempfile
from pathlib import Path
from bodytwin.framework.numerical_gate_contract_v1 import GateContract
SOURCE='src/bodytwin/cells/endocrine/spermatogenesis_sertoli.py'
PIN='584ae4a26749ffbdf49aa972dd85dfe9408692862a083bd0f38b9561ffa9c8e4'
OVERALL='overall_pass_strict_all'


def assess(record,names):
    gates=record.get('gates')
    if not isinstance(gates,dict) or set(gates)!=set(names) or any(not isinstance(v,dict) or 'pass' not in v for v in gates.values()):
        raise ValueError('Exact nested gate inventory required')
    return GateContract(names,names,'overall').assess({'gates':{k:gates[k]['pass'] for k in names},'overall':record.get(OVERALL)})


def differences(a,b,path=''):
    if isinstance(a,dict) and isinstance(b,dict) and a.keys()==b.keys():
        return [p for k in a for p in differences(a[k],b[k],path+'/'+k)]
    if isinstance(a,list) and isinstance(b,list) and len(a)==len(b):
        return [p for i,(x,y) in enumerate(zip(a,b)) for p in differences(x,y,path+'/'+str(i))]
    return [] if a==b else [path]


def controls():
    bad=[{}, {'gates':{'a':True},OVERALL:True}, {'gates':{'a':{}},OVERALL:True},
         {'gates':{'a':{'pass':1}},OVERALL:True}, {'gates':{'a':{'pass':True}},OVERALL:1},
         {'gates':{'a':{'pass':False}},OVERALL:True}, {'gates':{'a':{'pass':True},'extra':{'pass':True}},OVERALL:True}]
    refused=0
    for record in bad:
        try:assess(record,('a',))
        except ValueError:refused+=1
    return refused==7 and not assess({'gates':{'a':{'pass':False}},OVERALL:False},('a',)).accepted


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    source=Path(__file__).resolve().parents[2]/SOURCE;raw=source.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PIN:raise ValueError('Frozen source changed')
    keys=[]
    for node in ast.parse(raw).body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Subscript) and isinstance(target.value,ast.Name) and target.value.id=='gates':keys.append(ast.literal_eval(target.slice))
    if len(keys)!=21 or len(set(keys))!=21:raise ValueError('Nested source inventory changed')
    names=tuple(sorted(keys));rows=[];records=[];parsed=[]
    with tempfile.TemporaryDirectory(prefix='nested-gates-') as directory:
        for leg in (0,1):
            out=Path(directory)/str(leg);out.mkdir();env=os.environ.copy();env.update(BODYTWIN_OUT=str(out),CUDA_VISIBLE_DEVICES='',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            try:code=subprocess.run([sys.executable,str(source)],cwd=out,env=env,capture_output=True,timeout=200).returncode
            except subprocess.TimeoutExpired:code=124
            files=[out/'spermatogenesis_sertoli'/f'spermatogenesis_sertoli_{suffix}.json' for suffix in ('results','evidence')]
            pair=[p.read_bytes() if p.exists() else b'' for p in files];record=json.loads(pair[0]) if pair[0] else {};accepted=False
            try:accepted=assess(record,names).accepted
            except ValueError:pass
            rows.append(dict(source_hash=PIN,returncode=code,output_bytes=len(pair[0]),output_hash=hashlib.sha256(pair[0]).hexdigest(),alias_exact=pair[0]==pair[1] and bool(pair[0]),accepted=accepted,
                             passed=sum(v.get('pass') is True for v in record.get('gates',{}).values()),finite=record.get('_meta',{}).get('nan_inf_scan_clean') is True))
            records.append(pair[0]);parsed.append(record)
    changed=differences(*parsed)
    gates=dict(full_repeat=records[0]==records[1] and bool(records[0]),required=all(r['accepted'] for r in rows),finite=all(r['finite'] for r in rows),aliases=all(r['alias_exact'] for r in rows),completed=all(r['returncode']==0 for r in rows),controls=controls())
    report=dict(rows=rows,gates=gates,changed_leaf_paths=changed,timestamp_only_difference=changed==['/_meta/generated_utc'],scope='Full record identity remains required; timestamp-only diagnosis does not promote numerical or biological validity.')
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
