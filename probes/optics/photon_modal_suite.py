"""Bounded remote photon determinism matrix; no scheduler or local CUDA context."""
from pathlib import Path
import hashlib
import json
import modal

REFERENCE_SHA='656a2caccfb6ad12ab154138799821e9a5dfb92021cf4a5c68e1a7cbdc076ba8'
SOURCES=('src/bodytwin/geometry/optics/ocular_photon_v1.py',
         'src/bodytwin/geometry/optics/tissue_photon_mc_v6.py',
         'probes/optics/ocular_transport_probe.py',
         'probes/optics/photon_cross_capture.py')
FAMILIES={'ocular':('ocular_transport_probe.py','ocular_transport_10000'),
          'scattering':('photon_cross_capture.py','photon_cross_capture')}
app=modal.App('photon-determinism-suite')
image=modal.Image.from_registry('nvidia/cuda:12.9.1-devel-ubuntu24.04',add_python='3.12').pip_install(
    'numpy==2.5.3','scipy==1.18.1','trimesh==5.1.0','warp-lang==1.13.0')


def capture(sources,reference):
    import os
    import subprocess
    import tempfile
    if hashlib.sha256(reference).hexdigest()!=REFERENCE_SHA:
        raise ValueError('Refraction oracle hash mismatch')
    inventory=subprocess.run(['nvidia-smi','--query-gpu=name','--format=csv,noheader'],capture_output=True,text=True,check=True,timeout=15).stdout.strip()
    with tempfile.TemporaryDirectory() as folder:
        root=Path(folder);(root/'reports').mkdir()
        for name,data in sources.items():
            if name not in SOURCES:raise ValueError('Source outside fixed manifest')
            target=root/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
        (root/'reference.py').write_bytes(reference)
        env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OCULAR_PHOTONS='10000',LENS_REFERENCE_SOURCE=str(root/'reference.py'))
        exits={}
        for family,(script,_) in FAMILIES.items():
            child=subprocess.run(['python',str(root/'probes/optics'/script)],env=env,capture_output=True,text=True,timeout=240)
            exits[family]=child.returncode
            if child.returncode:print(child.stderr)
        return dict(inventory=inventory,exits=exits,source_hashes={k:hashlib.sha256(v).hexdigest() for k,v in sources.items()},reference_sha256=hashlib.sha256(reference).hexdigest(),outputs={p.name:p.read_bytes() for p in (root/'reports').iterdir()})


CAPTURES={sku:app.function(image=image,gpu=sku,cpu=8,memory=16384,timeout=540,startup_timeout=60,retries=0,name='capture_'+sku)(capture) for sku in ('L4','A10G','H100')}


def validate_arrays(family,report,data):
    keys=('absorption','terminal','counters','paths','exits','history','retinal_histogram') if family=='ocular' else ('absorption','terminal','counters')
    arms=(532,650) if family=='ocular' else (0,1)
    expected={f'{arm}_{leg}_{key}' for arm in arms for leg in (0,1) for key in keys}
    if family=='ocular':
        expected|={f'{arm}_reference_{key}' for arm in arms for key in ('positions','directions')}
        expected|={f'{arm}_direct_{key}' for arm in arms for key in ('probability','attenuation')}
    if set(data)!=expected or len(report['rows'])!=4:raise ValueError('Missing or unexpected array inventory')
    identities=set()
    for row in report['rows']:
        arm=row['wavelength'] if family=='ocular' else int(row['reflect']);leg=row['leg']
        identities.add((arm,leg))
        if set(row['hashes'])!=set(keys):raise ValueError('Incomplete hash manifest')
        for key in keys:
            array=data[f'{arm}_{leg}_{key}']
            if hashlib.sha256(array.tobytes()).hexdigest()!=row['hashes'][key]:raise ValueError('Array hash mismatch')
        if data[f'{arm}_{leg}_terminal'].shape!=(10000,4):raise ValueError('Wrong packet count')
    if identities!={(a,l) for a in arms for l in (0,1)}:raise ValueError('Missing arm or repeat')
    return True


def matrix(folder,expected):
    import numpy as np
    captures={};arrays={};reports={};rows=[];repeats=[]
    for backend in CAPTURES:
        sub=folder/backend;captures[backend]=json.loads((sub/'capture.json').read_text())
        for family,(_,tag) in FAMILIES.items():
            reports[backend,family]=json.loads((sub/(tag+'.json')).read_text())
            with np.load(sub/(tag+'_arrays.npz')) as archive:arrays[backend,family]={k:archive[k] for k in archive.files}
            validate_arrays(family,reports[backend,family],arrays[backend,family])
    def compare(a,b):
        if a.shape!=b.shape or a.dtype!=b.dtype:return dict(exact=False,differences=None,max_abs=None)
        return dict(exact=a.tobytes()==b.tobytes(),differences=int(np.count_nonzero(a!=b)),max_abs=float(np.max(np.abs(a.astype(float)-b.astype(float)))) if a.size else 0.)
    for family in FAMILIES:
        baseline=arrays['L4',family]
        for backend in CAPTURES:
            data=arrays[backend,family]
            if set(data)!=set(baseline):raise ValueError('Incomplete array capture')
            for key,a in data.items():
                if backend!='L4':rows.append(dict(family=family,backend=backend,array=key,**compare(baseline[key],a)))
                if '_0_' in key:repeats.append(dict(family=family,backend=backend,array=key,**compare(a,data[key.replace('_0_','_1_',1)])))
    artifact_names={'capture.json'}|{tag+suffix for _,tag in FAMILIES.values() for suffix in ('.json','_arrays.npz')}
    gates=dict(artifact_integrity=all((folder/backend).is_dir() and {p.name for p in (folder/backend).iterdir()}==artifact_names and all(p.stat().st_size>0 for p in (folder/backend).iterdir()) for backend in CAPTURES),architectures=all(token in captures[b]['inventory'] for b,token in [('L4','L4'),('A10G','A10'),('H100','H100')]),
               sources=all(c['source_hashes']==expected and c['reference_sha256']==REFERENCE_SHA and all(v==0 for v in c['exits'].values()) for c in captures.values()),
               original_gates=all(all(r['gates'].values()) for r in reports.values()),
               same_backend=all(r['exact'] for r in repeats),cross_backend=all(r['exact'] for r in rows))
    return dict(gates=gates,rows=rows,repeats=repeats,source_hashes=expected,reference_sha256=REFERENCE_SHA,
                runtime=dict(python='3.12',numpy='2.5.3',scipy='1.18.1',trimesh='5.1.0',warp='1.13.0',cuda='12.9.1'),
                scope='Two synthetic photon families,10000 packets per arm/leg,three remote architectures. No local capture, scheduler, throughput or anatomical claim.')


@app.local_entrypoint()
def main(lens_reference:str,output_dir:str):
    root=Path(__file__).resolve().parents[2]
    reference=Path(lens_reference).read_bytes()
    if hashlib.sha256(reference).hexdigest()!=REFERENCE_SHA:raise ValueError('Refraction oracle hash mismatch')
    folder=Path(output_dir);folder.mkdir(parents=True,exist_ok=False)
    sources={n:(root/n).read_bytes() for n in SOURCES};expected={n:hashlib.sha256(data).hexdigest() for n,data in sources.items()}
    for backend,function in CAPTURES.items():
        result=function.remote(sources,reference);sub=folder/backend;sub.mkdir()
        for name,data in result.pop('outputs').items():(sub/name).write_bytes(data)
        (sub/'capture.json').write_text(json.dumps(result,indent=2)+'\n')
        if any(result['exits'].values()):raise RuntimeError('Capture failed; partial artifacts preserved')
    first=matrix(folder,expected);second=matrix(folder,expected)
    a=json.dumps(first,sort_keys=True,indent=2)+'\n';b=json.dumps(second,sort_keys=True,indent=2)+'\n'
    (folder/'matrix.json').write_text(a);(folder/'matrix_repeat.json').write_text(b)
    print(json.dumps(dict(gates=first['gates'],matrix_bytes_exact=a==b)))
    if a!=b or not all(first['gates'].values()):raise RuntimeError('Strict matrix gate failed; artifacts preserved')
