"""Audit complete frozen ocular arrays across explicitly captured CUDA backends."""
import hashlib
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]


def compare(a,b):
    if a.shape!=b.shape or a.dtype!=b.dtype:
        return dict(byte_identical=False,shape_dtype_match=False,differing_elements=None,max_abs_difference=None)
    return dict(byte_identical=a.tobytes()==b.tobytes(),shape_dtype_match=True,
                differing_elements=int(np.count_nonzero(a!=b)),
                max_abs_difference=float(np.max(np.abs(a.astype(np.float64)-b.astype(np.float64)))) if a.size else 0.)


def main():
    folder=ROOT/'reports/ocular_cross';captures={};reports={};arrays={}
    for backend in ('L4','A10G','H100'):
        sub=folder/backend;captures[backend]=json.loads((sub/'capture.json').read_text());reports[backend]=json.loads((sub/'ocular_transport_10000.json').read_text())
        with np.load(sub/'ocular_transport_10000_arrays.npz') as archive:arrays[backend]={k:archive[k] for k in archive.files}
    expected={n:hashlib.sha256((ROOT/n.replace('examples/optics', 'probes/optics')).read_bytes().replace(b'probes/optics', b'examples/optics')).hexdigest() for n in ('src/bodytwin/geometry/optics/ocular_photon_v1.py','examples/optics/ocular_transport_probe.py')}
    rows=[]
    for backend in ('A10G','H100'):
        for key,a in arrays['L4'].items():
            b=arrays[backend].get(key)
            rows.append(dict(reference='L4',backend=backend,array=key,**(compare(a,b) if b is not None else dict(byte_identical=False,shape_dtype_match=False,differing_elements=None,max_abs_difference=None))))
    repeat=[]
    for backend,data in arrays.items():
        for key,a in data.items():
            if '_0_' in key:
                repeat.append(dict(backend=backend,array=key,**compare(a,data[key.replace('_0_','_1_',1)])))
    gates=dict(capture_valid=all(c['exit']==0 and c['source_hashes']==expected for c in captures.values()) and
               all(r['reference_sha256']=='656a2caccfb6ad12ab154138799821e9a5dfb92021cf4a5c68e1a7cbdc076ba8' and r['transport_sha256']==expected['src/bodytwin/geometry/optics/ocular_photon_v1.py'] and r['photons_per_arm']==10000 for r in reports.values()) and
               all(token in captures[b]['inventory'] for b,token in [('L4','L4'),('A10G','A10'),('H100','H100')]),
               original_gates=all(all(r['gates'].values()) for r in reports.values()),
               same_backend_repeat=all(r['byte_identical'] for r in repeat),
               cross_backend_exact=all(r['byte_identical'] for r in rows) and all(set(a)==set(arrays['L4']) for a in arrays.values()))
    result=dict(gates=gates,rows=rows,repeat_rows=repeat,source_hashes=expected,
                scope='Two 10000-packet runs per wavelength and remote backend; exact-byte gate; local backend missing; no throughput, anatomy or clinical claim.')
    (ROOT/'reports/ocular_cross_backend.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(gates=gates,failed_rows=[r for r in rows if not r['byte_identical']])))
    return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
