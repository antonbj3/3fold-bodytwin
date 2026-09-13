"""Full-ledger cross-backend scattering audit over frozen captures."""
import hashlib
import json
from pathlib import Path
import numpy as np
from ocular_cross_backend_probe import compare
ROOT=Path(__file__).resolve().parents[2]


def main():
    captures={};reports={};arrays={}
    for backend in ('L4','A10G','H100'):
        sub=ROOT/'reports/photon_cross'/backend
        captures[backend]=json.loads((sub/'capture.json').read_text());reports[backend]=json.loads((sub/'photon_cross_capture.json').read_text())
        with np.load(sub/'photon_cross_capture_arrays.npz') as archive:arrays[backend]={k:archive[k] for k in archive.files}
    expected={n:hashlib.sha256((ROOT/n.replace('examples/optics', 'probes/optics')).read_bytes().replace(b'probes/optics', b'examples/optics')).hexdigest() for n in ('src/bodytwin/geometry/optics/tissue_photon_mc_v6.py','examples/optics/photon_cross_capture.py')}
    rows=[];repeat=[]
    for backend in ('A10G','H100'):
        for key,a in arrays['L4'].items():
            rows.append(dict(reference='L4',backend=backend,array=key,**compare(a,arrays[backend][key])))
    for backend,data in arrays.items():
        for reflect in (0,1):
            for key in ('absorption','terminal','counters'):
                repeat.append(dict(backend=backend,reflect=reflect,array=key,**compare(data[f'{reflect}_0_{key}'],data[f'{reflect}_1_{key}'])))
    gates=dict(capture_valid=all(c['exit']==0 and c['source_hashes']==expected for c in captures.values()) and all(r['source_sha256']==expected['src/bodytwin/geometry/optics/tissue_photon_mc_v6.py'] and r['photons']==10000 for r in reports.values()) and all(token in captures[b]['inventory'] for b,token in [('L4','L4'),('A10G','A10'),('H100','H100')]),
               physical_gates=all(all(r['gates'].values()) for r in reports.values()),same_backend_repeat=all(r['byte_identical'] for r in repeat),cross_backend_exact=all(r['byte_identical'] for r in rows) and all(set(a)==set(arrays['L4']) for a in arrays.values()))
    result=dict(gates=gates,rows=rows,repeat_rows=repeat,source_hashes=expected,scope='Frozen tissue scattering and reflection, 10000 packets per arm/leg; no local capture or expanded physical/statistical claim.')
    (ROOT/'reports/photon_cross_backend.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(gates=gates,failed_rows=[r for r in rows if not r['byte_identical']])));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
