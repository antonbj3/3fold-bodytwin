"""Exercise matrix integrity rejection against complete committed photon captures."""
import copy
import json
from pathlib import Path
import numpy as np
from photon_modal_suite import validate_arrays
ROOT=Path(__file__).resolve().parents[2]


def leg():
    rows=[]
    for family,folder,tag in [('ocular','ocular_cross','ocular_transport_10000'),('scattering','photon_cross','photon_cross_capture')]:
        root=ROOT/'reports'/folder/'L4';report=json.loads((root/(tag+'.json')).read_text())
        with np.load(root/(tag+'_arrays.npz')) as archive:data={k:archive[k] for k in archive.files}
        valid=validate_arrays(family,report,data);key=next(k for k in data if k.endswith('_terminal'))
        modified=dict(data);modified[key]=data[key].copy();modified[key][0,0]+=1
        missing=dict(data);missing.pop(key)
        duplicate=copy.deepcopy(report);duplicate['rows'][1]=copy.deepcopy(duplicate['rows'][0])
        controls={}
        for name,r,a in [('modified',report,modified),('missing',report,missing),('duplicate',duplicate,data)]:
            try:validate_arrays(family,r,a);controls[name]=False
            except ValueError:controls[name]=True
        rows.append(dict(family=family,valid=valid,rejections=controls))
    return rows


def main():
    first=leg();second=leg();gates=dict(valid=all(r['valid'] for r in first),tamper=all(all(r['rejections'].values()) for r in first),repeat=json.dumps(first,sort_keys=True)==json.dumps(second,sort_keys=True))
    report=dict(gates=gates,rows=first,rejected_controls=sum(sum(r['rejections'].values()) for r in first))
    (ROOT/'reports/photon_suite_integrity.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
