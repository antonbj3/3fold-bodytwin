"""Preregistered exact synthetic tetrahedral material-boundary checks."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.tetra_interfaces_v1 import tetra_interfaces
from tetra_interface_mechanism import fixtures


def leg():
    rows=[]
    for name,(v,t,labels) in fixtures().items():
        try:mesh=tetra_interfaces(v,t,labels)
        except ValueError as exc:
            rows.append(dict(case=name,rejected=True,reason=str(exc)));continue
        regions=[]
        for region in np.unique(labels):
            faces=np.concatenate((mesh['faces'][mesh['back']==region],mesh['faces'][mesh['front']==region][:,::-1]))
            p=mesh['vertices'][faces]
            six=float(np.einsum('ij,ij->i',p[:,0],np.cross(p[:,1],p[:,2])).sum())
            c=v[t[labels==region]]
            expected=float(np.abs(np.einsum('ij,ij->i',c[:,1]-c[:,0],np.cross(c[:,2]-c[:,0],c[:,3]-c[:,0]))).sum())
            edges=faces[:,[0,1,1,2,2,0]].reshape(-1,2)
            unique,inverse,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
            balance=np.zeros(len(unique),np.int64)
            np.add.at(balance,inverse,np.where(edges[:,0]<edges[:,1],1,-1))
            regions.append(dict(region=int(region),six_volume=six,expected_six_volume=expected,
                                bad_incidence=int(np.count_nonzero(count!=2)),unbalanced_edges=int(np.count_nonzero(balance))))
        rows.append(dict(case=name,rejected=False,triangles=len(mesh['faces']),regions=regions,
                         hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in mesh.items()}))
    return rows


if __name__=='__main__':
    a,b=leg(),leg();good=[r for r in a if not r['rejected']];by_name={r['case']:r for r in a}
    gates=dict(full_repeat=a==b,
               triangle_counts=all(not by_name[k]['rejected'] and by_name[k]['triangles']==n for k,n in dict(single=4,same=6,layers=7,inverted=7).items()),
               exact_volumes=all(r['six_volume']==r['expected_six_volume'] for x in good for r in x['regions']),
               oriented_manifold=all(r['bad_incidence']==r['unbalanced_edges']==0 for x in good for r in x['regions']),
               inversion_invariant=by_name['layers'].get('hashes')==by_name['inverted'].get('hashes'),
               invalid_rejected=all(by_name[k]['rejected'] for k in ('duplicate','triple')))
    report=dict(rows=a,gates=gates,scope='Synthetic conforming-cell control only; arbitrary nonadjacent overlaps not checked.')
    (ROOT/'reports/tetra_interface_probe.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));raise SystemExit(0 if all(gates.values()) else 2)
