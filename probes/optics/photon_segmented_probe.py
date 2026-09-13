"""Compare segmented candidate enumeration on captured complete photon queries."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import warp as wp
from bodytwin.geometry.optics import tissue_photon_paths_v1 as baseline
from bodytwin.geometry.optics import tissue_photon_segmented_v1 as candidate

@wp.kernel
def reference(mesh:wp.uint64,q:wp.array2d(dtype=wp.float64),out:wp.array(dtype=int)):
    i=wp.tid()
    out[i]=baseline._candidate_hit(mesh,wp.vec3d(q[i,0],q[i,1],q[i,2]),wp.vec3d(q[i,3],q[i,4],q[i,5]),q[i,6])

@wp.kernel
def segmented(mesh:wp.uint64,q:wp.array2d(dtype=wp.float64),out:wp.array(dtype=int)):
    i=wp.tid()
    out[i]=candidate._candidate_hit(mesh,wp.vec3d(q[i,0],q[i,1],q[i,2]),wp.vec3d(q[i,3],q[i,4],q[i,5]),q[i,6])


def run(fixture):
    data=np.load(fixture)
    v,f,q,expected=(data[k] for k in ('vertices','faces','query','expected'))
    before={k:hashlib.sha256(data[k].tobytes()).hexdigest() for k in data.files}
    mesh=wp.Mesh(points=wp.array(v,dtype=wp.vec3),indices=wp.array(f.ravel(),dtype=int))
    queries=wp.array(q,dtype=wp.float64);out=wp.zeros(len(q),dtype=int)
    rows=[];captures=[]
    for name,kernel in (('baseline',reference),('candidate',segmented)):
        def launch():wp.launch(kernel,len(q),inputs=[mesh.id,queries,out])
        launch();wp.synchronize()
        legs=[];arrays=[]
        for leg in range(2):
            wp.synchronize();start=time.perf_counter()
            for _ in range(2):launch()
            wp.synchronize();legs.append((time.perf_counter()-start)*1000/2)
            arrays.append(out.numpy())
        rows.append(dict(name=name,milliseconds=legs,mismatches=[int(np.count_nonzero(a!=expected)) for a in arrays],hashes=[hashlib.sha256(a.tobytes()).hexdigest() for a in arrays]))
        captures.append(arrays)
    gates=dict(query_parity=all(n==0 for row in rows for n in row['mismatches']),
               full_repeat=all(pair[0].tobytes()==pair[1].tobytes() for pair in captures),
               input_immutable=before=={k:hashlib.sha256(data[k].tobytes()).hexdigest() for k in data.files} and queries.numpy().tobytes()==q.tobytes(),
               faster_both=max(rows[1]['milliseconds'])<min(rows[0]['milliseconds']))
    return dict(rows=rows,gates=gates,queries=len(q),overall_pass=all(gates.values()))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fixture',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    report=run(a.fixture);a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));raise SystemExit(0 if report['overall_pass'] else 2)
