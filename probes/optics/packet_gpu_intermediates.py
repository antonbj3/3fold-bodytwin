"""Observe GPU packet arithmetic stages beside the frozen attenuation kernel."""
import hashlib,json
from pathlib import Path
import numpy as np
import warp as wp
ROOT=Path(__file__).resolve().parents[2];PACKET=2**30


@wp.kernel
def stages(lengths:wp.array(dtype=wp.float64),values:wp.array2d(dtype=wp.float64),integers:wp.array(dtype=wp.int64)):
    i=wp.tid()
    e=wp.exp(-lengths[i]);weight=wp.float64(1073741824.)*e
    shifted=weight+wp.float64(.5)
    values[i,0]=e;values[i,1]=weight;values[i,2]=shifted
    integers[i]=wp.int64(wp.floor(shifted))


def main():
    with np.load(ROOT/'reports/packet_rounding_mechanism_arrays.npz') as a:
        lengths=a['lengths'];cpu=a['terminal'][:,1];continuous=a['continuous']
    with np.load(ROOT/'reports/packet_rounding_gpu_arrays.npz') as a:frozen=a['terminal_0'][:,1]
    rows=[];arrays={}
    for leg in (0,1):
        inp=wp.array(lengths,dtype=wp.float64,device='cuda:0')
        values=wp.zeros((len(lengths),3),dtype=wp.float64,device='cuda:0');integers=wp.zeros(len(lengths),dtype=wp.int64,device='cuda:0')
        wp.launch(stages,len(lengths),inputs=[inp,values,integers],device='cuda:0')
        v=values.numpy();i=integers.numpy();difference=i!=cpu
        ulp=np.abs(v[:,1].copy().view(np.int64)-continuous.view(np.int64))
        rows.append(dict(value_hash=hashlib.sha256(v.tobytes()).hexdigest(),integer_hash=hashlib.sha256(i.tobytes()).hexdigest(),
                         frozen_integer_differences=int(np.count_nonzero(i!=frozen)),cpu_integer_differences=int(difference.sum()),
                         continuous_differences=int(np.count_nonzero(v[:,1]!=continuous)),maximum_ulp_difference=int(ulp.max()),
                         scaling_exact=bool(np.array_equal(v[:,0]*PACKET,v[:,1])),
                         all_integer_differences_have_continuous_difference=bool(np.all(v[difference,1]!=continuous[difference]))))
        arrays[f'values_{leg}']=v;arrays[f'integers_{leg}']=i
    gates=dict(full_repeat=rows[0]==rows[1],frozen_equivalence=all(r['frozen_integer_differences']==0 for r in rows),
               scaling=all(r['scaling_exact'] for r in rows),localisation=all(r['all_integer_differences_have_continuous_difference'] for r in rows))
    report=dict(rows=rows,gates=gates,scope='Captured isolated arithmetic stages only; no compiler-instruction cause or universal transcendental bound.')
    (ROOT/'reports/packet_gpu_intermediates.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/packet_gpu_intermediates_arrays.npz',**arrays)
    print(json.dumps(report));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
