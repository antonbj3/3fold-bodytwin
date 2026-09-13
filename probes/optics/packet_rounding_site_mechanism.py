"""Localise binary64 intermediate rounding on frozen packet-boundary inputs."""
from decimal import Decimal, localcontext
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];PACKET=2**30


def main():
    with np.load(ROOT/'reports/packet_rounding_mechanism_arrays.npz') as a:
        lengths=a['lengths'];cpu=a['terminal'][:,1];continuous=a['continuous']
    with np.load(ROOT/'reports/packet_rounding_decimal_arrays.npz') as a:reference=a['escaped_112_0']
    rows=[];arrays={}
    for leg in (0,1):
        rounded_exp=[];rounded_weight=[]
        with localcontext() as context:
            context.prec=112
            for length in lengths:
                e=(-Decimal.from_float(float(length))).exp()
                rounded_exp.append(float(e)*PACKET)
                rounded_weight.append(float(e*PACKET))
        values=dict(rounded_exp=np.array(rounded_exp),rounded_weight=np.array(rounded_weight),cpu_exp=continuous)
        table=[]
        for name,v in values.items():
            integers=np.floor(v+.5).astype(np.int64)
            arrays[f'{name}_values_{leg}']=v;arrays[f'{name}_integers_{leg}']=integers
            table.append(dict(route=name,decimal_integer_differences=int(np.count_nonzero(integers!=reference)),cpu_integer_differences=int(np.count_nonzero(integers!=cpu)),cpu_continuous_differences=int(np.count_nonzero(v!=continuous)),hash=hashlib.sha256(integers.tobytes()).hexdigest()))
        rows.append(table)
    gates=dict(full_repeat=rows[0]==rows[1] and all(arrays[k].tobytes()==arrays[k[:-1]+'1'].tobytes() for k in arrays if k.endswith('_0')),
               scaling_routes_exact=arrays['rounded_exp_values_0'].tobytes()==arrays['rounded_weight_values_0'].tobytes(),
               cpu_reconstructed=rows[0][2]['cpu_integer_differences']==0)
    report=dict(rows=rows,gates=gates,scope='Binary64 intermediate observer on saved scalar inputs; no GPU instruction attribution or baseline replacement.')
    (ROOT/'reports/packet_rounding_site.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/packet_rounding_site_arrays.npz',**arrays)
    print(json.dumps(report));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
