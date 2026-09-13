"""Independent precision-stability observer for saved packet-rounding inputs."""
from decimal import Decimal, localcontext, ROUND_HALF_UP
import hashlib
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
PACKET=2**30


def main():
    with np.load(ROOT/'reports/packet_rounding_mechanism_arrays.npz') as a:
        lengths=a['lengths'];cpu=a['terminal'][:,1]
    with np.load(ROOT/'reports/packet_rounding_gpu_arrays.npz') as a:gpu=a['terminal_0'][:,1]
    rows=[];arrays={}
    for precision in (80,112):
        for leg in (0,1):
            rounded=[];margin=[]
            with localcontext() as ctx:
                ctx.prec=precision
                for length in lengths:
                    value=Decimal(PACKET)*(-Decimal.from_float(float(length))).exp()
                    integer=value.to_integral_value(rounding=ROUND_HALF_UP)
                    half=value.to_integral_value(rounding='ROUND_FLOOR')+Decimal('.5')
                    rounded.append(int(integer));margin.append(abs(value-half))
            escaped=np.array(rounded,np.int64)
            arrays[f'escaped_{precision}_{leg}']=escaped
            rows.append(dict(precision=precision,leg=leg,hash=hashlib.sha256(escaped.tobytes()).hexdigest(),
                             cpu_differences=int(np.count_nonzero(escaped!=cpu)),gpu_differences=int(np.count_nonzero(escaped!=gpu)),
                             max_cpu_difference=int(np.max(np.abs(escaped-cpu))),max_gpu_difference=int(np.max(np.abs(escaped-gpu))),
                             escaped=int(escaped.sum()),absorbed=int((PACKET-escaped).sum()),minimum_half_margin=str(min(margin))))
    gates=dict(repeat=rows[0]['hash']==rows[1]['hash'] and rows[2]['hash']==rows[3]['hash'],
               precision_stability=rows[0]['hash']==rows[2]['hash'],
               energy=all(r['escaped']+r['absorbed']==len(lengths)*PACKET for r in rows),
               difference_bound=all(max(r['max_cpu_difference'],r['max_gpu_difference'])<=1 for r in rows))
    report=dict(rows=rows,gates=gates,scope='Precision stability on exact represented lengths; no baseline changes or arbitrary-input rounding proof.')
    (ROOT/'reports/packet_rounding_decimal.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/packet_rounding_decimal_arrays.npz',**arrays)
    print(json.dumps(report));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
