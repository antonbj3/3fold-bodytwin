"""Represented-input rounding boundary table before a cross-reference claim."""
import hashlib
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];PACKET=1<<30


def inputs():
    k=np.random.default_rng(20260913).integers(1,PACKET-1,2048,dtype=np.int64)
    central=-np.log((k.astype(float)+.5)/PACKET)
    lengths=np.column_stack((np.nextafter(central,-np.inf),central,np.nextafter(central,np.inf))).ravel()
    targets=np.repeat(k.astype(float)+.5,3)
    return lengths,targets


def measure():
    lengths,targets=inputs();continuous=np.array([float(PACKET)*np.exp(-float(length)) for length in lengths]);escaped=np.floor(continuous+.5).astype(np.int64)
    terminal=np.column_stack((PACKET-escaped,escaped,np.zeros(len(escaped),np.int64)))
    delta=continuous-targets
    arrays=dict(lengths=lengths,targets=targets,continuous=continuous,terminal=terminal)
    return dict(intervals=len(lengths),maximum_half_offset=float(np.max(np.abs(delta))),below_half=int(np.count_nonzero(delta<0)),exact_half=int(np.count_nonzero(delta==0)),above_half=int(np.count_nonzero(delta>0)),absorbed=int(terminal[:,0].sum()),escaped=int(escaped.sum()),hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in arrays.items()}),arrays


def main():
    first,arrays=measure();second,again=measure();gates=dict(repeat=first==second and all(a.tobytes()==again[k].tobytes() for k,a in arrays.items()),energy=first['absorbed']+first['escaped']==first['intervals']*PACKET,near_boundary=first['maximum_half_offset']<=1e-6)
    report=dict(measurement=first,gates=gates,scope='Isolated represented interval rounding inputs; no full geometric path or clinical claim.')
    (ROOT/'reports/packet_rounding_mechanism.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/packet_rounding_mechanism_arrays.npz',**arrays);print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
