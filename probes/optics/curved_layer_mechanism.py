"""Measure curved nested surface discretization before a labelled layer fixture."""
import hashlib
import json
from pathlib import Path
import numpy as np
import trimesh
ROOT=Path(__file__).resolve().parents[2]


def measure():
    rows=[]
    for sections in (32,64,128):
        surfaces=[trimesh.creation.cylinder(radius=r,height=h,sections=sections) for r,h in ((30.,120.),(27.,114.))]
        ideal=[np.pi*30**2*120,np.pi*27**2*114]
        actual=[float(m.volume) for m in surfaces]
        errors=[abs(actual[1]-ideal[1])/ideal[1],abs((actual[0]-actual[1])-(ideal[0]-ideal[1]))/(ideal[0]-ideal[1])]
        rows.append(dict(sections=sections,vertices=[len(m.vertices) for m in surfaces],faces=[len(m.faces) for m in surfaces],
                         hashes=[hashlib.sha256(m.vertices.tobytes()+m.faces.tobytes()).hexdigest() for m in surfaces],
                         closed_oriented=all(m.is_watertight and m.is_winding_consistent and m.volume>0 for m in surfaces),
                         volumes=actual,shell_volume=actual[0]-actual[1],region_relative_errors=errors,
                         radial_sagitta=30*(1-np.cos(np.pi/sections)),
                         conservative_clearance=min(30*np.cos(np.pi/sections)-27,3.),
                         candidate_volume_gate=all(e<=.002 for e in errors)))
    return rows


if __name__=='__main__':
    a,b=measure(),measure()
    gates=dict(full_repeat=a==b,closed_oriented=all(r['closed_oriented'] for r in a),
               nested_positive=all(r['shell_volume']>0 and r['conservative_clearance']>0 for r in a))
    report=dict(rows=a,gates=gates,scope='Synthetic nested geometry only; no anatomy/transport certificate.')
    (ROOT/'reports/curved_layer_mechanism.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));raise SystemExit(0 if all(gates.values()) else 1)
