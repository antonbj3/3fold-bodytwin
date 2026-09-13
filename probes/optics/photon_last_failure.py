"""Certify observational equivalence while classifying failed photon paths."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh
import warp as wp
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from tissue_photon_last_failure import simulate


def sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def main():
    wp.init();m=trimesh.creation.box(extents=[60,60,60]);m.apply_translation([30,30,30])
    base=json.loads((ROOT/'reports/photon_cube_v2_l4.json').read_text())['legs'][0]['hashes']
    rows=[]
    for _ in range(2):
        out=simulate(m.vertices,m.faces,[[0,0,1,1],[.005,1,.01,1.37]],np.zeros(12,dtype=int),np.ones(12,dtype=int),
            source=[29,29,.0001],direction=[0,0,1],initial_region=1,origin=[0,0,0],shape=[60,60,60])
        terminal=out['terminal'];canonical=terminal[:,:4].copy();canonical[:,3]=(canonical[:,3]!=0)
        actual={'absorption':sha(out['absorption']),'terminal':sha(canonical),'counters':sha(out['counters'])}
        detail=[]
        for tid in np.flatnonzero(terminal[:,3]):
            t=terminal[tid]
            detail.append({'photon':int(tid),'reason':int(t[3]),'position_mm':(t[4:7]/1e9).tolist(),
                           'direction':(t[7:10]/1e9).tolist(),'ray_distance_mm':float(t[10]/1e9),'boundary_hit':int(t[11]),
                           'residual_integer':int(t[2])})
        rows.append({'sites':detail,'canonical_hashes':actual,'diagnostic_sha256':sha(terminal),
                     'baseline_identical':actual==base});print(json.dumps(rows[-1]),flush=True)
    gates={'two_complete_results_identical':rows[0]==rows[1],'frozen_paths_unchanged':all(r['baseline_identical'] for r in rows)}
    report={'legs':rows,'gates':gates,'reason_map':{'1':'region_label','2':'ray_distance','3':'grid_exit','4':'negative_segment','5':'segment_limit'}}
    (ROOT/'reports/photon_last_failure.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True);return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
