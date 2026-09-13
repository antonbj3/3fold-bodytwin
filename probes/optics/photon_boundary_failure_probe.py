"""Nonperturbing boundary/cap mechanism measurement against frozen v3 hashes."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh
import warp as wp
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from tissue_photon_boundary_failure import simulate


def sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def quantiles(a):return np.quantile(a,[0,.25,.5,.75,1]).tolist() if len(a) else []


def main():
    wp.init();baseline=json.loads((ROOT/'reports/photon_boundary_l4.json').read_text())
    cases=[]
    for case in ('cube60b','slab'):
        base=next(c for c in baseline['cases'] if c['case']==case)['legs'][0]['hashes']
        side=100 if case=='slab' else 60;xy=50 if case=='slab' else 29
        medium=[.01,10,.9,1] if case=='slab' else [.005,1,.01,1.37]
        mesh=trimesh.creation.box(extents=[side]*3);mesh.apply_translation([side/2]*3)
        rows=[]
        for _ in range(2):
            out=simulate(mesh.vertices,mesh.faces,[[0,0,1,1],medium],np.zeros(12,dtype=int),np.ones(12,dtype=int),
                source=[xy,xy,.0001],direction=[0,0,1],initial_region=1,origin=[0,0,0],shape=[side]*3,reflect=case=='cube60b')
            terminal=out['terminal'];canonical=terminal[:,:4].copy();canonical[:,3]=(canonical[:,3]!=0)
            hashes={'absorption':sha(out['absorption']),'terminal':sha(canonical),'counters':sha(out['counters'])}
            groups=[]
            for code in range(6):
                selected=(terminal[:,3]==code)&(terminal[:,2]>0)
                a=terminal[selected];ids=np.flatnonzero(selected)
                samples=[]
                for tid,t in zip(ids[:8],a[:8]):
                    samples.append({'photon':int(tid),'residual_integer':int(t[2]),'position':(t[4:7]/1e9).tolist(),
                        'direction':(t[7:10]/1e9).tolist(),'distance':float(t[10]/1e9),'hit':int(t[11]),
                        'events':int(t[12]),'hits':int(t[13]),'zero_losses':int(t[14]),'trailing_zero_losses':int(t[15])})
                groups.append({'code':code,'count':len(a),'residual_integer_sum':int(a[:,2].sum()),
                    'weight_quantiles':quantiles(a[:,2]),'event_quantiles':quantiles(a[:,12]),
                    'zero_loss_quantiles':quantiles(a[:,14]),'trailing_zero_quantiles':quantiles(a[:,15]),
                    'distance_quantiles':quantiles(a[:,10]/1e9),'known_hit_count':int(a[:,11].sum()),'samples':samples})
            rows.append({'groups':groups,'canonical_hashes':hashes,'baseline_identical':hashes==base,'diagnostic_sha256':sha(terminal)})
        gates={'frozen_paths_unchanged':all(r['baseline_identical'] for r in rows),'complete_repeat_exact':rows[0]==rows[1]}
        cases.append({'case':case,'legs':rows,'gates':gates})
        print(json.dumps({'case':case,'gates':gates,'groups':rows[0]['groups']}),flush=True)
    report={'cases':cases,'reason_map':{'0':'event_cap','1':'region_label','2':'ray_distance','3':'grid_exit','4':'negative_segment','5':'segment_limit'},
            'all_pass':all(all(c['gates'].values()) for c in cases)}
    (ROOT/'reports/photon_boundary_failure.json').write_text(json.dumps(report,indent=2)+'\n')
    return 0 if report['all_pass'] else 1


if __name__=='__main__':raise SystemExit(main())
