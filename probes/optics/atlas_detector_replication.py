"""Independent-seed evaluation of frozen detector centers and selected IDs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_detector_information import matrices, score
from atlas_oxygen_bounds import observe


def run(original, captures):
    with np.load(original, allow_pickle=False) as a:
        centers, _ = candidates(a['exits'], a['terminal'], a['source'])
    selected = (5, 7, 9); reference = (0, 4, 8)
    rows = []
    for filename in captures:
        with np.load(filename, allow_pickle=False) as a:
            out = {k:a[k] for k in a.files if k != 'source'}
        paths = out['paths']; terminal = out['terminal']; exits = out['exits']
        masks = (np.linalg.norm(exits[:,None,:3] - centers[None,:,:], axis=2) <= 5) & (terminal[:,1] > 0)[:,None]
        jacobian = matrices(paths, masks)
        chosen = score(jacobian, selected); baseline = score(jacobian, reference)
        tail = observe(paths, masks[:,selected[0]], terminal[:,2] > 0, .8)[0]['unfinished_upper']
        rows.append(dict(hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()},
                         jacobian_hash=hashlib.sha256(jacobian.tobytes()).hexdigest(),
                         selected_ids=list(selected), reference_ids=list(reference),
                         selected_counts=masks[:,selected].sum(axis=0).tolist(),
                         reference_counts=masks[:,reference].sum(axis=0).tolist(),
                         selected_sigma_min=chosen, reference_sigma_min=baseline, gain=chosen/baseline,
                         absorbed=int(terminal[:,0].sum()),escaped=int(terminal[:,1].sum()),
                         residual=int(terminal[:,2].sum()),leaks=int(terminal[:,3].sum()),caps=int(out['counters'][2]),
                         photons=len(paths),unfinished_intensity_upper_per_detector=tail))
    gates=dict(full_repeat=rows[0]==rows[1],
               energy_and_leaks=all(r['absorbed']+r['escaped']+r['residual']==r['photons']*2**30 and r['leaks']==0 for r in rows),
               support=all(min(r['selected_counts']+r['reference_counts'])>=20 for r in rows),
               positive_rank=all(np.isfinite([r['selected_sigma_min'],r['reference_sigma_min']]).all() and min(r['selected_sigma_min'],r['reference_sigma_min'])>0 for r in rows),
               independent_gain=all(r['gain']>1 for r in rows))
    return dict(rows=rows,gates=gates,center_hash=hashlib.sha256(centers.tobytes()).hexdigest(),
                scope='New-seed conditional information at previously frozen centers and IDs; capped proposal, geometry and clinical acceptance remain separate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',required=True)
    p.add_argument('captures',nargs=2)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();r=run(a.original,a.captures)
    a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
    raise SystemExit(0 if all(r['gates'].values()) else 2)
