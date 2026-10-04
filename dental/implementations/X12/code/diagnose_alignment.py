from dental_release.paths import expand as _release_expand
import json, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from archive_io import nifti, read_member
from pair_atlas import transformations, transform, ROOT, FDI
sys.path.insert(0, _release_expand('@DENTAL_INPUT_ROOT@/workspace/cells/geometry'))
from tf2_io import ZIP, ROOT as TFROOT, read_mha_bytes
start = time.time()
lookup = {r['name']: r for r in json.load(open(ROOT / 'raw/local_records.json'))}
out = []
with zipfile.ZipFile(ZIP) as z:
    for pid in ['P1', 'P2', 'P4', 'P5', 'P9']:
        case = 'ToothFairy2P_' + pid[1:].zfill(3)
        (t, sp, h) = read_mha_bytes(z.read(TFROOT + '/labelsTr/' + case + '.mha'))
        (nimg, _) = nifti(lookup['Pulpy3D/' + pid + '/gt_ian.nii.gz'])
        n = np.asanyarray(nimg.dataobj) > 0
        tn = (t == 3) | (t == 4)
        (pp, _) = nifti(lookup['Pulpy3D/' + pid + '/gt_instance.nii.gz'])
        p = np.asanyarray(pp.dataobj)
        (semimg, _) = nifti(lookup['Pulpy3D/' + pid + '/gt_pulp_mandible.nii.gz'])
        sem = np.asanyarray(semimg.dataobj) > 0
        nerve = []
        semantic = []
        for (perm, flip) in transformations(n.shape, t.shape):
            nn = transform(n, perm, flip)
            dice = 2 * np.sum(nn & tn) / max(1, int(nn.sum()) + int(tn.sum()))
            nerve.append({'dice': float(dice), 'perm': perm, 'flip': flip})
        for (perm, flip) in transformations(p.shape, sem.shape):
            pp = transform(p > 0, perm, flip)
            if np.array_equal(pp, sem):
                semantic.append({'perm': perm, 'flip': flip})
        nerve.sort(key=lambda r: r['dice'], reverse=True)
        best = nerve[0]
        mapped = transform(p, best['perm'], best['flip'])
        coords = np.argwhere(mapped > 0)
        values = mapped[tuple(coords.T)]
        containment = float(np.mean(t[tuple(coords.T)] == values))
        out.append({'case': pid, 'nerve_best': nerve[:2], 'instance_to_semantic_exact': semantic, 'same_IAN_rule_pulp_containment': containment})
        print(json.dumps(out[-1]), flush=True)
        (ROOT / 'raw/R3_alignment_diagnosis.json').write_text(json.dumps(out, indent=2) + '\n')
print('wall_s', time.time() - start)
