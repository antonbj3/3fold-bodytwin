"""Bounded one-time import of two matched same-tooth variants, without source edits."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import shutil, numpy as np
from common import R, read, dump, sha
src = Path(_release_expand('@DENTAL_WORK_ROOT@/LANE_X1_CROWN_LOOP/L005_lo_M1_k6'))
manifest = []
for (vid, d) in [('medium', 'M1'), ('thick', 'M2')]:
    for (name, target) in [('model.npz', f'inputs/geometry/{d}_model.npz'), ('crown.stl', f'exports/{d}/crown.stl')]:
        p = src / vid / name
        q = R / target
        q.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, q)
        manifest.append(dict(source=str(p), local=target, source_sha256=sha(p), sha256=sha(q), bytes=q.stat().st_size))
    a = np.load(src / vid / 'case_h0.12.npz')
    dump(f'inputs/geometry/{d}_grid.json', dict(origin=a['origin'].tolist(), h=float(a['h']), shape=a['shape'].tolist(), z_m=float(a['z_m']), margin_shift_mm=-0.2 if vid == 'medium' else 0.2))
dump('MATCHED_SOURCE_MANIFEST.json', manifest)
print('Matched additional MB', sum((m['bytes'] for m in manifest)) / 1000000.0)
