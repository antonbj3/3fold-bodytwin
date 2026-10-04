"""Invoked in a namespace with selected public inputs, model, runtime and output only."""
from dental_release.paths import expand as _release_expand
import sys, json, time, resource, hashlib
from pathlib import Path
import numpy as np, torch
sys.path.insert(0, '/model')
from toothcraft_plugin import Plugin
start = time.perf_counter()
p = Plugin('/model', '/weights')
out = Path('/output')
rows = []
for (n, src) in enumerate(sorted(Path('/inputs').glob('*.npz'))):
    with np.load(src, allow_pickle=False) as a:
        t = dict(a)
    tic = time.perf_counter()
    d = p.generate({'tsdf': t['tsdf'], 'seed': 123})
    dest = out / (src.stem + '.npz')
    np.savez_compressed(dest, tsdf=d['tsdf'])
    rows.append(dict(key=src.stem, status=d['status'], seconds=time.perf_counter() - tic, sha256=hashlib.sha256(dest.read_bytes()).hexdigest(), finite=bool(np.isfinite(d['tsdf']).all()), min=float(d['tsdf'].min()), max=float(d['tsdf'].max()), seed=123))
    (out / 'RECEIPT.json').write_text(json.dumps(dict(rows=rows, seconds=time.perf_counter() - start, torch=torch.__version__, cuda=torch.version.cuda, peak_vram_bytes=torch.cuda.max_memory_allocated(), max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, weights_loaded_strict=True, actual_diffusion_transitions=p.diff.num_timesteps, private_reference_unmounted=not Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/GENCAD_V4')).exists()), indent=2) + '\n')
    print(src.stem, 'saved', flush=True)
