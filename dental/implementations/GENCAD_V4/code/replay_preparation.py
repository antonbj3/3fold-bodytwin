"""Recompute X1B preparations from the local source-triangle snapshot."""
from common import *
sys.path.insert(0, str(DATA / 'participant_code'))
sys.path.insert(0, str(DATA / 'participant_code/legacy/vendor'))
sys.path.insert(0, str(ROOT / 'code/vendor'))
from legacy.geometry import height
from scipy import ndimage
from score_panel import verified_npz
import crown_fit_geometry as G
import crown_design_geometry as CDG

def run():
    from integrity import verify
    verify()
    start = time.perf_counter()
    wf = read(ROOT / 'FROZEN_WHOLE_INPUTS.json')
    q = read(ROOT / 'PREREG_R3.json')['parameters']
    rows = []
    scratch = DATA / 'replay_scratch'
    scratch.mkdir(exist_ok=True)
    CDG.SCRATCH = str(scratch)
    for r in wf['records']:
        if r['status'] != 'PREPARED':
            continue
        key = r['key']
        a = verified_npz(DATA / 'whole_private' / key / 'reference.npz', wf['private_files'][key + '/reference.npz'])
        p = verified_npz(DATA / 'whole_inputs' / key / 'preparation.npz', wf['public_files'][key + '/preparation.npz'])
        h = float(p['step'])
        origin = p['origin']
        shape = p['preparation'].shape
        g = G.Grid(origin, origin + h * (np.array(shape) - 1), h)
        (xs, ys, zs) = g.axes()
        (xx, yy) = np.meshgrid(xs, ys, indexing='ij')
        roof = height(a['source_triangles'], np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)
        base = float(a['source_triangles'][:, :, 2].min() - q['basal_extension_mm'])
        occ = np.isfinite(roof)[:, :, None] & (zs[None, None, :] >= base) & (zs[None, None, :] <= roof[:, :, None])
        phi = (ndimage.distance_transform_edt(~occ, sampling=h) - ndimage.distance_transform_edt(occ, sampling=h)).astype(np.float32)
        (prep, _) = CDG.preparation(dict(grid=g, phi_T=phi, z_m=float(p['margin_z'])), q['occlusal_reduction_mm'], q['axial_reduction_mm'], q['spacer_mm'], cache=False)
        cavity = G.design_cavity(prep, g, np.full(prep.shape, q['spacer_mm']), levels=(0.0,))
        pe = float(np.max(abs(prep - p['preparation'])))
        ce = float(np.max(abs(cavity - p['cavity'])))
        source = float(np.max(abs(phi - a['phi_source'])))
        rows.append(dict(key=key, source_field_error_mm=source, preparation_error_mm=pe, cavity_error_mm=ce, pass_exact=source == pe == ce == 0.0))
        for f in scratch.glob('*.npz'):
            f.unlink()
    scratch.rmdir()
    out = dict(rows=rows, all_exact=all((r['pass_exact'] for r in rows)), seconds=time.perf_counter() - start, scope='Independent replay from retained local source triangles, same source-derived virtual envelope and X1B operations. No original archive dependency.')
    dump(ROOT / 'raw/PREPARATION_REPLAY.json', out)
    if not out['all_exact']:
        raise ValueError('Preparation replay differs')
    print('preparation replay', len(rows), 'exact')
if __name__ == '__main__':
    run()
