"""Measured arch reader and rule-based virtual preparation. No physical preparation claimed."""
from common import *
sys.path.insert(0, str(DATA / 'participant_code'))
sys.path.insert(0, str(DATA / 'participant_code/legacy/vendor'))
sys.path.insert(0, str(ROOT / 'code/vendor'))
from source_ops import pair, make_site
from legacy.geometry import height, crop
import crown_fit_geometry as G
import crown_design_geometry as CDG
from scipy import ndimage
import resource

def run():
    begin = time.perf_counter()
    pr = read(ROOT / 'PREREG_R3.json')
    q = pr['parameters']
    records = []
    for row in pr['selection']:
        (arches, sources) = pair(row)
        for fam in pr['families']:
            key = row['case_key'] + '_' + fam
            st = time.perf_counter()
            rec = dict(key=key, case_key=row['case_key'], family=fam, split=row['split'], dataset=row['dataset'], status='FAILED')
            try:
                site = make_site(arches, fam, compute_obstacles=False)
                low = arches['lower']
                local = (low['v'] - site['base']) @ site['R']
                tri = local[low['f'][low['owner'] == site['source_fdi']]]
                v = tri.reshape(-1, 3)
                margin = float(np.quantile(v[:, 2], q['margin_quantile']))
                base = float(v[:, 2].min() - q['basal_extension_mm'])
                h = q['grid_mm']
                lo = v.min(0) - q['padding_mm']
                lo[2] = base - q['padding_mm']
                hi = v.max(0) + q['padding_mm']
                grid = G.Grid(lo, hi, h)
                if np.prod(grid.shape) > q['max_grid_nodes']:
                    raise ValueError('frozen grid-size limit')
                (xs, ys, zs) = grid.axes()
                (xx, yy) = np.meshgrid(xs, ys, indexing='ij')
                xy = np.c_[xx.ravel(), yy.ravel()]
                roof = height(tri, xy).reshape(xx.shape)
                occupied = np.isfinite(roof)[:, :, None] & (zs[None, None, :] >= base) & (zs[None, None, :] <= roof[:, :, None])
                if not occupied.any():
                    raise ValueError('empty source-supported insertion envelope')
                phi = (ndimage.distance_transform_edt(~occupied, sampling=h) - ndimage.distance_transform_edt(occupied, sampling=h)).astype(np.float32)
                scratch = DATA / 'whole_work' / key
                scratch.mkdir(parents=True, exist_ok=True)
                CDG.SCRATCH = str(scratch)
                (prep, _) = CDG.preparation(dict(grid=grid, phi_T=phi, z_m=margin), q['occlusal_reduction_mm'], q['axial_reduction_mm'], q['spacer_mm'], cache=False)
                cavity = G.design_cavity(prep, grid, np.full(prep.shape, q['spacer_mm']), levels=(0.0,))
                pub = DATA / 'whole_inputs' / key
                pub.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(pub / 'preparation.npz', preparation=prep, cavity=cavity, origin=grid.origin, step=h, margin_z=margin, source_base=site['base'], source_R=site['R'])
                z = zs[None, None, :]
                srcocc = np.clip(0.5 - phi / h, 0, 1) * (z >= margin)
                prepocc = np.clip(0.5 - prep / h, 0, 1) * (z >= margin)
                rec.update(status='PREPARED', virtual_source_volume_mm3=float(srcocc.sum() * h ** 3), virtual_removed_volume_mm3=float(np.maximum(srcocc - prepocc, 0).sum() * h ** 3), grid_shape=grid.shape, margin_z_mm=margin, source_fdi=site['source_fdi'])
                private = DATA / 'whole_private' / key
                private.mkdir(parents=True, exist_ok=True)
                adj = {}
                fdi = site['source_fdi']
                for (label, number) in [('mesial', fdi - 1 if fdi % 10 > 1 else 41), ('distal', fdi + 1)]:
                    adj[label] = local[low['f'][low['owner'] == number]]
                up = (arches['upper']['v'] - site['base']) @ site['R']
                op = crop(up[arches['upper']['f']], v[:, :2], pad=1.0)
                np.savez_compressed(private / 'reference.npz', source_triangles=tri, antagonist_triangles=op, mesial_triangles=adj['mesial'], distal_triangles=adj['distal'], phi_source=phi)
                dump(private / 'PROVENANCE.json', dict(sources=sources, source_role='Entire X11 target surface retained as private source; not a measured preparation', source_fdi=site['source_fdi'], reference_uses_source='All source triangles above the preregistered virtual margin; no virtual closure is called observed'))
                for p in scratch.glob('*.npz'):
                    p.unlink()
                scratch.rmdir()
            except (ValueError, RuntimeError, KeyError) as e:
                rec['reason'] = str(e)
            rec['seconds'] = time.perf_counter() - st
            records.append(rec)
            dump(ROOT / 'raw/WHOLE_PREPARATION_CHECKPOINT.json', records)
            print('whole prep', key, rec['status'], round(rec['seconds'], 1), rec.get('reason', ''), flush=True)
    dump(DATA / 'whole_inputs/RECORDS.json', records)
    freeze(ROOT / 'FROZEN_WHOLE_INPUTS.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R3.json'), public_files=inventory(DATA / 'whole_inputs'), private_files=inventory(DATA / 'whole_private'), records=records, seconds=time.perf_counter() - begin, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, geometry_scope='Virtual finite insertion envelope and X1B preparation; not measured prepared anatomy'))
    state('WHOLE_INPUTS_FROZEN', 'Actual X1B preparation operation complete', 'Generate complete virtual shells without mounting original target surfaces')
if __name__ == '__main__':
    run()
