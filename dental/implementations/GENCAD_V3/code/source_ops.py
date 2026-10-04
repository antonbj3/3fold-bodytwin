from dental_release.paths import expand as _release_expand
import hashlib, zipfile
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from legacy.geometry import *
from legacy.vendor.obstacle import constraints
from util import read, sha
ZIPS = {d: Path(_release_expand('@DENTAL_DATA_ROOT@/geometry')) / d / f for (d, f) in [('Bits2Bites', 'Bits2Bites_v01.zip'), ('Bite2Text', 'Bite2Text.zip')]}
FDI = {'molar_crown': 36, 'premolar_crown': 34, 'anterior_crown': 31, 'bridge3': 35, 'implant_crown': 46, 'inlay': 36, 'onlay': 46, 'veneer': 31, 'lattice_onlay': 37}

def pair(row):
    meta = read(row['metadata'])
    out = {}
    sources = []
    if sha(row['metadata']) != row['metadata_sha256']:
        raise ValueError('label metadata drift')
    with zipfile.ZipFile(ZIPS[row['dataset']]) as z:
        for jaw in ['upper', 'lower']:
            suffix = '/' + row['case'] + '/' + jaw + '.stl' if row['dataset'] == 'Bits2Bites' else row['case'] + '/ios/ios_' + jaw + '.stl'
            matches = [s for s in z.namelist() if s.endswith(suffix)]
            if len(matches) != 1:
                raise ValueError('Missing/nonunique geometry member ' + suffix)
            member = matches[0]
            blob = z.read(member)
            info = meta['arches'][jaw]
            if hashlib.sha256(blob).hexdigest() != info['input_sha256']:
                raise ValueError('scan hash mismatch')
            (v, f) = parse_stl(blob)
            lp = Path(row['metadata']).parent / (jaw + '_labels.npz')
            if sha(lp) != info['labels']['sha256']:
                raise ValueError('label hash mismatch')
            labs = np.load(lp, allow_pickle=False)['labels']
            if len(labs) != len(v):
                raise ValueError('label indexing')
            fl = labs[f]
            owner = np.where(fl[:, 0] == fl[:, 1], fl[:, 0], np.where(fl[:, 0] == fl[:, 2], fl[:, 0], np.where(fl[:, 1] == fl[:, 2], fl[:, 1], 0)))
            out[jaw] = dict(v=v, f=f, labels=labs, owner=owner, meta=info)
            sources.append(dict(member=member, zip_path=str(ZIPS[row['dataset']]), bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(), label_sha256=sha(lp)))
    return (out, sources)

def make_site(pair, family, compute_obstacles=True):
    arch = pair['lower']
    fdi = FDI[family]
    pp = arch['v'][arch['labels'] == fdi]
    if len(pp) < 100:
        raise ValueError('missing predicted tooth ' + str(fdi))
    fr = arch['meta']['frame']
    occ = unit(fr['occlusal_unit'])
    right = unit(fr['right_unit'])
    x = unit(right - occ * (right @ occ))
    y = unit(np.cross(occ, x))
    R = np.stack([x, y, occ], axis=1)
    base = pp.mean(0)
    if family == 'veneer':
        lm = arch['meta']['landmarks']['teeth'].get(str(fdi), {})
        normal = unit(lm['buccal_unit'])
        x = unit(np.cross(occ, normal))
        y = unit(np.cross(normal, x))
        R = np.stack([x, y, normal], axis=1)
    elif family == 'bridge3':
        a = arch['v'][arch['labels'] == 34]
        b = arch['v'][arch['labels'] == 36]
        if min(len(a), len(b)) < 100:
            raise ValueError('bridge abutments missing')
        x = unit(b.mean(0) - a.mean(0))
        x = unit(x - occ * (x @ occ))
        y = unit(np.cross(occ, x))
        R = np.stack([x, y, occ], axis=1)
        pp = arch['v'][np.isin(arch['labels'], [34, 35, 36])]
        base = pp.mean(0)
    local = (pp - base) @ R
    band = local[local[:, 2] <= np.quantile(local[:, 2], 0.35)] if family != 'veneer' else local
    (lo, hi) = np.quantile(band[:, :2], [0.02, 0.98], axis=0)
    center = (hi + lo) / 2
    half = (hi - lo) / 2
    if np.min(half) < 0.7 or np.max(half) > 18:
        raise ValueError('invalid predicted footprint')
    if family == 'inlay':
        half *= 0.48
    if family == 'onlay':
        half *= 0.82
    if family == 'veneer':
        half *= 0.8
    base = base + R[:, :2] @ center
    (xy, faces, uv) = grid(half, family)
    target_labels = [34, 35, 36] if family == 'bridge3' else [fdi]
    lv = (arch['v'] - base) @ R
    tri = lv[arch['f']]
    target_tri = crop(tri[np.isin(arch['owner'], target_labels)], xy)
    native = height(target_tri, xy)
    if np.isfinite(native).mean() < 0.25:
        raise ValueError('too little native projected reference coverage')
    ni = ~np.isin(arch['labels'], target_labels + [0])
    neighbors = lv[ni]
    if len(neighbors) < 100:
        raise ValueError('missing adjacent geometry')
    dd = np.linalg.norm(neighbors[:, :2], axis=1)
    near = neighbors[dd < max(half) * 2 + 5]
    if len(near) < 30:
        near = neighbors[np.argsort(dd)[:200]]
    ztop = float(np.quantile(near[:, 2], 0.85))
    U = ((pair['upper']['v'] - base) @ R)[pair['upper']['f']]
    U = crop(U, xy)
    if family == 'veneer':
        U = np.empty((0, 3, 3))
    if len(U) and compute_obstacles:
        (A, b) = constraints(U, xy, faces, 0.0)
    else:
        A = csr_matrix((0, len(xy)))
        b = np.empty(0)
    if A.shape[0] and (not np.isfinite(b).all() or np.min(A.data) < -1e-07 or np.max(np.abs(np.asarray(A.sum(1)).ravel() - 1)) > 1e-07):
        raise ValueError('invalid continuous obstacle rows')
    ceiling = height(U, xy, True)
    return dict(xy=xy, faces=faces, uv=uv, A=A, b=b, ceiling=ceiling, reference=native, target_tri=target_tri, antagonist=U, weights=area_weights(xy, faces), neighbor_points=near[np.linspace(0, len(near) - 1, min(len(near), 1500), dtype=int)], ztop=ztop, half=half, base=base, R=R, source_fdi=fdi)

def wall(family, level):
    if family == 'veneer':
        return (0.6, 'IPS_e_max_CAD_veneer')
    if family == 'inlay':
        return (1.0, 'IPS_e_max_CAD_inlay')
    if family == 'onlay':
        return (1.5, 'IPS_e_max_CAD_partial')
    if family == 'lattice_onlay':
        return (0.5, 'engineering_lattice_no_clinical_IFU')
    if family == 'bridge3':
        return (0.5 if level in ['easy', 'normal'] else 1.0, 'KATANA_HTML_PLUS' if level in ['easy', 'normal'] else 'KATANA_STML')
    anterior = family == 'anterior_crown'
    return ((0.4 if anterior else 0.5) if level in ['easy', 'normal'] else 1.0, 'KATANA_HTML_PLUS' if level in ['easy', 'normal'] else 'IPS_e_max_CAD_adhesive_1mm')
