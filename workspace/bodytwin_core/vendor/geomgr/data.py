"""geomgr.data: build the frozen inputs (LOCAL ONLY: needs trimesh, P1 cache, VSD, Keast).

  python3 -m geomgr.data --out <dir>     -> femur_pop.npz, tibia_pop.npz, vsd_obs.npz, registry_femur_r.json,
                                           registry_tibia_r.json, data_manifest.json

Sources (read-only): results/P1/cache/{registered,landmarks,template}.npz (P1 registrations/anchors),
VSD ManualLandmarks + Bones via results/P1/p1_common (load_vsd_raters, load_vsd_femur),
Keast 2023 tibia SSM (ShapeModels/tibia/tibiaShapeModel.mat + Data/participant-characteristics.csv).
INTERN: attachments are TLEM 2.0-derived (license unclear) -> license='internal'.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE.parents[1]
P1 = RES / 'P1'
KEAST = Path('external_mount')

from . import core as C  # noqa: E402
from .identity import Entity, Registry  # noqa: E402


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _p1():
    if str(P1) not in sys.path:
        sys.path.insert(0, str(P1))
    import p1_common  # noqa: F401  (read-only import: functions only)
    return p1_common


# ------------------------------------------------------------------ femur
def femur_anchor_lists(F, tri, bc):
    return [list(map(int, F[t])) for t in tri], [list(map(float, b / b.sum())) for b in np.clip(bc, 0, None)]


def consensus_factory(TV, F, tpos):
    import trimesh
    tm = trimesh.Trimesh(TV, F, process=False)

    def consensus(idx):
        c, _, tri = trimesh.proximity.closest_point(tm, tpos[idx].mean(0))
        bc = trimesh.triangles.points_to_barycentric(TV[F[tri]], c)
        return tri, bc
    return consensus


def build_femur(out):
    import trimesh
    pc = _p1()
    R = np.load(P1 / 'cache/registered.npz')
    keys = [str(k) for k in R['keys']]
    Xp, F = R['Xp'], R['F'].astype(np.int64)
    LMK = np.load(P1 / 'cache/landmarks.npz')
    TV = np.load(P1 / 'cache/template.npz')['V']
    tri, bc = LMK['tri'], LMK['bc']
    idx_imp = [i for i, k in enumerate(keys) if k.startswith('imp:')]
    idx_vsd = [i for i, k in enumerate(keys) if k.startswith('vsd:')]
    idx_tlem = [i for i, k in enumerate(keys) if k.startswith('tlem:')]
    order = idx_imp + idx_vsd + idx_tlem
    S = Xp[order]
    names = [keys[i] for i in order]
    groups = [k[4:7] if k.startswith('imp:') else (k[4:] if k.startswith('vsd:') else 'TLEM') for k in names]
    dataset = ['imperial' if k.startswith('imp:') else ('vsd' if k.startswith('vsd:') else 'tlem') for k in names]

    # --- VSD raters (mirrored to right like P1), raw CT vertices
    subj, side, arr = pc.load_vsd_raters()          # (19,22,20,3)
    meanL = np.nanmean(arr, axis=2)
    nS, nL, nR = arr.shape[0], arr.shape[1], arr.shape[2]
    assert [f'vsd:{s}' for s in subj] == [names[i] for i in range(70, 89)], 'VSD order mismatch'

    def pts(Ld):
        d = {n: Ld[i] for i, n in enumerate(C.LM_NAMES)}
        d['HC6'] = C.sphere_fit(np.array([d[k] for k in C.HEAD6]))[0]
        return d
    obs_lm = np.full((nS, nR, 23, 3), np.nan)
    obs_feat = np.full((nS, nR, len(C.FEMUR_FEATURES)), np.nan)
    feat_mean20 = np.zeros((nS, len(C.FEMUR_FEATURES)))
    dev_anat = np.full((nS, nR, 23, 3), np.nan)
    frames = np.zeros((nS, 3, 3))
    for s, sid in enumerate(subj):
        Vraw = np.asarray(pc.load_vsd_femur(sid, side[s]).vertices)
        L0 = pts(meanL[s])
        _, Ax = C.knee_frame(L0)
        frames[s] = Ax
        feat_mean20[s] = C.femur_features_from_landmarks(L0, Vraw)
        P0 = np.array([L0[n] for n in C.PT_NAMES])
        for r in range(nR):
            if np.isnan(arr[s, :, r]).any():
                continue
            Lr = pts(arr[s, :, r])
            Pr = np.array([Lr[n] for n in C.PT_NAMES])
            obs_lm[s, r] = Pr
            obs_feat[s, r] = C.femur_features_from_landmarks(Lr, Vraw)
            dev_anat[s, r] = (Pr - P0) @ Ax.T
    # --- template positions of each VSD subject's mean landmarks (P1 p1_landmarks recipe)
    tpos = np.zeros((nS, nL, 3))
    for s, sid in enumerate(subj):
        V = S[names.index(f'vsd:{sid}')]
        m = trimesh.Trimesh(V, F, process=False)
        c, _, tr = trimesh.proximity.closest_point(m, meanL[s])
        b = trimesh.triangles.points_to_barycentric(V[F[tr]], c)
        tpos[s] = pc.bary_points(TV, F, tr, b)
    consensus = consensus_factory(TV, F, tpos)
    tri_all, bc_all = consensus(list(range(nS)))
    # sanity: consensus of all == P1 anchors
    p_ours = pc.bary_points(TV, F, tri_all, bc_all)
    p_p1 = pc.bary_points(TV, F, tri, bc)
    cons_vs_p1 = float(np.abs(p_ours - p_p1).max())
    loo_tri = np.zeros((nS, nL), int)
    loo_bc = np.zeros((nS, nL, 3))
    ndef = np.full((nS, nS, 23, 3), np.nan)       # [s, q]: residual of q with anchors built without s and q (anat)
    for s in range(nS):
        t_, b_ = consensus([i for i in range(nS) if i != s])
        loo_tri[s], loo_bc[s] = t_, b_
        for q in range(nS):
            if q == s:
                continue
            t2, b2 = consensus([i for i in range(nS) if i not in (s, q)])
            Vq = S[names.index(f'vsd:{subj[q]}')]
            A = pts(pc.bary_points(Vq, F, t2, b2))
            L0 = pts(meanL[q])
            res = np.array([L0[n] - A[n] for n in C.PT_NAMES])
            ndef[s, q] = res @ frames[q].T
    # all-19 definition noise (for the deployed service): consensus without q only
    ndef_all = np.zeros((nS, 23, 3))
    for q in range(nS):
        t2, b2 = consensus([i for i in range(nS) if i != q])
        Vq = S[names.index(f'vsd:{subj[q]}')]
        A = pts(pc.bary_points(Vq, F, t2, b2))
        L0 = pts(meanL[q])
        ndef_all[q] = np.array([L0[n] - A[n] for n in C.PT_NAMES]) @ frames[q].T
    att_idx, att_w = femur_anchor_lists(F, LMK['att_tri'], LMK['att_bc'])
    lm_idx, lm_w = femur_anchor_lists(F, tri, bc)
    vsd_meta = []
    for sid in subj:
        import scipy.io as sio
        Mm = sio.loadmat(pc.VSD / 'Bones' / f'{sid}.mat', simplify_cells=True)['M']
        vsd_meta.append([float(Mm['height']), float(Mm['weight']), float(Mm['age'])])
    np.savez_compressed(
        out / 'femur_pop.npz', S=S, F=F, names=np.array(names), groups=np.array(groups), dataset=np.array(dataset),
        template_V=TV, lm_tri=tri, lm_bc=bc, loo_tri=loo_tri, loo_bc=loo_bc,
        att_names=LMK['att_names'], att_muscle=LMK['att_muscle'], att_tri=LMK['att_tri'], att_bc=LMK['att_bc'])
    np.savez_compressed(
        out / 'vsd_obs.npz', subj=np.array(subj), side=np.array(side), obs_lm=obs_lm, obs_feat=obs_feat,
        feat_mean20=feat_mean20, dev_anat=dev_anat, frames=frames, ndef=ndef, ndef_all=ndef_all,
        meta_h_w_age=np.array(vsd_meta), feat_names=np.array(C.FEMUR_FEATURES), pt_names=np.array(C.PT_NAMES))
    reg = femur_registry(F, TV, lm_idx, lm_w, att_idx, att_w, [str(a) for a in LMK['att_names']])
    (out / 'registry_femur_r.json').write_text(json.dumps(reg.to_json(), indent=0))
    return dict(consensus_all_vs_P1_anchor_max_mm=cons_vs_p1, n_imperial=len(idx_imp), n_vsd=len(idx_vsd),
                n_tlem=len(idx_tlem), n_rating_sets=int(np.isfinite(obs_lm[:, :, 0, 0]).sum()),
                registry_manifest=reg.manifest_hash(), n_entities=len(reg.entities))


def femur_registry(F, TV, lm_idx, lm_w, att_idx, att_w, att_names):
    part = 'femur_r'
    reg = Registry(part, len(TV), F, version_note='P1 template C34RFE decimated (5002 v); anchors = VSD 19-rater consensus (P1)')
    for n, ix, w in zip(C.LM_NAMES, lm_idx, lm_w):
        reg.add(Entity(part, 'landmark', n, indices=ix, weights=w, source='atlas',
                       rule='VSD ManualLandmarks consensus of 19 subjects projected on template (results/P1/p1_landmarks.py)'))
    seen = {}
    for n, ix, w in zip(att_names, att_idx, att_w):
        base = n.replace(' ', '_')
        seen[base] = seen.get(base, 0) + 1
        nm = base if seen[base] == 1 else f'{base}#{seen[base]}'
        reg.add(Entity(part, 'attachment', nm, indices=ix, weights=w, source='atlas_internal', license='internal',
                       rule='TLEM 2.0 Table A3 point projected on template (results/P1/p1_landmarks.py)'))
    lmb = lambda n: f'{part}/landmark/{n}'  # noqa: E731
    reg.add(Entity(part, 'derived', 'HJC', rule='sphere_center', deps=[lmb(k) for k in C.HEAD6], source='derived',
                   note='HC6: least-squares sphere through the 6 head landmarks'))
    reg.add(Entity(part, 'derived', 'KJC', rule='midpoint', deps=[lmb('MEC'), lmb('LEC')], source='derived'))
    # regions on the template (derived rules, priority order makes them disjoint)
    L = {n: np.asarray(w) @ TV[np.asarray(ix)] for n, ix, w in zip(C.LM_NAMES, lm_idx, lm_w)}
    hc, hr = C.sphere_fit(np.array([L[k] for k in C.HEAD6]))
    kc = 0.5 * (L['MEC'] + L['LEC'])
    mech = C.unit(hc - kc)
    Lm = np.linalg.norm(hc - kc)
    h = (TV - kc) @ mech / Lm
    xml = (TV - kc) @ C.unit(L['LEC'] - L['MEC'])
    neck_c = np.mean([L[k] for k in ('SNI', 'ANI', 'INI', 'PNI')], 0)
    taken = np.zeros(len(TV), bool)
    rules = [('head', np.linalg.norm(TV - hc, axis=1) <= hr + 2.0, f'|v - HJC| <= r_head + 2 mm (r={hr:.2f})'),
             ('greater_trochanter', np.linalg.norm(TV - L['SGT'], axis=1) <= 20.0, '|v - SGT| <= 20 mm'),
             ('lesser_trochanter', np.linalg.norm(TV - L['LT'], axis=1) <= 12.0, '|v - LT| <= 12 mm'),
             ('neck', np.linalg.norm(TV - neck_c, axis=1) <= 20.0, '|v - neck centre| <= 20 mm'),
             ('medial_condyle', (h <= 0.08) & (xml < 0), 'height <= 0.08 L_mech above KJC, medial of KJC'),
             ('lateral_condyle', (h <= 0.08) & (xml >= 0), 'height <= 0.08 L_mech above KJC, lateral of KJC'),
             ('shaft', (h >= 0.30) & (h <= 0.70), '0.30 <= height/L_mech <= 0.70')]
    for nm, mask, rule in rules:
        m = mask & ~taken
        taken |= m
        reg.add(Entity(part, 'region', nm, indices=np.flatnonzero(m).tolist(), rule=rule, source='derived',
                       note='vertex set on template; synthetic region definition (no anatomical facit)'))
    reg.add(Entity(part, 'frame', 'ISB', rule='isb_femur_right',
                   deps=[f'{part}/derived/HJC', f'{part}/derived/KJC', lmb('MEC'), lmb('LEC')], source='derived',
                   note='Wu et al. 2002 femur frame from HJC, KJC, epicondyles'))
    reg.add(Entity(part, 'frame', 'knee', rule='knee_p1', deps=[lmb('MEC'), lmb('LEC'), lmb('PMC'), lmb('PLC')],
                   source='derived', note='results/P1/p1_common.knee_frame'))
    return reg


# ------------------------------------------------------------------ tibia (Keast 2023)
def build_tibia(out):
    import scipy.io as sio
    import trimesh
    from scipy.spatial import cKDTree
    d = sio.loadmat(KEAST / 'ShapeModels/tibia/tibiaShapeModel.mat', simplify_cells=True)['tibiaShapeModel']
    Nd = d['nodes'].reshape(30, 3500, 3)
    F = (d['F'] - 1).astype(np.int64)
    regdir = KEAST / 'ShapeModels/tibia/registered'
    files = sorted(os.listdir(regdir))
    trees = [cKDTree(Nd[j]) for j in range(30)]
    case = [None] * 30
    match = []
    for fn in files:
        m = trimesh.load(regdir / fn, process=False)
        dd = [trees[j].query(m.vertices)[0].max() for j in range(30)]
        j = int(np.argmin(dd))
        case[j] = fn.split('-')[0]
        match.append(float(dd[j]))
    assert None not in case and len(set(case)) == 30
    cov = {}
    with open(KEAST / 'Data/participant-characteristics.csv') as f:
        for row in csv.DictReader(f):
            cov[row['deidentified_record_number']] = row
    H = np.array([float(cov[c]['Height']) for c in case])
    W = np.array([float(cov[c]['Weight']) for c in case])
    A = np.array([float(cov[c]['Age']) if cov[c]['Age'] else np.nan for c in case])
    mesh = trimesh.Trimesh(d['meanPoints'], F, process=False)
    reg, rules = tibia_registry(d['meanPoints'], F)
    np.savez_compressed(out / 'tibia_pop.npz', S=Nd, F=F, case=np.array(case), height_cm=H, weight_kg=W, age_yr=A,
                        mean_points=d['meanPoints'])
    (out / 'registry_tibia_r.json').write_text(json.dumps(reg.to_json(), indent=0))
    return dict(n=30, row_match_max_mm=float(max(match)), mean_watertight=bool(mesh.is_watertight),
                mean_volume_mm3=float(mesh.volume), age_missing=int(np.isnan(A).sum()), rules=rules,
                registry_manifest=reg.manifest_hash(), n_entities=len(reg.entities))


def tibia_registry(V, F):
    part = 'tibia_r'
    V = np.asarray(V, float)
    c = V.mean(0)
    _, _, vt = np.linalg.svd(V - c)
    ax = vt[0]
    h = (V - c) @ ax
    # proximal end = wider cross-section (plateau)
    top = V[h > np.percentile(h, 95)]
    bot = V[h < np.percentile(h, 5)]
    def spread(P):
        Q = P - P.mean(0)
        Q = Q - np.outer(Q @ ax, ax)
        return float(np.sqrt((Q ** 2).sum(1).mean()))
    if spread(top) < spread(bot):
        ax, h = -ax, -h
    i_prox = int(np.argmax(h))
    i_mm = int(np.argmin(h))
    distal = V[h < np.percentile(h, 6)]
    ml = V[i_mm] - distal.mean(0)
    ml = ml - (ml @ ax) * ax
    ml = ml / np.linalg.norm(ml)                        # medial direction
    prox = h > np.percentile(h, 90)
    pc_ = V[prox].mean(0)
    s_ml = (V - pc_) @ ml
    med = prox & (s_ml > 8.0)
    lat = prox & (s_ml < -8.0)
    i_mtp = int(np.flatnonzero(med)[np.argmax(h[med])])
    i_ltp = int(np.flatnonzero(lat)[np.argmax(h[lat])])
    reg = Registry(part, len(V), F, version_note='Keast 2023 tibia SSM topology (3500 points)')
    for nm, i, rule in (('TPROX', i_prox, 'most proximal vertex along the mean long axis'),
                        ('MM', i_mm, 'most distal vertex (medial malleolus tip)'),
                        ('MTP', i_mtp, 'most proximal vertex of the medial plateau (>8 mm medial of proximal centroid)'),
                        ('LTP', i_ltp, 'most proximal vertex of the lateral plateau (>8 mm lateral)')):
        reg.add(Entity(part, 'landmark', nm, indices=[i], weights=[1.0], source='derived', rule=rule,
                       note='geometric rule on the Keast mean shape; vertex anchor'))
    L = (h.max() - h.min())
    hh = (h - h.min()) / L
    taken = np.zeros(len(V), bool)
    for nm, mask, rule in (('plateau', hh >= 0.88, 'top 12 % of length'), ('distal', hh <= 0.10, 'bottom 10 % of length'),
                           ('shaft', (hh >= 0.30) & (hh <= 0.70), '30-70 % of length')):
        m = mask & ~taken
        taken |= m
        reg.add(Entity(part, 'region', nm, indices=np.flatnonzero(m).tolist(), rule=rule, source='derived'))
    lb = lambda n: f'{part}/landmark/{n}'  # noqa: E731
    reg.add(Entity(part, 'derived', 'PLAT_C', rule='midpoint', deps=[lb('MTP'), lb('LTP')], source='derived'))
    reg.add(Entity(part, 'frame', 'tibia', rule='tibia_right', deps=[lb('MTP'), lb('LTP'), lb('MM')], source='derived',
                   note='origin plateau midpoint, y to proximal from MM, z lateral'))
    return reg, dict(i_prox=i_prox, i_mm=i_mm, i_mtp=i_mtp, i_ltp=i_ltp, length_mm=float(L))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    man = dict(femur=build_femur(out), tibia=build_tibia(out))
    man['files'] = {p.name: sha(p) for p in sorted(out.glob('*')) if p.suffix in ('.npz', '.json') and p.name != 'data_manifest.json'}
    man['sources'] = {'P1/cache/registered.npz': sha(P1 / 'cache/registered.npz'),
                      'P1/cache/landmarks.npz': sha(P1 / 'cache/landmarks.npz'),
                      'P1/cache/template.npz': sha(P1 / 'cache/template.npz'),
                      'keast tibiaShapeModel.mat': sha(KEAST / 'ShapeModels/tibia/tibiaShapeModel.mat')}
    (out / 'data_manifest.json').write_text(json.dumps(man, indent=1))
    print(json.dumps({k: v for k, v in man.items() if k != 'files'}, indent=1, default=str))


if __name__ == '__main__':
    main()
