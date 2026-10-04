"crown_design_geometry.py — DENT-DESIGN-CROWN geometri: verklig molar -> SDF -> syntetisk preparation -> krona/die.\n\nRe-uses K2 : s Geometric Core (cells/manufacturing/crown_fit_geometry.py , dental workspace 2026 -09 -23 :\ndental frame, exact redistans A143, make_preparation, morphological opening/ closure). No copied; import via path.\n\nFall (dokumenterat, syntetiskt): OpenMandible basmodell (git e1f8cef, DOI 10.1016/j.dental.2021.01.009),\nTooth_L6 = FDI 36 (lower left first molar) with real neighbors L5 ( FDI 35 ) and L7 ( FDI 37 ) in the same frame.\nAntagonist: OpenMandible lacks an upper jaw -> antagonist represented by a sphere (radie R_ant ) inserted in\nthe crown's Central Fossa (ANTAGANDE, see REPORT §Fall). Data license for STL UNKNOWN -> is processed local , is not disseminated.\n\nDesign variables bearing geometry (x):\n  t_occ minimum occlusive crown wall thickness [mm ] (preparation = normal-erosion of the anatomy with t_occ + s )\n  t_ax axial reduction = chamfer depth = axial wall thickness at margin [mm ]\n  r_f occlusal fissure rounding [mm ]: morphological closure of the anatomic outer shape (fills fissures,\n          Adds only material where a bullet with radius r_f does not reach -> can never penetrate an antagonist-\n          sphere with radius >= r_f ; see g_occl )\n  s cement gap (set-up, internal) [mm ); preparation reduced t + s , intaglio = prep dilated s\n  r_tool  ball end mill radius at the sintered scale [mm]: outer surface concavities < r_tool cannot be milled -> closure\nFixed (hot/preparation parameters): convergence 3° per wall (6° total), line angle radius 0,5 mm ( min ( 0,5 , t_occ ));\nchamfer holes 0,9 · t_ax , margin plane z_m = lowest point of enamel + 1,0 mm (as K2 ).\n\nEnheter: mm. SDF negative inside. Large fields -> SCRATCH ( INTENSO ).\n"
from dental_release.paths import expand as _release_expand
import os
import sys
import json
import time
import hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'manufacturing'))
import crown_fit_geometry as G
SCRATCH = os.environ.get('CROWN_SCRATCH', _release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/datasets/dental_3fold_extra/scratch_CROWN'))
OM_STL = G.OM_STL
Z_BOT_BELOW_MARGIN = 3.0

def arch_frame(R, c):
    "Buckal and mesial direction for L6 in the dental frame from the entire dental centroids of the arc ( OpenMandible - frame ).\n    Buckal = away from the centre of the arc (projected on the tooth xy plane); mesial = against L5 ."
    import trimesh
    cents = {}
    for side in 'LR':
        for i in range(1, 8):
            p = os.path.join(OM_STL, f'Tooth_{side}{i}_Enamel.stl')
            if os.path.exists(p):
                cents[f'{side}{i}'] = trimesh.load(p).vertices.mean(0)
    P = np.array(list(cents.values()))
    arch_c = P.mean(0)
    to_tf = lambda v: R @ v
    b = to_tf(cents['L6'] - arch_c)
    b[2] = 0
    b /= np.linalg.norm(b)
    m = to_tf(cents['L5'] - cents['L6'])
    m[2] = 0
    m -= m @ b * b
    m /= np.linalg.norm(m)
    return {'buccal': b, 'mesial': m, 'centroids_model': {k: v.tolist() for (k, v) in cents.items()}}

def neighbours_in_frame(R, c):
    """Grannarnas emalj (L5, L7) som trianglar i tandramen -> proximal kontroll."""
    import trimesh
    out = {}
    for t in ('L5', 'L7'):
        m = trimesh.load(os.path.join(OM_STL, f'Tooth_{t}_Enamel.stl'))
        out[t] = ((m.vertices - c) @ R.T, m.faces)
    return out

def build_case(h, cache=True):
    "The SDF of the tooth frame on grid with site for die down to z_m - 3 mm . Cache i SCRATCH."
    os.makedirs(SCRATCH, exist_ok=True)
    fn = os.path.join(SCRATCH, f'case_h{h}.npz')
    if cache and os.path.exists(fn):
        d = np.load(fn)
        grid = G.Grid(d['origin'], d['origin'] + d['h'] * (d['shape'] - 1), float(d['h']))
        return dict(grid=grid, phi_T=d['phi_T'], z_m=float(d['z_m']), R=d['R'], c=d['c'], buccal=d['buccal'], mesial=d['mesial'])
    ms = G.load_tooth('L6')
    (R, c) = G.tooth_frame(ms)
    Ve = (ms['Enamel'].vertices - c) @ R.T
    z_m = float(Ve[:, 2].min() + 1.0)
    V_all = np.vstack([(ms[p].vertices - c) @ R.T for p in ('Enamel', 'Dentin')])
    top = V_all[V_all[:, 2] > z_m - Z_BOT_BELOW_MARGIN - 0.5]
    M = 2.0
    lo = [top[:, 0].min() - M, top[:, 1].min() - M, z_m - Z_BOT_BELOW_MARGIN - 0.5 + 0.5 * h]
    hi = [top[:, 0].max() + M, top[:, 1].max() + M, Ve[:, 2].max() + 1.5]
    grid = G.Grid(lo, hi, h)
    (phi_T, _) = G.tooth_sdf(ms, R, c, grid, levels=(0.0,))
    af = arch_frame(R, c)
    np.savez(fn, phi_T=phi_T, origin=grid.origin, h=grid.h, shape=np.array(grid.shape), z_m=z_m, R=R, c=c, buccal=af['buccal'], mesial=af['mesial'])
    return dict(grid=grid, phi_T=phi_T, z_m=z_m, R=R, c=c, buccal=af['buccal'], mesial=af['mesial'])

def preparation(case, t_occ, t_ax, s=0.0, cache=True):
    "Synthetic full crown repair for crown with minimum wall thickness t_occ (occlusal, normal to the surface) and\n    Chamfer depth t_ax, with cement gap s included in the reduction (CAD: thickness = crown material, gap added inside).\n    Based on K2 make_preparation ( cells/manufacturing/crown_fit_geometry.py ) but the occlusive reduction is a\n    exactly normal-erosion of the anatomy with ( t_occ + s ) instead of an axial displacement (the displacement gave\n    normal thickness t*cos(slope): 0.64 mm minimum at 1.0 mm in the pilot). Axial: the contour at z_m displaced inwards ( t_ax + s );\n    3° per wall ; preparation = average of axial and occlusive terms (the most reduced wins), line angles\n    opens (radie min (0,5 , t_occ )), closes chamfer hole ( 0,9 · t_ax ). Under z_m untouched tooth."
    (grid, h) = (case['grid'], case['grid'].h)
    key = f'prepE_h{h}_o{t_occ:.3f}_a{t_ax:.3f}_s{s:.3f}.npz'
    fn = os.path.join(SCRATCH, key)
    if cache and os.path.exists(fn):
        d = np.load(fn)
        return (d['prep'], d['margin'])
    (phi_T, z_m) = (case['phi_T'], case['z_m'])
    (X, Y, Z) = grid.mesh()
    (ro, ra) = (t_occ + s, t_ax + s)
    (phi2, zk) = G.slice_sdf2d(phi_T, grid, z_m)
    lat = phi2[:, :, None] + ra + np.maximum(Z - z_m, 0.0) * np.tan(np.radians(3.0))
    del X, Y
    occ = G.redistance_lowmem(phi_T, grid, levels=(0.0, -ro)) + ro
    above = np.maximum(lat, occ)
    P1 = np.minimum(above, Z - z_m).astype(np.float32)
    lar = min(0.5, t_occ)
    P1 = G.redistance_lowmem(P1, grid, levels=(0.0, -lar))
    P1 = G.opening(P1, grid, lar)
    prep = np.where(Z > z_m + 1e-09, P1, phi_T).astype(np.float32)
    cr = 0.9 * t_ax
    prep = G.redistance_lowmem(prep, grid, levels=(0.0, cr))
    closed = G.closing(prep, grid, cr, nxt=(0.0, 0.03, 0.06, 0.1))
    band = Z < z_m + 2 * cr + 0.3
    prep = G.redistance_lowmem(np.where(band, closed, prep).astype(np.float32), grid, levels=(0.0, 0.03, 0.06, 0.1))
    g2 = G.Grid(grid.origin[:2], grid.origin[:2] + grid.h * (np.array(grid.shape[:2]) - 1), grid.h)
    (Pz, Nz) = G.zero_crossings(phi2, g2)
    margin = np.column_stack([Pz, np.full(len(Pz), zk)])
    np.savez(fn, prep=prep, margin=margin)
    return (prep, margin)

def outer_form(case, r_f=0.0, r_tool=0.0, cache=True):
    "The crown's outer shape: anatomical surface (OpenMandible L6 ) closed with r = max(r_f , r_tool ).\n    r_tool : the ball cutter cannot milling outer concaves with curvature radius < r (material remains) -> closure ."
    grid = case['grid']
    r = max(r_f, r_tool)
    if r <= 0:
        return case['phi_T']
    fn = os.path.join(SCRATCH, f'outer_h{grid.h}_r{r:.3f}.npz')
    if cache and os.path.exists(fn):
        return np.load(fn)['phi']
    phi = G.redistance_lowmem(case['phi_T'], grid, levels=(0.0, r))
    phi = G.closing(phi, grid, r)
    np.savez(fn, phi=phi)
    return phi

def crown_and_die(case, prep, phi_O, s):
    " Crown  = outer shape   s  ( A143 -redistans).\n    Die = intlio (above margin, the cement is counted as die-material)  tooth  under margin, kapad  z_m  -  3   mm ."
    (grid, z_m) = (case['grid'], case['z_m'])
    Z = grid.axes()[2][None, None, :]
    phi_I = G.dilate(G.redistance_lowmem(prep, grid, levels=(0.0, s)), grid, s) if s > 0 else prep
    crown = G.redistance_lowmem(np.maximum(np.maximum(phi_O, -phi_I), z_m - Z).astype(np.float32), grid)
    die_full = np.where(Z > z_m, phi_I, case['phi_T'])
    die = G.redistance_lowmem(np.maximum(die_full, z_m - Z_BOT_BELOW_MARGIN - Z).astype(np.float32), grid)
    return (crown, die, phi_I)

def volume(phi, grid, zmin=None):
    "Volume [ mm ^ 3 ] from SDF with subvoxel coating clip( 0,5 - phi/h, 0 , 1 )."
    occ = np.clip(0.5 - phi / grid.h, 0, 1)
    if zmin is not None:
        Z = grid.axes()[2][None, None, :]
        occ = occ * (Z >= zmin)
    return float(occ.sum() * grid.h ** 3)

def thickness_stats(crown, phi_I, phi_O, grid, z_m, t_ax=1.0):
    "Wall thickness = distance from the inlet surface points to the outer surface (- phi_O, exactly close to the surface).\n    Occlusal : Intakelio-normal n_z >= 0,5 and z > z_m + 2 · t_ax + 0,5 (above the chamfer). Axial: n_z < 0,5 and\n    z > z_m + 1,0 (the chamfer's thin outlet against the margin is excluded and is separately declared as ‘margin_zone')."
    (V, F) = G.surface_mesh(phi_I, grid)
    V = V[V[:, 2] > z_m + 0.05]
    V = V[grid.sample(phi_O, V, cval=1.0) < 0]
    t = -grid.sample(phi_O, V, cval=0.0)
    g = np.gradient(phi_I, grid.h)
    n = np.stack([grid.sample(gi, V, cval=0.0) for gi in g], 1)
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    occl = (n[:, 2] >= 0.5) & (V[:, 2] > z_m + 2 * t_ax + 0.5)
    ax = (n[:, 2] < 0.5) & (V[:, 2] > z_m + 1.0)
    mz = V[:, 2] <= z_m + 1.0
    q = lambda a, p: float(np.percentile(a, p)) if len(a) else None
    return {'occl_min': q(t[occl], 0.1), 'occl_p5': q(t[occl], 5), 'occl_median': q(t[occl], 50), 'axial_min': q(t[ax], 0.1), 'axial_p5': q(t[ax], 5), 'axial_median': q(t[ax], 50), 'margin_zone_p50': q(t[mz], 50), 'n_occl': int(occl.sum()), 'n_axial': int(ax.sum())}

def surface_for_fe(phi, grid, target_edge, out_stl, smooth_iter=0):
    "Zero level -> marching cubes -> pymeshlab isotropic ommeshning (target edge target_edge mm ) -> STL .\n    Returns (V, F) and quality measures. Control : waterproof, 2 -diversity."
    import pymeshlab
    import trimesh
    (V, F) = G.surface_mesh(phi, grid)
    ms = pymeshlab.MeshSet()
    ms.add_mesh(pymeshlab.Mesh(V, F))
    ms.meshing_remove_duplicate_vertices()
    ms.meshing_isotropic_explicit_remeshing(targetlen=pymeshlab.PureValue(target_edge), iterations=6, featuredeg=60.0, adaptive=False)
    m = ms.current_mesh()
    (V2, F2) = (m.vertex_matrix(), m.face_matrix())
    tm = trimesh.Trimesh(V2, F2, process=True)
    tm.fix_normals()
    info = {'n_faces': int(len(tm.faces)), 'watertight': bool(tm.is_watertight), 'volume_mm3': float(tm.volume), 'edge_mean_mm': float(tm.edges_unique_length.mean())}
    dv = grid.sample(phi, np.asarray(tm.vertices), cval=1.0)
    info['vertex_sdf_abs_p99_um'] = float(np.percentile(np.abs(dv), 99) * 1000)
    tm.export(out_stl)
    return (np.asarray(tm.vertices), np.asarray(tm.faces), info)

def proximal_and_occlusal_checks(case, crown, grid, R_ant=2.5, d_fossa=None):
    "g conditions that geometry can test without a registered antagonist:\n    - Proximal: minimum distance crown -> neighbor enamel (L5, L7) [mm]; < 0 = penetration (ANTAGANDE : OpenMandible -\n      the arc is in contact; the original tooth gives the reference value)\n    - occlusal : antagonist sphere (radie R_ant ) set in the crown's fossa; penetrating the sphere set in\n      The original form (before design change) = interference [mm ]."
    nb = neighbours_in_frame(case['R'], case['c'])
    out = {}
    for (t, (V, F)) in nb.items():
        d = grid.sample(crown, V, cval=10.0)
        out[f'prox_min_mm_{t}'] = float(d.min())
    return out
