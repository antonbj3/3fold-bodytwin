from common import *
from mechanics import *
from scipy.ndimage import map_coordinates
from scipy.optimize import minimize, LinearConstraint
import time, trimesh, zipfile, xml.etree.ElementTree as ET

def three_mf(mesh, path):
    ns = 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
    ET.register_namespace('', ns)
    r = ET.Element('{' + ns + '}model', unit='millimeter')
    res = ET.SubElement(r, '{' + ns + '}resources')
    o = ET.SubElement(res, '{' + ns + '}object', id='1', type='model')
    m = ET.SubElement(o, '{' + ns + '}mesh')
    v = ET.SubElement(m, '{' + ns + '}vertices')
    f = ET.SubElement(m, '{' + ns + '}triangles')
    for p in mesh.vertices:
        ET.SubElement(v, '{' + ns + '}vertex', dict(zip(['x', 'y', 'z'], [format(float(z), '.17g') for z in p])))
    for p in mesh.faces:
        ET.SubElement(f, '{' + ns + '}triangle', dict(zip(['v1', 'v2', 'v3'], map(str, p))))
    ET.SubElement(ET.SubElement(r, '{' + ns + '}build'), '{' + ns + '}item', objectid='1')
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for (name, data) in [('[Content_Types].xml', b'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>'), ('_rels/.rels', b'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>'), ('3D/3dmodel.model', ET.tostring(r))]:
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 3, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, data)

def run():
    start = time.perf_counter()
    a = np.load(X1 / 'inputs/geometry/D1_model.npz')
    grid = read(X1 / 'inputs/geometry/D1_grid.json')
    V = a['cV']
    T = a['cT']
    f = boundary(T)
    bn = np.unique(f)
    query = ((V - np.asarray(grid['origin'])) / grid['h']).T
    cav = map_coordinates(a['cavity'], query, order=1, mode='nearest')
    outer = map_coordinates(a['outer'], query, order=1, mode='nearest')
    zm = grid['z_m']
    fixed = bn[(abs(cav[bn]) < 0.2) | (V[bn, 2] < zm + 0.2)]
    fixed = np.unique(fixed)
    top = (V[:, 2] > zm + 1.0) & (outer < 0.25)
    weight = np.clip((V[:, 2] - zm - 1.0) / 2.0, 0, 1) * top
    weight[fixed] = 0
    xy = V[:, :2]
    q = (xy - xy.mean(0)) / np.maximum(np.max(abs(xy - xy.mean(0)), axis=0), 1e-08)
    (x, y) = q.T
    basis = np.c_[np.ones(len(V)), x, y, np.exp(-4 * ((x - 0.35) ** 2 + (y - 0.2) ** 2))]
    modes = np.zeros((len(V), 3, 4))
    modes[:, 2, :] = weight[:, None] * basis
    saved = np.load(X1 / 'inputs/fe/L005_lo_M1_k6_thin_3Y.npz')
    center = np.asarray(saved['lc__axial'])[0, :3]
    topnodes = bn[(outer[bn] < 0.2) & (V[bn, 2] > zm + 1.0)]
    ww = np.exp(-np.sum((V[topnodes, :2] - center[:2]) ** 2, axis=1) / (2 * 0.6 ** 2))
    ww /= ww.sum()
    load = np.zeros_like(V)
    load[topnodes, 2] = -ww
    model = Elastic(V, T, fixed, load, modes=modes)
    zero = np.zeros(4)
    base = model.evaluate(zero)
    (inv, B, v0, sign) = geometry(V, T)
    dv = []
    for j in range(4):
        dP = np.concatenate([np.zeros((len(T), 4, 1)), modes[T, :, j]], axis=2)
        dv.append(v0 * np.einsum('eij,eji->e', inv, dP))
    dv = np.array(dv).T
    volume_grad = dv.sum(0)
    M = modes[:, 2, :]
    constraints = [LinearConstraint(M, -0.2, 0.2), LinearConstraint(dv, -0.5 * v0, np.inf), LinearConstraint(volume_grad[None, :], -np.inf, 0.0)]
    history = []

    def fun(p):
        r = model.evaluate(p)
        return (r['J'], r['gradient'])

    def cb(p):
        r = model.evaluate(p)
        history.append(dict(J=r['J'], volume_mm3=r['volume_mm3']))
    sol = minimize(fun, zero, jac=True, method='SLSQP', constraints=constraints, bounds=[(-0.2, 0.2)] * 4, callback=cb, options=dict(maxiter=25, ftol=1e-08))
    final = model.evaluate(sol.x)
    gc = gradient_check(model, sol.x)
    model.evaluate(sol.x)
    ratio = float(np.exp(base['J'] - final['J']))
    new = final['vertices']
    identity = float(np.max(abs(new[fixed] - V[fixed])))
    vnew = geometry(new, T)[2]
    dest = ROOT / 'exports/whole_STS_D1'
    dest.mkdir(parents=True, exist_ok=True)
    meshes = []
    for (name, verts) in [('reference', V), ('optimized', new)]:
        mesh = trimesh.Trimesh(verts, f, process=True)
        if mesh.volume < 0:
            mesh.invert()
        mesh.export(dest / (name + '.stl'))
        three_mf(mesh, dest / (name + '.3mf'))
        loaded = trimesh.load_mesh(dest / (name + '.stl'), process=True)
        err = abs(loaded.volume / mesh.volume - 1)
        meshes.append(dict(name=name, watertight=bool(loaded.is_watertight), winding_consistent=bool(loaded.is_winding_consistent), volume_mm3=float(loaded.volume), roundtrip_volume_relative_error=err, pass_gate=bool(loaded.is_watertight and loaded.volume > 0 and (err <= 1e-05))))
    np.savez_compressed(DATA / 'WHOLE_RAW.npz', original_vertices=V, optimized_vertices=new, tetra=T, faces=f, fixed_nodes=fixed, base_principal=base['positive_principal'], optimized_principal=final['positive_principal'], base_stress=base['stress_MPa_per_N'], optimized_stress=final['stress_MPa_per_N'], load=load, modes=modes)
    out = dict(round='R3', claim_type='capability', source_geometry=str(X1 / 'inputs/geometry/D1_model.npz'), source_sha256=sha(X1 / 'inputs/geometry/D1_model.npz'), coefficients_mm=sol.x, optimizer_success=bool(sol.success), optimizer_message=str(sol.message), iterations=int(sol.nit), FEM_queries=model.calls, model_force_quantile_ratio=ratio, calibrated_force05_N=None, base_volume_mm3=base['volume_mm3'], optimized_volume_mm3=final['volume_mm3'], intaglio_fixed_coordinate_identity_error_mm=identity, min_signed_tetra_volume_ratio=float(np.min(vnew / v0)), equilibrium_relative_residual=final['equilibrium_relative_residual'], adjoint_check=gc, exports=meshes, history=history, wall_full_surface_certificate='UNKNOWN', insertion='Inherited X1 cavity coordinates retained exactly; exterior deformed, complete trajectory not re-certified', cement='Intaglio nodes identical; physical seating film UNKNOWN', milling='UNKNOWN_FULL_CROWN_TOOL_AND_SHAFT_ACCESS', antagonist='UNKNOWN_NO_REGISTERED_STS_ANTAGONIST', material_volume_constraint='Fixed, not substance reduction', substance_removal_mm3=next((r['removed_mm3'] for r in read(X1 / 'inputs/X1_ROWS.json') if r['id'] == 'L005_lo_M1_k6_thin_3Y')), substance_removal_scope='Inherited virtual preparation volume; same preparation for both designs, not measured clinical removal', quantile_scope='Conditional rigid-support fixed-load independent-volume-flaw model with borrowed m; no absolute calibration', external_referent=read(ROOT / 'PREREG_R3.json')['external_referent'], gates=dict(adjoint=gc['pass_gate'], equilibrium=final['equilibrium_relative_residual'] <= 1e-08, volume=final['volume_mm3'] <= base['volume_mm3'] * (1 + 1e-08), intaglio_identity=identity == 0, min_tetra_volume=float(np.min(vnew / v0)) >= 0.5 - 1e-08, exports=all((m['pass_gate'] for m in meshes)), all_physical_constraints=False), seconds=time.perf_counter() - start)
    dump(ROOT / 'raw/R3.json', out)
    freeze(ROOT / 'FROZEN_WHOLE_PREDICTIONS.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R3.json'), predicted_ratio=ratio, ratio_operational_rejection_window=[ratio * np.exp(-0.3), ratio * np.exp(0.3)], calibrated_force05_N=None, exports=inventory(dest), raw_sha256=sha(DATA / 'WHOLE_RAW.npz'), measurement_status='NO_PHYSICAL_MEASUREMENT', physical_constraints='UNRESOLVED_WALL_MILLING_ANTAGONIST', resolution='PER_TOOTH', time_scale='HANDOVER'))
    state('R3_WHOLE_PREDICTIONS_FROZEN', out['gates'], 'Verify exports and create locked lab comparison; no physical approval claim')
    (ROOT / 'HANDOFF_R3.md').write_text('R3 whole STS crown exported with locally differentiated sparse elasticity and volume/intaglio constraints. Published force scale still cannot validate it, full tool/shaft, antagonist and whole wall certificates are UNKNOWN. Frozen matched force ratio allows same-batch lab rejection without fitting after measurement. Next construction: acquire a rigid/compliant support contrast and fracture origin for the optimized/reference pair, or add registered whole-crown antagonist/CAM certificates.\n')
    print(clean({k: out[k] for k in ['model_force_quantile_ratio', 'gates', 'FEM_queries', 'seconds']}))
if __name__ == '__main__':
    run()
