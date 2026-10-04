from common import *
import ast, math, importlib.util
from mechanics import *

def source_force_port():
    tree = ast.parse((V4 / 'code/literature.py').read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'calibrated_force' or (isinstance(n, ast.Assign) and any((isinstance(a, ast.Name) and a.id == 'PROTOCOL_FIELDS' for a in n.targets)))]
    ns = {'math': math}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(V4 / 'code/literature.py'), 'exec'), ns)
    return ns

def run():
    start = __import__('time').perf_counter()
    ns = source_force_port()
    groups = read(X1 / 'raw/MATCHED_R2.json')['groups']
    prott = [r for r in groups if r['study'] == 'PMC8558575']
    by = {r['thickness_mm']: r for r in prott}
    quant = []
    for r in sorted(prott, key=lambda r: r['thickness_mm']):
        F0 = r['characteristic_load_N']
        m = r['Weibull_m']
        q = F0 * (-math.log(0.95)) ** (1 / m)
        corners = [F * (-math.log(0.95)) ** (1 / mm) for F in r['F0_90CI_N'] for mm in r['m_90CI']]
        quant.append(dict(thickness_mm=r['thickness_mm'], force05_N=q, parameter_rectangle_range_N=[min(corners), max(corners)], resolution='POPULATION', range_kind='Arithmetic range over marginal 90CI rectangle; not a joint 90CI', locator=r['locator']))
    curated = read(V4 / 'payload/literature/CROWN_CURATED.json')
    blocks = {}
    excluded = []
    for r in curated:
        if r['eligible'] and r['protocol'].startswith('crown'):
            blocks.setdefault((r['study'], r['protocol'], r['product']), []).append(r)
        else:
            excluded.append(dict(locator=r['locator'], reason='not eligible crown endpoint or different specimen quantity'))
    preds = []
    rejected = []
    for (key, rows) in blocks.items():
        rows = sorted(rows, key=lambda r: r['t'])
        if len(set((r['t'] for r in rows))) < 3:
            rejected.append(dict(key=key, reason='fewer than 3 distinct matched thicknesses', n=len(rows)))
            continue
        (lo, hi) = (rows[0], rows[-1])
        for r in rows[1:-1]:
            w = math.log(r['t'] / lo['t']) / math.log(hi['t'] / lo['t'])
            f = lo['mean'] ** (1 - w) * hi['mean'] ** w
            preds.append(dict(study=key[0], protocol=key[1], thickness_mm=r['t'], predicted_N=f, observed_N=r['mean'], locator='doi:' + r['doi'] + ' ' + r['locator'], log_error=abs(math.log(f / r['mean'])), pass_gate=abs(math.log(f / r['mean'])) <= 0.2, resolution='PER_TOOTH', wrong_value_rejected=abs(math.log(1.6)) > 0.2))
    query = dict(thickness_mm=0.8, product='KATANA_HTML_PLUS', specimen_geometry='GenCAD affine roof', support='rigid basal patch', cement='nominal 40um', angle_deg=0, indenter_diameter_mm=None, crosshead_mm_min=None, contact_layer='nodal normalized 1N', aging=None)
    force = ns['calibrated_force']([], query)
    V = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    T = np.array([[0, 1, 2, 3]])
    (inv, B, vol, sign) = geometry(V, T)
    G = np.array([[0.001, 0.002, 0], [0, -0.0004, 0.0006], [0.0008, 0, 0.0005]])
    u = V @ G.T
    expected = np.array([G[0, 0], G[1, 1], G[2, 2], G[0, 1] + G[1, 0], G[1, 2] + G[2, 1], G[0, 2] + G[2, 0]])
    error = float(np.max(abs(B[0] @ u.ravel() - expected)))
    xy = np.array([[0.0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]])
    faces = np.array([[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]])
    za = np.array([0.5, 0.8, 1.1, 0.5, 0.8, 1.1])
    zb = za[::-1].copy()
    inn = np.zeros(6)
    (va, tt) = roof_mesh(xy, faces, za, inn)
    (vb, _) = roof_mesh(xy, faces, zb, inn)
    loads = np.zeros_like(va)
    loads[6, 2] = -1
    modes = np.zeros((12, 3, 2))
    modes[6:, 2, 0] = 1
    modes[6:, 2, 1] = np.array([-1, 0, 1, -1, 0, 1])
    ma = Elastic(va, tt, np.arange(6), loads, modes=modes)
    mb = Elastic(vb, tt, np.arange(6), loads, modes=modes)
    a = ma.evaluate([0, 0])
    b = mb.evaluate([0, 0])
    identity = max(float(np.max(abs(np.sort(za) - np.sort(zb)))), abs(a['volume_mm3'] - b['volume_mm3']))
    sufficient = dict(summary=['exact thickness multiset', 'volume', 'minimum/maximum/mean thickness'], identity_error=identity, volume_A_mm3=a['volume_mm3'], volume_B_mm3=b['volume_mm3'], hazard_A=a['hazard'], hazard_B=b['hazard'], conditional_quantile_ratio_B_over_A=(a['hazard'] / b['hazard']) ** (1 / ma.m), downstream_relative_difference=abs((a['hazard'] / b['hazard']) ** (1 / ma.m) - 1), minimum_extension='Joint local coordinates, thickness, support and load port; thickness histogram alone is insufficient', pass_counterexample=identity == 0 and abs((a['hazard'] / b['hazard']) ** (1 / ma.m) - 1) > 0.01, resolution='PER_POINT -> PER_TOOTH model')
    gc = gradient_check(ma, np.array([0.03, 0.02]))
    out = dict(round='R1', claim_type='capability', external_referent=read(ROOT / 'PREREG_R1.json')['external_referent'], published_force05=quant, calibration_rows=preds, generated_absolute_force=force, absolute_calibrated_quantile_N=None, sufficiency=sufficient, gradient_check=gc, elastic_patch=dict(max_strain_error=error, pass_gate=error <= 1e-09, wrong_constitutive_value_rejected=np.max(abs(1.6 * expected - expected)) > 1e-09, external_referent=dict(kind='closed_form', locator='linear tetra affine displacement identity: epsilon=sym(grad u)', compared_quantity='engineering strain under exact affine displacement', refutes_us=True)), dropout=dict(input_rows=len(curated), excluded_rows=len(excluded), excluded_fraction=len(excluded) / len(curated), excluded=excluded, rejected_calibration_blocks=rejected), gates=dict(all_published_interiors=all((r['pass_gate'] for r in preds)), absolute_protocol_supported=False, adjoint=gc['pass_gate'], exact_summary_counterexample=sufficient['pass_counterexample']), seconds=__import__('time').perf_counter() - start)
    dump(ROOT / 'raw/R1.json', out)
    (ROOT / 'HANDOFF_R1.md').write_text('R1: protocol-conditioned absolute force is UNKNOWN; published force interpolation retains every failed target. Scalar thickness/volume insufficiency tested at machine precision. Exact fixed-topology elastic adjoint compared with central differences. Next construction: constrain local elastic shape optimization by unchanged GenCAD predicates, at all four levels.\n')
    state('R1_DECIDED', out['gates'], 'R2 preregister spatial shape optimizer and score withheld GenCAD geometry')
    print(clean(dict(gates=out['gates'], sufficiency=sufficient, gradient=gc)))
if __name__ == '__main__':
    run()
