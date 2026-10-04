from common import *
from native_port import model, deck_contract, q_matrix
from support_coverage import interpolate10, assemble, reduce_solver, kkt_solver

def run():
    t0 = time.perf_counter()
    p = DATA / 'NATIVE_SOLVER/solve_WarnNodeMissTiedContact.nam'
    if not p.exists():
        raise ValueError('Actual native solver omission artifact missing')
    missed = np.array([int(l.strip()) - 1 for l in p.read_text().splitlines() if l.strip() and (not l.strip().startswith('*'))])
    (a, g, m) = model()
    port = np.load(DATA / 'NATIVE_PORT.npz')
    slaves = port['slaves']
    active = ~np.isin(slaves, missed)
    corners = slaves < len(a['cV'])
    chosen = active & corners
    usable = active & (port['distances_mm'] <= 0.06)
    native_distance = port['distances_mm'][active]
    nc = len(a['cV'])
    V = np.vstack([a['cV'], a['dV']])
    T = np.vstack([a['cT'], a['dT'] + nc])
    ne = len(a['cT'])
    fixed = nc + m.fix_nodes[m.fix_nodes < len(a['dV'])]
    E = np.r_[np.full(ne, 201400.0), np.full(len(a['dT']), 18000.0)]
    nu = np.r_[np.full(ne, 0.3), np.full(len(a['dT']), 0.31)]
    (K, *_) = assemble(V, T, E, nu)
    ss = slaves[chosen]
    mt = port['master_nodes'][chosen, :3] + nc
    w = port['barycentric'][chosen]
    (Q, map_) = q_matrix(len(V), ss, mt, w)
    (reduced, q) = reduce_solver(K, Q, map_[fixed])
    full = kkt_solver(K, ss, mt, w, fixed)
    deck = deck_contract()
    Lc = interpolate10(a['cV'], m.Vc, a['cT'], m.Tc)
    Ld = interpolate10(a['dV'], m.Vd, a['dT'], m.Td)
    rows = []
    dis = []
    base = read(R / 'raw/R2_SUPPORT_COVERAGE.json')['contrasts_P1_conditional']['native_all']
    for (i, f10) in enumerate(deck['loads']):
        f = np.vstack([Lc.T @ f10[:m.nc], Ld.T @ f10[m.nc:]])
        (u, res) = reduced(f)
        (uk, rk) = full(f)
        err = float(np.linalg.norm(uk - u) / np.linalg.norm(u))
        energy = float(0.5 * np.dot(f.ravel(), u))
        rows.append({'angle_deg': [0, 30][i], 'strain_energy_N_mm': energy, 'relative_change_from_all_declared_P1': energy / base[i]['strain_energy_N_mm_per_unit_load'] - 1, 'reduced_residual': res, 'same_Q_KKT_relative_displacement_error': err, 'KKT_residual': rk, 'same_Q_gate': err <= 1e-07 and max(res, rk) <= 1e-08})
        dis.append(u)
    native_set_receipt = {'declared': slaves.tolist(), 'omitted': missed.tolist(), 'effective': slaves[active].tolist(), 'usable': slaves[usable].tolist()}

    def accepts(receipt):
        return len(receipt['omitted']) == 0 and len(receipt['effective']) == len(receipt['declared']) and (len(receipt['usable']) == len(receipt['declared']))
    positive = {'declared': [1, 2], 'omitted': [], 'effective': [1, 2], 'usable': [1, 2]}
    fault = {'declared': [1, 2], 'omitted': [2], 'effective': [1], 'usable': [1]}
    fault2 = {'declared': [1, 2], 'omitted': [], 'effective': [1, 2], 'usable': [1]}
    out = {'round': 'R3_COMPILED_SUPPORT', 'claim_type': 'capability', 'native_solver_warning_sha256': sha(p), 'counts': {'declared_nodes': len(slaves), 'compiled_omitted_nodes': len(missed), 'compiled_omitted_fraction': len(missed) / len(slaves), 'compiled_active_nodes': int(active.sum()), 'effective_distance_qualified_nodes': int(usable.sum()), 'effective_over_006_nodes': int((active & ~usable).sum()), 'P1_active_corner_slaves': len(ss)}, 'effective_projection_mm': {'max': float(native_distance.max()), 'p95': float(np.quantile(native_distance, 0.95)), 'resolution': 'PER_POINT'}, 'P1_diagnostics': rows, 'receipt_controls': {'valid_positive': accepts(positive), 'omitted_fault_rejected': not accepts(fault), 'distance_coverage_fault_rejected': not accepts(fault2), 'actual_source_receipt_accepted': accepts(native_set_receipt)}, 'gate': {'declared_native_coverage': len(missed) == 0, 'effective_projection': float(native_distance.max()) <= 0.06, 'same_Q_full_FE': all((x['same_Q_gate'] for x in rows)), 'optimization_allowed': False}, 'native_full_coefficient_parity': 'UNKNOWN_NOT_EXTRACTED', 'cost': cost(t0)}
    dump(R / 'raw/R3_COMPILED_SUPPORT.json', out)
    dump(R / 'raw/NATIVE_SUPPORT_RECEIPT.json', native_set_receipt)
    np.savez_compressed(DATA / 'COMPILED_SUPPORT.npz', active_mask=active, usable_mask=usable, omitted_nodes=missed, P1_active_corner_slaves=ss, P1_displacements=np.array(dis))
    state('R3_COMPILED_SUPPORT_DECIDED', out['gate'], 'Freeze negative predictions, bind source facits, spatial failure export and clean replay')
    print(json.dumps(clean(out), indent=2))
if __name__ == '__main__':
    run()
