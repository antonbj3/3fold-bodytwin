from common import *
import csv, time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import trimesh

def run():
    A = read(ROOT / 'raw/A.json')
    B = read(ROOT / 'raw/B_FIT.json')
    BB = read(ROOT / 'raw/B_BOUNDS.json')
    C = read(ROOT / 'raw/C2.json')
    D = read(ROOT / 'raw/D.json')
    V = read(ROOT / 'raw/VERIFICATION.json')
    z = np.load(DATA / 'B_FIELDS.npz')
    gr = np.load(DATA / 'B_GRIDS.npz')
    reg = np.load(DATA / 'C_ADVERSE_REGION.npz')
    dest = ROOT / 'exports'
    dest.mkdir(exist_ok=True)
    exports = []
    for (name, key) in [('reference', 'vertices_reference'), ('optimized', 'vertices_optimized')]:
        mesh = trimesh.Trimesh(z[key], z['faces'], process=False)
        mesh.fix_normals()
        path = dest / (name + '_research_only.stl')
        mesh.export(path)
        loaded = trimesh.load_mesh(path, process=True)
        exports.append(dict(name=name, path=str(path), sha256=sha(path), watertight=bool(loaded.is_watertight), winding_consistent=bool(loaded.is_winding_consistent), positive_volume=bool(loaded.volume > 0), volume_mm3=float(loaded.volume), relative_volume_serialization_error=float(abs(loaded.volume / mesh.volume - 1)), units='mm externally declared; STL has no unit field', qualification='RESEARCH_GEOMETRY_ONLY; full wall/insertion/cement/CAM/antagonist NOT_CERTIFIED'))
    if not (ROOT / 'FROZEN_EXPORTS.json').exists():
        freeze(ROOT / 'FROZEN_EXPORTS.json', dict(exports=exports))
    else:
        old = read(ROOT / 'FROZEN_EXPORTS.json')['exports']
        assert all((x['sha256'] == y['sha256'] for (x, y) in zip(exports, old)))
    with (ROOT / 'PER_TOOTH.csv').open('w') as f:
        fields = ['key', 'family', 'eligible', 'force_target_N', 'measured_bite_pose', 'full_milling_certificate', 'strength_Pareto', 'force_Pareto', 'substance_Pareto', 'reason']
        w = csv.DictWriter(f, fields)
        w.writeheader()
        for row in A['cohort']:
            r = {k: row.get(k) for k in fields}
            r.update(strength_Pareto='UNKNOWN', force_Pareto='UNKNOWN', substance_Pareto='UNKNOWN')
            w.writerow(r)
    cost = dict(A_seconds=A['seconds'], B_fit_seconds=B['seconds'], B_interval_seconds=BB['seconds'], C_initial_seconds=read(ROOT / 'raw/C.json')['seconds'], C2_seconds=C['seconds'], D_seconds=D['seconds'], verification_seconds=V['seconds'], preparation_and_source_review_seconds=None, physical_acquisition_seconds=None, agent_tokens=None, peak_measured_process_MiB=max(B['peak_rss_MiB'], V['peak_rss_MiB']), B_factorizations=sum(B['factorizations'].values()), D_cells=sum((r['cells_evaluated'] for r in D['supports'].values())), B_cells=sum((r['continuous']['cells_evaluated'] for r in BB['contrasts']['optimized'].values())), own_data_bytes=sum((p.stat().st_size for p in DATA.glob('*') if p.is_file())), threads_current=1, initial_C_HiGHS_threads='UNRECORDED; preserved limitation, C2 explicit1')
    cost['own_data_bytes_including_replays'] = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    cost['resource_gate_waits'] = {name: (ROOT / 'raw' / name).read_text().count('[ heavy_run ] wait') for name in ['B_RUN.log', 'B_BOUNDS_RUN.log', 'D_RUN.log', 'VERIFY_RUN.log', 'RUN_ALL_FAILED_SCRIPT_EDIT.log', 'REFIT_RUN.log', 'RUN_ALL.log'] if (ROOT / 'raw' / name).exists()}
    cost['resource_gate_logged_sleep_seconds'] = 30 * sum(cost['resource_gate_waits'].values())
    cost['verification_initial_seconds'] = read(ROOT / 'raw/VERIFICATION_INITIAL.json')['seconds']
    failed_launcher = ROOT / 'raw/RUN_ALL_FAILED_SCRIPT_EDIT.log'
    if failed_launcher.exists():
        measurements = [json.loads(line) for line in failed_launcher.read_text().splitlines() if line.startswith('{')]
        cost['verification_in_failed_launcher_seconds'] = measurements[-1]['seconds'] if measurements else None
    if (ROOT / 'raw/REFIT_REPLAY.json').exists():
        cost['isolated_optimizer_replay'] = read(ROOT / 'raw/REFIT_REPLAY.json')
        replay_fit = read(ROOT / 'sources/HISTORICAL_REFIT_FIT.json')
        cost['peak_measured_process_MiB'] = max(cost['peak_measured_process_MiB'], replay_fit['peak_rss_MiB'])
        cost['optimizer_factorizations_including_replay'] = cost['B_factorizations'] + sum(replay_fit['factorizations'].values())
    artifacts = [dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(DATA.glob('*')) if p.is_file()]
    result = dict(claim_type='capability', outcome='NO_CERTIFIED_INVERSE_CROWN_BREAKTHROUGH; CONDITIONAL_RELATIVE_ORDER_OPERATOR_AND_PROCESS_ROBUSTNESS_THRESHOLD', review_state='PENDING_INDEPENDENT_REVIEW', external_referent=dict(kind='independent_measurement', locator=['doi:10.4047/jap.2021.13.5.269 Table1, Methods, Fig5', 'doi:10.3390/ma17020365 Table2, Methods', 'doi:10.5281/zenodo.10597292', 'https://osf.io/xctdy/'], compared_quantity='Published fitted crown-force distribution and protocol parameters; local independently acquired tooth geometry. No measured optimized/reference fracture pair exists.', refutes_us=True), facit_scope='Published data anchors scenarios and rejects universal protocol transfer; it does NOT validate our generated relative fracture gain', cohort=dict(total=18, by_family={s: dict(total=6, full_certified_gain=0, physical_Pareto='UNKNOWN') for s in ['anterior', 'premolar', 'molar']}, excluded=18, excluded_fraction=1, reason='Missing verified bite pose/force targets/full CAM; predecessor geometry offset failure'), separate_specimen=dict(dataset='STS-Tooth3D', source=X1 / 'inputs/geometry/D1_model.npz', type='molar', patient_cohort_substitution=False, preparation_changed=False, substance_saving_mm3=0.0, substance_scope='Same virtual preparation, not measured retained dentin', reference='anatomical-shape ceramic crown'), B_fit=B, B_continuous_Weyl=BB['contrasts'], C_density_threshold=C['worst'], C_initial_control_failures_preserved=True, D_dual_certificate=D, certified_physical_gain=False, absolute_fracture_quantile_N=None, hard_conditions=dict(insertion_Farkas='MISSING_FOR_FINAL_SHAPE', minimum_wall_material='NOT_CERTIFIED', cement_gap='UNKNOWN_PHYSICAL; unchanged intaglio only', milling='UNKNOWN_FULL_TOOL_SHAFT_HOLDER', antagonist_nonpenetration='UNKNOWN_REGISTERED_BITE', tooth_force_targets='UNKNOWN', FE_continuum_error='MISSING', floating_point_enclosure='MISSING'), cost=cost, artifacts=artifacts, verification=V, exports=exports, edges=[dict(producer='load-basis local stress tensors', consumer='continuous relative hazard', resolution='PER_POINT', time_scale='SIMULTANEOUS'), dict(producer='spatial flaw/process intensity', consumer='relative hazard', resolution='PER_POINT', time_scale='SIMULTANEOUS'), dict(producer='manufactured geometry and film', consumer='loaded fracture specimen', resolution='PER_SURFACE_REGION', time_scale='HANDOVER')], phenomenological_debts=['common volume flaw intensity and exponent: replace with fracture-origin/process-matched spatial measurement', 'support stiffness/tie and fixed load patch: replace with measured support/contact/compliance', 'material wall/CAM and cement: replace with same-object continuous geometric/metrology certificates'])
    result['sufficiency'] = A['sufficiency']
    result['sufficiency']['external_referent'] = dict(kind='our_own_fixture', locator=str(ROOT / 'raw/A.json'), compared_quantity='Exact two-region summary versus spatially weighted hazard; mathematical counterexample only', refutes_us=True)
    result['protocol_deviations'] = read(ROOT / 'PROTOCOL_DEVIATIONS.json')
    result['parameter_scope'] = 'm is a sensitivity range over published whole-crown fitted confidence endpoints; angle0..30deg interpolates distinct published protocols as our scenario, not a measured patient load interval'
    result['full_multiobjective_spec'] = 'DERIVATION.md; all18 physical Pareto fronts UNKNOWN; separate specimen only tests relative mechanics at fixed preparation'
    result['reproduction'] = dict(command='./run_all.sh', fresh_optimizer_command='./run_all.sh --refit', fit_replay=cost.get('isolated_optimizer_replay'))
    save(ROOT / 'results.json', result)
    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=(14, 4.7))
    ax = fig.add_subplot(131)
    mindex = int(np.argmin(abs(gr['ms'] - 6)))
    for (label, color) in [('R3', '#b24d34'), ('R5', '#777777'), ('optimized', '#126a80')]:
        ax.plot(gr['angles'], 100 * (gr[label + '__18000'][:, mindex] - 1), label=label, color=color)
    ax.axhline(0, c='black', lw=0.7)
    ax.axhline(2, c='#555555', ls='--', label='frozen +2% gate')
    ax.set(xlabel='Load angle (degree)', ylabel='Relative fracture-quantile change (%)', title='Conditional FE; 18 GPa support\nfixed patch, m ≈ 6')
    ax.legend(fontsize=8)
    ax = fig.add_subplot(132)
    for (i, (name, r)) in enumerate(D['supports'].items()):
        lo = 100 * (r['minimum_force_ratio_lower'] - 1)
        hi = 100 * (r['minimum_force_ratio_sample_upper'] - 1)
        ax.plot([lo, hi], [i, i], c='#126a80', lw=4)
        ax.scatter([lo, hi], [i, i], s=18, c='#126a80')
    ax.axvline(0, c='black', lw=0.7)
    ax.axvline(2, c='#555555', ls='--')
    ax.set_yticks([0, 1, 2], ['8.6 GPa', '18 GPa', 'Rigid'])
    ax.set(xlabel='Worst ratio change over angle/m box (%)', title='Analytic discrete-model enclosure\nnot a physical certificate')
    ax = fig.add_subplot(133, projection='3d')
    cent = z['vertices_optimized'][z['tetra']].mean(1)
    mask = reg['adverse_mask']
    ix = np.arange(0, len(cent), 12)
    ax.scatter(*cent[ix].T, c=np.where(mask[ix], '#b24d34', '#bbbbbb'), s=2, alpha=0.6)
    ax.set(xlabel='x (mm)', ylabel='y (mm)', zlabel='z (mm)', title=f"Adverse flaw-density region\nreversal at κ = {C['worst']['kappa']:.4f}")
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/inverse_crown.png', dpi=170)
    plt.close(fig)
    lines = ['# Inverse crown: conditional mechanics resolved, physical design still open', '', 'No crown in the frozen six-jaw/18-tooth cohort can be claimed manufacturable and robustly stronger. The physical Pareto front remains UNKNOWN for all six anterior, six premolar and six molar sites. All source array hashes were checked; no failed case was replaced.', '', 'A separate whole STS molar now has a runnable minimax shape optimizer and an analytic continuous angle/modulus ordering bound. This is a conditional discrete FE capability, not empirical crown validation.', '', '| Support scenario | Worst sampled force ratio | Analytic lower bound | Useful >1.02 gain |', '|---|---:|---:|---|']
    for (s, r) in D['supports'].items():
        lines.append(f'| {s} MPa'.replace('rigid MPa', 'Rigid') + f" | {BB['contrasts']['optimized'][s]['sampled_min']:.8f} | {r['minimum_force_ratio_lower']:.8f} | FAIL |")
    lines += ['', f"All ratios: PER_TOOTH; continuous m∈[2,15.6], angle∈[0,30]° in one plane, fixed Gaussian load patch, homogeneous volume flaws. Support stiffness is a finite scenario set, not an enclosed continuum. Same virtual preparation and crown volume; zero claimed substance saving. The worst sampled gain is {100 * (min((x['sampled_min'] for x in BB['contrasts']['optimized'].values())) - 1):.4f}%, below the frozen 2% criterion.", '', f"The changed flaw representation finds ranking reversal at common spatial density contrast κ={C['worst']['kappa']:.8f}, PER_TOOTH in one tested scenario. This is a derived process-uniformity requirement, not measured flaw heterogeneity. The exact summary witness has identity error 0 and hazards17/8 after spatial weighting. Minimum extension: joint local stress and flaw measure.", '', 'External referents: Prott2021 Table1/Methods/Fig5 (doi:10.4047/jap.2021.13.5.269), Chen2024 Table2/Methods (doi:10.3390/ma17020365), STS geometry (doi:10.5281/zenodo.10597292), and original Teeth3DS geometry (https://osf.io/xctdy/). Published crown-force fits anchor sensitivity parameters; none is a matched optimized-crown fracture measurement. All10 Prott groups are retained as source context, 0/10 qualify as direct validation of this design.', '', 'B’s initial Weyl bounds failed the frozen 0.002 log-width criterion. D changes the bound to a joint convex dual construction, tightening the range substantially. D’s stricter0.0002 target still fails for both compliant supports and passes only for rigid support. Exact arithmetic proves the stated inequality; directed floating arithmetic and continuum FE error are missing. The small positive lower bounds therefore do not prove a physical gain.', '', f"Controls: inherited R3/R4/R5 shapes were evaluated on the same information. R3’s worst ratio drops to {min((x['sampled_min'] for x in BB['contrasts']['R3'].values())):.6f}; R4 also reverses. Independent direct loading, adjoint differences and bounded LPs pass after the preserved C numerical failure was repaired by objective scaling. No optimizer superiority or global optimality is claimed.", '', f"Resource record: B {B['seconds']:.1f}s, {cost['B_factorizations']} factorizations, measured peak {cost['peak_measured_process_MiB']:.1f}MiB. Own data {cost['own_data_bytes'] / 1000000.0:.1f}MB. Historical/source-reading cost is unmeasured; the initial HiGHS thread count was unrecorded and C2 explicitly uses one thread.", '', 'Missing hard conditions: final insertion/Farkas certificate, material wall bound, spatial seated cement film, complete cutter/shaft/holder access, registered antagonist nonpenetration and measured tooth-force targets. The surface tie retains the predecessor’s failed projection gate. Relative strength does not erase these missing prerequisites.', '', 'Run `./run_all.sh`. Read README_DEMO.md, DERIVATION.md and LAB_PROTOCOL.md. Research-only STL files are not released manufacturing designs. Status PENDING_INDEPENDENT_REVIEW.', '', '![Conditional inverse-crown results](figures/inverse_crown.png)', '']
    (ROOT / 'RESULTS.md').write_text('\n'.join(lines))
    if not (ROOT / 'FROZEN_PREDICTIONS.json').exists():
        freeze(ROOT / 'FROZEN_PREDICTIONS.json', dict(claim_type='capability', geometry_exports=exports, field_sha256=B['field_sha256'], relative_quantile_bounds={k: {kk: vv for (kk, vv) in val.items() if kk in ['minimum_force_ratio_lower', 'minimum_force_ratio_sample_upper']} for (k, val) in D['supports'].items()}, spatial_reversal_threshold=C['worst'], physical_fracture_measurement='NOT_RUN', individual_pair_failure_order='UNKNOWN; quantile ordering does not determine an individual random pair', physical_certification='NOT_ADMITTED', absolute_force_N=None))
    print(json.dumps(dict(report='results.json', exports=len(exports), figure='figures/inverse_crown.png', verification=V['pass_gate'])))
if __name__ == '__main__':
    run()
