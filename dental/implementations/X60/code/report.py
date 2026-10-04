from common import *
import collections, csv, time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def run():
    r1 = read(ROOT / 'raw/R1.json')
    r2 = read(ROOT / 'raw/R2.json')
    r3 = read(ROOT / 'raw/R3.json')
    r4 = read(ROOT / 'raw/R4.json')
    r5 = read(ROOT / 'raw/R5.json')
    ver = read(ROOT / 'raw/VERIFICATION.json')
    rows = read(ROOT / 'raw/SCORED_ROWS.json')
    lab = read(ROOT / 'raw/LAB_CONTROLS.json')
    summaries = []
    for name in ['constraint_optimizer', 'X1B_morphology', 'volume_1.0', 'volume_1.1']:
        rr = [r for r in rows if r['participant'] == name]
        passrows = [r for r in rr if r['L1'] == 'PASS']
        ratios = [r['optimization']['model_force_quantile_ratio'] for r in passrows if 'optimization' in r]
        summaries.append(dict(participant=name, tasks=len(rr), statuses=dict(collections.Counter((r['L1'] for r in rr))), quality_Pareto=sum((r.get('quality_Pareto', False) for r in rr)), conditional_model_gain_gt2pct=sum((r.get('meaningful_conditional_gain', False) for r in rr)), median_conditional_force_ratio=float(np.median(ratios)) if ratios else None, optimizer_nonconvergence=sum((r.get('optimization', {}).get('optimizer_success') == False for r in rr)), resolution='PER_TOOTH'))
    costs = {tag: read(DATA / 'predictions' / tag / 'COST.json') for tag in ['volume_1.0', 'volume_1.1']}
    interior = r1['calibration_rows']
    out = dict(claim_type='capability', status='SPATIAL_OPTIMIZER_AND_FROZEN_SPECIMENS_DELIVERED_ABSOLUTE_CALIBRATION_AND_MANUFACTURABILITY_UNKNOWN', review_state='PENDING_INDEPENDENT_REVIEW', capability='Executable local crown/bridge roof and whole-crown shape adjoints, GenCAD v4 spatial scoring, frozen paired laboratory contrasts with explicit calibration refusal', external_referent=dict(kind='independent_measurement', locator=sorted({r['locator'] for r in interior}), compared_quantity='Protocol-specific held-out interior arithmetic mean crown fracture force N; no optimized-crown endpoint has been measured', refutes_us=True), secondary_referents=[dict(kind='published_dataset', locator='doi:10.5281/zenodo.10597292', compared_quantity='STS source tooth geometry', refutes_us=True), dict(kind='published_dataset', locator='https://ditto.ing.unimore.it/bits2bites/', compared_quantity='Registered source roof and antagonist proximity masks', refutes_us=True)], results=dict(published_calibration_pass=sum((r['pass_gate'] for r in interior)), published_calibration_targets=len(interior), absolute_calibrated_force05_N=None, summaries=summaries, sufficiency=r1['sufficiency'], whole_rigid_ratio=r3['model_force_quantile_ratio'], whole_nearest_tie=r4['contrasts'], whole_surface_tie=r5['contrasts'], surface_tie_interface=r5['interface_projection_mm'], whole_final_gates=r5['gates'], physical_Pareto='UNKNOWN: uncalibrated force, clinical substance removal, full CAM and antagonist are missing'), dropout=dict(source_cases=64, selected_cases=12, selection_excluded_cases=52, selection_is_not_observation_dropout=True, requested_tasks=192, source_unavailable_tasks=r2['source_unavailable_tasks'], geometric_abstentions_per_candidate=61, geometric_abstention_fraction=61 / 192, abstention_reason='locked wall/film/continuous antagonist and bridge section inequalities infeasible; exact inherited morphology LP reports abstention', literature=r1['dropout']), controls=dict(R1_adjoint=r1['gradient_check'], R2_adjoint=r2['all_adjoint_checks'], R2_original_fault_gate=r2['all_faults_rejected'], valid_coordinate_faults=len(ver['valid_coordinate_faults']), valid_coordinate_all_rejected=ver['all_faults_rejected'], exports=ver['all_exports_pass'], export_faults=ver['all_export_faults_rejected'], lab_faults=all((r['positive_fixture'] and r['wrong_force_rejected'] and r['missing_support_refused'] for r in lab)), equal_information_reference='Frozen GenCAD v4 constraint_optimizer / X1B_morphology actually scored; no method novelty asserted'), cost=dict(generator_seconds={k: v['wall_seconds'] for (k, v) in costs.items()}, generator_peak_RSS_MiB={k: v['max_rss_MiB'] for (k, v) in costs.items()}, query_counts={k: sum((d['FEM_queries'] for d in read(DATA / 'predictions' / k / 'DIAGNOSTICS.json') if 'FEM_queries' in d)) for k in costs}, source_preparation_fit_acquisition_seconds='UNKNOWN_INHERITED; no new acquisition', scoring_seconds=r2['seconds'], whole_fit_and_validation_seconds={k: r['seconds'] for (k, r) in [('R3', r3), ('R4', r4), ('R5', r5)]}, discovery='Implementation/search effort is not zero; session elapsed bounded by prereg/report timestamps, detailed token attribution UNKNOWN', first_prereg_utc=read(ROOT / 'PREREG_R1.json')['frozen_utc'], report_utc=now(), physical_manufacture_measurement='NOT_RUN', fallback='61/192 geometric abstentions, one nonconverged SLSQP per roof budget retained; physical UNKNOWN returned', thread_limit=4, GPU=False, data_bytes=budget()), resolution_contract=dict(elastic_coordinates='PER_POINT', force_law_parameters='POPULATION', shape_quantile_ratio='PER_TOOTH conditional model', cement='PER_POINT nominal spacing only; physical film UNKNOWN', substance_removal='PER_TOOTH virtual source preparation; actual removal UNKNOWN', time_scale_simulation='SIMULTANEOUS', design_to_lab='HANDOVER'), phenomenological_debts=[dict(quantity='3Y E/nu', measurement='Batch modulus under hydration/aging protocol'), dict(quantity='Force-law m used as volume-flaw m', measurement='Matched crown individual forces, fracture origin and flaw-region characterization'), dict(quantity='Fixed Gaussian load footprint', measurement='Registered antagonist load/contact patch under fixed force'), dict(quantity='18GPa perfect-tied die', measurement='Support compliance, film and seating under matched geometry'), dict(quantity='Whole interface', measurement='Prepared/assembled microCT or surface metrology; source-tie p95 currently fails half-grid-step gate')], rigorous_linear_sensitivity_enclosure='MISSING; exact local adjoint for fixed topology/load/tie only', artifacts=dict(data_root=str(DATA), data=inventory(DATA), frozen_prediction_sha256=sha(ROOT / 'FROZEN_PREDICTIONS.json'), per_task='raw/SCORED_ROWS.json', figure='figures/crown_optimizer.png', exports=inventory(ROOT / 'exports')))
    if (ROOT / 'raw/REPLAY_VERIFICATION.json').exists():
        out['fresh_replay'] = read(ROOT / 'raw/REPLAY_VERIFICATION.json')
        out['cost']['fresh_replay_seconds'] = out['fresh_replay']['seconds']
        out['cost']['fresh_replay_peak_RSS_MiB'] = out['fresh_replay']['peak_children_RSS_MiB']
    if (ROOT / 'raw/R1_QUANTITY_ORIGINS.json').exists():
        out['primary_quantity_origins'] = read(ROOT / 'raw/R1_QUANTITY_ORIGINS.json')
    if (ROOT / 'raw/SOURCE_MODEL_DIFFERENCES.json').exists():
        out['source_model_differences'] = read(ROOT / 'raw/SOURCE_MODEL_DIFFERENCES.json')
    if (ROOT / 'raw/SUMMARY_IDENTITY_DETAIL.json').exists():
        out['summary_identity_detail'] = read(ROOT / 'raw/SUMMARY_IDENTITY_DETAIL.json')
    for debt in out['phenomenological_debts']:
        debt['resolution'] = 'PHENOMENOLOGICAL'
    out['cost']['validation_FEM_query_count'] = 'Not separately aggregated; explicit FD probes and fresh-replay complete wall/memory cost recorded'
    out['external_referent']['compared_quantity'] = 'Chen measured arithmetic mean force N; Prott reported force-law parameters and explicitly Weibull-derived arithmetic means. No optimized-crown fracture force measurement exists.'
    dump(ROOT / 'results.json', out)
    with (ROOT / 'raw/PER_TASK.csv').open('w') as f:
        keys = ['task_id', 'split', 'family', 'level', 'participant', 'L1', 'anatomy_rmse_mm', 'contact_symdiff_mm2', 'negative_gap_area_mm2', 'digital_roof_volume_mm3', 'virtual_roof_removal_quadrature_mm3', 'conditional_model_force05_N', 'calibrated_fracture_force05_N', 'quality_Pareto']
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    (fig, axs) = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ax = axs[0, 0]
    for (i, s) in enumerate(summaries):
        ax.bar(i, s['statuses'].get('PASS', 0), color=['gray', 'silver', 'teal', 'darkcyan'][i])
    ax.set_xticks(range(4), ['Morphology LP', 'X1B closing', 'Adjoint V≤V₀', 'Adjoint V≤1.1V₀'], rotation=12)
    ax.set_ylim(0, 192)
    ax.set_ylabel('PASS / 192 tasks')
    ax.set_title('Unchanged GenCAD roof predicates')
    ax = axs[0, 1]
    rr = [r for r in rows if r['participant'] == 'volume_1.0' and r['L1'] == 'PASS']
    ax.scatter([r['anatomy_rmse_mm'] for r in rr], [r['optimization']['model_force_quantile_ratio'] for r in rr], s=16, c=['teal' if r.get('quality_Pareto') else 'gray' for r in rr], alpha=0.7)
    ax.axhline(1.02, color='red', ls='--')
    ax.set_xlabel('Source-roof RMSE (mm)')
    ax.set_ylabel('Conditional model force ratio')
    ax.set_title('PER_TOOTH; uncalibrated rigid basal model')
    ax = axs[1, 0]
    xx = np.arange(3)
    ax.bar(xx, [r3['model_force_quantile_ratio'], r4['contrasts']['0']['robust_optimized_ratio'], r5['contrasts']['0']['optimized_ratio']], color=['teal', 'orange', 'purple'])
    ax.axhline(1.02, color='red', ls='--')
    ax.axhline(1, color='gray')
    ax.set_ylim(0.97, 1.08)
    ax.set_xticks(xx, ['Rigid, one load', '18GPa, node tie', '18GPa, surface tie'], rotation=12)
    ax.set_ylabel('Optimized / reference force ratio')
    ax.set_title('Whole crown: support changes the design result')
    ax = axs[1, 1]
    vals = [r['log_error'] for r in interior]
    ax.bar(np.arange(len(vals)), vals, color=['teal' if r['pass_gate'] else 'red' for r in interior])
    ax.axhline(0.2, color='red', ls='--')
    ax.set_ylabel('|log(predicted / measured force)|')
    ax.set_xlabel('Published protocol-specific interior targets')
    ax.set_title('External fracture calibration: every target retained')
    (ROOT / 'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT / 'figures/crown_optimizer.png', dpi=180)
    fig.savefig(ROOT / 'figures/crown_optimizer.pdf')
    plt.close(fig)
    table = '\n'.join((f"| {s['participant']} | {s['statuses'].get('PASS', 0)}/192 | {s['conditional_model_gain_gt2pct']} | {s['quality_Pareto']} |" for s in summaries))
    text = f"# Spatial crown optimization with an externally bounded force claim\n\nA callable elastic shape optimizer now produces frozen local crown/bridge geometries and whole STS crown specimens. The published force reference does **not** support an absolute calibrated fraction for these new geometries: {sum((r['pass_gate'] for r in interior))}/{len(interior)} protocol-specific interpolation targets pass the inherited ±0.20 log gate. The unmatched product, shape, support/contact and flaw mechanism keep generated force05 **UNKNOWN**.\n\n| Participant, PER_TOOTH restricted roof | Geometry PASS | Conditional gain >2% | Digital quality Pareto |\n|---|---:|---:|---:|\n{table}\n\nPanel: 12 metadata-hash-selected cases (8 prior test, 4 auxiliary), four crown/bridge families and all four unchanged easy/normal/hard/boundary requirements. All 192 tasks and 1,920 participant rows are retained. 61/192 candidate tasks (31.77%) abstain for locked geometric infeasibility, with no source-site dropout. The remaining 52/64 panel cases are a declared resource selection, not failed observations. Both candidate volume budgets reach identical solutions; extra allowed volume is inactive. One SLSQP iteration-limit result per budget remains visible despite geometric PASS. No global optimum or algorithm superiority is claimed.\n\nThe exact summary counterexample has identical thickness multiset and volume (identity error **0.0**, machine precision), yet model force05 differs **{100 * r1['sufficiency']['downstream_relative_difference']:.2f}%**. The minimum extension retains thickness location jointly with support and load. A mean or histogram cannot carry the mechanics edge. Adjoint/FD and equilibrium checks pass, but a rigorous enclosure over the design box is **MISSING**.\n\nWhole-crown succession, PER_TOOTH conditional model: R3 rigid support gives ratio **{r3['model_force_quantile_ratio']:.6f}** at unchanged volume and exactly unchanged intaglio nodes. Evaluating that frozen design on a compliant 18GPa die gives **{r4['contrasts']['0']['R3_frozen_ratio_on_compliant_support']:.6f}** axially. R4 two-load minimax gives **{r4['contrasts']['0']['robust_optimized_ratio']:.6f}**, failing the fixed >1.02 gate. R5 surface-barycentric tie gives **{r5['contrasts']['0']['optimized_ratio']:.6f}** for both loads, also failing. Its interface projection p95 is **{r5['interface_projection_mm']['p95']:.6f} mm**, exceeding the preregistered 0.06mm half-grid-step gate. The nearer surface tie is therefore not a certified physical port.\n\nSubstance removal: the GenCAD preparation is fixed across participants, so only the same virtual roof-envelope removal is reported per task. Whole STS virtual preparation removal is inherited and unchanged. Actual retained dentin and clinical removal are UNKNOWN. The only Pareto flag is v4's unchanged three-coordinate digital quality Pareto; physical strength/removal/CAM Pareto is UNKNOWN.\n\nAll 262 exported candidate roof shells and the complete STS specimens are closed digital meshes. Roof shells omit complete axial anatomy. Whole STS exports retain an unvalidated rigid or perfectly tied model; full wall, CAM tool/shaft, registered antagonist, physical film and manufacture are UNKNOWN. Files are explicitly research specimens, not certified manufacturable clinical crowns.\n\nControls: the original R2 antagonist injections reached the 100mm schema limit and did not exercise that predicate; their failed gate is preserved. The subsequent 37 valid-coordinate or structural faults all refuse the intended bad input. STL/3MF unit, missing-face and wrong-volume faults reject; locked lab comparison rejects wrong forces and missing support. Frozen v4 morphology and X1B controls were actually scored on the same tasks. No clinical recommendation or independent scientific acceptance.\n\nSources: Prott et al., [Table1](https://doi.org/10.4047/jap.2021.13.5.269); Chen et al., [Table2](https://doi.org/10.3390/ma17020365); full per-target locators and quantity origins are in raw/R1.json. Prott's 5% force quantiles are derived from published fitted force-law parameters, not independent measured quantiles; marginal-CI corner ranges are not joint prediction intervals. X23's own LS2 rod population was inspected but not transferred to zirconia.\n\n![Crown optimization and its force-calibration limit](figures/crown_optimizer.png)\n\nRun ./run_all.sh to verify locked inputs, regenerate all constructions in a separate replay directory, compare frozen geometry/contrasts and inject rejecting faults. Next construction: reproduce the native X1 interface classification and surface tie, quantify C3D4/C3D10 and measured support/contact discrepancy, then freeze a matched reference/optimized batch before metrology and fracture. Status PENDING_INDEPENDENT_REVIEW.\n"
    (ROOT / 'RESULTS.md').write_text(text)
    return out
if __name__ == '__main__':
    run()
