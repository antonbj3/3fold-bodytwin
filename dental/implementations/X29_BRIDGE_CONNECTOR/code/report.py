from dental_release.paths import expand as _release_expand
import os, collections, time
import numpy as np
os.environ.setdefault('MPLCONFIGDIR', str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'raw/mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import R, read, dump, sha, now, state

def run():
    r1 = read('raw/R1_FIT.json')
    r2 = read('raw/R2_COMPARISON.json')
    r3 = read('raw/R3_COMPARISON.json')
    fe = read('raw/BRIDGE_FE.json')
    r4 = read('raw/R4_PAIRED.json')
    r5 = read('raw/R5_HAZARD.json')
    v = read('VERIFICATION.json')
    (fig, axs) = plt.subplots(2, 2, figsize=(12, 8.2), layout='constrained')
    ax = axs[0, 0]
    colors = {'3Y': '#2563a6', '4Y': '#e68a19', '5Y': '#348b60'}
    for (i, g) in enumerate(r2['groups']):
        (lo, hi) = g['mean_CI95_N']
        ax.errorbar(g['mean_N'], i, xerr=[[g['mean_N'] - lo], [hi - g['mean_N']]], fmt='o', color=colors[g['material']], capsize=4)
    ax.set_yticks(range(6), [f"{g['material']}  {('4×2.25' if i % 2 == 0 else '3×3')} mm" for (i, g) in enumerate(r2['groups'])])
    ax.invert_yaxis()
    ax.set_ylim(6.2, -0.6)
    ax.axvline(1500, color='#a72929', ls='--')
    ax.set(xlabel='System fracture load N; 95% CI of group mean', title='All 9 mm²: unaged end-supported KATANA')
    ax.text(0.02, 0.03, 'Source: PMC10471503 Table 2; n=6/group', transform=ax.transAxes, fontsize=9)
    ax = axs[0, 1]
    for (i, g) in enumerate(r3['groups']):
        (lo, hi) = g['mean_CI95_N']
        ax.errorbar(g['mean_N'], i, xerr=[[g['mean_N'] - lo], [hi - g['mean_N']]], fmt='o', color='#7044a0', capsize=4)
    ax.set_yticks([0, 1], ['Axial aged group', '30° aged/test group'])
    ax.set_ylim(1.6, -0.6)
    ax.axvline(1500, color='#a72929', ls='--')
    ax.set(xlabel='System fracture load N; 95% CI of group mean', title='Aged 3Y cantilever: 21.1 mm² (minimum 12)')
    ax.text(0.02, 0.08, 'F30/F0 = 0.390 [0.323, 0.470]\nSource: PMC10788321 Tables 1–2; n=8/group', transform=ax.transAxes, fontsize=9)
    ax = axs[1, 0]
    base = R / 'raw/fe/h4_b2.25_s0.18'
    g = np.load(base / 'geometry.npz')
    s = np.load(base / 'axial_stress.npz')['s1']
    cf = g['cf']
    sel = np.arange(len(cf))[::5]
    im = ax.scatter(cf[sel, 0], cf[sel, 2], c=np.minimum(np.maximum(s[sel], 0), 0.4), s=4, cmap='magma', rasterized=True)
    ax.axvline(-5.5, color='gray', ls=':')
    ax.axvline(5.5, color='gray', ls=':')
    ax.set(xlabel='x mm', ylabel='z mm', title='Fresh three-unit elastic field; axial1 N, surface samples')
    fig.colorbar(im, ax=ax, label='Tensile stress MPa/N (unvalidated absolute values)')
    ax = axs[1, 1]
    for (name, c) in [('Connector', '#2563a6'), ('Contact', '#e68a19')]:
        key = 'q_' + name.lower()
        for (i, lc) in enumerate(['axial', 'offaxis30']):
            vals = [x[key] for x in r4['contrasts'] if x['load_case'] == lc]
            x0 = i + (-0.1 if name == 'Connector' else 0.1)
            ax.plot([x0, x0], [min(vals), max(vals)], color=c, lw=4, label=name if i == 0 else None)
            ax.scatter([x0], [vals[-1]], color=c)
    ax.axhline(1, color='gray', ls='--')
    ax.set_xticks([0, 1], ['Axial', '30°'])
    ax.set(ylim=(0.8, 1.85), ylabel='Conditional force ratio:16 / 9 mm²', title='Area 9 → 16 mm²: mode-specific proxy ratios')
    ax.legend()
    ax.text(0.02, 0.03, 'Lines: two numerical meshes, not confidence intervals\nOverall paired ratio5% gate FAIL; no calibrated force prediction', transform=ax.transAxes, fontsize=8)
    fig.suptitle('Published fracture loads and conditional bridge-response queries', fontsize=13)
    fig.savefig(R / 'bridge_connector.png', dpi=180)
    fig.savefig(R / 'bridge_connector.pdf')
    plt.close(fig)
    screening = read('raw/SCREENING.json')
    relevant = [s for s in screening if s['relevant']]
    reasons = collections.Counter((s['reason'] for s in screening if not s['retained_system']))
    allruns = fe['runs'] + r4['new_runs']
    artifacts = []
    for p in sorted(list((R / 'raw/fe').glob('*/*')) + list((R / 'raw/replay').glob('*'))):
        if p.is_file() and p.suffix in ['.npz', '.frd', '.inp', '.stl', '.json']:
            artifacts.append(dict(path=str(p.resolve()), sha256=sha(p), bytes=p.stat().st_size))
    source2 = R / 'raw/sources/PMC10788321.xml'
    out = dict(lane='X29-bridge-connector', created_utc=now(), claim_type='information_link', achieved='Two published studies now answer setup-conditioned laboratory mean margins; area alone cannot determine fracture mode or load-direction/history capacity.', overall_status='INFORMATION_LINK_DEMONSTRATED; CONNECTOR_DIMENSIONING_PHYSICALLY_UNVALIDATED', external_referent=r3['external_referent'], external_referents=[r2['external_referent'], r3['external_referent']], source2=dict(path=str(source2), sha256=sha(source2), license='CC BY4.0; verify source permissions retained'), per_material=dict(z3Y='KATANA9mm2 unaged mean margins exceed1500N; LavaPlus21.1mm2 aged cantilever fails1500N even axial. Product/regime-specific, no transferable capacity law.', z4Y='KATANA9mm2 4x2.25 lower95%mean CI1572.8N>1500N, below16mm2 rule; excess for this benchchallenge only.3x3 CI crosses1500N.', z5Y='9mm2 meanCI crosses1500N; UTML molarspan outside premolar-only guide. Excess/deficit unresolved.', lithium_disilicate='UNKNOWN_NO_LOCAL_TARGET_MEASUREMENT'), R1=r1, R2=r2, R3=r3, R4=r4, R5=r5, gates=dict(source_controls=v['all_normal_pass'], injected_faults=v['all_faults_rejected'], area_study_transfer=False, bridge_surface_geometry=all((r['surface_p99_error_mm'] <= 0.06 for r in allruns)), absolute_stress_refinement=fe['refinement_gate_pass'], paired_ratio=r4['paired_mesh_gate_pass'], hazard_ratio=r5['all_mesh_gate_pass']), exclusions=dict(connector_query_sources=len(screening), dental_fracture_relevant=len(relevant), R1_retained_system_studies=1, R1_rejected_sources=len(screening) - 1, R1_rejected_fraction=(len(screening) - 1) / len(screening), R1_relevant_rejected=len(relevant) - 1, R1_relevant_rejected_fraction=(len(relevant) - 1) / len(relevant), uncensored_connector_studies=0, uncensored_connector_group_exclusion_fraction=1.0, reasons=dict(reasons), R3_separately_retained_cantilever_studies=1, R3_scope='diagnostic information link; not pooled into end-supported regression'), cost=dict(fresh_initial_load_cases=12, fresh_initial_GEOMETRY_solve_seconds=sum((r['seconds'] for r in allruns)), max_solver_rss_bytes=max((r['max_solver_rss_bytes'] for r in allruns)), threads=4, prep_fit_validation_discovery='See COMMANDS.md. Source scan includes exclusions and two failed solver attempts; no paidagent/GPU/download.', physical_tests_performed=0, query_cost='six source margins; one matched directional contrast; cachedsurfaceintegrals; no fitted strength'), data_directory=_release_expand('@DENTAL_WORK_ROOT@/X29-bridge-connector'), large_artifacts=artifacts, frozen_predictions_sha256=sha(R / 'FROZEN_PREDICTIONS.json'), quantity_resolution=dict(published_force_means_SD_CI_and_contrasts='POPULATION; individual endpoints PER_TOOTH unavailable', connector_area_height_width='PER_SURFACE_REGION', stress_tensor_field='PER_POINT', regionalstresspercentiles_hazard_and_ratios='PER_SURFACE_REGION', applied1000N_required1500N_and_unknownWeibullm='PHENOMENOLOGICAL', surface_exportdeviation='PER_SURFACE_REGION'), phenomenological_debts=['1500N/1000N laboratory challenge replaces no patientload', 'unmeasuredsamebatch mode strengths', 'unknownWeibullm/sigma0', 'contactfracture law', 'cementtie/supportcompliance', 'asbuiltconnectorradius'], graph_status='PENDING_INDEPENDENT_REVIEW; missing appropriate connector fracture target; rank/packet STALE_INPUT; no graph mutation')
    if (R / 'raw/FRESH_REPLAY.json').exists():
        out['fresh_replay'] = read('raw/FRESH_REPLAY.json')
    if (R / 'COST_LEDGER.json').exists():
        out['cost']['detailed_ledger'] = read('COST_LEDGER.json')
    dump('results.json', out)
    return out
if __name__ == '__main__':
    run()
