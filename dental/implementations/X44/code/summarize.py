from dental_release.paths import expand as _release_expand
from common import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def run():
    st = time.perf_counter()
    reports = {n: json.loads((L / f'raw/{n}.json').read_text()) for n in ['K09', 'K13', 'K50_R1', 'K50_R2']}
    k9 = reports['K09']
    k13 = reports['K13']
    a = reports['K50_R1']
    b = reports['K50_R2']
    inv = json.loads((L / 'CHAIN_INVENTORY.json').read_text())
    q = np.load(L / 'raw/K09_PROBES.npz')
    (fig, axs) = plt.subplots(1, 3, figsize=(13, 3.8), dpi=160)
    colors = {1: '#2682a2', 2: '#d27f31', 3: '#6f8071'}
    for r in (3, 2, 1):
        m = q['region'] == r
        axs[0].scatter(q['points_mm'][m, 0], q['points_mm'][m, 1], s=3, alpha=0.7, c=colors[r], label={1: 'tooth', 2: 'PDL', 3: 'bone'}[r])
    axs[0].set(title='Complete source mandible (sampled probes)', xlabel='x [mm]', ylabel='y [mm]')
    axs[0].axis('equal')
    axs[0].legend(fontsize=8)
    source = json.loads((ROOT / 'results/F2_pdl_nonlinear/jepsen2023_fig2b.json').read_text())['points']
    x = np.array([r['deflection_mm'] for r in source])
    f = np.array([r['F_mean_N'] for r in source])
    axs[1].plot(x[:11], f[:11], 'o-', label='published loading')
    axs[1].plot(x[11:], f[11:], 's-', label='published unloading')
    xx = np.linspace(0, 0.2, 100)
    axs[1].plot(xx, a['posterior']['k_mean_N_per_mm'] * xx, '--', label='R1 linear fit')
    axs[1].set(title='Phase is measured information', xlabel='crown deflection [mm]', ylabel='force magnitude [N]')
    axs[1].legend(fontsize=8)
    vals = [a['loading_pass_count'], b['held_pass_count']]
    tot = [a['loading_count'], b['held_count']]
    axs[2].bar(['R1 loading', 'R2 phase-aware'], np.array(vals) / tot, color=['#b85140', '#2682a2'])
    axs[2].set_ylim(0, 1.3)
    axs[2].set(ylabel='held-point fraction in tolerance', title='Observation fit; material promotion refused')
    for (j, (v, n)) in enumerate(zip(vals, tot)):
        axs[2].text(j, v / n + 0.04, f'{v}/{n}', ha='center')
    axs[2].text(0.5, 1.15, 'Historical replay, not blind validation', ha='center', fontsize=8)
    fig.tight_layout()
    fig.savefig(L / 'figures/missing_chains.png')
    fig.savefig(L / 'figures/missing_chains.pdf')
    plt.close(fig)
    manifest = []
    for p in sorted(D.iterdir()):
        if p.is_file():
            manifest.append(dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)))
    sources = [ROOT / 'cells/physics/fe_reference_io.py', ROOT / 'results/F1_pdl_constitutive/runs.jsonl', ROOT / 'results/F1_pdl_constitutive/REPORT.md', ROOT / 'results/F2_pdl_nonlinear/jepsen2023_fig2b.json', ROOT / 'results/F2_pdl_nonlinear/digitize_jepsen.py', Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/shared-storage/scratch/dental_F2/lit/jepsen_fig2.jpg')), Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Open-Full-Jaw/LICENSE')), Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/bodytwin/scripts/physics_exp/viscoelastic.py'))]
    inputs = [dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)) for p in sources]
    write(L / 'SOURCE_MANIFEST.json', dict(source_inputs=inputs, data_files=manifest, source_public_dataset='https://github.com/diku-dk/Open-Full-Jaw', source_measurement='https://pmc.ncbi.nlm.nih.gov/articles/PMC9889448/', no_remote_dataset_fetched=True))
    code = [dict(path=str(p.relative_to(L)), sha256=sha(p)) for p in sorted((L / 'code').glob('*.py'))] + [dict(path='run_all.sh', sha256=sha(L / 'run_all.sh'))]
    result = dict(lane='X44-missing-chains', schema='missing-chains-runnable-v1', claim_type='capability', status='PENDING_INDEPENDENT_REVIEW', capability='Published full-mandible source region/body query, guarded regional constitutive dispatch, and protocol-preserving population observation/calibration alarm', selected_chains=['K09', 'K13', 'K50'], external_referent=dict(kind='published_dataset', locator='https://github.com/diku-dk/Open-Full-Jaw', compared_quantity='Source native material volumes/interfaces and region identity; measurement subclaims use separate primary force-deflection referents below', refutes_us=True), chains=reports, inventory=inv, graph_edges=[dict(source='Open-Full-Jaw native mesh', target='K09', quantity='body/region/SDF', unit='mm', resolution='PER_POINT', timescale='SIMULTANEOUS', status='PENDING_INDEPENDENT_REVIEW'), dict(source='K09 body/region', target='K13 law dispatch', quantity='source body and material identity', unit='dimensionless', resolution='PER_POINT', timescale='SIMULTANEOUS', status='PENDING_INDEPENDENT_REVIEW'), dict(source='Jepsen Fig2b', target='K50 observation state', quantity='group force-deflection/phase', unit='N,mm', resolution='POPULATION', timescale='HANDOVER', status='PENDING_INDEPENDENT_REVIEW'), dict(source='K50 alarm', target='K13 validity group', quantity='protocol/model-form validity', unit='dimensionless', resolution='POPULATION', timescale='HANDOVER', status='PROPOSED_SCOPE_GUARD_NOT_PATIENT_MATERIAL_EDGE')], numerical_replays_are_not_independent_evidence=True, validation_counts=dict(K09_native_elements=517202, K09_owner_probes=960, K09_material_sign_queries=2880, K09_direct_distance_queries=9, K13_invalid_parameter_calls=8, K50_R1_held_loading=6, K50_R1_unloading=11, K50_R2_held=10), rejection_summary=dict(chains_selected=3, chains_deferred=9, chain_deferred_fraction=0.75, PDL_parameter_defaults_rejected=2, PDL_parameter_defaults_candidates=3, unknown_actual_PDL_element_fraction=k13['coverage']['unknown_fraction'], R1_loading_rejected=3, R1_loading_candidates=6, R1_unloading_rejected=9, R1_unloading_candidates=11, R2_material_promotions_rejected=1, R2_material_promotions_candidates=1), cost=dict(stage_reports={k: r['cost'] for (k, r) in reports.items()}, full_cold_cost_includes_inherited_FE='F1 source runs contain seconds for all preserved fit/solve/refinement attempts; none are new observations. No speedup claim. Historical human acquisition/model discovery cost UNKNOWN.', inherited_F1_logged_solver_seconds=sum((float(json.loads(x).get('seconds', 0)) for x in (ROOT / 'results/F1_pdl_constitutive/runs.jsonl').read_text().splitlines())), scientific_acquisition_cost='UNKNOWN; no new lab measurement', preparation='Source extraction and hashes included in K09 stage; report/source searches unmeasured separately', questions=0, fallback='UNKNOWN_NOT_ACQUIRED_SIGNED_PAIRED_TOOTH_OR_MANUFACTURING_DATA'), data_manifest=manifest, input_manifest=inputs, code_manifest=code, full_physical_chain_closed=False, limitations=['K09 exact source geometry of one full mandible; PDL generated in source pipeline, not patient histology', 'K13 whole-tooth/bone are scenario closures; upper-incisor population E conditional; actual mandibular PDL remains unknown', 'K50 observation phase branches are group curves, not fitted work-conjugate constitutive tissue law', 'Source split historically visible; not blinded independent validation', 'No physical fabrication, patient load, clinical advice, source graph or external checkout edits'], next_construction='A signed same-tooth, time/phase-resolved force-displacement port with initial prestress and device channels, independent geometry plus axial/lateral directions; freeze before acquisition')
    write(L / 'results.json', result)
    write(L / 'GRAPH_COVERAGE_PROPOSAL.json', dict(review_state='PENDING_INDEPENDENT_REVIEW', existing_target='DENT-MAT-PDL-CONSTITUTIVE', missing_scoped_ports=[dict(chain='K09', proposal='DENT-FIELD-NATIVE-MULTIREGION-OWNER', resolution='PER_POINT', timescale='SIMULTANEOUS', scope=_release_expand('native @DENTAL_CASE_ID@ computational source geometry, no edited-field claim')), dict(chain='K50', proposal='DENT-OBS-PHASED-FORCE-INVALIDATION', resolution='POPULATION', timescale='HANDOVER', scope='published incisor mean protocol, material promotion guard only')]))
    state('READY_FOR_INDEPENDENT_REVIEW', dict(K09=k9['outcome'], K13=k13['outcome'], K50_R1=a['outcome'], K50_R2=b['outcome']), 'Acquire signed same-tooth data; preserve source native partition; no global clinical/material promotion')
    print(json.dumps(dict(selected=['K09', 'K13', 'K50'], data_bytes=sum((x['bytes'] for x in manifest)), physical_chain_closed=False, results_sha256=sha(L / 'results.json'))))
if __name__ == '__main__':
    run()
