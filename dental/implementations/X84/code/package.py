from dental_release.paths import expand as _release_expand
from pathlib import Path
import hashlib
P = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, v):
    Path(p).write_text(json.dumps(v, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
import json, datetime, resource, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    tick = time.perf_counter()
    r = [json.load(open(P / f'rounds/R{i}/results.json')) for i in [1, 2, 3, 4]]
    (g, c, m, exact) = r
    (fig, ax) = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    pts = []
    for p in g['pairs']:
        u = np.array(p['upper_point_cbct_mm'])
        l = np.array(p['lower_point_cbct_mm'])
        ax[0, 0].plot([u[0], l[0]], [u[2], l[2]], color='0.7', lw=0.6)
        pts.extend([u, l])
        ax[0, 0].scatter(u[0], u[2], c='#3b76b3', s=15)
        ax[0, 0].scatter(l[0], l[2], c='#d78835', s=15)
    ax[0, 0].set(xlabel='CBCT x (conditional mm)', ylabel='CBCT z (conditional mm)', title='Demo1 candidate pose: closed bite UNKNOWN\nMinimum projected gap 4.229 mm')
    h = [e['height_edit_mm'] * 1000 for e in g['edits']]
    gap = [e['target_gap_after_virtual_first_touch_mm'] for e in g['edits']]
    ax[0, 1].plot(h, gap, 'o-', c='#343d54')
    ax[0, 1].axhline(0, color='0.4')
    ax[0, 1].set(xlabel='Virtual FDI36 height edit (micrometres)', ylabel='Target gap after virtual first touch (mm)', title='Demo1 axial scenario: contact still absent\nPhysical height decision UNKNOWN')
    x = [v['applied_N'] for v in c['rows']]
    y = [v['corrected_N'] for v in c['rows']]
    ax[1, 0].scatter(x, y, c=['#c54141' if not v['pass_gate'] else '#35954e' for v in c['rows']], s=22)
    xx = np.array([0, 400])
    ax[1, 0].plot(xx, xx, 'k-', lw=0.7)
    ax[1, 0].fill_between(xx, 0.9 * xx, 1.1 * xx, color='gray', alpha=0.15)
    ax[1, 0].set(xlabel='Externally applied load (N)', ylabel='Published gain corrected readout (N)', title='Hattori1994 Table1: 27/48 outside10%\nOther subject and historical sensor')
    vals = [v['reported_force_N'] for v in m['external_rows']]
    ax[1, 1].plot(range(1, 7), vals, 'o-', c='#3b76b3')
    ax[1, 1].axhline(100, color='#c54141', ls='--', label='Fixed research load100N')
    ax[1, 1].set(xlabel='Clench acquisition', ylabel='Reported FDI36 force (N)', title='Hattori1994 Table2: 96–177N\nReference subject; Demo1 force UNKNOWN')
    ax[1, 1].legend()
    fig.suptitle('X84: same-patient geometry is usable; patient force still needs a measurement', fontsize=13)
    fig.savefig(P / 'figures/patient_load.png', dpi=150)
    plt.close(fig)
    result = dict(schema='X84-patient-bound-load-v1', claim_type='capability', lane='X84-patient-bound-load', patient_id=g['patient_id'], review_state='PENDING_INDEPENDENT_REVIEW', model='UNKNOWN_PHYSICAL; FIXED_AXIAL_PROJECTION + OBSERVATION_BOUND_STATIC_PORT', outcome='SAME_PATIENT_CONTACT_WITNESSES_AND_STATIC_MEASUREMENT_CONSUMER_EXECUTABLE; K16_PHYSICAL_LOAD_STILL_MISSING', became_possible='27 source-tooth surfaces yield27 continuous axial pair witnesses in patient CBCT frame; static calibrated contact intervals can be consumed without fitting support stiffness; no Demo1 force or clinical height decision established', external_referent=m['external_referent'], external_referents=[g.get('external_referent', dict(kind='published_dataset', locator='https://doi.org/10.5281/zenodo.8027553', compared_quantity='Demo1 released source-paired IOS crown coordinates, not contact/force truth', refutes_us=True)), c['external_referent'], m['external_referent']], key_results=dict(candidate_minimum_gap=dict(value=g['source_pose_minimum_gap_mm'], unit='mm conditional', resolution='PER_POINT'), virtual_target36_gap=dict(value=g['target36_first_touch_gap_mm'], unit='mm conditional', resolution='PER_TOOTH'), target_height_edits=dict(value=g['edits'], unit='mm conditional', resolution='PER_TOOTH'), patient_force=dict(value=[0, None], unit='N', upper_endpoint='UNBOUNDED_WITHOUT_TOTAL', resolution='PER_TOOTH', status='UNKNOWN'), historical_gain_failures=dict(value=c['n_fail'], denominator=len(c['rows']), resolution='PER_TOOTH', relative_tolerance=0.1), reference_subject_FDI36=dict(value=m['external_force36_observed_range_N'], unit='N printed', resolution='PER_TOOTH', sensor_error='UNKNOWN'), fixed100N_display_exceedances=dict(value=m['fixed100N_display_exceedances'], denominator=6, resolution='PER_TOOTH')), rounds=[dict(round='R' + str(i), claim_type=v['claim_type'], result_file=str(P / f'rounds/R{i}/results.json'), sha256=sha(P / f'rounds/R{i}/results.json')) for (i, v) in zip([1, 2, 3, 4], r)], sufficiency=[exact['force_summary_test'], exact['wrench_summary_test']], balance_only_predecessor_sufficiency=[g['sufficiency'], m['sufficiency']], attrition=dict(geometry_source_triangles=g['source_triangles'], degenerate_projection_rejections=sum((v['rejected_degenerate_axial_projection'] for v in g['dropout'])), projected_pairs_retained=len(g['pairs']), possible_tooth_pair_count=g['all_upper_lower_pairs'], no_projection_overlap_rejections=g['all_upper_lower_pairs'] - len(g['pairs']), pair_dropout_fraction=g['pair_projection_dropout_fraction'], reference_force_rows=m['attrition']), reference_benchmark_design='Retrospective after reading source tables; published gain was fitted on the same Table1. Residual check refutes a uniform bound; it is not heldout calibration validation.', negative_results=['No acquired closed bite in Demo1;4.229mm candidate separation is a geometric observation, not physical pose truth', 'Uniform historical gain outside10% on27/48 reference loads', 'All84 published tooth force entries rejected for Demo1 identity binding', 'Patient force upper endpoint unbounded without measured total', 'Same total gives73N different FDI36 force in positive diagonal contact-compliance model; same total and FDI36 force gives102.203Nmm moment difference on real geometry', 'Printed force zeros are censored observations, not certified zero forces', 'No published-to-Demo1 same-subject force identity found; source search is not proof of global absence'], controls=dict(R1_all_pass=g['controls']['all_pass'], R2_injected_errors_all_rejected=c['all_injected_errors_rejected'], R3_all_pass=m['all_control_gates_pass'], R3_fault_count=len(m['faults']), R4_exact_equilibrium_all_pass=exact['all_pass'], same_information_comparison='LP/direct inverse parity; no algorithm advantage claimed'), enclosure='Conditional affine digital geometry and analytic box+sum force bounds; formal floating-point enclosure MISSING. Physical pose/scale/sensor bounds UNKNOWN. No unbounded linear sensitivity is presented as a physical certificate.', cost=dict(rounds=[v['cost'] for v in r], documentation_figure_wall_s=time.perf_counter() - tick, discovery_and_reading='UNKNOWN', upstream_fit='UNKNOWN_REUSED', physical_measurement='NOT_RUN', RAM_peak_MiB=max((v['cost']['maxrss_MiB'] for v in r)), threads_max=4, GPU=False, large_arrays=[], intermediate_limit_bytes=3000000000), missing_measurement_file=str(P / 'NEXT_MEASUREMENT.json'), edge_contracts=[dict(producer='X69 acquired source-paired tooth geometry', consumer='X77/X78 FDI36 crown static load boundary', resolution='PER_POINT', timescale='SIMULTANEOUS', status='CANDIDATE_COORDINATE_LINK; physical loaded bite UNKNOWN'), dict(producer='Prospective Demo1 normal contact-force observation', consumer='FDI36 static force query and crown load port', resolution='PER_POINT', timescale='SIMULTANEOUS', status='MISSING_MEASUREMENT'), dict(producer='Prospective Demo1 height-response measurement', consumer='X82 height-force interval operation on this patient', resolution='PER_TOOTH', timescale='SIMULTANEOUS', status=_release_expand('NOT_RUN; historical @DENTAL_CASE_ID@ data not transferred'))], source_manifest=dict(file=str(P / 'raw/SOURCE_MANIFEST.json'), sha256=sha(P / 'raw/SOURCE_MANIFEST.json')), claims_outside_scope=['patient-specific lifetime', 'sensor-free force inferred from inserted sensor', 'clinical crown-height prescription', 'tangential/full3D force', 'measured pressure field on Demo1'])
    write(P / 'results.json', result)
    if not (P / 'FROZEN_PREDICTIONS_R3.json').exists():
        fp = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(P / 'PREREG_R3.json'), patient_id=g['patient_id'], static_query='UNKNOWN until matched calibration+loaded bite acquisition', height_edits_mm=[-0.05, 0, 0.05], fixed_research_load_N=100, physical_acquisition='NOT_RUN', reference_benchmark_already_read=True, scope='Prospective Demo1 prediction frozen before any physical acquisition; does not claim to predate historical published measurement')
        write(P / 'FROZEN_PREDICTIONS_R3.json', fp)
        (P / 'FROZEN_PREDICTIONS_R3.json.sha256').write_text(sha(P / 'FROZEN_PREDICTIONS_R3.json') + '\n')
    print(result['outcome'])
if __name__ == '__main__':
    main()
