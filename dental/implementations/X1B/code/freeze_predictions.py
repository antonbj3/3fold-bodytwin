"""Freeze final operational contrast predictions before any physical test."""
import math
from scipy.stats import norm
from common import R, read, dump, sha, now

def run():
    r2 = read('raw/MATCHED_R2.json')
    r1 = read('FROZEN_PREDICTIONS_R1_DIAGNOSTIC.json')
    g = read('PREREG_ANGLE_DOUBLE_RATIO_R3.json')['metrics']
    targets = []
    for d in r2['family']:
        if d['design'] == 'D1':
            continue
        q = d['FE_conditional_ratios']['30'] / d['FE_conditional_ratios']['0']
        delta = abs(math.log(q))
        z = norm.ppf(0.975) + norm.ppf(0.8)
        n = math.ceil(4 * g['force_CV_planning'] ** 2 * z * z / delta ** 2)
        margin = max(0, delta - 0.3)
        operational_n = math.ceil(4 * g['force_CV_planning'] ** 2 * z * z / margin ** 2) if margin > 0 else None
        targets.append(dict(design=d['design'], reference_design='D1', Q_FE=q, Q_separable_control=1.0, acceptance_window_Q=[q * math.exp(-0.3), q * math.exp(0.3)], control_acceptance_window_Q=[math.exp(-0.3), math.exp(0.3)], window_kind='operational log-error tolerance; NOT a confidence/prediction interval', n_min_per_group=12, estimated_n_for_80percent_power_against_Q_equal_1=n, estimated_n_for_80percent_power_rejecting_control_tolerance=operational_n, planning_CV=0.2, planning_independent_groups=True, discriminability=dict(absolute_log_effect=delta, windows_overlap=q * math.exp(0.3) >= math.exp(-0.3), maximum_CV_for_n12_80percent_power_against_control_tolerance=margin * math.sqrt(12) / (2 * z), primary=operational_n is not None and operational_n <= 12), claim='Prospective falsifiable elastic/contact contrast. Clinical validity and fracture mechanism UNKNOWN.'))
    for d in r2['family']:
        if d['design'] == 'D1':
            continue
        r1['designs'].append(dict(design=d['design'], material='3Y', geometry_parameters={'id': d['source_id'], 't_occ': d['nominal_thickness_mm']}, gap=d['gap'], load_predictions={a: dict(point_N=f, interval_N=None, validity='UNKNOWN', target_n=12, interval_kind='No validated absolute interval available; use paired contrast only') for (a, f) in d['FE_absolute_diagnostic_N'].items()}, conditional_ratio_prediction=d['literature_conditional_prediction']))
    f = dict(schema='X1b-final-v2', frozen_utc=now(), measurement_status='NO_PHYSICAL_MEASUREMENT', lab_protocol_sha256=sha(R / 'LAB_PROTOCOL.md'), lab_template_sha256=sha(R / 'LAB_MEASUREMENTS_TEMPLATE.csv'), prereg_sha256={p.name: sha(p) for p in sorted(R.glob('PREREG_*.json'))}, code_sha256={p.name: sha(p) for p in sorted((R / 'code').glob('*.py')) if p.name not in ('prepare_package.py', 'prepare_matched.py')}, external_table_sha256=sha(R / 'raw/LITERATURE_GROUPS.json'), designs=r1['designs'], matched_family=r2['family'], double_ratio_predictions=targets, calibrated_absolute_prediction_status='UNKNOWN: R1 unsupported study transfer and R2 wide bands remain failed', note='D1/D2/D3 are legacy absolute diagnostics. D1/M1/M2 is the new same-tooth laboratory family. Sharp R3 windows are operational rejection tolerances, not calibrated statistical intervals.', do_not_refit='Never change this artifact after measurement; new information requires a new version before held-out tests')
    p = R / 'FROZEN_PREDICTIONS.json'
    if p.exists():
        old = read(p.name)
        f['frozen_utc'] = old['frozen_utc']
        if f != old:
            raise ValueError('Final frozen prediction drift')
    else:
        dump(p.name, f)
    p.with_suffix('.sha256').write_text(sha(p) + '\n')
    dump('raw/R3_DOUBLE_RATIO.json', dict(targets=targets, same_information_control='separable model computes Q=1 identically; changes to material/thickness exponent cannot change it', physical_double_ratio_observed=False))
    return targets
if __name__ == '__main__':
    print(run())
