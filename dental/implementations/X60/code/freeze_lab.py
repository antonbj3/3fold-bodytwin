from common import *
from compare_lab import compare

def run():
    freeze(ROOT / 'PREREG_LAB.json', dict(claim_type='capability', prospective=True, observable='same-batch optimized/reference arithmetic mean fracture-load ratio at 0 and30deg; common Weibull shape is a closure linking ratio to quantiles', minimum_numeric_endpoints_per_group=6, bootstrap_replicates=2000, bootstrap_seed=20261003, fixed_log_error_tolerance=0.3, reject_on_missing_or_mismatched=['material_batch', 'die_batch', 'cement_batch', 'footprint_id'], required_origin=['ceramic_bulk', 'intaglio'], absolute_force05='UNKNOWN_REQUIRES_INDIVIDUAL_FORCE_LAW_FIT_AND_VALIDATED_TRANSFER', full_cost='Physical manufacture, metrology and fracture experiments NOT PERFORMED; n6 is usability minimum, not a power guarantee'))
    freeze(ROOT / 'FROZEN_PREDICTIONS.json', dict(claim_type='capability', surface_tie_prediction_sha256=sha(ROOT / 'FROZEN_SURFACE_TIE_PREDICTIONS.json'), roof_prediction_sha256={tag: sha(ROOT / f'FROZEN_PREDICTIONS_{tag}.json') for tag in ['volume_1.0', 'volume_1.1']}, whole_prediction_sha256=sha(ROOT / 'FROZEN_WHOLE_PREDICTIONS.json'), robust_prediction_sha256=sha(ROOT / 'FROZEN_ROBUST_PREDICTIONS.json'), roof_exports_sha256=sha(ROOT / 'FROZEN_ROOF_EXPORTS.json'), lab_prereg_sha256=sha(ROOT / 'PREREG_LAB.json'), measurement='NO_PHYSICAL_MEASUREMENT', calibrated_fracture_force05_N=None, physical_manufacture='UNKNOWN', resolution='PER_POINT -> PER_TOOTH; source law POPULATION', time_scale='HANDOVER'))
    (ROOT / 'measurements_template.csv').write_text('design,angle_deg,force_N,material_batch,die_batch,cement_batch,footprint_id,fracture_origin\nreference,0,,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN\nsurface_tie_optimized,0,,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN\nreference,30,,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN\nsurface_tie_optimized,30,,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN,UNKNOWN\n')
    frozen = read(ROOT / 'FROZEN_SURFACE_TIE_PREDICTIONS.json')
    ctrl = []
    for angle in [0, 30]:
        ratio = frozen['contrasts'][str(angle)]['optimized_ratio']
        rows = []
        for (design, f) in [('reference', 1000.0), ('surface_tie_optimized', ratio * 1000.0)]:
            rows.extend((dict(design=design, angle_deg=str(angle), force_N=str(f), material_batch='fixture', die_batch='fixture', cement_batch='fixture', footprint_id='fixture', fracture_origin='ceramic_bulk') for _ in range(6)))
        good = compare(rows, frozen)
        bad = [dict(r, force_N=str(float(r['force_N']) * 1.6)) if r['design'] == 'surface_tie_optimized' else r for r in rows]
        wrong = compare(bad, frozen)
        missing = [dict(r, die_batch='UNKNOWN') for r in rows]
        unk = compare(missing, frozen)
        get = lambda rr: next((r for r in rr if r['angle_deg'] == angle))
        ctrl.append(dict(angle_deg=angle, positive_fixture=get(good)['status'] == 'WITHIN_OPERATIONAL_WINDOW', wrong_force_rejected=get(wrong)['status'] == 'REFUTED_FIXED_PREDICTION', missing_support_refused=get(unk)['status'] == 'UNKNOWN', external_referent_kind='our_own_fixture', fixture_is_not_physical_evidence=True))
    dump(ROOT / 'raw/LAB_CONTROLS.json', ctrl)
if __name__ == '__main__':
    run()
