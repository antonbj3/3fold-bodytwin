from common import *
from integrity import verify
from quality import contact_map

def run():
    verify()
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    faces = np.array([[0, 1, 2]])
    missing = contact_map(xy, faces, np.full(3, np.nan), np.full(3, np.nan))
    empty = contact_map(xy, faces, np.full(3, 0.25), np.full(3, 0.25))
    observed = contact_map(xy, faces, np.full(3, 0.0625), np.full(3, 0.25))
    checks = dict(missing_not_zero=missing['contact_symdiff_mm2'] is None and missing['negative_gap_area_mm2'] is None, observed_empty_is_zero=empty['contact_symdiff_mm2'] == 0.0 and empty['empty_contact_union'] is True, observed_area_correct=observed['contact_symdiff_mm2'] == 0.5, injected_missing_perfect_score_rejected=not missing['contact_symdiff_mm2'] == 0.0)
    dump(ROOT / 'rounds/R7.json', dict(round='R7', claim_type='capability', decision='MISSING_SUPPORT_REFUSED' if all(checks.values()) else 'FAIL', checks=checks, missing=missing, observed_empty=empty, observed_nonempty=observed, external_referent=read(ROOT / 'PREREG_R7.json')['external_referent'], scope='Absent observation differs from observed zero. Same thresholds and frozen participant outputs; retrospective repair.'))
    if not all(checks.values()):
        raise ValueError('Support gate failed')
if __name__ == '__main__':
    run()
