from common import *
import numpy as np
from screening import rates

def main():
    pairs = []
    y = np.r_[np.ones(10, int), np.zeros(90, int)]
    a = np.zeros(100, int)
    b = np.r_[np.ones(20, int), np.zeros(80, int)]
    qa = rates(y, a)
    qb = rates(y, b)
    sa = float(np.mean(y == a))
    sb = float(np.mean(y == b))
    pairs.append(dict(summary='Overall accuracy', summary_A=sa, summary_B=sb, identity_error=abs(sa - sb), state_A_confusion=qa, state_B_confusion=qb, downstream='Sensitivity', downstream_A=qa['sensitivity'], downstream_B=qb['sensitivity'], downstream_difference=abs(qa['sensitivity'] - qb['sensitivity']), minimum_extension='Class-conditioned 2x2 confusion matrix; accuracy alone cannot identify TPR at fixed Sp.'))
    y = np.r_[np.ones(10, int), np.zeros(10, int)]
    a = y.copy()
    a[:2] = 0
    b = y.copy()
    b[-2:] = 1
    qa = rates(y, a)
    qb = rates(y, b)
    sa = float(np.mean(y != a))
    sb = float(np.mean(y != b))
    pairs.append(dict(summary='Total report disagreement', summary_A=sa, summary_B=sb, identity_error=abs(sa - sb), state_A_confusion=qa, state_B_confusion=qb, downstream='Positive-report agreement', downstream_A=qa['sensitivity'], downstream_B=qb['sensitivity'], downstream_difference=abs(qa['sensitivity'] - qb['sensitivity']), minimum_extension='Separate positive and negative report-pair agreement. 7-16% total disagreement is not a class-conditioned sensitivity ceiling.'))
    a = np.array([-4.0, 0.0, 4.0, 8.0])
    b = np.array([2.0, 2.0, 2.0, 2.0])
    sa = float(np.median(a))
    sb = float(np.median(b))
    pairs.append(dict(summary='Median signed tooth margin (mm)', summary_A=sa, summary_B=sb, identity_error=abs(sa - sb), state_A=a.tolist(), state_B=b.tolist(), downstream='Any local negative sign', downstream_A=bool(np.any(a < 0)), downstream_B=bool(np.any(b < 0)), downstream_difference=1, minimum_extension='Minimum signed margin for global any-negative; per-tooth signed margin for localization. Geometry signs not externally calibrated diagnoses.'))
    assert all((r['summary_A'] == r['summary_B'] and r['identity_error'] == 0 and (r['downstream_difference'] > 0) for r in pairs))
    write(P / 'raw/SUMMARY_SUFFICIENCY.json', dict(pairs=pairs, identity_level='EXACT_BINARY64', resolution='PER_ARCH for rates; PER_TOOTH for local signed-margin states', external_referent=read(P / 'PREREG_SUMMARY_SUFFICIENCY.json')['external_referent'], rigorous_statement='Explicit finite vectors exhibit different downstream values at identical summaries; not an affine physical sensitivity approximation.', injection=dict(corrupted_identical_summary=2.01, identity_check_rejects=2.01 != 2.0)))
if __name__ == '__main__':
    main()
