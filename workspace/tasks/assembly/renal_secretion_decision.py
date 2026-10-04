"""Decide whether a creatinine rise means the kidney got worse, and quantify what was blocked.

THE DECISION. A drug raises serum creatinine. Standard practice reads that as a fall in glomerular
filtration and may reduce or stop a dose. The alternative is that filtration is untouched and only
the tubular SECRETION of creatinine is blocked, which changes the marker and not the organ. Those
two readings differ in what should happen to the patient, so separating them is a decision, and it
is a decision outside the eye, made from quantities that are measured in the clinic.

THE EXTERNAL FACIT, one study, both arms in the same subjects (PMID 32989831,
DOI 10.1002/jcph.1750, "Tucatinib Inhibits Renal Transporters OCT2 and MATE Without Impacting Renal
Function in Healthy Subjects"), carried in results/LANE_READ_CREATIVE/CONSUMABLE.json as READ-A001
through READ-A013:

    iohexol clearance, a filtration marker that is not secreted   129.8 -> 129.0 mL/min
    metformin renal clearance, mostly secreted                    30.0 -> 17.6 L/h
    in-vitro IC50, metformin as substrate   OCT2 14.7, MATE1 0.34, MATE2-K 0.135 uM
    in-vitro IC50, creatinine as substrate  OCT2 0.107, MATE1 0.0855 uM

WHAT THIS FILE COMPUTES RATHER THAN RESTATES. The published ratio 0.57 is a ratio of two clearances
and says nothing about which process moved. Splitting each clearance into filtration and secretion
turns the same three measured numbers into the blocked fraction of secretion, and that fraction then
implies, through the in-vitro IC50s, the unbound inhibitor concentration each transporter would need
in order to be the route responsible. Those implied concentrations differ by two orders of magnitude
between MATE and OCT2, so the route is decidable from numbers that are already in hand.

THE EQUALLY INFORMED CONTROL is creatinine-based estimated GFR. It sees the same patient and the
same creatinine and has no transporter information; it is what practice uses. The facit says it is
wrong here, because iohexol moved by 0.6 percent.

ASSUMPTION THAT CARRIES THE RESULT, stated because it is the first thing to attack: metformin is
taken as negligibly protein bound and not reabsorbed, so its filtration clearance equals GFR. If
plasma protein binding were material, the filtration term would be smaller and the blocked secretion
fraction would be LARGER, not smaller, so the direction of the verdict survives the assumption
failing. PENDING_INDEPENDENT_REVIEW. No clinical recommendation is implied.
"""
from __future__ import annotations

L_PER_H_TO_ML_PER_MIN = 1000.0 / 60.0

GFR_BASELINE_ML_MIN = 129.8          # READ-A002, iohexol
GFR_ON_DRUG_ML_MIN = 129.0           # READ-A003, iohexol
CL_METFORMIN_BASELINE_L_H = 30.0     # READ-A009
CL_METFORMIN_ON_DRUG_L_H = 17.6      # READ-A009
IC50_METFORMIN_UM = {'OCT2': 14.7, 'MATE1': 0.34, 'MATE2-K': 0.135}
IC50_CREATININE_UM = {'OCT2': 0.107, 'MATE1': 0.0855}
# The published geometric least-squares-mean ratio and its 90 percent interval, which the first
# version of this file did not use: it worked from the arithmetic means alone and therefore gave the
# blocked fraction without an interval. A swarm job (BT-FW48-AUTO-99bbf4a78c1398) quoted the interval
# from the same DOI while doing something else, and it is the input that turns the verdict below into
# a range instead of a point. Propagating it is the honest version of the same decision.
LSM_RATIO = 0.571
LSM_RATIO_90CI = (0.516, 0.632)


def split_renal_clearance(cl_total_ml_min: float, gfr_ml_min: float,
                          unbound_fraction: float = 1.0) -> dict:
    """Filtration is GFR times the unbound fraction; whatever is left over is secretion."""
    filtration = gfr_ml_min * unbound_fraction
    return {'filtration_ml_min': filtration,
            'secretion_ml_min': cl_total_ml_min - filtration,
            'secreted_share': (cl_total_ml_min - filtration) / cl_total_ml_min}


def implied_unbound_inhibitor_uM(blocked_fraction: float, ic50_uM: float) -> float:
    """Invert one-site inhibition: a blocked fraction f needs I = IC50 * f / (1 - f)."""
    if not 0.0 < blocked_fraction < 1.0:
        raise ValueError('blocked fraction must lie strictly between 0 and 1')
    return ic50_uM * blocked_fraction / (1.0 - blocked_fraction)


def blocked_fraction_at(i_uM: float, ic50_uM: float) -> float:
    return i_uM / (i_uM + ic50_uM)


def blocked_fraction_from_ratio(ratio: float) -> float:
    """Blocked secretion implied by a total-clearance ratio, with filtration held at its measured
    values in each arm. Filtration is measured separately, so only the secretion term moves."""
    cl_base = CL_METFORMIN_BASELINE_L_H * L_PER_H_TO_ML_PER_MIN
    sec_base = cl_base - GFR_BASELINE_ML_MIN
    sec_drug = ratio * cl_base - GFR_ON_DRUG_ML_MIN
    return 1.0 - sec_drug / sec_base


def main() -> None:
    base = split_renal_clearance(CL_METFORMIN_BASELINE_L_H * L_PER_H_TO_ML_PER_MIN,
                                 GFR_BASELINE_ML_MIN)
    drug = split_renal_clearance(CL_METFORMIN_ON_DRUG_L_H * L_PER_H_TO_ML_PER_MIN,
                                 GFR_ON_DRUG_ML_MIN)
    blocked = 1.0 - drug['secretion_ml_min'] / base['secretion_ml_min']
    filt_change = GFR_ON_DRUG_ML_MIN / GFR_BASELINE_ML_MIN - 1.0
    ratio = CL_METFORMIN_ON_DRUG_L_H / CL_METFORMIN_BASELINE_L_H

    print('DECISION: is the kidney worse, or is only the marker blocked?')
    print(f'  filtration, measured by a marker that is not secreted : {filt_change * 100:+.2f} %')
    print(f'  total metformin renal clearance                       : {(ratio - 1) * 100:+.2f} %')
    print(f'  secreted share of metformin clearance at baseline     : {base["secreted_share"]:.4f}')
    print(f'  tubular secretion of metformin, blocked fraction      : {blocked:.4f}')
    print('  verdict: secretion is blocked, filtration is not -> a creatinine rise under this drug')
    print('           is a marker artefact unless filtration is measured independently.')

    print('\nWHICH ROUTE carries it: unbound inhibitor concentration each transporter would need')
    for name, ic50 in IC50_METFORMIN_UM.items():
        print(f'  {name:<8} IC50 {ic50:>6.3f} uM -> needs I = {implied_unbound_inhibitor_uM(blocked, ic50):.4f} uM')
    i_mate2k = implied_unbound_inhibitor_uM(blocked, IC50_METFORMIN_UM['MATE2-K'])
    spread = implied_unbound_inhibitor_uM(blocked, IC50_METFORMIN_UM['OCT2']) / i_mate2k
    print(f'  the OCT2 route demands {spread:.0f} times the concentration the MATE2-K route does,')
    print('  so the two hypotheses are separated by the numbers already measured, not by opinion.')

    print('\nCONSEQUENCE for creatinine, at the concentration the MATE2-K route implies')
    for name, ic50 in IC50_CREATININE_UM.items():
        f_cr = blocked_fraction_at(i_mate2k, ic50)
        print(f'  {name:<8} creatinine IC50 {ic50:.4f} uM -> secretion blocked {f_cr * 100:.1f} %'
              f'  ({f_cr / blocked:.2f} times the block on metformin)')
    print('  creatinine is blocked HARDER than metformin, which is why the marker moves while')
    print('  filtration does not. The direction is a prediction of the same three measured numbers.')

    print('\nCONTROL, equally informed: creatinine-based eGFR sees the same patient and the same')
    print('  creatinine and has no transporter term, so it reports impairment. The facit is the')
    print(f'  iohexol clearance, which moved {filt_change * 100:+.2f} %.')

    print('\nTHE VERDICT AS A RANGE, from the published 90 percent interval on the clearance ratio')
    lo_f = blocked_fraction_from_ratio(LSM_RATIO_90CI[1])   # a HIGHER ratio means LESS blocked
    hi_f = blocked_fraction_from_ratio(LSM_RATIO_90CI[0])
    print(f'  blocked secretion fraction : {lo_f:.4f} to {hi_f:.4f}, '
          f'point estimate {blocked_fraction_from_ratio(LSM_RATIO):.4f}')
    print(f'    the arithmetic-mean figure above, {blocked:.4f}, sits inside that interval')
    i_lo = implied_unbound_inhibitor_uM(lo_f, IC50_METFORMIN_UM['MATE2-K'])
    i_hi = implied_unbound_inhibitor_uM(hi_f, IC50_METFORMIN_UM['MATE2-K'])
    print(f'  implied unbound concentration on the MATE2-K route : {i_lo:.4f} to {i_hi:.4f} uM')
    print(f'  the OCT2 route needs {implied_unbound_inhibitor_uM(hi_f, IC50_METFORMIN_UM["OCT2"]):.1f} uM '
          f'even at the interval end most favourable to it, so the route verdict survives the interval.')

    print('\nFALSIFIER, stated before the next measurement: if the measured unbound plasma')
    print(f'  concentration of the inhibitor is below {i_mate2k:.4f} uM, then not even the most')
    print('  sensitive transporter in this set can account for the blocked fraction, and a route')
    print('  outside OCT2/MATE1/MATE2-K is required. That single measurement decides it.')

    import json, os
    out = {
        'filtration_change_fraction': filt_change,
        'total_clearance_ratio': ratio,
        'baseline_secreted_share': base['secreted_share'],
        'secretion_blocked_fraction': blocked,
        'implied_unbound_inhibitor_uM': {k: implied_unbound_inhibitor_uM(blocked, v)
                                         for k, v in IC50_METFORMIN_UM.items()},
        'creatinine_block_at_mate2k_concentration': {k: blocked_fraction_at(i_mate2k, v)
                                                     for k, v in IC50_CREATININE_UM.items()},
        'external_facit': {'pmid': '32989831', 'doi': '10.1002/jcph.1750'},
        'published_lsm_ratio': LSM_RATIO,
        'published_lsm_ratio_90ci': list(LSM_RATIO_90CI),
        'blocked_fraction_90ci': [blocked_fraction_from_ratio(LSM_RATIO_90CI[1]),
                                  blocked_fraction_from_ratio(LSM_RATIO_90CI[0])],
        'blocked_fraction_point_from_lsm': blocked_fraction_from_ratio(LSM_RATIO),
        'implied_unbound_inhibitor_uM_mate2k_90ci': [
            implied_unbound_inhibitor_uM(blocked_fraction_from_ratio(LSM_RATIO_90CI[1]),
                                         IC50_METFORMIN_UM['MATE2-K']),
            implied_unbound_inhibitor_uM(blocked_fraction_from_ratio(LSM_RATIO_90CI[0]),
                                         IC50_METFORMIN_UM['MATE2-K'])],
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    d = 'results/ASSEMBLY_RENAL_SECRETION'
    os.makedirs(d, exist_ok=True)
    with open(d + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('\nwrote', d + '/decision.json')


if __name__ == '__main__':
    main()
