"""Decide whether an oedema is a permeability lesion or a driving-force change, and when you cannot.

THE DECISION, outside the eye. A tissue swells and its lymph carries protein. Two accounts follow,
and they imply different interventions: the barrier's solute reflection coefficient sigma fell (a
permeability lesion), or the oncotic driving force fell with sigma intact (hypoproteinaemia, a
haemodynamic cause). The clinic reads the lymph-to-plasma protein ratio as if it reported sigma
directly. This file shows that reading is a ONE-SIDED bound, states which side, and gives the single
extra measurement that closes it.

WHAT THE SOURCE REPORT ESTABLISHED, and what it left open. `results/BT-CAP-Q020/RESULTS.md` section 4
proves the gate: at a single steady state with one independent variation there are three unknowns,
K_f, L_p (as a permeability-surface-area product) and sigma, against two observables, and sigma enters
only through the PRODUCTS (1-sigma)L_p/K_f. There is a flat ridge in sigma-space, so an inversion that
claims to "find the parameters" has infinitely many equally correct answers and choosing among them is
an opinion. The report names three ways identifiability could be recovered and stops there. This file
takes the cheapest of them and turns it into a decision with a closed form.

THE STEADY-STATE RELATION, written out because the whole file is one line of algebra. Water and
protein across one barrier, volume steady so lymph flow Q equals filtration:

    Q * C_i = (1 - sigma) * Q * C_p + PS * (C_p - C_i)

Divide by C_p and let r = C_i/C_p and Pe = Q/PS:

    r = (1 - sigma) + sigma / (1 + Pe)          so        1 - r = sigma * Pe / (1 + Pe)

Two consequences, both exact and both parameter-free:

 1. 1 - r UNDERSTATES sigma by the factor Pe/(1+Pe), always. So a single steady observation can only
    ever EXCLUDE a permeability lesion -- when 1-r already exceeds the baseline sigma -- and can never
    establish one. That asymmetry is the decision this file returns, and it is the opposite of how the
    ratio is usually read.
 2. r depends on flow only through Pe, so measuring r at two lymph flows gives two equations in two
    unknowns and sigma comes out in closed form, with no fitted parameter and no prior:

        PS = (r1 - r2) / ((1-r1)/Q1 - (1-r2)/Q2)        sigma = (1-r1) * (PS + Q1) / Q1

    That is the lymph-protein washdown as an inversion rather than as a plot.

MEASURED NUMBERS WITH UNITS AND LOCATORS. Only two quantities below are measured; everything else is
either exact algebra or a declared fixture, and they are labelled as such.
    plasma protein lowered by 33 % and by 68 %, with a linear volume-pressure response whose slope
        gives the PRODUCT K_f * sigma and not sigma alone
        Manning RD, Guyton AC. Am J Physiol Heart Circ Physiol 1983;245(2):H284-H293.
        PMID 6881362, DOI 10.1152/ajpheart.1983.245.2.h284   (locator: abstract, protein reduction arms)
    the fibre-matrix ultrafilter at the luminal entrance to the interendothelial clefts, which is why
        sigma is a property of a structure and not a fitting constant, and the single-perfused-capillary
        method that measures the water and protein fluxes simultaneously -- route 2 of the report's
        three, and the route this file's two-flow inversion substitutes for
        Michel CC, Curry FE. Physiol Rev 1999;79(3):703-761.
        PMID 10390517, DOI 10.1152/physrev.1999.79.3.703   (locator: section on the fibre-matrix model)
    the interstitial pressure-volume curve is NOT linear, which is why the oedema plateaus and why no
        volume reading can be converted to a pressure reading without that curve
        Guyton AC. Circ Res 1965;16(5):452-460. PMID 14289154, DOI 10.1161/01.RES.16.5.452
        A convention trap the source report caught and this file repeats so it is not re-made: the
        string 16.5.452 in that DOI is volume, issue and first page. It is not a year. Reading it as a
        year dates the work to the 1990s and mis-sorts any bibliography keyed on it.
    the interstitial matrix and the lymphatic pump as the elements that set transcapillary balance,
        i.e. that Q in the relation above is pumped and can be varied, which is what makes the second
        measurement physically available
        Granger HJ. Microvasc Res 1979;18(3):209-216. PMID 386049, DOI 10.1016/0026-2862(79)90029-3

THE EQUALLY INFORMED CONTROL is the single-observation reading sigma = 1 - r. It sees the same patient,
the same lymph sample and the same plasma sample, and it is what practice uses. It is not "another
method": it is this method with the second flow missing. The file reports its exact bias, which is
1/(1+Pe) in relative terms, so the control is scored rather than dismissed.

FALSIFIER, written before any measurement. If r is found to be INDEPENDENT of lymph flow over a flow
range of three-fold or more, then Pe >> 1 throughout, the naive reading is already correct, the second
measurement buys nothing and this file is pointless. The same data kill it in a second way: if the
two-flow inversion returns a negative PS or a sigma outside [0, 1], the one-barrier steady-state
relation above is the wrong model for that tissue and no amount of extra flows repairs it. Both
outcomes are detected and reported by name, not swallowed.

WHAT THIS FILE WILL NOT DO. It does not report a sigma for any real tissue, because no report in the
harvest carries a measured lymph-to-plasma ratio at two flows with its locator. The verdict it returns
on real data is therefore CANNOT DECIDE plus the one-sided bound, and that is the honest output.

review_state: PENDING_INDEPENDENT_REVIEW. No clinical recommendation is implied.
"""
from __future__ import annotations

import json
import os

PROTEIN_REDUCTION_ARMS_PCT = (33.0, 68.0)   # PMID 6881362


def ratio_from_parameters(sigma: float, PS: float, Q: float) -> float:
    """Forward model: lymph-to-plasma protein ratio r at lymph flow Q. Dimensionless out.

    Q and PS carry the same unit (volume per time); only their ratio enters, which is why the
    inversion below is insensitive to whether that unit is mL/min or mL/h as long as it is the SAME
    unit in both. Mixing them is the cheapest way to get a wrong sigma that looks plausible.
    """
    if PS <= 0.0 or Q <= 0.0:
        raise ValueError('PS and Q must be positive and in the same volume-per-time unit')
    return (1.0 - sigma) + sigma / (1.0 + Q / PS)


def naive_sigma(r: float) -> float:
    """The equally informed control: read 1 - r as sigma. A LOWER bound, never an estimate."""
    return 1.0 - r


def relative_bias_of_naive(Q: float, PS: float) -> float:
    """Exact relative shortfall of the control, 1/(1+Pe). Negative means underestimate."""
    return -1.0 / (1.0 + Q / PS)


def invert_two_flows(r1: float, Q1: float, r2: float, Q2: float) -> dict:
    """Closed-form sigma and PS from the protein ratio at two lymph flows.

    The denominator is a difference of two quotients, so this is exactly the shape the cancellation
    screen exists for: if the two flows are close, both numerator and denominator go to zero together
    and the quotient is set by measurement noise rather than by physiology. `conditioning` below is
    that risk made numerical, and `repair_flow_ratio_for` inverts it into the flow separation the
    measurement must achieve.
    """
    a = (1.0 - r1) / Q1
    b = (1.0 - r2) / Q2
    den = a - b
    if den == 0.0:
        return {'verdict': 'CANNOT DECIDE', 'why': 'the two flows gave the same quotient; '
                                                   'no information about sigma in this pair'}
    PS = (r1 - r2) / den
    if PS <= 0.0:
        return {'verdict': 'MODEL REJECTED', 'why': f'implied PS = {PS:.6g} is not positive; the '
                                                    'one-barrier steady-state relation is wrong here',
                'PS': PS}
    sigma = (1.0 - r1) * (PS + Q1) / Q1
    if not 0.0 <= sigma <= 1.0:
        return {'verdict': 'MODEL REJECTED', 'why': f'implied sigma = {sigma:.6g} is outside [0, 1]',
                'sigma': sigma, 'PS': PS}
    return {'verdict': 'IDENTIFIED', 'sigma': sigma, 'PS': PS,
            'Pe_at_Q1': Q1 / PS, 'Pe_at_Q2': Q2 / PS,
            'naive_sigma_at_Q1': naive_sigma(r1),
            'naive_relative_bias_at_Q1': relative_bias_of_naive(Q1, PS)}


def conditioning(r1: float, Q1: float, r2: float, Q2: float, dr: float = 0.01) -> dict:
    """How far sigma moves when each measured ratio is perturbed by `dr`, dimensionless in both.

    This is the number that decides whether the measurement is worth making, so it is reported with
    the verdict and not as an appendix.
    """
    base = invert_two_flows(r1, Q1, r2, Q2)
    if base['verdict'] != 'IDENTIFIED':
        return {'verdict': base['verdict']}
    out = {}
    for label, pert in (('r1', (dr, 0.0)), ('r2', (0.0, dr))):
        up = invert_two_flows(r1 + pert[0], Q1, r2 + pert[1], Q2)
        out[f'dsigma_per_{label}'] = ((up.get('sigma', float('nan')) - base['sigma']) / dr
                                     if up['verdict'] == 'IDENTIFIED' else float('nan'))
    out['sigma'] = base['sigma']
    return out


def one_sided_exclusion(r: float, sigma_baseline: float) -> dict:
    """The decision available from ONE observation, which is the asymmetric one.

    1 - r is a lower bound on sigma. So a measured ratio can refute a permeability lesion but cannot
    demonstrate one. Returns one choice out of three, never a parameter.
    """
    lower = naive_sigma(r)
    if lower > sigma_baseline:
        return {'choice': 'PERMEABILITY LESION EXCLUDED',
                'why': f'lower bound on sigma is {lower:.4f}, above the baseline {sigma_baseline:.4f};'
                       ' the barrier cannot have become more permeable',
                'sigma_lower_bound': lower}
    return {'choice': 'CANNOT DECIDE FROM ONE FLOW',
            'why': f'lower bound on sigma is {lower:.4f}, which is consistent with the baseline '
                   f'{sigma_baseline:.4f} and with any smaller sigma; measure r at a second lymph '
                   'flow',
            'sigma_lower_bound': lower}


def repair_flow_ratio_for(target_dsigma: float, sigma: float, PS: float, Q1: float,
                          dr: float = 0.01, max_ratio: float = 64.0) -> float:
    """Smallest second-to-first flow ratio at which a `dr` error in r moves sigma by less than
    `target_dsigma`. The repair response: scale the flow separation toward its ideal and read whether
    the headline survives."""
    ratio = 1.25
    while ratio <= max_ratio:
        Q2 = Q1 * ratio
        c = conditioning(ratio_from_parameters(sigma, PS, Q1), Q1,
                         ratio_from_parameters(sigma, PS, Q2), Q2, dr=dr)
        worst = max(abs(c.get('dsigma_per_r1', float('inf')) * dr),
                    abs(c.get('dsigma_per_r2', float('inf')) * dr))
        if worst < target_dsigma:
            return ratio
        ratio *= 1.25
    return float('nan')


def main() -> None:
    print('DECISION: permeability lesion, or driving-force change? And can one sample decide it?')
    print('  relation: r = (1-sigma) + sigma/(1+Pe),  Pe = Q/PS,  r = lymph/plasma protein ratio')

    print('\nEXACTNESS CHECK on the inversion, declared fixture, not a measurement.')
    print('  sigma and PS are chosen, r is computed forward at two flows, then recovered.')
    for sigma, PS, Q1, Q2 in ((0.90, 0.05, 0.02, 0.10),
                              (0.70, 0.20, 0.05, 0.40),
                              (0.99, 0.01, 0.10, 0.50),
                              (0.50, 1.00, 0.20, 2.00)):
        r1 = ratio_from_parameters(sigma, PS, Q1)
        r2 = ratio_from_parameters(sigma, PS, Q2)
        got = invert_two_flows(r1, Q1, r2, Q2)
        err_s = abs(got['sigma'] - sigma)
        err_p = abs(got['PS'] - PS)
        print(f'  sigma {sigma:.2f} PS {PS:.3g} (same unit as Q) -> r1 {r1:.4f} r2 {r2:.4f} '
              f'-> recovered sigma {got["sigma"]:.6f} PS {got["PS"]:.6g}  '
              f'errors {err_s:.1e} / {err_p:.1e}')
        assert err_s < 1e-9 and err_p < 1e-9

    print('\nTHE CONTROL, SCORED. Single-observation reading sigma = 1 - r, same sample, no second flow.')
    sigma, PS = 0.90, 0.05
    rows = []
    for Q in (0.005, 0.02, 0.05, 0.15, 0.50, 2.00):
        r = ratio_from_parameters(sigma, PS, Q)
        nb = relative_bias_of_naive(Q, PS)
        rows.append({'Q_same_unit_as_PS': Q, 'Pe': Q / PS, 'r': r,
                     'naive_sigma': naive_sigma(r), 'relative_bias': nb})
        print(f'  Pe {Q / PS:7.2f}  r {r:.4f}  control reads sigma {naive_sigma(r):.4f} '
              f'against true {sigma:.2f}   relative bias {nb * 100:+.1f} %')
    print('  the bias is exactly -1/(1+Pe): at a resting lymph flow with Pe below 1 the control')
    print('  understates sigma by more than half, and it understates it in ONE direction always.')

    print('\nTHE DECISION AVAILABLE FROM ONE SAMPLE, which is one-sided')
    for r in (0.05, 0.20, 0.743):
        d = one_sided_exclusion(r, sigma_baseline=0.90)
        print(f'  r = {r:.3f} -> {d["choice"]}: {d["why"]}')
    print('  so a high protein ratio is uninformative about sigma, and a LOW one is decisive.')
    print('  That is the reverse of reading the ratio as a permeability index.')

    print('\nCANCELLATION SCREEN. The inversion denominator is (1-r1)/Q1 - (1-r2)/Q2, a difference of')
    print('  quotients, so two nearly equal flows make it a ratio of two small numbers and the answer')
    print('  is then set by measurement noise. Measured as the sigma shift per 0.01 absolute error in r:')
    sweep = []
    for ratio in (1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0):
        Q1 = 0.02
        Q2 = Q1 * ratio
        c = conditioning(ratio_from_parameters(sigma, PS, Q1), Q1,
                         ratio_from_parameters(sigma, PS, Q2), Q2)
        worst = max(abs(c['dsigma_per_r1']), abs(c['dsigma_per_r2'])) * 0.01
        if worst != worst:
            sweep.append({'flow_ratio': ratio, 'sigma_shift_per_0.01_in_r': None,
                          'note': 'a 0.01 error in r drives the implied sigma outside [0,1]; '
                                  'the inversion is inadmissible at this flow separation'})
            print(f'  flow ratio {ratio:5.1f}x -> a 0.01 error in r pushes sigma outside [0,1]: '
                  'inversion inadmissible')
            continue
        sweep.append({'flow_ratio': ratio, 'sigma_shift_per_0.01_in_r': worst})
        print(f'  flow ratio {ratio:5.1f}x -> sigma moves {worst:.4f} per 0.01 error in r')
    need = repair_flow_ratio_for(0.05, sigma, PS, 0.02)
    print(f'  REPAIR RESPONSE: scaling the flow separation toward its ideal, a flow ratio of '
          f'{need:.2f}x or more')
    print('  holds the sigma error under 0.05 for a 0.01 error in r. The headline -- that the control')
    print('  is a one-sided bound with bias -1/(1+Pe) -- is unaffected by the conditioning, because it')
    print('  is an inequality on a single reading and uses no difference at all. The conditioning')
    print('  constrains only the two-flow inversion, and it is reported as a protocol requirement.')

    print('\nVERDICT ON REAL TISSUE: CANNOT DECIDE. No report in this harvest carries a measured')
    print('  lymph-to-plasma protein ratio at two lymph flows with a locator, so no sigma is reported')
    print('  for any tissue here. What is delivered is the protocol requirement and the one-sided')
    print(f'  bound, plus the measured anchor that the slope of a {PROTEIN_REDUCTION_ARMS_PCT[0]:.0f} % '
          f'and {PROTEIN_REDUCTION_ARMS_PCT[1]:.0f} % plasma-protein')
    print('  reduction gives the PRODUCT K_f*sigma (PMID 6881362) and therefore cannot close the gate.')

    out = {
        'relation': 'r = (1-sigma) + sigma/(1+Pe), Pe = Q/PS, r dimensionless, Q and PS same volume/time unit',
        'control_is': 'single-observation reading sigma = 1 - r',
        'control_exact_relative_bias': '-1/(1+Pe), always an underestimate of sigma',
        'one_sided_result': 'a single steady observation can exclude a permeability lesion, never establish one',
        'inversion_closed_form': {'PS': '(r1-r2)/((1-r1)/Q1-(1-r2)/Q2)',
                                  'sigma': '(1-r1)*(PS+Q1)/Q1'},
        'fixture_roundtrip_max_abs_error': 1e-9,
        'control_scored_rows': rows,
        'conditioning_sweep': sweep,
        'required_flow_ratio_for_sigma_error_below_0.05': need,
        'measured_anchors': [
            {'value': 33.0, 'unit': '% reduction in plasma protein', 'pmid': '6881362',
             'doi': '10.1152/ajpheart.1983.245.2.h284', 'locator': 'abstract, protein reduction arm',
             'identifies': 'the product K_f*sigma, not sigma'},
            {'value': 68.0, 'unit': '% reduction in plasma protein', 'pmid': '6881362',
             'doi': '10.1152/ajpheart.1983.245.2.h284', 'locator': 'abstract, protein reduction arm',
             'identifies': 'the product K_f*sigma, not sigma'},
        ],
        'structural_anchors': [
            {'pmid': '10390517', 'doi': '10.1152/physrev.1999.79.3.703',
             'locator': 'fibre-matrix model section', 'carries': 'sigma as a property of a structure; '
                                                                 'simultaneous water and solute flux method'},
            {'pmid': '14289154', 'doi': '10.1161/01.RES.16.5.452',
             'locator': 'Circ Res 1965;16(5):452-460',
             'carries': 'interstitial pressure-volume curve is nonlinear',
             'convention_trap': '16.5.452 is volume/issue/page, not a year'},
            {'pmid': '386049', 'doi': '10.1016/0026-2862(79)90029-3',
             'locator': 'Microvasc Res 1979;18(3):209-216',
             'carries': 'lymphatic pump sets Q, which is what makes the second flow available'},
        ],
        'source_report': 'results/BT-CAP-Q020/RESULTS.md section 4',
        'verdict_on_real_tissue': 'CANNOT DECIDE; no measured two-flow protein ratio with a locator',
        'falsifier': ('r independent of lymph flow over a three-fold range kills the file; a negative '
                      'implied PS or a sigma outside [0,1] rejects the one-barrier model itself'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    d = 'results/ASSEMBLY_REFLECTION_COEFFICIENT'
    os.makedirs(d, exist_ok=True)
    with open(d + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('\nwrote', d + '/decision.json')


if __name__ == '__main__':
    main()
