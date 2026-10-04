"""Decide whether a patterned ligand source can beat a uniform bath, before building the actuator.

THE DECISION, outside the eye. A localised source is to be placed in tissue to activate a membrane
receptor more effectively than a uniform bath of the SAME TOTAL DOSE. The design has a geometry, a
patch spacing in micrometres, and a schedule, a pulse duration in seconds. The question that decides
whether the actuator is worth building is whether the receptor can still tell the patches apart: a
receptor that wanders across several patches during the pulse averages them, and the patterned source
then delivers exactly what the bath delivers. One number comes out of this file: the minimum patch
spacing in micrometres, for a stated pulse duration, below which the design IS the bath.

WHY A NECESSARY CONDITION IS THE HONEST OUTPUT. The source report establishes that the full inverse
problem here is not identifiable: over a fixed run the ligand diffusivity, the receptor diffusivity and
the duration enter mainly through products and diffusive length scales, so many parameter sets map to
nearly the same peak time and activation ratio, and picking one is drawing a random member off a level
set. `results/BT-CAP-Q080/RESULTS.md` says so in its own section on uncertainty. This file therefore
does not invert anything. It returns a GATE that any admissible design must pass, computed from a
length and a time that were both measured, and it states that passing the gate is necessary and not
sufficient.

MEASURED NUMBERS, EACH WITH UNIT AND LOCATOR.
    nested confinement compartments, two scales      210 nm and 730 nm
    residence time before a hop, matched to each     about 45 ms and about 760 ms
        Suzuki K, Ritchie K, Kajikawa E, Fujiwara T, Kusumi A. Biophys J 2005;88(5):3659-3680.
        PMID 15681644, DOI 10.1529/biophysj.104.048538
        (locator: single-molecule hop-diffusion result, nested compartment sizes and residence times)
    receptors redistribute into ligand-induced microdomains, and complex size scales with receptor
    DENSITY -- which is why free ligand concentration is not the only state variable
        Calebiro D et al. PNAS 2013;110(2):743-748.
        PMID 23267088, DOI 10.1073/pnas.1205798110
    equilibrium affinity of the tracer-adjacent antagonist, pKb 6.5 to 6.7, i.e. 200 to 316 nM
        Briddon SJ et al. PNAS 2004;101(13):4673-4678.
        PMID 15070776, DOI 10.1073/pnas.0400420101
        PROVENANCE CAVEAT, carried from the source report rather than hidden: that is an equilibrium
        affinity for a COMPETITIVE ANTAGONIST accumulated from cAMP and inositol-phosphate assays. It
        is not an on-rate/off-rate pair for a tracer. It is used below only to state the concentration
        SCALE a patch must exceed to activate anything, and no ratio in this file depends on it.

THE QUANTITY THIS FILE COMPUTES RATHER THAN RESTATES. The two measured length-and-time pairs are
independent, and each one implies a macroscopic diffusion coefficient through the two-dimensional
relation <r^2> = 4 D t, so D = L^2 / (4 tau). Computing both is a cross-check on the hop picture that
the source report does not run: the two nested scales agree on D to within a factor of 1.40, which is
why the averaging length below is quoted as a window rather than a point. From that window the
receptor's own averaging length over a pulse is sqrt(4 D T), and that is the gate.

A UNIT TRAP, named because this project has been bitten by this class. The source cell carries a
receptor diffusivity in um^2/s and a concentration in uM, and a sibling cell carries a crista length in
um. `um` and `uM` differ by case and by DIMENSION -- a length against a concentration -- and a
case-folding step in the upstream job counted them as a shared dimension. Nothing here converts between
them. The two nanometre compartment sizes are converted to micrometres explicitly, in one place.

THE EQUALLY INFORMED CONTROL is the uniform bath at equal total dose: the path that already works,
needs no actuator, and sees the same receptor diffusivity, the same pulse duration and the same dose.
This file scores it exactly rather than dismissing it: when the proposed spacing falls below the
averaging length, the patterned design and the bath are the same stimulus as far as the receptor
population is concerned, so the control is not merely competitive, it is identical, and the whole cost
of the actuator buys nothing.

FALSIFIER, written before any measurement. If a measured activation-integral ratio comes out at or
below 1 for a design whose spacing EXCEEDS the computed averaging length, the gate is not the binding
constraint and this file is wrong about which quantity decides. Second falsifier, cheaper: if a third
independent length-and-time pair for the same receptor implies a macroscopic D outside the 0.175 to
0.245 um^2/s window by more than a factor of two, the averaging length is wrong, and every spacing this
file prints moves as the square root of that error -- the file prints that sensitivity so the
falsification is quantitative rather than rhetorical.

review_state: PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import math
import os

COMPARTMENT_NM = (210.0, 730.0)       # PMID 15681644
RESIDENCE_S = (0.045, 0.760)          # PMID 15681644, matched to the compartment sizes in order
PKB_RANGE = (6.5, 6.7)                # PMID 15070776
NM_PER_UM = 1000.0


def macroscopic_D_um2_s(compartment_nm: float, residence_s: float, dimension: int = 2) -> float:
    """Macroscopic diffusion coefficient implied by one hop scale. um^2/s out.

    <r^2> = 2*d*D*t in d dimensions, so D = L^2 / (2*d*tau). The membrane case is d = 2, giving
    L^2/(4 tau). `dimension` is exposed because the convention is the first thing an independent
    reader will challenge, and because the AGREEMENT RATIO between the two scales is invariant to it,
    which is the point of the screen below.
    """
    L_um = compartment_nm / NM_PER_UM
    return L_um ** 2 / (2.0 * dimension * residence_s)


def averaging_length_um(D_um2_s: float, pulse_duration_s: float, dimension: int = 2) -> float:
    """Root-mean-square displacement of a receptor over the pulse. um out.

    This is a single product under a square root. There is no difference of two terms anywhere in it,
    so the cancellation screen has nothing to find here: the headline cannot arise from two terms
    extinguishing each other. The screen is run instead on the D cross-check, which IS a ratio.
    """
    return math.sqrt(2.0 * dimension * D_um2_s * pulse_duration_s)


def kd_nM_from_pkb(pkb: float) -> float:
    return 10.0 ** (9.0 - pkb)


def decide(patch_spacing_um: float, pulse_duration_s: float,
           patch_peak_concentration_nM: float | None = None) -> dict:
    """ONE CHOICE: GO, NO-GO, or NO-GO ON DOSE. Plus the margin that produced it."""
    D_lo = min(macroscopic_D_um2_s(c, t) for c, t in zip(COMPARTMENT_NM, RESIDENCE_S))
    D_hi = max(macroscopic_D_um2_s(c, t) for c, t in zip(COMPARTMENT_NM, RESIDENCE_S))
    L_lo = averaging_length_um(D_lo, pulse_duration_s)
    L_hi = averaging_length_um(D_hi, pulse_duration_s)
    kd_hi_nM = kd_nM_from_pkb(PKB_RANGE[0])      # the weaker affinity, i.e. the higher Kd
    out = {
        'patch_spacing_um': patch_spacing_um,
        'pulse_duration_s': pulse_duration_s,
        'macroscopic_D_window_um2_s': [D_lo, D_hi],
        'averaging_length_window_um': [L_lo, L_hi],
        'minimum_admissible_spacing_um': L_hi,
        'margin': patch_spacing_um / L_hi,
        'kd_upper_nM': kd_hi_nM,
    }
    if patch_spacing_um < L_lo:
        out['choice'] = 'NO-GO'
        out['why'] = (f'the receptor crosses {L_lo:.2f} to {L_hi:.2f} um during a {pulse_duration_s:g} s '
                      f'pulse, more than the {patch_spacing_um:g} um spacing, so the population averages '
                      'the patches: this design IS the uniform bath and the actuator buys nothing')
    elif patch_spacing_um < L_hi:
        out['choice'] = 'CANNOT DECIDE'
        out['why'] = (f'{patch_spacing_um:g} um lies inside the averaging-length window '
                      f'{L_lo:.2f} to {L_hi:.2f} um, which is as wide as the two measured hop scales '
                      'disagree; a third hop measurement is required before this design is called '
                      'either way')
    elif patch_peak_concentration_nM is not None and patch_peak_concentration_nM < kd_hi_nM:
        out['choice'] = 'NO-GO ON DOSE'
        out['why'] = (f'spacing clears the gate with margin {out["margin"]:.2f}, but a peak of '
                      f'{patch_peak_concentration_nM:g} nM is below the {kd_hi_nM:.0f} nM affinity '
                      'scale, so the patch does not occupy the receptor whatever its geometry')
    else:
        out['choice'] = 'GO, NECESSARY CONDITION ONLY'
        out['why'] = (f'{patch_spacing_um:g} um exceeds the averaging length {L_hi:.2f} um by a factor '
                      f'{out["margin"]:.2f}, so the receptor population can still resolve the pattern. '
                      'This is necessary and NOT sufficient: the activation advantage over the bath '
                      'also needs the ligand field to hold a gradient, which no number here settles')
    return out


def main() -> None:
    print('DECISION: can a patterned ligand source beat a uniform bath of equal total dose?')
    print('  gate: can the receptor still resolve the pattern during the pulse?')

    print('\nCROSS-CHECK the source report does not run: two measured hop scales, one macroscopic D')
    Ds = []
    for c_nm, t_s in zip(COMPARTMENT_NM, RESIDENCE_S):
        D = macroscopic_D_um2_s(c_nm, t_s)
        Ds.append(D)
        print(f'  compartment {c_nm:5.0f} nm, residence {t_s * 1000:6.1f} ms -> D = {D:.4f} um^2/s')
    agreement = max(Ds) / min(Ds)
    print(f'  the two independent scales agree on D to within a factor {agreement:.2f}, which is why the')
    print('  averaging length below is a WINDOW and not a point. Two nested scales measured in one')
    print('  membrane landing within 1.4x of one macroscopic D is a nontrivial check on the hop picture.')

    print('\nCANCELLATION SCREEN, with the repair response. The agreement ratio is a quotient of two')
    print('  computed quantities, so it is the place where a convention could flatter the result. The')
    print('  convention is the dimension in <r^2> = 2 d D t. Scale it toward its correct value for a')
    print('  membrane, d = 2, and away from it, and read the ratio:')
    screen = []
    for d in (1, 2, 3):
        ds = [macroscopic_D_um2_s(c, t, dimension=d) for c, t in zip(COMPARTMENT_NM, RESIDENCE_S)]
        screen.append({'dimension': d, 'D_um2_s': ds, 'agreement_ratio': max(ds) / min(ds)})
        print(f'  d = {d} -> D = {ds[0]:.4f} and {ds[1]:.4f} um^2/s, agreement ratio '
              f'{max(ds) / min(ds):.4f}')
    print('  the ratio is identical across all three: the convention cancels exactly, so the agreement')
    print('  is a property of the two measurements and not of the factor chosen. The headline does not')
    print('  get worse under the repair. The absolute D does move, and the gate below therefore moves')
    print('  as its square root, which is reported rather than buried.')

    print('\nTHE GATE, as a table of minimum admissible spacings')
    rows = []
    for T in (0.1, 1.0, 10.0, 60.0, 300.0):
        lo = averaging_length_um(min(Ds), T)
        hi = averaging_length_um(max(Ds), T)
        rows.append({'pulse_duration_s': T, 'averaging_length_um': [lo, hi],
                     'minimum_admissible_spacing_um': hi})
        print(f'  pulse {T:6.1f} s -> receptor averages over {lo:6.3f} to {hi:6.3f} um; spacing must '
              f'exceed {hi:6.3f} um')
    print(f'  for comparison the receptor\'s own confinement compartments are '
          f'{COMPARTMENT_NM[0] / NM_PER_UM:.3f} and {COMPARTMENT_NM[1] / NM_PER_UM:.3f} um, so any')
    print('  pulse longer than about a second already carries the receptor across many compartments.')

    print(f'\nTHE DOSE SCALE, from pKb {PKB_RANGE[0]} to {PKB_RANGE[1]}: Kd between '
          f'{kd_nM_from_pkb(PKB_RANGE[1]):.0f} and {kd_nM_from_pkb(PKB_RANGE[0]):.0f} nM.')
    print('  A patch peak below that does not occupy the receptor whatever the geometry. Used as a')
    print('  threshold only; no ratio in this file depends on it, because its provenance is an')
    print('  equilibrium antagonist affinity and not a tracer rate pair.')

    print('\nWORKED VERDICTS, one choice each')
    verdicts = []
    for spacing, T, conc in ((1.0, 10.0, None), (5.0, 10.0, None), (3.0, 10.0, None),
                             (20.0, 10.0, 50.0), (20.0, 10.0, 1000.0), (0.5, 0.1, None)):
        d = decide(spacing, T, conc)
        verdicts.append({'inputs': {'patch_spacing_um': spacing, 'pulse_duration_s': T,
                                    'patch_peak_concentration_nM': conc}, 'output': d})
        print(f'  spacing {spacing:5.1f} um, pulse {T:5.1f} s, peak '
              f'{"n/a" if conc is None else f"{conc:g} nM":>8}  -> {d["choice"]}')
        print(f'      {d["why"]}')

    print('\nCONTROL, equally informed: the uniform bath at equal total dose. Scored, not dismissed --')
    print('  in every NO-GO row above the control delivers the SAME stimulus to the receptor population,')
    print('  so its score equals the patterned design\'s and the actuator cost is pure loss. The file')
    print('  returns GO only where the control is strictly distinguishable, and even then only as a')
    print('  necessary condition.')

    out = {
        'macroscopic_D_from_each_hop_scale_um2_s': dict(zip([f'{c:.0f}nm' for c in COMPARTMENT_NM], Ds)),
        'agreement_ratio_between_the_two_scales': agreement,
        'dimension_convention_screen': screen,
        'gate_table': rows,
        'kd_window_nM': [kd_nM_from_pkb(PKB_RANGE[1]), kd_nM_from_pkb(PKB_RANGE[0])],
        'worked_verdicts': verdicts,
        'control': 'uniform bath at equal total dose; identical stimulus whenever spacing < averaging length',
        'measured_anchors': [
            {'value': list(COMPARTMENT_NM), 'unit': 'nm', 'quantity': 'nested confinement compartments',
             'pmid': '15681644', 'doi': '10.1529/biophysj.104.048538',
             'locator': 'single-molecule hop-diffusion result'},
            {'value': list(RESIDENCE_S), 'unit': 's', 'quantity': 'residence time before a hop',
             'pmid': '15681644', 'doi': '10.1529/biophysj.104.048538',
             'locator': 'same result, matched to the compartment sizes'},
            {'value': list(PKB_RANGE), 'unit': 'pKb', 'quantity': 'equilibrium antagonist affinity',
             'pmid': '15070776', 'doi': '10.1073/pnas.0400420101',
             'locator': 'cAMP and inositol-phosphate assays',
             'caveat': 'competitive antagonist equilibrium affinity, not a tracer on/off pair'},
            {'value': None, 'unit': None,
             'quantity': 'complex size scales with receptor density, so concentration is not the only state',
             'pmid': '23267088', 'doi': '10.1073/pnas.1205798110', 'locator': 'single-molecule arm'},
        ],
        'source_report': 'results/BT-CAP-Q080/RESULTS.md',
        'scope': 'necessary condition on geometry only; sufficiency needs a ligand-field gradient this file does not compute',
        'falsifier': ('an activation-integral ratio at or below 1 for a spacing above the averaging '
                      'length refutes the gate; a third hop pair implying a D outside 0.175-0.245 '
                      'um^2/s by more than 2x moves every spacing as its square root'),
        'unit_trap_recorded': 'um is a length and uM is a concentration; no conversion is made between them',
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    d = 'results/ASSEMBLY_FOCAL_DELIVERY'
    os.makedirs(d, exist_ok=True)
    with open(d + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('\nwrote', d + '/decision.json')


if __name__ == '__main__':
    main()
