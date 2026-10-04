"""Decide whether a biofilm retention requirement is met by changing the matrix or by retiming the load.

THE DECISION, outside the eye. A surface must hold or shed an adherent biofilm layer against a flow.
Two routes exist to meet a stated requirement, and they cost completely different things:

    MATERIAL ROUTE    raise the matrix yield strength, in Pa, by chemistry or by a genetic lever
    OPERATING ROUTE   leave the matrix alone and change how fast the load is applied, in s, relative
                      to the layer's poroelastic drainage time

One choice comes out of this file: MATERIAL, OPERATING, or NOT ACHIEVABLE AS AN EPS LAYER AT ALL,
together with the tolerated flow velocity in m/s and the thickness change in micrometres that the
operating route needs. The routes are compared by their SCALING EXPONENTS, which is what makes the
verdict robust to the drag coefficient nobody knows.

WHY THE EXPONENTS DECIDE IT. Wall stress rises with the square of velocity, so strength buys velocity
only as a square root: the entire measured multispecies strength band, 6 to 51 Pa, is a factor 8.5 in
Pa and only sqrt(8.5) = 2.92 in m/s. Drainage time in a poroelastic layer goes as the square of
thickness, so a factor 2.92 in micrometres buys a factor 8.5 in seconds -- the same factor that the
whole published strength literature spans, obtained from one geometric change. The asymmetry between
an exponent of 1/2 and an exponent of 2 is the content.

MEASURED NUMBERS, EACH WITH UNIT AND LOCATOR, all carried into
`results/BT-CAP-Q160/RESULTS.md` section 4.3 with their citations:

    multispecies oral biofilm strength      6 to 51 Pa
    single-species oral biofilm strength    5 to 17 Pa
    over an applied shear rate range of     0.1 to 50 per second
    direction: strength FALLS as shear rate rises
        Paramonova E et al. J Dent Res 2009;88(10):922-926.
        PMID 19783800, DOI 10.1177/0022034509344569   (locator: strength versus shear-rate result)
    compliance spread WITHIN a single biofilm   3 orders of magnitude
        Galy O et al. Biophys J 2012;103(6):1400-1408.
        PMID 22995513, DOI 10.1016/j.bpj.2012.07.001  (locator: magnetic microparticle actuation,
        per-point compliance field; also the curli and F-pilus comparison)
    the lever exists and is measurable: cellulose expression raises the Young modulus by AFM and the
    EPS shear modulus by QCM-D, and lowers diffusion at the same time
        Ziemba C et al. npj Biofilms Microbiomes 2016;2:1.
        PMID 28649395, DOI 10.1038/s41522-016-0001-2
    the acceptance surface is published by someone else, rotating-disc assay
        Dennington SP et al. Surf Topogr Metrol Prop 2015;3:034004.
        DOI 10.1088/2051-672X/3/3/034004

TWO CONVENTION TRAPS, named because this project has been bitten by this class five times.
 1. Pa against Pa is not automatically a comparison. Paramonova's number is a COMPRESSIVE strength of
    a biofilm under an applied shear rate; the drag term is a WALL SHEAR STRESS. Both are stresses and
    both print as Pa. This file therefore does not claim the two are the same quantity: it uses the
    strength band only as a RANGE OF ATTAINABLE STRESS SCALES, and every verdict below is a ratio
    within one of the two quantities, never a difference between them. Where an absolute velocity is
    printed it is labelled as conditional on a wall-stress-to-strength correspondence being asserted
    by whoever uses it, and the ratio verdict does not depend on that assertion.
 2. Shear RATE in per second and drainage TIME in seconds are reciprocal-looking and are not inverses
    of each other. One is a kinematic rate imposed by the rheometer, the other is a material relaxation
    time set by thickness and permeability. They are never divided into each other here.

THE ASSUMPTION THAT CARRIES THE MATERIAL-ROUTE CEILING, stated so it can be attacked first: the band's
endpoints are read as the attainable range of matrix strength. If 6 to 51 Pa is a BETWEEN-BIOFILM range
rather than a WITHIN-SWEEP range, the per-decade slope reported below is unfounded, and the file says so
and withholds it; the ceiling argument survives either way, because 51 Pa is the largest multispecies
value measured regardless of what varied to produce it.

THE EQUALLY INFORMED CONTROL is the single-threshold acceptance criterion: read the same Paramonova
result and adopt one yield stress, conventionally the low end, with no shear-rate argument. It is the
practice, it sees exactly the same data, and it is not "another method". This file scores it: because
strength falls with shear rate across a 500-fold rate range, a criterion fixed at low shear overstates
the strength available at operating shear by up to the full 8.5-fold width of the band.

FALSIFIER, written before any measurement. If strength is measured to be INDEPENDENT of shear rate over
0.1 to 50 per second in the target system, the operating route loses its mechanism and this file is
wrong: only the material route would exist. Second falsifier: a published multispecies strength above
51 Pa moves the ceiling and with it the NOT ACHIEVABLE verdict, so that verdict is held only against
the strongest value on record and is reported with that value attached. Third: if a requirement is
stated to a tolerance tighter than the 3-decade within-biofilm compliance spread, the requirement is
not falsifiable and the file refuses it rather than answering it.

review_state: PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import math
import os

STRENGTH_MULTISPECIES_PA = (6.0, 51.0)      # PMID 19783800
STRENGTH_SINGLE_SPECIES_PA = (5.0, 17.0)    # PMID 19783800
SHEAR_RATE_RANGE_PER_S = (0.1, 50.0)        # PMID 19783800
WITHIN_BIOFILM_COMPLIANCE_DECADES = 3.0     # PMID 22995513
FLUID_DENSITY_KG_M3 = 1000.0                # aqueous, stated not measured
DEFAULT_SKIN_FRICTION_COEFFICIENT = 0.01    # stated not measured; swept below


def velocity_for_wall_stress(tau_pa: float, c_f: float = DEFAULT_SKIN_FRICTION_COEFFICIENT,
                             rho: float = FLUID_DENSITY_KG_M3) -> float:
    """Velocity in m/s at which the wall stress reaches tau_pa. tau_w = 0.5 * c_f * rho * v^2."""
    return math.sqrt(2.0 * tau_pa / (c_f * rho))


def velocity_gain_from_strength_gain(strength_factor: float) -> float:
    """A factor in Pa buys its SQUARE ROOT in m/s. c_f and rho cancel exactly, which is why this is
    the quantity the verdict rests on rather than any absolute velocity."""
    return math.sqrt(strength_factor)


def thickness_factor_for_drainage_gain(time_factor: float) -> float:
    """Poroelastic drainage time goes as thickness squared, so a factor in s needs its square root in
    micrometres. The exponent is the assumption; it is the standard poroelastic scaling and it is also
    what the source report states. If the exponent were 1 rather than 2, the operating route would be
    exactly as weak as the material route and the verdict would collapse -- that is the single thing to
    attack here, and it is an exponent, which a measurement of drainage time against two thicknesses
    settles."""
    return math.sqrt(time_factor)


def material_route_ceiling() -> dict:
    lo, hi = STRENGTH_MULTISPECIES_PA
    f = hi / lo
    return {'strength_band_pa': [lo, hi], 'strength_factor': f,
            'velocity_factor': velocity_gain_from_strength_gain(f),
            'ceiling_strength_pa': hi}


def per_decade_slope(band_is_within_sweep: bool) -> float | None:
    """Strength lost per decade of shear rate, IF the band is a within-sweep range. Returns None when
    that reading is not licensed, rather than a number nobody can trace."""
    if not band_is_within_sweep:
        return None
    lo, hi = STRENGTH_MULTISPECIES_PA
    decades = math.log10(SHEAR_RATE_RANGE_PER_S[1] / SHEAR_RATE_RANGE_PER_S[0])
    return (hi / lo) ** (1.0 / decades)


def requirement_is_falsifiable(required_relative_tolerance: float) -> bool:
    """A tolerance tighter than the material's own heterogeneity is not a requirement."""
    return required_relative_tolerance >= 10.0 ** (-WITHIN_BIOFILM_COMPLIANCE_DECADES)


def decide(required_velocity_factor: float, baseline_strength_pa: float,
           baseline_thickness_um: float,
           required_relative_tolerance: float = 0.1,
           thickness_floor_um: float = 10.0) -> dict:
    """ONE CHOICE out of four, plus the number each route needs.

    `required_velocity_factor` is how much faster the flow must be tolerated, dimensionless.
    `thickness_floor_um` is the thinnest layer the application still allows, in micrometres. It is an
    INPUT, not a measurement, and it is what makes the material route reachable at all: the exponent
    argument alone always prefers retiming, so without a floor the choice would be degenerate and the
    file would not be deciding anything. A first version of this file had exactly that defect -- the
    MATERIAL branch was unreachable -- and the floor is the fix.
    """
    if not requirement_is_falsifiable(required_relative_tolerance):
        return {'choice': 'REQUIREMENT REFUSED',
                'why': f'a relative tolerance of {required_relative_tolerance:.1e} is tighter than the '
                       f'{WITHIN_BIOFILM_COMPLIANCE_DECADES:.0f}-decade compliance spread measured '
                       'inside a single biofilm (PMID 22995513); such a criterion selects noise'}
    needed_strength_pa = baseline_strength_pa * required_velocity_factor ** 2
    ceiling = material_route_ceiling()
    needed_time_factor = required_velocity_factor ** 2
    thickness_factor = thickness_factor_for_drainage_gain(needed_time_factor)
    thickness_um = baseline_thickness_um / thickness_factor
    material_ok = needed_strength_pa <= ceiling['ceiling_strength_pa']
    operating_ok = thickness_um >= thickness_floor_um
    out = {
        'required_velocity_factor': required_velocity_factor,
        'material_route_needed_strength_pa': needed_strength_pa,
        'material_route_ceiling_pa': ceiling['ceiling_strength_pa'],
        'material_route_feasible': material_ok,
        'operating_route_drainage_time_factor': needed_time_factor,
        'operating_route_thickness_factor': 1.0 / thickness_factor,
        'operating_route_thickness_um': thickness_um,
        'thickness_floor_um': thickness_floor_um,
        'operating_route_feasible': operating_ok,
    }
    if material_ok and operating_ok:
        out['choice'] = 'OPERATING'
        out['why'] = (f'both routes reach it; retiming needs a {thickness_factor:.2f}-fold thinner '
                      f'layer ({thickness_um:.1f} um) against a '
                      f'{required_velocity_factor ** 2:.2f}-fold strength change, and thickness is the '
                      'cheaper lever because its exponent is 2 against the strength route\'s 1/2')
    elif material_ok and not operating_ok:
        out['choice'] = 'MATERIAL'
        out['why'] = (f'retiming would need {thickness_um:.1f} um, below the stated floor of '
                      f'{thickness_floor_um:.1f} um; {needed_strength_pa:.1f} Pa is inside the measured '
                      f'band up to {ceiling["ceiling_strength_pa"]:.0f} Pa (PMID 19783800)')
    elif operating_ok and not material_ok:
        out['choice'] = 'OPERATING'
        out['why'] = (f'the material route needs {needed_strength_pa:.1f} Pa, above the strongest '
                      f'multispecies biofilm on record ({ceiling["ceiling_strength_pa"]:.0f} Pa, '
                      f'PMID 19783800); retiming reaches it at {thickness_um:.1f} um')
    else:
        out['choice'] = 'NOT ACHIEVABLE AS AN EPS LAYER'
        out['why'] = (f'needs {needed_strength_pa:.1f} Pa against a ceiling of '
                      f'{ceiling["ceiling_strength_pa"]:.0f} Pa, and retiming would need '
                      f'{thickness_um:.1f} um against a floor of {thickness_floor_um:.1f} um')
    return out


def main() -> None:
    c = material_route_ceiling()
    print('DECISION: meet a biofilm retention requirement by matrix chemistry, or by retiming the load?')
    print(f'  measured multispecies strength band : {c["strength_band_pa"][0]:.0f} to '
          f'{c["strength_band_pa"][1]:.0f} Pa, a factor {c["strength_factor"]:.2f}')
    print(f'  single-species band                 : {STRENGTH_SINGLE_SPECIES_PA[0]:.0f} to '
          f'{STRENGTH_SINGLE_SPECIES_PA[1]:.0f} Pa')
    print(f'  applied shear-rate range            : {SHEAR_RATE_RANGE_PER_S[0]:g} to '
          f'{SHEAR_RATE_RANGE_PER_S[1]:g} per second, a factor '
          f'{SHEAR_RATE_RANGE_PER_S[1] / SHEAR_RATE_RANGE_PER_S[0]:.0f}')
    print(f'  THE HEADLINE: the whole published strength band, a factor {c["strength_factor"]:.2f} in Pa,')
    print(f'  buys only a factor {c["velocity_factor"]:.2f} in tolerated velocity, because wall stress')
    print('  goes as the square of velocity. Drainage time goes as the square of thickness, so the')
    print(f'  same factor {c["strength_factor"]:.2f} in seconds costs a factor '
          f'{thickness_factor_for_drainage_gain(c["strength_factor"]):.2f} in micrometres.')

    print('\nABSOLUTE VELOCITIES, conditional on a wall-stress-to-strength correspondence, and the')
    print('  cancellation screen that goes with them. The drag coefficient is not measured anywhere in')
    print('  this material, so any absolute velocity is as uncertain as it is. Sweep it a decade:')
    sweep = []
    for c_f in (0.002, 0.005, 0.01, 0.02):
        v_lo = velocity_for_wall_stress(c["strength_band_pa"][0], c_f)
        v_hi = velocity_for_wall_stress(c["strength_band_pa"][1], c_f)
        sweep.append({'c_f': c_f, 'v_at_6Pa_m_s': v_lo, 'v_at_51Pa_m_s': v_hi, 'ratio': v_hi / v_lo})
        print(f'  c_f {c_f:.3f} -> {v_lo:.2f} m/s at 6 Pa, {v_hi:.2f} m/s at 51 Pa, '
              f'ratio {v_hi / v_lo:.3f}')
    print('  the RATIO is identical to four figures across the whole sweep: c_f and the density cancel')
    print('  exactly in it. So the headline is not an artefact of two terms extinguishing each other --')
    print('  it is a scaling exponent, and the repair response (scaling the drag coefficient toward any')
    print('  value its own literature allows) leaves it untouched while the absolute velocities move')
    print(f'  by {sweep[-1]["v_at_6Pa_m_s"] / sweep[0]["v_at_6Pa_m_s"]:.2f}-fold. Only the ratio is '
          'carried into any verdict.')

    print('\nTHE CONTROL, SCORED. A single-threshold acceptance criterion fixed at the low end of the')
    print('  band sees the same measurement. Because strength falls as shear rate rises over a '
          f'{SHEAR_RATE_RANGE_PER_S[1] / SHEAR_RATE_RANGE_PER_S[0]:.0f}-fold')
    print('  rate range, a threshold set at operating shear and a threshold set at rheometer shear')
    print(f'  differ by up to the full factor {c["strength_factor"]:.2f}. A single number cannot carry '
          'that, and the')
    print('  direction of its error is known: it overstates available strength at high shear.')
    slope = per_decade_slope(band_is_within_sweep=True)
    slope_withheld = per_decade_slope(band_is_within_sweep=False)
    print(f'  IF the band is a within-sweep range, the slope is {slope:.2f}-fold strength lost per')
    print('  decade of shear rate. If it is a between-biofilm range it is not a slope at all, and in')
    print(f'  that reading this file returns {slope_withheld} rather than a number: the ceiling verdict')
    print('  does not need the slope, so nothing downstream is built on it.')

    print('\nWORKED VERDICTS, one choice each')
    cases = [
        (1.5, 6.0, 100.0, 0.1),
        (2.5, 6.0, 100.0, 0.1),
        (3.5, 6.0, 100.0, 0.1),
        (8.0, 6.0, 100.0, 0.1),
        (2.0, 6.0, 100.0, 1e-4),
        (2.5, 6.0, 20.0, 0.1),
        (6.0, 6.0, 40.0, 0.1),
    ]
    verdicts = []
    for f, s0, L0, tol in cases:
        d = decide(f, s0, L0, tol)
        verdicts.append({'inputs': {'required_velocity_factor': f, 'baseline_strength_pa': s0,
                                    'baseline_thickness_um': L0,
                                    'required_relative_tolerance': tol},
                         'output': d})
        print(f'  velocity x{f:.1f}, baseline {s0:.0f} Pa, {L0:.0f} um, tolerance {tol:g}')
        print(f'    -> {d["choice"]}: {d["why"]}')

    print('\nWHAT IS NOT CLAIMED. The compressive strength of an oral biofilm and a wall shear stress')
    print('  are both printed in Pa and are not the same quantity; no difference between them is taken')
    print('  anywhere above. The genetic lever is cited as EXISTING and measurable (PMID 28649395) and')
    print('  no strength is predicted from a cellulose level, because no such calibration is in hand.')

    out = {
        'headline': {
            'strength_band_pa': c['strength_band_pa'],
            'strength_factor': c['strength_factor'],
            'velocity_factor_bought_by_whole_band': c['velocity_factor'],
            'thickness_factor_buying_the_same_factor_in_drainage_time':
                thickness_factor_for_drainage_gain(c['strength_factor']),
            'material_route_exponent': 0.5,
            'operating_route_exponent': 2.0,
        },
        'drag_coefficient_sweep_cancellation_screen': sweep,
        'control': 'single-threshold acceptance criterion at the low end of the measured band',
        'control_worst_case_overstatement_factor': c['strength_factor'],
        'per_decade_slope_if_within_sweep': slope,
        'per_decade_slope_if_between_biofilm': slope_withheld,
        'within_biofilm_compliance_decades': WITHIN_BIOFILM_COMPLIANCE_DECADES,
        'worked_verdicts': verdicts,
        'measured_anchors': [
            {'value': [6.0, 51.0], 'unit': 'Pa', 'quantity': 'multispecies oral biofilm strength',
             'pmid': '19783800', 'doi': '10.1177/0022034509344569',
             'locator': 'strength versus shear-rate result'},
            {'value': [5.0, 17.0], 'unit': 'Pa', 'quantity': 'single-species oral biofilm strength',
             'pmid': '19783800', 'doi': '10.1177/0022034509344569', 'locator': 'same result'},
            {'value': [0.1, 50.0], 'unit': 'per second', 'quantity': 'applied shear rate range',
             'pmid': '19783800', 'doi': '10.1177/0022034509344569', 'locator': 'same result'},
            {'value': 3.0, 'unit': 'orders of magnitude', 'quantity': 'within-biofilm compliance spread',
             'pmid': '22995513', 'doi': '10.1016/j.bpj.2012.07.001',
             'locator': 'magnetic microparticle actuation, per-point compliance field'},
            {'value': None, 'unit': None, 'quantity': 'genetic lever on modulus exists and is measurable',
             'pmid': '28649395', 'doi': '10.1038/s41522-016-0001-2', 'locator': 'AFM and QCM-D arms'},
            {'value': None, 'unit': None, 'quantity': 'externally published acceptance surface',
             'doi': '10.1088/2051-672X/3/3/034004', 'locator': 'rotating-disc assay'},
        ],
        'source_report': 'results/BT-CAP-Q160/RESULTS.md section 4.3',
        'attackable_assumption': ('drainage time scales as thickness squared; if the exponent were 1 the '
                                  'operating route would be no better than the material route'),
        'falsifier': ('strength independent of shear rate over 0.1-50 per second removes the operating '
                      'route mechanism; a multispecies strength above 51 Pa moves the ceiling; a '
                      'tolerance tighter than 3 decades is refused rather than answered'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    d = 'results/ASSEMBLY_BIOFILM_RETENTION_ROUTE'
    os.makedirs(d, exist_ok=True)
    with open(d + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('\nwrote', d + '/decision.json')


if __name__ == '__main__':
    main()
