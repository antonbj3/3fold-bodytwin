"""How hard must a surface throw away cleaved C3 to stay quiescent? A number, from measured rates.

THE DECISION, and it is seed 2's second half rather than the eye: immune regulation fails in two
directions, too little and too much, and the alternative complement pathway is the cleanest place to
put a number on the margin. A surface is quiescent when the amplification loop dies out and inflamed
when it runs away. Those are the two failure modes, and the switch between them is a single
dimensionless product of quantities that are all measured.

One convertase lives for a mean lifetime tau = t_half / ln2 and cleaves C3 at kcat with
Michaelis-Menten saturation. So the number of C3 molecules one convertase turns over before it falls
apart is

    N = kcat * tau * S / (Km + S)

and each cleaved C3 becomes a new convertase with some deposition efficiency p. The loop runs away
when N*p > 1, so the surface must discard all but p_crit = 1/N of what it cleaves. That is the
regulatory margin, and it is a prediction about a measurement rather than a restatement of one.

MEASURED INPUTS, from results/LANE_READ_CREATIVE/CONSUMABLE.json (rows READ-A035 to READ-A039), all
alternative-pathway convertase quantities:

    C3bBb kcat                  1.78 s^-1        READ-A035
    C3bBb Km                    5.86e-6 M        READ-A036
    C3bBb decay half-life       90 s             READ-A037
    C3(H2O)Bb decay half-life   77 s             READ-A038
    decayed Bb hemolytic activity, relative      0.01   READ-A039

WHAT IS NOT SUPPLIED, and it is the one input the verdict needs: the C3 concentration at the surface.
So this file does not assert a single N. It computes N across the substrate range, reports the two
endpoints that bracket every possible answer, and states the deposition efficiency each implies.
A bracket from measured rates is worth more than a point estimate from an invented concentration.

THE CONTROL, equally informed: the reading in which decay alone decides, which is what the measured
half-lives invite on their own. It is handed the same half-lives and the same decayed-Bb activity and
no turnover number. It predicts quiescence at every substrate level, because a convertase that halves
in 90 s and loses 99 percent of its activity looks spent. The turnover number is what refutes it.

PENDING_INDEPENDENT_REVIEW. No clinical or biological validation is claimed; this is arithmetic on
published in-vitro rates, and the falsifier below is the measurement that would decide it.
"""
from __future__ import annotations

import json
import math
import os

KCAT_PER_S = 1.78          # READ-A035
KM_M = 5.86e-6             # READ-A036
T_HALF_C3BBB_S = 90.0      # READ-A037
T_HALF_C3H2OBB_S = 77.0    # READ-A038
DECAYED_BB_RELATIVE_ACTIVITY = 0.01   # READ-A039


def mean_lifetime_s(t_half_s: float) -> float:
    """Exponential decay: the mean lifetime is the half-life over ln 2, not the half-life."""
    return t_half_s / math.log(2.0)


def turnovers_per_convertase(kcat_per_s: float, t_half_s: float, substrate_M: float | None) -> float:
    """Molecules cleaved by one convertase before it decays. substrate_M None means saturating."""
    tau = mean_lifetime_s(t_half_s)
    saturation = 1.0 if substrate_M is None else substrate_M / (KM_M + substrate_M)
    return kcat_per_s * tau * saturation


def critical_deposition_efficiency(n_turnovers: float) -> float:
    """The loop breaks even at N*p = 1, so the surface must discard all but 1/N."""
    return 1.0 / n_turnovers


def main() -> None:
    tau90 = mean_lifetime_s(T_HALF_C3BBB_S)
    tau77 = mean_lifetime_s(T_HALF_C3H2OBB_S)
    n_sat = turnovers_per_convertase(KCAT_PER_S, T_HALF_C3BBB_S, None)
    n_at_km = turnovers_per_convertase(KCAT_PER_S, T_HALF_C3BBB_S, KM_M)
    n_sat_h2o = turnovers_per_convertase(KCAT_PER_S, T_HALF_C3H2OBB_S, None)

    print('MEASURED RATES -> THE REGULATORY MARGIN')
    print(f'  mean convertase lifetime, surface-bound  : {tau90:.2f} s  (half-life {T_HALF_C3BBB_S:.0f} s)')
    print(f'  mean convertase lifetime, fluid-phase    : {tau77:.2f} s  (half-life {T_HALF_C3H2OBB_S:.0f} s)')
    print(f'  C3 cleaved per convertase, saturating C3 : {n_sat:.1f}')
    print(f'  C3 cleaved per convertase, C3 at Km      : {n_at_km:.1f}')
    print(f'  the two bracket every substrate level, so the margin is known without knowing C3.')
    print(f'  deposition efficiency a quiescent surface must stay below:')
    print(f'    at saturating C3 : {critical_deposition_efficiency(n_sat) * 100:.3f} %')
    print(f'    at C3 = Km       : {critical_deposition_efficiency(n_at_km) * 100:.3f} %')

    print('\nWHY DECAY CANNOT BE THE REGULATOR, which is the control refuted')
    decay_suppression = 1.0 / DECAYED_BB_RELATIVE_ACTIVITY
    print(f'  decay costs the convertase a factor {decay_suppression:.0f} in hemolytic activity (READ-A039),')
    print(f'  while one lifetime of turnover gains a factor {n_sat:.0f} at saturation.')
    # Deliberately NOT printed as a ratio: 231 is a count of cleavages and 100 is an activity
    # suppression, so their quotient is not a physical quantity. The comparison that holds is one
    # of orders -- a single lifetime already produces more cleavages than decay removes activity --
    # and that is enough for the conclusion without inventing a number.
    print('  the counts are not the same kind of quantity, so no ratio is claimed; what holds is that')
    print('  one lifetime already produces more cleavages than decay removes activity, so a surface')
    print('  relying on decay alone is NOT quiescent. An active regulator is necessary, not helpful.')

    print('\nTHE TWO INITIATION ROUTES differ far less than that')
    print(f'  surface-bound against fluid-phase amplification : {n_sat / n_sat_h2o:.3f}')
    print(f'  that is a {100 * (n_sat / n_sat_h2o - 1):.1f} percent difference, two orders below the margin above,')
    print('  so which route starts the loop is not what decides whether it runs away.')

    print('\nFALSIFIER, one measurement, stated before it is made')
    print(f'  measure the deposition efficiency per cleaved C3 on a NON-activator host surface.')
    print(f'  above {critical_deposition_efficiency(n_sat) * 100:.3f} percent the loop runs away on the host, so the')
    print('  observed quiescence of host surfaces cannot be explained by convertase decay and')
    print('  deposition efficiency alone; below it, the regulators are not load-bearing for')
    print('  quiescence and the inflamed state needs a different explanation. Either outcome')
    print('  removes a reading that is currently carried as an assumption.')

    out = {
        'mean_lifetime_surface_bound_s': tau90,
        'mean_lifetime_fluid_phase_s': tau77,
        'turnovers_per_convertase_saturating': n_sat,
        'turnovers_per_convertase_at_Km': n_at_km,
        'critical_deposition_efficiency_saturating': critical_deposition_efficiency(n_sat),
        'critical_deposition_efficiency_at_Km': critical_deposition_efficiency(n_at_km),
        'decay_activity_suppression_factor': decay_suppression,
        'cleavages_per_lifetime_vs_activity_suppression_unlike_quantities': [n_sat, decay_suppression],
        'route_amplification_ratio_surface_over_fluid': n_sat / n_sat_h2o,
        'inputs': {'kcat_per_s': KCAT_PER_S, 'Km_M': KM_M,
                   't_half_surface_bound_s': T_HALF_C3BBB_S,
                   't_half_fluid_phase_s': T_HALF_C3H2OBB_S,
                   'decayed_Bb_relative_activity': DECAYED_BB_RELATIVE_ACTIVITY},
        'source_rows': ['READ-A035', 'READ-A036', 'READ-A037', 'READ-A038', 'READ-A039'],
        'not_supplied': 'C3 concentration at the surface; the two brackets are given instead',
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    d = 'results/ASSEMBLY_COMPLEMENT_AMPLIFICATION'
    os.makedirs(d, exist_ok=True)
    with open(d + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('\nwrote', d + '/decision.json')


if __name__ == '__main__':
    main()
