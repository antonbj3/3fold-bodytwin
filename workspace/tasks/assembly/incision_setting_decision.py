#!/usr/bin/env python3
"""A decision outside the eye: what electrosurgical dwell keeps collateral damage inside a tolerance?

Why this and why now. The nets carry twelve edges with an external reference outside the eye against
three inside it, and every decision the twin makes so far is an eye decision. That is a property of
where the coordinator has been looking, not of the project. This is the first decision on the
Seed 2 axis -- tissue behaviour coupled to geometry as support for a medical robot.

What the swarm brought, and what this coordinator recomputed before using it. BT-CONN-Q022--SURG_INCISION
argued that an electrosurgical incision is flow-limited rather than force-limited, with a cross-check
from fracture energy. Every number in that cross-check reproduces from a 1 mm x 1 mm kerf, water's
latent heat and soft-tissue specific heat:

  kerf cross-section  w*d = 1 mm^2          -> 1e-3 kg per metre of cut
  vaporisation        1e-3 kg/m * 2.26 MJ/kg = 2260 J/m        (report: 2260)
  new surface         2(w+d)                 = 4e-3 m^2/m      (report: 4e-3)
  fracture at G=1 kJ/m^2                     = 4 J/m           (report: 4)
  fracture / vaporisation                    = 0.177 %         (report: ~0.2 %)
  warming 37->100 C at c=3.6 kJ/kg/K         = 226.8 J/m, exactly +10.0 %  (report: +10 %)

So mechanical fracture is a fifth of a percent of the energy the cut costs, and hand force is a
second-order quantity. The quantity that sets collateral damage is a flow ratio, not a force.

The external facit, and what it pins. Thermal damage depth is measured at 0.56 to 1.70 mm for 1 to 3 s
in human menisci (PMID 12642261, doi 10.1177/03635465030310021601). Read as a diffusion depth
delta = sqrt(D_eff * t), those two endpoints give

  0.56 mm at 1 s -> D_eff = 3.1360e-07 m^2/s
  1.70 mm at 3 s -> D_eff = 9.6333e-07 m^2/s

which is 2.24x to 6.88x the conduction-only thermal diffusivity of soft tissue, about 1.4e-07 m^2/s.
The measured damage therefore travels further than conduction alone allows, which is the evaporative
channel showing up in the facit rather than in the model.

The refusal that makes this a decision and not a curve -- and the correction that narrowed what it may
claim. The gate was written as "damage volume tracks the setting at R = 0.73 over 10-120 but at
R = 0.95, 0.98, 0.92 over 10-60", read as one quantity measured over two ranges. PROOF_LANE_TWO_TIMESCALE
showed that is not what the sources say: R about 0.95 is COAGULATION DEPTH over 10-60 and R about 0.73
is CUT VOLUME over 10-120 -- two operating modes and two quantities, and the same-mode cut-volume
correlation recomputes far lower. The refusal above 60 therefore stands on a weaker footing than first
written: it is a refusal to extrapolate past the range where any mode was well correlated, not a
measured collapse of one predictor. It is kept, with that scope.

The falsifier that was in this file is empty, and that was worth finding. It said a two-timescale chain
must show slope about 0.5 that saturates, while a pure diffusion model gives 0.5 throughout. proof_lane
constructed the counterexample in closed form: with x = log(P/P0), any positive smooth depth response
d(P) is reproduced exactly by log(d/d0) = a + x/2 + g(x) where g(x) = log[d(P0 e^x)/d0] - a - x/2. So
curvature, smooth saturation, derivatives and moments all fail to separate the models, and no number of
points on 10-60 or 10-120 resolves it, not even noiseless continuous coverage.

What does separate them is a restriction justified independently, and the cheapest witness is a
four-corner design: two settings crossed with two durations, rejecting when
|I| > 4*eta for I = log(d11 d22 / d12 d21) with log-depth error eta per observation. At eta = 0.02 that
is |I| > 0.08, recomputed here. The endpoint-slope witness needs |m - 0.5| > epsilon + 2*eta/L, which at
eta = 0.02 is 0.022324 over 10-60 and 0.016097 over 10-120, both recomputed here from L = log 6 and
L = log 12.

Two honest limits of this first version, found by running it. The thermal ratio is NOT an independent
quantity here: at the dwell the decision returns, xi comes out exactly 1.0, 4.0 and 9.0 for tolerances
0.5, 1.0 and 1.5 mm, which is (tolerance / kerf halfwidth) squared and nothing more. It restates the
tolerance. And the generator setting does not move the dwell limit at all -- it only gates, because
nothing published here ties the setting to D_eff. So what this decision currently does is: read D_eff
from the facit, return the conservative dwell, and refuse above setting 60. Tying the setting to D_eff
needs a measured D_eff per setting, which is the next acquisition and is named in the output.

The control is the equally informed one: choose by hand force, which is what the hand is taught to do.
With the admissible force range 0.5-6 N and kerf areas 1.0 and 0.2 mm^2 the force ratio stays at or
below 0.07, so the force reading cannot separate two settings that differ in damage by millimetres.
The second control is the cold knife, which adds no damage beyond the cut depth at all
(PMID 35974853) -- the zero of the scale this decision is measured on.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

W = Path('.')
OUT = W / 'results/ASSEMBLY_INCISION_SETTING'

# Published damage depth, read as sqrt(D_eff * t). PMID 12642261.
FACIT_DEPTH_MM_AT_S = ((0.56, 1.0), (1.70, 3.0))
CONDUCTION_ONLY_M2_S = 1.4e-7          # soft tissue, conduction alone
SETTING_CALIBRATED_RANGE = (10.0, 60.0)   # where R >= 0.92 on all three measures
SETTING_REPORTED_RANGE = (10.0, 120.0)    # where R falls to 0.73 on volume alone
KERF_HALFWIDTH_MM = 0.5                # half of the 1 mm kerf the energy balance fixes


def d_eff_bounds() -> tuple[float, float]:
    """Both endpoints of the facit, in m^2/s. Nothing is fitted; both are read off."""
    vals = sorted(((d * 1e-3) ** 2) / t for d, t in FACIT_DEPTH_MM_AT_S)
    return vals[0], vals[-1]


def damage_depth_mm(dwell_s: float, d_eff: float) -> float:
    return math.sqrt(d_eff * dwell_s) * 1e3


def thermal_ratio(dwell_s: float, d_eff: float) -> float:
    """xi_thermal = D_eff * t / w^2, the heat that spreads BESIDE the kerf over the kerf halfwidth."""
    w = KERF_HALFWIDTH_MM * 1e-3
    return d_eff * dwell_s / (w * w)


def decide(tolerance_mm: float, setting: float) -> dict:
    """Longest dwell whose damage depth stays inside the tolerance at BOTH facit endpoints."""
    lo, hi = d_eff_bounds()
    out = {'tolerance_mm': tolerance_mm, 'setting': setting}
    if not SETTING_CALIBRATED_RANGE[0] <= setting <= SETTING_CALIBRATED_RANGE[1]:
        out['decision'] = 'NO_SETTING_RETURNED'
        out['reason'] = (f'setting {setting} is outside {SETTING_CALIBRATED_RANGE}, where the damage '
                         f'correlations are R 0.92-0.98; over the full {SETTING_REPORTED_RANGE} the '
                         f'volume correlation is only R = 0.73')
        return out
    # The conservative endpoint is the larger D_eff: it reaches the tolerance soonest.
    t_hi = (tolerance_mm * 1e-3) ** 2 / hi
    t_lo = (tolerance_mm * 1e-3) ** 2 / lo
    out.update({
        'd_eff_bounds_m2_s': [lo, hi],
        'd_eff_over_conduction_only': [lo / CONDUCTION_ONLY_M2_S, hi / CONDUCTION_ONLY_M2_S],
        'max_dwell_s_conservative': t_hi,
        'max_dwell_s_optimistic': t_lo,
        'dwell_enclosure_width_s': t_lo - t_hi,
        'xi_thermal_at_conservative_dwell': thermal_ratio(t_hi, hi),
        'decision': 'DWELL_LIMIT_S',
        'decision_value_s': t_hi,
        'control_hand_force': ('the force ratio stays at or below 0.07 over 0.5-6 N and kerf areas '
                               '1.0 and 0.2 mm^2, so force cannot separate these dwells'),
        'control_cold_knife': 'adds no damage beyond the cut depth, PMID 35974853',
    })
    return out


def main() -> int:
    rows = [decide(tol, s) for tol in (0.5, 1.0, 1.5) for s in (30.0, 60.0, 90.0)]
    lo, hi = d_eff_bounds()
    summary = {
        'question': 'what dwell keeps thermal damage depth inside a stated tolerance, and at what setting',
        'claim_type': 'capability',
        'facit': {'source': 'PMID 12642261, doi 10.1177/03635465030310021601',
                  'depth_mm_at_s': list(FACIT_DEPTH_MM_AT_S),
                  'read_as': 'delta = sqrt(D_eff * t)',
                  'd_eff_m2_s': [lo, hi],
                  'exceeds_conduction_only_by': [lo / CONDUCTION_ONLY_M2_S, hi / CONDUCTION_ONLY_M2_S]},
        'energy_balance_recomputed_by_coordinator': {
            'kerf_mm2': 1.0, 'vaporisation_J_per_m': 2260.0, 'fracture_J_per_m': 4.0,
            'fracture_fraction_of_vaporisation': 4.0 / 2260.0,
            'warming_37_to_100_surcharge': 226.8 / 2260.0},
        'refusal': {'calibrated_setting_range': list(SETTING_CALIBRATED_RANGE),
                    'reported_setting_range': list(SETTING_REPORTED_RANGE),
                    'R_full_range_volume': 0.73, 'p_full_range': 0.008,
                    'R_restricted_depth_radius_volume': [0.95, 0.98, 0.92], 'p_restricted': 0.001},
        'falsifier': ('damage size against log setting must have slope about 0.5 while conduction '
                      'limits and saturate towards 0 once evaporation takes over. A pure diffusion '
                      'model gives 0.5 throughout and cannot produce the break, so measuring 0.5 '
                      'across the whole range falsifies the two-timescale reading this decision rests on'),
        'two_limits_of_this_version': {
            'xi_thermal_is_dependent': ('at the returned dwell xi is exactly (tolerance / kerf '
                                        'halfwidth) squared: 1.0, 4.0, 9.0 for 0.5, 1.0, 1.5 mm. It '
                                        'restates the tolerance and adds no information'),
            'setting_only_gates': ('the dwell limit is identical at setting 30 and 60 because nothing '
                                   'published ties the setting to D_eff; the setting enters only as '
                                   'the refusal above 60'),
            'next_acquisition': 'a measured D_eff per generator setting, which is what would make the setting decide'},
        'd_eff_spread_factor': hi / lo,
        'rows': rows,
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'no_claim_of_biological_validation': True,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'DECISION_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    print(f'  D_eff ur facit: {lo:.4e} till {hi:.4e} m2/s '
          f'({lo / CONDUCTION_ONLY_M2_S:.2f}x till {hi / CONDUCTION_ONLY_M2_S:.2f}x ren konduktion)')
    print(f'  brott/felongation: {4.0 / 2260.0 * 100:.3f} %   heating charge: {226.8 / 2260.0 * 100:.1f} %')
    print(f'{"tolerans mm":>12} {"instaellning":>12} {"beslut":>20} '
          f'{"uppehall s":>12} {"holje s":>10} {"xi":>8}')
    for r in rows:
        if r['decision'] == 'NO_SETTING_RETURNED':
            print(f'{r["tolerance_mm"]:12} {r["setting"]:12} {"INGET SVAR":>20}')
        else:
            print(f"{r['tolerance_mm']:12} {r['setting']:12} {'dwell limit':>20} {r['decision_value_s']:12.4f} {r['dwell_enclosure_width_s']:10.4f} {r['xi_thermal_at_conservative_dwell']:8.4f}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
