#!/usr/bin/env python3
"""Extend the surgical decision from one variable to three: sphere, cylinder and axis.

Why. The power decision used nine fields of a 54-column clinical dataset and the 20 spherical cases.
The other 69 patients received a TORIC lens, where the decision has three parts — how much power, how
much cylinder, and at what axis — and they were excluded for want of a vector prediction.

Where the value is. A toric plan needs the total corneal astigmatism, anterior plus posterior. Standard
calculators ESTIMATE the posterior contribution from a population regression, because it is usually not
measured. This dataset measures it: anterior radii with their axis and posterior radii with their axis,
from a tomographer. So the twin can use the measured posterior cornea where practice uses an estimate,
and the two can be scored against the same patients' measured postoperative refraction. That is an
information-link claim in its strict form: the control is current practice WITHOUT the new information.

Method. Each corneal surface is a sphero-cylinder, and sphero-cylinders at different axes add in power
vector form (M, J0, J45), which is exact for thin-element combination:
    M = (S1 + S2)/2 ... with cylinder C at axis A: J0 = -(C/2)cos(2A), J45 = -(C/2)sin(2A).
The combined cornea is converted back to principal meridians, each meridian is carried through the same
thick paraxial chain the power decision uses, and the two meridional refractions become the predicted
sphere, cylinder and axis. The implanted toric cylinder is applied at its label axis.

Honest limits, inherited and new. The chain keeps the preoperative geometry frozen after surgery and a
modelled lens design rather than the implanted serial-number lens. The surgically induced astigmatism of
the incision is not modelled. Toric lens rotation after implantation is not modelled and is a known
clinical error source. Nothing here is a clinical recommendation.
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path

R15 = Path('results/LANE_EYE_OPTICAL_TWIN/r15')
OUT = Path('results/ASSEMBLY_TORIC_DECISION')
sys.path.insert(0, str(R15))
sys.path.insert(0, str(R15.parent))
from hydration_r4 import material, H0                      # our own stromal index

N_AQUEOUS = 1.336
VERTEX_M = 0.012
# Population posterior-astigmatism estimate standing in for a measurement, which is what a calculator
# without tomography does: posterior cornea adds against-the-rule cylinder of about 0.3 D at 90 deg.
#
# THE ADVANTAGE THIS CHAIN REPORTS IS NOT ROBUST TO THIS PAIR, measured 2026-10-04. The stand-in
# (-0.30 D at 90 deg) has power-vector J0 = -0.15, which is against-the-rule and is what the
# literature describes. Rotating it to axis 0 while keeping the sign gives J0 = +0.15, and on the
# same 69 eyes the population arm then falls from 0.6209 to 0.3636 D while the measured arm stays at
# 0.2808 D. So the benefit of measuring the posterior cornea is 0.3401 D under one convention and
# 0.0828 D under the other, and the data PREFER the convention that shrinks it. The lane
# LANE_CORNEA_SHAPE reached the same conclusion independently with its own field model (0.5874 ->
# 0.3541 D) and labels the axis-90 variant MISORIENTED.
#
# Nothing here settles which orientation a calculator should use; that needs a cited nomogram rather
# than a constant in this file. Until it has one, the headline must be quoted WITH the sensitivity,
# which main() now prints, and the 0.34 D figure must not be used alone.
POP_POSTERIOR_CYL_D = 0.30
POP_POSTERIOR_AXIS_DEG = 90.0
POP_POSTERIOR_AXIS_ALTERNATIVE_DEG = 0.0


def to_vector(sphere: float, cyl: float, axis_deg: float) -> tuple[float, float, float]:
    a = math.radians(axis_deg)
    return (sphere + cyl / 2.0, -(cyl / 2.0) * math.cos(2 * a), -(cyl / 2.0) * math.sin(2 * a))


def from_vector(M: float, J0: float, J45: float) -> tuple[float, float, float]:
    c = 2.0 * math.hypot(J0, J45)
    # 3/10 22:30, found by a Sol connection hunt and verified here: atan2(-J45, -J0) negates both
    # components, which rotates the double-angle vector by 180 degrees and therefore the AXIS by
    # 90. The round trip was wrong in every case: 0 came back as 90, 45 as 135, 30 as 120. With
    # J0 = -(C/2)cos(2A) and J45 = -(C/2)sin(2A) and C negative, (J0, J45) = |C/2|(cos 2A, sin 2A),
    # so 2A = atan2(J45, J0) with no negation. The sphere and cylinder always round-tripped
    # correctly, which is why the residual never showed it: the lens is aligned to this same
    # returned axis, so a common 90-degree rotation cancels in the residual and only an ABSOLUTE
    # axis reported out of here was wrong.
    axis = 0.5 * math.degrees(math.atan2(J45, J0)) % 180.0
    return (M + c / 2.0, -c, axis)                           # minus-cylinder convention


def surface_power(r1_mm: float, r2_mm: float, axis_deg: float, n_from: float, n_to: float):
    """A refracting surface with two principal radii, as a sphero-cylinder in power vector form.

    The steeper meridian of a surface lies 90 degrees from the flat one; the reported axis is the flat
    meridian, so the cylinder sits on the steep one.
    """
    p1 = (n_to - n_from) / (r1_mm / 1000.0)
    p2 = (n_to - n_from) / (r2_mm / 1000.0)
    sphere, cyl = p1, p2 - p1
    return to_vector(sphere, cyl, axis_deg)


def meridional_refraction(row: dict, cornea_sphere: float, cornea_cyl: float, cornea_axis: float,
                          iol_sphere: float, iol_cyl: float, iol_axis: float):
    """Carry each principal meridian through the chain that is ALREADY VERIFIED, not a new one.

    My first version wrote its own vergence recursion and produced a residual cylinder of 5.2 D, where
    real residual astigmatism after a toric lens is under 1 D. That version is withdrawn: the absolute
    scale was wrong, so no relative comparison resting on it could be reported. The fix is not to debug
    my arithmetic but to reuse `clinical_anchor.predict`, which carries the same chain and was scored
    against 20 patients with a mean absolute error of 0.47 D.

    `predict` takes the mean anterior curvature from R1 and R2 and the posterior from rpmean, so a
    meridian is obtained by handing it a row whose radii ARE that meridian's radii. The surface powers
    are converted back to radii with our own stromal index, which keeps the chain's own convention.
    """
    from clinical_anchor import predict
    nc = material(H0)['n_stroma']

    def project(s, c, a, theta_deg):
        M, J0, J45 = to_vector(s, c, a)
        th = math.radians(theta_deg)
        return M + J0 * math.cos(2 * th) + J45 * math.sin(2 * th)

    # THREE meridians, not two. With only 0 and 90 the J45 component is unrecoverable and every
    # predicted axis snaps to 0 or 90, which is what the first run did while real axes are oblique
    # (one patient measured 25 deg). Sampling 0, 60 and 120 determines (M, J0, J45) exactly.
    out = []
    for meridian in (0.0, 60.0, 120.0):
        # The combined cornea's power in this meridian, expressed as an equivalent anterior radius with
        # the posterior surface held at its mean, so the chain sees one consistent meridian.
        k = project(cornea_sphere, cornea_cyl, cornea_axis, meridian)
        r_post = row['IOLM_rpmean']
        k_post = (N_AQUEOUS - nc) / (r_post / 1000.0)
        k_ant = k - k_post
        r_ant_mm = (nc - 1.0) / k_ant * 1000.0
        p = project(iol_sphere, iol_cyl, iol_axis, meridian)
        r = dict(row)
        r['IOLM_R1'] = r_ant_mm
        r['IOLM_R2'] = r_ant_mm
        r['IOLP'] = p
        out.append(predict(r, 'thick')[0]['refraction_infinity_D'])
    # Least squares on P(th) = M + J0*cos(2th) + J45*sin(2th) for the three sampled meridians.
    import numpy as _np
    A = _np.array([[1.0, math.cos(2 * math.radians(t)), math.sin(2 * math.radians(t))]
                   for t in (0.0, 60.0, 120.0)])
    M, J0, J45 = _np.linalg.solve(A, _np.array(out))
    return from_vector(float(M), float(J0), float(J45))


def main() -> int:
    import openpyxl
    path = next(R15.glob('*.xlsx'))
    rows = list(openpyxl.load_workbook(path, data_only=True).active.values)
    headers = [str(h) for h in rows[1]]
    need = ['IOLM_CCT', 'IOLM_AL', 'IOLM_R1', 'IOLM_R2', 'IOLM_AXIS', 'IOLM_PR1', 'IOLM_PR2',
            'IOLM_PRAXIS', 'IOLP', 'IOLT', 'CASIA_POR_AQD post', 'sph', 'cyl', 'Achse']
    recs = []
    for n, raw in enumerate(rows[2:], 3):
        d = dict(zip(headers, raw))
        if not all(isinstance(d.get(k), (int, float)) for k in need if k not in ('Achse',)):
            continue
        if not d['IOLT']:
            continue                                          # toric rows only
        d['sheet_row'] = n
        recs.append(d)
    if not recs:
        print('no toric rows with all required fields'); return 1

    nc = material(H0)['n_stroma']
    results = []
    for d in recs:
        ant = surface_power(d['IOLM_R1'], d['IOLM_R2'], d['IOLM_AXIS'], 1.0, nc)
        post_measured = surface_power(d['IOLM_PR1'], d['IOLM_PR2'], d['IOLM_PRAXIS'], nc, N_AQUEOUS)
        post_estimated = to_vector(
            (N_AQUEOUS - nc) / (((d['IOLM_PR1'] + d['IOLM_PR2']) / 2) / 1000.0),
            -POP_POSTERIOR_CYL_D, POP_POSTERIOR_AXIS_DEG)

        measured_sph, measured_cyl, measured_axis = d['sph'], d['cyl'], d.get('Achse')
        row = dict(sheet_row=d['sheet_row'], implanted_sphere_D=d['IOLP'],
                   implanted_cyl_D=d['IOLT'], measured_sph_D=measured_sph,
                   measured_cyl_D=measured_cyl, measured_axis_deg=measured_axis)
        post_estimated_alt = to_vector(
            (N_AQUEOUS - nc) / (((d['IOLM_PR1'] + d['IOLM_PR2']) / 2) / 1000.0),
            -POP_POSTERIOR_CYL_D, POP_POSTERIOR_AXIS_ALTERNATIVE_DEG)
        for name, post in (('measured_posterior', post_measured), ('population_estimate', post_estimated),
                           ('population_estimate_alt', post_estimated_alt)):
            cs, cc, ca = from_vector(*[a + b for a, b in zip(ant, post)])
            # A toric lens is implanted ALIGNED TO THE STEEP CORNEAL MERIDIAN, so its cylinder axis
            # is the combined cornea's axis and not zero. With axis 0 the lens cylinder adds to the
            # cornea's instead of cancelling it, which is what inflated the first residual to 5.2 D and
            # the second to 2.4 D. The label cylinder is specified at the lens plane; the chain applies
            # it there, which is where predict() places the lens.
            # The sign of the lens cylinder relative to the cornea is a notation question, and getting
            # it wrong makes the lens ADD astigmatism instead of cancelling it: the first runs predicted
            # residuals of 1.6 to 2.5 D against corneal cylinders of 0.9 to 1.5 D, i.e. larger than the
            # defect being corrected, which is physically impossible for a correcting lens. Both
            # conventions are evaluated and the one whose residual is smaller than the corneal cylinder
            # is the physical one; the choice is recorded per patient rather than assumed.
            cands = {}
            for tag, lens_cyl in (('minus_at_corneal_axis', -abs(d['IOLT'])),
                                  ('plus_at_corneal_axis', abs(d['IOLT']))):
                cands[tag] = meridional_refraction(d, cs, cc, ca, d['IOLP'], lens_cyl, ca)
            tag = min(cands, key=lambda t: abs(cands[t][1]))
            pred = cands[tag]
            row[name] = dict(lens_cyl_convention=tag,
                             corneal_sphere_D=round(cs, 4), corneal_cyl_D=round(cc, 4),
                             corneal_axis_deg=round(ca, 2),
                             predicted_sph_D=round(pred[0], 4), predicted_cyl_D=round(pred[1], 4),
                             predicted_axis_deg=round(pred[2], 2))
        results.append(row)

    def cyl_error(name):
        return [abs(r[name]['predicted_cyl_D'] - r['measured_cyl_D']) for r in results]

    def vector_error(name):
        """Length of the difference between predicted and measured cylinder as a double-angle vector.

        Magnitude error alone scores |predicted cyl| against |measured cyl| and is blind to the axis,
        so a prediction with the right amount of cylinder at the wrong meridian scores perfectly.
        LANE_CORNEA_SHAPE r24 recomputed the same 69 eyes with both metrics and the ARM ORDERING
        CHANGES between them, so reporting one without the other picks a winner by choice of metric.
        The vector form is (J0, J45) = (-(C/2)cos2A, -(C/2)sin2A) and the difference is reported as a
        cylinder, i.e. twice the vector length, which is the usual convention.
        """
        out = []
        for r in results:
            ax_m = r['measured_axis_deg']
            if ax_m is None:
                continue
            _, j0p, j45p = to_vector(0.0, r[name]['predicted_cyl_D'], r[name]['predicted_axis_deg'])
            _, j0m, j45m = to_vector(0.0, r['measured_cyl_D'], float(ax_m))
            out.append(2.0 * math.hypot(j0p - j0m, j45p - j45m))
        return out

    summary = dict(
        toric_patients=len(results),
        corneal_cyl_measured_vs_estimated_mean_abs_difference_D=round(statistics.mean(
            abs(r['measured_posterior']['corneal_cyl_D'] - r['population_estimate']['corneal_cyl_D'])
            for r in results), 4),
        residual_cyl_mae_with_measured_posterior_D=round(statistics.mean(cyl_error('measured_posterior')), 4),
        residual_cyl_mae_with_population_estimate_D=round(statistics.mean(cyl_error('population_estimate')), 4),
        measured_posterior_closer_count=sum(
            1 for a, b in zip(cyl_error('measured_posterior'), cyl_error('population_estimate')) if a < b),
        population_estimate_closer_count=sum(
            1 for a, b in zip(cyl_error('measured_posterior'), cyl_error('population_estimate')) if b < a),
        # The same comparison with the stand-in rotated to the alternative axis. It is in the summary
        # rather than in a note because the advantage is a quarter of its headline size under the
        # other convention, and a reader who sees only the headline is being misled by omission.
        residual_cyl_mae_with_population_estimate_alternative_axis_D=round(
            statistics.mean(cyl_error('population_estimate_alt')), 4),
        measured_advantage_D=round(statistics.mean(cyl_error('population_estimate'))
                                   - statistics.mean(cyl_error('measured_posterior')), 4),
        measured_advantage_alternative_axis_D=round(
            statistics.mean(cyl_error('population_estimate_alt'))
            - statistics.mean(cyl_error('measured_posterior')), 4),
        eyes_with_measured_axis=len(vector_error('measured_posterior')),
        residual_vector_mae_with_measured_posterior_D=round(
            statistics.mean(vector_error('measured_posterior')), 4),
        residual_vector_mae_with_population_estimate_D=round(
            statistics.mean(vector_error('population_estimate')), 4),
        residual_vector_mae_with_population_estimate_alternative_axis_D=round(
            statistics.mean(vector_error('population_estimate_alt')), 4),
        control='current practice: the same chain with a population estimate of posterior astigmatism',
        claim_type='information_link',
        not_modelled=['surgically induced astigmatism of the incision',
                      'toric lens rotation after implantation',
                      'postoperative corneal geometry change'],
        review_state='PENDING_INDEPENDENT_REVIEW',
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'TORIC_V1.json').write_text(json.dumps(dict(summary=summary, rows=results),
                                                  indent=1, ensure_ascii=False))
    for k, v in summary.items():
        if isinstance(v, (int, float)):
            print(f'  {k} = {v}')
    print(f'  written: {OUT / "TORIC_V1.json"}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
