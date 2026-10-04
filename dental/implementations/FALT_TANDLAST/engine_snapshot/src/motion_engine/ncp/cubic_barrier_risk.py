"""Exact, domain-conditional uniqueness gate for a frozen contact operator.

Updated cubic force: f(g)=2*(eta/g**2+e)*(w-g)_+**2/w.
eta=m_eff*c_ref**2 is an energy in joules, e is stiffness in N/m.
The source intentionally freezes stiffness in each local Newton block and
updates it at the next assembly. Its block is not the total derivative of
this updated force; the difference is a semi-implicit design choice.
Inputs must bound ALL admissible solutions in the declared gap domain, with
fixed affine gap Jacobians, SPD mass, and monotone bounded Coulomb/clamp
friction. This function does not prove existence or a global gap floor.
"""
from fractions import Fraction

def _q(value):
    if isinstance(value, bool):
        raise ValueError('boolean is not a numeric bound')
    try:
        return Fraction(value)
    except (TypeError, ValueError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError('finite rational-compatible bounds required') from exc

def cubic_barrier_step_risk(*, gap_lower, barrier_width_upper, gap_mass_energy_upper, elastic_stiffness_upper, timestep, mu_upper, tangent_gram_upper):
    """Return AT_MOST_ONE_ON_DECLARED_GAP_DOMAIN or UNKNOWN, exactly.

    SI units: gap/width m; eta J; elastic N/m; timestep s; mu dimensionless;
    tangent Gram kg^-1. Fractions are evaluated exactly; float inputs refer
    to their exact represented values and must already be conservative bounds.
    No numerical outward rounding of upstream data is supplied here.
    For many contacts, gap/width/eta/e bound every contact separately and the
    Gram bound is the operator norm of the whole tangential Gram matrix.
    """
    (delta, w, eta, e, h, mu, t2) = map(_q, (gap_lower, barrier_width_upper, gap_mass_energy_upper, elastic_stiffness_upper, timestep, mu_upper, tangent_gram_upper))
    if delta < 0 or w <= 0 or eta < 0 or (e < 0) or (h <= 0) or (mu < 0) or (t2 < 0):
        raise ValueError('invalid gap, coefficient, timestep or operator bound')
    out = {'status': 'UNKNOWN', 'scope': 'DECLARED_GAP_DOMAIN', 'existence_certified': False, 'all_solution_gap_floor_certified': False, 'geometry_frozen': True, 'time_accuracy_certified': False}
    if mu == 0 or t2 == 0:
        out.update(status='AT_MOST_ONE_ON_DECLARED_GAP_DOMAIN', risk_ratio_exact='0', curvature_upper_exact=None, reason='no tangential friction coupling; monotone normal law')
        return out
    if delta == 0:
        out['reason'] = 'positive all-solution gap floor required'
        return out
    k = Fraction(0) if delta >= w else 4 * eta * (w - delta) / delta ** 3 + 4 * e * (w - delta) / w
    gamma = mu ** 2 * t2 * h ** 2 * k / 4
    out.update(curvature_upper_exact=str(k), risk_ratio_exact=str(gamma), status='AT_MOST_ONE_ON_DECLARED_GAP_DOMAIN' if gamma < 1 else 'UNKNOWN', reason='strict two-solution bound' if gamma < 1 else 'bound does not exclude multiplicity')
    return out
