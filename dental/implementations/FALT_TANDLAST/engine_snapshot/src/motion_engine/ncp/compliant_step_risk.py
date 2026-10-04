"""Conditional global motion uniqueness for a compliant implicit Coulomb step.

Caller obligations: G=J M^-1 J^T, M positive definite, static contact geometry;
normal law u_n+s_i(lambda_n) perpendicular to lambda_n>=0; s_i monotone;
`normal_response_upper` bounds 1/s_i' on the complete range of ALL admissible
roots (not just the returned root); `tangent_gram_upper` bounds
||J_t M^-1/2||_2^2; mu_max bounds every friction coefficient.

The test uses maximal dissipation and normal complementarity, including stick
and open contacts. It guarantees uniqueness of motion IF bounds are valid.
It does not establish existence, unique forces, continuum accuracy, or a
physical rigid-limit selection. Unsupported/infinite bounds return UNKNOWN.
No matrix factorization or root solve. See docs/compliant_step_risk.md.
"""
from __future__ import annotations
import math

def compliant_step_risk(*, mu_max: float, tangent_gram_upper: float, normal_response_upper: float) -> dict:
    """Return a sufficient test; nonacceptance means UNKNOWN, not multiplicity.

    Inputs have units 1, 1/kg, kg respectively. Bounds are supplied by the
    constitutive model/assembly; a local tangent stiffness alone is insufficient.
    The fixed scalar work is per query, after the caller's bound setup/update.
    All positive arithmetic rounds upward to avoid accepting the equality case.
    """
    vals = (mu_max, tangent_gram_upper, normal_response_upper)
    if any((math.isnan(x) or x < 0 for x in vals)):
        raise ValueError('bounds must be nonnegative and not NaN')
    if not math.isfinite(mu_max):
        raise ValueError('mu_max must be finite')
    if any((x == 0 for x in vals)):
        gamma = 0.0
    elif any((math.isinf(x) for x in vals)):
        gamma = math.inf
    else:
        gamma = mu_max
        for factor in (mu_max, tangent_gram_upper, normal_response_upper):
            gamma = math.nextafter(gamma * factor, math.inf)
        gamma = math.nextafter(gamma / 4.0, math.inf)
    accepted = gamma < 1.0
    return {'verdict': 'UNIQUE_MOTION_IF_BOUNDS_VALID' if accepted else 'UNKNOWN', 'risk_ratio_upper': gamma, 'refine_or_verify': not accepted, 'existence_certified': False, 'force_uniqueness_certified': False, 'time_accuracy_certified': False}

def initial_hertz_step_risk(*, h: float, mu_max: float, tangent_gram_upper: float, normal_row_norm_upper: float, gravity_mass_norm_upper: float, hertz_coefficient_upper: float, damping_eta_upper: float) -> dict:
    """Constant-work first-step route from cached geometry/material bounds.

    Scope: initial rest, zero prestrain, static frozen J, nonnegative reference
    gaps, v_free=h*a, and F=K_H*delta^1.5*max(1+eta*delta_dot,0).
    Bounds have units 1/kg, 1/sqrt(kg), sqrt(kg)*m/s², N/m^1.5, s/m.
    No later-step history, moving-boundary work or temporal error is covered.
    """
    if not math.isfinite(h) or h <= 0:
        raise ValueError('h must be finite and strictly positive')
    vals = (normal_row_norm_upper, gravity_mass_norm_upper, hertz_coefficient_upper, damping_eta_upper)
    if any((math.isnan(x) or x < 0 for x in vals)):
        raise ValueError('geometry and constitutive bounds must be nonnegative and not NaN')

    def product(*terms):
        if any((x == 0 for x in terms)):
            return 0.0
        value = 1.0
        for term in terms:
            value = math.nextafter(value * term, math.inf)
        return value
    depth = product(h, h, normal_row_norm_upper, gravity_mass_norm_upper)
    root_depth = math.nextafter(math.sqrt(depth), math.inf) if depth else 0.0
    elastic = product(1.5, h)
    damping = product(4.0, damping_eta_upper, depth)
    bracket = math.nextafter(elastic + damping, math.inf)
    response = product(h, hertz_coefficient_upper, root_depth, bracket)
    result = compliant_step_risk(mu_max=mu_max, tangent_gram_upper=tangent_gram_upper, normal_response_upper=response)
    result.update(normal_response_upper_kg=response, penetration_upper_m=depth, scope='initial_rest_zero_prestrain_static_geometry')
    return result
