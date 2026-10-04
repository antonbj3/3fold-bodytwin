"""Anatomical transverse ribbon reduction. No continuum-shell or biological guarantee."""
import numpy as np

def cantilever_stiffness(E, w, h, L, control=False):
    if min(E, w, h, L) <= 0:
        raise ValueError('Positive beam dimensions required')
    if not control:
        return E * w * h ** 3 / (4 * L ** 3)
    ei = E * w * h ** 3 / 12
    ke = ei / L ** 3 * np.array([[12, 6 * L, -12, 6 * L], [6 * L, 4 * L ** 2, -6 * L, 2 * L ** 2], [-12, -6 * L, 12, -6 * L], [6 * L, 2 * L ** 2, -6 * L, 4 * L ** 2]])
    return float(ke[2, 2] - ke[2, 3] * ke[3, 2] / ke[3, 3])

def neighbours(fdi, teeth):
    (q, p) = divmod(fdi, 10)
    prev = fdi - 1 if p > 1 else (q + 1 if q % 2 == 1 else q - 1) * 10 + 1
    return [k for k in [prev, fdi + 1] if str(k) in teeth and k % 10 <= 8]

def response(arch, active, delta, E, h, control=False):
    teeth = arch['teeth']
    a = teeth[str(active)]
    pa = np.array(a['contact_mm'])
    ca = np.array(a['centroid_mm'])
    u = delta * np.array(a['buccal_unit'])
    out = {}
    edges = []

    def add(k, F, M):
        if k not in out:
            out[k] = np.zeros(6)
        out[k] += np.r_[F, M]
    for k in neighbours(active, teeth):
        b = teeth[str(k)]
        pb = np.array(b['contact_mm'])
        cb = np.array(b['centroid_mm'])
        dr = pa - pb
        L = float(np.linalg.norm(dr))
        if L <= 0:
            raise ValueError('Zero beam span')
        w = (a['crown_height_mm'] + b['crown_height_mm']) / 2
        t = dr / L
        P = np.eye(3) - np.outer(t, t)
        stiff = cantilever_stiffness(E, w, h, L, control)
        Fa = -stiff * P @ u
        Fb = -Fa
        root_couple = np.cross(dr, Fb)
        add(active, Fa, np.cross(pa - ca, Fa))
        add(k, Fb, np.cross(pb - cb, Fb) + root_couple)
        edges.append({'anchor_fdi': k, 'span_mm': L, 'width_mm': w, 'k_N_per_mm': stiff, 'root_couple_Nmm': root_couple.tolist(), 'force_on_active_N': Fa.tolist()})
    if not edges:
        return {'status': 'REFUSED_NO_IMMEDIATE_NEIGHBOUR', 'active_fdi': active}
    totalF = np.sum([v[:3] for v in out.values()], axis=0)
    totalM = np.sum([v[3:] + np.cross(teeth[str(k)]['centroid_mm'], v[:3]) for (k, v) in out.items()], axis=0)
    F = out[active][:3]
    return {'status': 'CONDITIONAL_RIBBON_MODEL', 'active_fdi': active, 'activation_mm': delta, 'E_MPa': E, 'thickness_mm': h, 'active_force_norm_N': float(np.linalg.norm(F)), 'active_buccolingual_force_N': float(-F @ np.array(a['buccal_unit'])), 'active_moment_norm_Nmm': float(np.linalg.norm(out[active][3:])), 'wrenches': {str(k): v.tolist() for (k, v) in out.items()}, 'moment_origin': 'per-tooth labelled crown centroid', 'edges': edges, 'force_equilibrium_N': float(np.linalg.norm(totalF)), 'global_torque_equilibrium_Nmm': float(np.linalg.norm(totalM)), 'neighbour_count': len(edges), 'validity': 'Small-strain transverse cantilever; rigid teeth; no contact-opening/PDL/CoR; physical force UNKNOWN'}

def window(interval):
    (lo, hi) = interval
    if hi < 0.5:
        return 'BELOW_PUBLISHED_BODILY_RANGE'
    if lo > 1:
        return 'ABOVE_PUBLISHED_BODILY_RANGE'
    if lo >= 0.5 and hi <= 1:
        return 'WITHIN_PUBLISHED_BODILY_RANGE'
    return 'STRADDLES_PUBLISHED_BODILY_RANGE'
