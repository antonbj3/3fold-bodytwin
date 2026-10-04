"""Conditional same-problem root proofs plus exact rectangular path coverage.

This composes the reviewed true-residual Krawczyk checker. It never promotes
instability, a failed uniqueness bound or one local box to a global root count.
The producer must supply outward enclosures of the identical complete residual
at each centre and its true Jacobian over each entire box. A problem_id is the
SHA-256 binding source, state, h, domain, full residual and geometry. Identity
agreement is checked here; correctness of that producer binding is a premise.

Geometry is defined here as filled axis-aligned zero-thickness rectangles,
translated affinely in generalized coordinates. Paths are linear from initial
state to every point of the endpoint box. No deforming or nonlinear path is
represented. A separating axis with strict same order at both endpoints
covers the whole path and all endpoint uncertainty, with exact rationals.
"""
from fractions import Fraction as Q
import re
from .krawczyk_step import certify_root_box

def _q(value):
    if isinstance(value, bool):
        raise ValueError('boolean coordinate')
    return Q(value)

def _rectangles(panels, lo, hi):
    n = len(lo)
    out = []
    if not panels:
        raise ValueError('nonempty complete rectangle scene required')
    for panel in panels:
        origin = list(map(_q, panel['origin']))
        size = list(map(_q, panel['size']))
        B = [list(map(_q, row)) for row in panel['translation']]
        if len(origin) != 3 or len(size) != 2 or any((v <= 0 for v in size)):
            raise ValueError('positive rectangular XY embedding required')
        if len(B) != 3 or any((len(row) != n for row in B)):
            raise ValueError('translation map dimensions')
        bounds = []
        for axis in range(3):
            low = origin[axis]
            high = origin[axis] + (size[axis] if axis < 2 else 0)
            for (j, b) in enumerate(B[axis]):
                low += b * (lo[j] if b >= 0 else hi[j])
                high += b * (hi[j] if b >= 0 else lo[j])
            bounds.append((low, high))
        out.append(bounds)
    return out

def translated_rectangle_path(*, initial, lower, upper, panels):
    """Exact sufficient path proof for every endpoint in one closed box."""
    try:
        (x0, lo, hi) = [list(map(_q, v)) for v in (initial, lower, upper)]
        if not x0 or len(x0) != len(lo) or len(lo) != len(hi):
            raise ValueError('path state dimensions')
        if any((a > b for (a, b) in zip(lo, hi))):
            raise ValueError('reversed endpoint box')
        start = _rectangles(panels, x0, x0)
        finish = _rectangles(panels, lo, hi)
        axes = []
        for i in range(len(panels)):
            for j in range(i):
                proof = None
                for axis in range(3):
                    for (a, b) in ((i, j), (j, i)):
                        if start[a][axis][1] < start[b][axis][0] and finish[a][axis][1] < finish[b][axis][0]:
                            proof = (axis, a, b)
                            break
                    if proof is not None:
                        break
                if proof is None:
                    return {'safe': False, 'reason': 'no common strict separating axis over the whole path'}
                axes.append(proof)
        return {'safe': True, 'axes': axes, 'scope': 'FILLED_XY_RECTANGLES_AFFINE_TRANSLATION_LINEAR_PATH'}
    except (KeyError, ValueError, TypeError, ZeroDivisionError, OverflowError) as exc:
        return {'safe': False, 'reason': str(exc)}

def certify_step_roots(*, problem_id, domain_lower, domain_upper, initial, boxes, panels=None):
    """Prove >=2 same-problem roots; preserve UNKNOWN for global uniqueness.

    Each box has problem_id, centre, radius, force_lower/upper,
    jacobian_lower/upper, preconditioner and jacobian_role. Its force data must
    refer to that centre. All coordinates use the same generalized ordering.
    MULTIPLE proves mathematical alternatives. accepted_step is always False:
    this function neither chooses a branch nor verifies integration accuracy.
    """
    result = {'verdict': 'UNKNOWN', 'global_root_count': 'UNKNOWN', 'root_count_lower': 0, 'collision_free_root_count_lower': 0, 'accepted_step': False, 'box_results': [], 'condition': 'TRUE_COMPLETE_RESIDUAL_ENCLOSURES_AND_PRODUCER_PROBLEM_BINDING_VALID'}
    try:
        if not isinstance(problem_id, str) or re.fullmatch('[0-9a-f]{64}', problem_id) is None:
            raise ValueError('SHA-256 problem binding required')
        (lo, hi, x0) = [list(map(_q, v)) for v in (domain_lower, domain_upper, initial)]
        if not lo or len(lo) != len(hi) or len(lo) != len(x0) or any((a >= b for (a, b) in zip(lo, hi))):
            raise ValueError('nonempty bounded domain dimensions')
        accepted = []
        for box in boxes:
            row = {'accepted': False}
            result['box_results'].append(row)
            if box['problem_id'] != problem_id:
                row['reason'] = 'different source/state/h/residual problem'
                continue
            (c, r) = [list(map(_q, box[k])) for k in ('center', 'radius')]
            if len(c) != len(lo) or len(r) != len(lo) or any((v <= 0 for v in r)):
                row['reason'] = 'box dimensions/radius'
                continue
            (blo, bhi) = ([v - d for (v, d) in zip(c, r)], [v + d for (v, d) in zip(c, r)])
            if any((a < d or b > e for (a, b, d, e) in zip(blo, bhi, lo, hi))):
                row['reason'] = 'box outside declared domain'
                continue
            cert = certify_root_box(residual_id=problem_id, radius=r, force_lower=box['force_lower'], force_upper=box['force_upper'], jacobian_lower=box['jacobian_lower'], jacobian_upper=box['jacobian_upper'], preconditioner=box['preconditioner'], jacobian_role=box.get('jacobian_role', 'MISSING'))
            row.update(root_verdict=cert.verdict, contraction_upper=str(cert.contraction_upper), inclusion_upper=str(cert.inclusion_upper))
            if cert.verdict != 'UNIQUE_IN_BOX_IF_ENCLOSURES_VALID':
                row['reason'] = cert.reason
                continue
            if any((not any((b < d or e < a for (a, b, d, e) in zip(blo, bhi, ol, oh))) for (ol, oh) in accepted)):
                row['reason'] = 'not disjoint from previously counted closed boxes'
                continue
            accepted.append((blo, bhi))
            row['accepted'] = True
            result['root_count_lower'] += 1
            path = translated_rectangle_path(initial=x0, lower=blo, upper=bhi, panels=panels)
            row['path'] = path
            if path['safe']:
                result['collision_free_root_count_lower'] += 1
        if result['root_count_lower'] >= 2:
            result['global_root_count'] = 'AT_LEAST_TWO_IF_ENCLOSURES_VALID'
            result['verdict'] = 'PROVEN_MULTIPLE_IF_ENCLOSURES_VALID'
            if result['collision_free_root_count_lower'] >= 2:
                result['verdict'] = 'PROVEN_MULTIPLE_WITH_FREE_PATHS_IF_ENCLOSURES_VALID'
        elif result['root_count_lower'] == 1:
            result['verdict'] = 'LOCAL_UNIQUE_ONLY_IF_ENCLOSURES_VALID'
        return result
    except (KeyError, ValueError, TypeError, ZeroDivisionError, OverflowError) as exc:
        result.update(verdict='UNKNOWN', global_root_count='UNKNOWN', accepted_step=False, root_count_lower=0, collision_free_root_count_lower=0, reason=str(exc))
        return result
