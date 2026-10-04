from fractions import Fraction as Q
from copy import deepcopy
import pytest
from field_engine.experimental import normal_cone as c
N = [[1, 0], [-1, 0], [0, 1], [0, -1]]

def certificates():
    out = []
    for axis in range(2):
        for sign in (-1, 1):
            normal_id = 2 * axis + (1 if sign == 1 else 0)
            out.append({'axis': axis, 'sign': sign, 'rows': [normal_id, 8], 'y': [1, 1]})
    return out

def test_no_nonzero_and_input_binding():
    result = c.verify_no_direction(N, certificates())
    assert result.status == 'NO_NONZERO_DIRECTION'
    assert result.checked_charts == 4
    changed = deepcopy(N)
    changed[1] = [1, 0]
    assert c.verify_no_direction(changed, certificates()).status == 'UNKNOWN'

def test_antipodal_normals_have_boundary_direction():
    assert c.verify_direction([[1, 0], [-1, 0]], [0, 1]).status == 'YES_BOUNDARY'
    assert c.verify_direction([[1, 0], [-1, 0]], [1, 0]).status == 'UNKNOWN'
    assert c.verify_direction([[1, 0], [0, 1]], [1, 1]).status == 'YES_STRICT'

def test_zero_and_empty_cannot_certify():
    assert c.verify_direction(N, [0, 0]).status == 'UNKNOWN'
    assert c.verify_no_direction(N, []).status == 'UNKNOWN'
    assert c.verify_direction([], [1, 0]).status == 'UNKNOWN'

def test_duplicate_missing_and_forged_sparse_rows():
    cc = certificates()
    cc[-1] = deepcopy(cc[0])
    assert c.verify_no_direction(N, cc).status == 'UNKNOWN'
    cc = certificates()
    cc[0]['y'][0] = 2
    assert c.verify_no_direction(N, cc).status == 'UNKNOWN'
    cc = certificates()
    cc[0]['rows'] = [8, 8]
    assert c.verify_no_direction(N, cc).status == 'UNKNOWN'

def test_row_envelope_requires_margin():
    assert c.verify_direction([[1, 0]], [1, 1], relaxation='1/2').status == 'YES_STRICT'
    assert c.verify_direction([[1, 0]], [1, 3], relaxation='1/2').status == 'UNKNOWN'
    assert c.verify_direction([[1, 0]], [1, 2], relaxation='1/2').status == 'YES_BOUNDARY'

def test_chart_bounds_and_relaxation():
    (A, b) = c.cone_chart([[1, 2], [3, 4]], 1, -1, relaxation='1/7')
    assert A[0] == (Q(-1), Q(-2))
    assert A[-1] == (Q(0), Q(1))
    assert b == (Q(1, 7), Q(1, 7), Q(1), Q(1), Q(1), Q(1), Q(-1))

@pytest.mark.parametrize('axis,sign', [(True, 1), (0, True), (2, 1), (0, 0)])
def test_invalid_chart(axis, sign):
    with pytest.raises(ValueError):
        c.cone_chart(N, axis, sign)

def test_approximate_inputs_are_not_exact():
    assert c.verify_direction([[0.1, 0]], [1, 0]).status == 'UNKNOWN'
