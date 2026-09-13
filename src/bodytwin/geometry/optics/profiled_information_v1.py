"""Local information after profiling explicitly supplied nuisance directions.

Inputs must already include noise whitening and desired parameter scales. The
function performs no calibration, data selection, file writes or input mutation.
Ranks are numerical SVD ranks, not certified identifiability over an interval.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class ProfiledInformation:
    projected: np.ndarray
    gram: np.ndarray
    singular_values: np.ndarray
    nuisance_rank: int
    interest_rank: int
    nuisance_cutoff: float
    interest_cutoff: float


def profile_information(interest, nuisance):
    """Profile nuisance columns; input shapes are (observations, parameters).

    Zero nuisance columns are permitted. Interest must have at least one column.
    Numerical-rank cutoffs use float64 epsilon, largest dimension and spectral
    norm. Interest rank uses the original interest norm so pure projection noise
    is not reported as identifiable signal. Returned raw arrays are not clipped.
    """
    interest = np.array(interest, dtype=np.float64, copy=True)
    nuisance = np.array(nuisance, dtype=np.float64, copy=True)
    if interest.ndim != 2 or nuisance.ndim != 2:
        raise ValueError('Jacobians must be two-dimensional')
    if not interest.shape[0] or not interest.shape[1] or nuisance.shape[0] != interest.shape[0]:
        raise ValueError('Nonempty interest and matching observation rows required')
    if not np.isfinite(interest).all() or not np.isfinite(nuisance).all():
        raise ValueError('Finite Jacobians required')
    epsilon = np.finfo(np.float64).eps
    nuisance_cutoff = 0.
    nuisance_rank = 0
    projected = interest.copy()
    if nuisance.shape[1]:
        u, singular, _ = np.linalg.svd(nuisance, full_matrices=False)
        nuisance_cutoff = float(epsilon*max(nuisance.shape)*singular[0])
        nuisance_rank = int(np.count_nonzero(singular > nuisance_cutoff))
        basis = u[:, :nuisance_rank]
        projected -= basis@(basis.T@interest)
    spectrum = np.linalg.svd(projected, compute_uv=False)
    spectrum = np.pad(spectrum, (0, interest.shape[1]-len(spectrum)))
    interest_cutoff = float(epsilon*max(interest.shape)*np.linalg.norm(interest, ord=2))
    interest_rank = int(np.count_nonzero(spectrum > interest_cutoff))
    gram = projected.T@projected
    for array in (projected, spectrum, gram):
        array.setflags(write=False)
    return ProfiledInformation(projected, gram, spectrum, nuisance_rank, interest_rank,
                               nuisance_cutoff, interest_cutoff)


def selftest():
    arrays = []
    rows = []
    for _ in range(2):
        interest = np.array([[3.,0.],[0.,2.],[0.,0.]])
        nuisance = np.array([[0.],[0.],[1.]])
        before = interest.tobytes(), nuisance.tobytes()
        base = profile_information(interest, nuisance)
        redundant = profile_information(interest, np.column_stack((nuisance, 2*nuisance, -nuisance)))
        full = profile_information(interest, np.eye(3))
        deficient = profile_information(np.array([[1.,1.],[0.,0.]]), np.zeros((2,0)))
        underdetermined = profile_information(np.ones((1,2)), np.zeros((1,0)))
        bad = [(np.ones(3), nuisance), (interest, np.ones(3)), (interest, np.ones((2,1))),
               (np.empty((0,2)), np.empty((0,1))), (np.full_like(interest,np.nan), nuisance), (interest, np.full_like(nuisance,np.inf))]
        refused = 0
        for left, right in bad:
            try:
                profile_information(left, right)
            except ValueError:
                refused += 1
        rows.append(dict(known=bool(np.max(np.abs(base.singular_values-[3.,2.]))<=1e-12),
                         redundant=bool(redundant.nuisance_rank==1 and np.max(np.abs(base.gram-redundant.gram))<=1e-12),
                         full_null=bool(np.linalg.norm(full.projected)<=1e-12 and full.interest_rank==0),
                         rank_loss=deficient.interest_rank==1 and underdetermined.interest_rank==1 and len(underdetermined.singular_values)==2,
                         unchanged=before==(interest.tobytes(),nuisance.tobytes()), refusals=refused==6))
        arrays.append(np.concatenate([np.concatenate((r.projected.ravel(),r.gram.ravel(),r.singular_values,
                       [r.nuisance_rank,r.interest_rank,r.nuisance_cutoff,r.interest_cutoff]))
                       for r in (base,redundant,full,deficient,underdetermined)]))
    return dict(**rows[0], full_repeat=rows[0]==rows[1] and arrays[0].tobytes()==arrays[1].tobytes())


if __name__ == '__main__':
    import json
    result = selftest()
    print(json.dumps(result))
    raise SystemExit(0 if all(result.values()) else 2)
