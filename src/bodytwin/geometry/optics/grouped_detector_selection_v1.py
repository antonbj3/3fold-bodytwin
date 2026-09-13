"""Enumerate explicitly feasible groups using profiled local information.

Callers supply already-whitened Jacobians and parameter scales. Shared nuisance
columns refer to the same parameters across groups; local nuisance columns are
independent for each selected group. Feasibility is explicit, never inferred.
"""
from dataclasses import dataclass
from numbers import Integral
import numpy as np
from bodytwin.geometry.optics.profiled_information_v1 import profile_information


@dataclass(frozen=True)
class GroupedSelection:
    indices: tuple
    candidates: tuple
    scores: np.ndarray
    ranks: np.ndarray


def select_groups(interest, shared_nuisance, local_nuisance, candidates):
    """Select among equally sized, explicitly feasible group-index tuples.

    Jacobians have shape (groups, observations_per_group, parameters). Empty
    shared/local nuisance dimensions are permitted. Invalid or duplicate candidate
    tuples are errors. The best candidate must resolve all interest columns under
    the profiler's explicit numerical-rank rule. No arrays are written to disk.
    """
    arrays = [np.array(value, dtype=np.float64, copy=True)
              for value in (interest, shared_nuisance, local_nuisance)]
    if any(value.ndim != 3 or not np.isfinite(value).all() for value in arrays):
        raise ValueError('Finite three-dimensional Jacobians required')
    interest, shared, local = arrays
    groups, observations, parameters = interest.shape
    if not groups or not observations or not parameters or any(value.shape[:2] != interest.shape[:2] for value in arrays[1:]):
        raise ValueError('Nonempty interest and matching group/observation dimensions required')
    canonical = []
    for candidate in candidates:
        candidate = tuple(candidate)
        if not candidate or any(not isinstance(i, Integral) or isinstance(i, (bool,np.bool_)) or i < 0 or i >= groups for i in candidate):
            raise ValueError('Invalid candidate indices')
        if len(set(candidate)) != len(candidate):
            raise ValueError('Candidate contains duplicate group')
        canonical.append(tuple(sorted(int(i) for i in candidate)))
    if not canonical or len(set(canonical)) != len(canonical) or len({len(t) for t in canonical}) != 1:
        raise ValueError('Distinct nonempty candidates of equal size required')
    canonical = tuple(sorted(canonical))
    scores, ranks = [], []
    for candidate in canonical:
        count = len(candidate)
        nuisance = np.zeros((count*observations, shared.shape[2]+count*local.shape[2]))
        nuisance[:, :shared.shape[2]] = shared[list(candidate)].reshape(count*observations, shared.shape[2])
        for slot, group in enumerate(candidate):
            nuisance[slot*observations:(slot+1)*observations,
                     shared.shape[2]+slot*local.shape[2]:shared.shape[2]+(slot+1)*local.shape[2]] = local[group]
        result = profile_information(interest[list(candidate)].reshape(count*observations, parameters), nuisance)
        ranks.append(result.interest_rank)
        scores.append(float(result.singular_values[-1]) if result.interest_rank == parameters else 0.)
    if not np.isfinite(scores).all() or max(scores) <= 0:
        raise ValueError('No candidate resolves all interest parameters')
    winner = max(range(len(canonical)), key=lambda i: scores[i])
    scores = np.array(scores); ranks = np.array(ranks, dtype=np.int64)
    scores.setflags(write=False); ranks.setflags(write=False)
    return GroupedSelection(canonical[winner], canonical, scores, ranks)


def selftest():
    import itertools
    interest = np.zeros((4,2,1)); interest[:,0,0] = [1,2,2,3]
    shared = np.zeros((4,2,0)); local = np.ones((4,2,1))
    candidates = list(itertools.combinations(range(4),2))
    before = [v.tobytes() for v in (interest,shared,local)]
    first = select_groups(interest,shared,local,candidates)
    second = select_groups(interest,shared,local,reversed(candidates))
    expected = np.array([np.sqrt(sum(interest[i,0,0]**2 for i in t)/2) for t in candidates])
    malformed = [[], [(0,0)], [(0,4)], [(0,1),(1,0)], [(0,),(1,2)], [(False,1)], [(0.,1)]]
    refused = 0
    for cases in malformed:
        try: select_groups(interest,shared,local,cases)
        except ValueError: refused += 1
    invalid_arrays = [(interest[:,0],shared,local), (interest,shared[:3],local),
                      (np.full_like(interest,np.nan),shared,local), (np.ones_like(interest),shared,local)]
    for args in invalid_arrays:
        try: select_groups(*args,candidates)
        except ValueError: refused += 1
    return dict(known=bool(np.max(np.abs(first.scores-expected))<=1e-12),
                tie=first.indices==second.indices==(1,3),
                complete_repeat=first.candidates==second.candidates and first.scores.tobytes()==second.scores.tobytes() and first.ranks.tobytes()==second.ranks.tobytes(),
                nonmutation=before==[v.tobytes() for v in (interest,shared,local)],
                refusals=refused==11)


if __name__ == '__main__':
    import json
    result=selftest(); print(json.dumps(result))
    raise SystemExit(0 if all(result.values()) else 2)
