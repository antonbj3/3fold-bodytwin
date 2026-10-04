"""Exact finite relational propagation, conditional on a complete oracle.

This module does not certify the caller's contact law or transition oracle.
Witnesses must be replayable; `complete=True` is an explicit oracle obligation.
UNKNOWN is sticky: a partial frontier is never promoted to complete reachability.
No probabilities, tolerance merging, entropy collapse or implicit material choice.
"""
from dataclasses import dataclass, field
from fractions import Fraction as Q
from collections import Counter
from time import perf_counter
import math

def _require_finite(value):
    if isinstance(value, (int, Q)):
        return
    try:
        finite = math.isfinite(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError('measurement numbers must be finite') from exc
    if not finite:
        raise ValueError('measurement numbers must be finite')

@dataclass(frozen=True)
class Child:
    state: tuple
    witness: tuple

@dataclass(frozen=True)
class Step:
    children: tuple
    complete: bool = False
    reason: str = 'ORACLE_COVERAGE_UNVERIFIED'

@dataclass
class Reachability:
    status: str
    layers: list
    edges: list = field(default_factory=list)
    stats: list = field(default_factory=list)
    reason: str = ''
    seconds: float = 0.0

def propagate(initial, transition, steps, *, verify, max_states=5000):
    """Merge exact FULL state identities, retaining all distinct predecessor edges.

    The caller must include every future-relevant memory/parameter in its state
    or immutable oracle context. Histories are DAG paths, not probabilities.
    On failure, `layers` contains only completed earlier layers; no branch is lost.
    """
    if steps < 0 or max_states < 1:
        raise ValueError('nonnegative steps and positive max_states required')
    start = perf_counter()
    current = dict.fromkeys(initial, 1)
    if len(current) > max_states:
        return Reachability('UNKNOWN', [current], reason='INITIAL_BUDGET')
    out = Reachability('COMPLETE', [current])
    for k in range(steps):
        (nxt, edges, generated, impossible) = (Counter(), [], 0, 0)
        (split_excess, branching_parents) = (0, 0)
        for (parent, multiplicity) in current.items():
            answer = transition(parent)
            if not answer.complete:
                (out.status, out.reason) = ('UNKNOWN', answer.reason)
                out.seconds = perf_counter() - start
                return out
            children = tuple(dict.fromkeys(answer.children))
            if not children:
                impossible += 1
            child_states = {child.state for child in children}
            split_excess += max(len(child_states) - 1, 0)
            branching_parents += int(len(child_states) > 1)
            seen_children = set()
            for child in children:
                if not verify(parent, child):
                    (out.status, out.reason) = ('UNKNOWN', 'INVALID_WITNESS')
                    out.seconds = perf_counter() - start
                    return out
                edges.append((parent, child.state, child.witness))
                if child.state not in seen_children:
                    nxt[child.state] += multiplicity
                    generated += 1
                    seen_children.add(child.state)
                if len(nxt) > max_states:
                    (out.status, out.reason) = ('UNKNOWN', 'FRONTIER_BUDGET')
                    out.seconds = perf_counter() - start
                    return out
        out.stats.append(dict(step=k + 1, states=len(nxt), generated=generated, merges=generated - len(nxt), impossible_parents=impossible, split_excess=split_excess, branching_parents=branching_parents, histories=sum(nxt.values()), witness_edges=len(edges)))
        out.edges.append(edges)
        out.layers.append(dict(nxt))
        current = dict(nxt)
    out.seconds = perf_counter() - start
    return out

def observe(states, readout, value, error):
    """Bounded-error set intersection. An empty answer is inconsistent evidence."""
    _require_finite(value)
    _require_finite(error)
    if error < 0:
        raise ValueError('negative measurement error')
    compatible = []
    for state in states:
        center = readout(state)
        _require_finite(center)
        if abs(center - value) <= error:
            compatible.append(state)
    return tuple(compatible)

def rank_measurements(states, sensors):
    """Minimax survivor count, then disjoint gap/resolution/cost.

    Sensors are explicit dicts: id, readout, error, resolution, cost, unit, scope.
    Resolution enters the disjoint-band score only; error must already include any
    quantization uncertainty. Closed bands touching at one endpoint overlap.
    Uniform counting is a design objective, not a probability or entropy claim.
    """
    states = tuple(states)
    if not states:
        return []
    ranked = []
    for m in sensors:
        for key in ('error', 'resolution', 'cost'):
            _require_finite(m[key])
        if m['error'] < 0 or m['resolution'] <= 0 or m['cost'] <= 0 or (not m['unit']) or (not m['scope']):
            raise ValueError('invalid sensor contract')
        centers = [m['readout'](s) for s in states]
        for center in centers:
            _require_finite(center)
        bands = [(v - m['error'], v + m['error']) for v in centers]
        for (lo, hi) in bands:
            _require_finite(lo)
            _require_finite(hi)
        worst = max((sum((lo <= x <= hi for (lo, hi) in bands)) for (x, _) in bands))
        gaps = [max(c - b, a - d) for (i, (a, b)) in enumerate(bands) for (c, d) in bands[i + 1:]]
        for gap_value in gaps:
            _require_finite(gap_value)
        gap = min(gaps) if gaps else Q(0)
        score = gap / m['resolution'] / m['cost'] if gap > 0 else Q(0)
        _require_finite(score)
        ranked.append(dict(id=m['id'], worst_survivors=worst, guaranteed_removed=len(states) - worst, min_gap=gap, farkas_score=score, unit=m['unit'], scope=m['scope'], bands=bands))
    return sorted(ranked, key=lambda a: (a['worst_survivors'], -a['farkas_score'], a['id']))
