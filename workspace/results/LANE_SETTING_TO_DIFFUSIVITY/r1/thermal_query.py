#!/usr/bin/env python3
"""Local heat queries for an explicitly synthetic half-space contract.

This is classical heat-kernel superposition. It does not model RF deposition,
vapor transport or histological injury. Native setting queries remain unknown.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scipy.optimize import brentq


def step_response(depth_m: float, age_s: float, alpha: float, k: float) -> float:
    """Temperature rise per unit planar heat flux, in K/(W/m^2)."""
    if depth_m < 0 or alpha <= 0 or k <= 0:
        raise ValueError("Invalid depth or material coefficient")
    if age_s <= 0:
        return 0.0
    z = depth_m / (2 * math.sqrt(alpha * age_s))
    return 2 * math.sqrt(alpha * age_s) / k * (
        math.exp(-z * z) / math.sqrt(math.pi) - z * math.erfc(z)
    )


def temperature(depth_m: float, time_s: float, intervals: list[dict],
                alpha: float, k: float, initial_C: float) -> float:
    """Evaluate one target without a spatial state or stored temperature history."""
    terms = []
    previous_end = 0.0
    for row in intervals:
        start, end, flux = row['start_s'], row['end_s'], row['flux_W_m2']
        if start < previous_end or end <= start or flux < 0:
            raise ValueError("Source intervals must be ordered, disjoint and nonnegative")
        previous_end = end
        terms.append(flux * (step_response(depth_m, time_s - start, alpha, k)
                             - step_response(depth_m, time_s - end, alpha, k)))
    return initial_C + math.fsum(terms)


def threshold_depth(flux: float, time_s: float, alpha: float, k: float,
                    initial_C: float, threshold_C: float) -> float:
    """Return a synthetic isotherm depth in m, never a damage depth."""
    if threshold_C <= initial_C:
        raise ValueError("Threshold must exceed initial temperature")
    rise = threshold_C - initial_C
    if flux * step_response(0, time_s, alpha, k) <= rise:
        return 0.0
    upper = 16 * math.sqrt(alpha * time_s)
    return brentq(lambda x: flux * step_response(x, time_s, alpha, k) - rise,
                  0.0, upper, xtol=1e-15)


def threshold_time(flux: float, depth_m: float, alpha: float, k: float,
                   initial_C: float, threshold_C: float,
                   horizon_s: float = 64.0) -> float | None:
    """First threshold crossing under continuous synthetic planar heating."""
    rise = threshold_C - initial_C
    if flux <= 0 or rise <= 0:
        raise ValueError("Positive flux and temperature margin required")
    if flux * step_response(depth_m, horizon_s, alpha, k) < rise:
        return None
    return brentq(lambda t: flux * step_response(depth_m, t, alpha, k) - rise,
                  0, horizon_s, xtol=1e-13)


def query(request: dict) -> dict:
    base = {'review_state': 'PENDING_INDEPENDENT_REVIEW',
            'scientific_admission': False, 'native_dwell_s': None,
            'damage_depth_mm': None, 'physical_error_bound_K': None}
    if request.get('contract') != 'SYNTHETIC_PLANAR_HALFSPACE_V1':
        return dict(base, status='UNKNOWN',
                    reason='Matched deposition, tissue/device geometry, injury law and joint error contract are unacquired.')
    required = ('alpha_m2_s', 'k_W_m_K', 'initial_C', 'depth_m', 'time_s', 'intervals')
    if any(key not in request for key in required):
        raise ValueError("Incomplete synthetic query")
    value = temperature(request['depth_m'], request['time_s'], request['intervals'],
                        request['alpha_m2_s'], request['k_W_m_K'], request['initial_C'])
    return dict(base, status='CONDITIONAL_MODEL_VALUE', temperature_C=value,
                source_functional_K=value-request['initial_C'],
                physical_validity='Synthetic constant-property planar conduction only; no RF or damage validation.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('request', type=Path)
    args = parser.parse_args()
    print(json.dumps(query(json.loads(args.request.read_text())), indent=2))
