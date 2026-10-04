"""Conditional inverse hard-tissue budget and heat-distance consumer ports.

Inputs are annotation-derived total hard tissue, not separated dentin. Results
exclude/retain a declared uniform erosion construction; no clinical feasibility.
"""
import json
import numpy as np
from scipy.special import erfc
from pathlib import Path

def remaining_lower(distance_mm, voxel_radius_mm, depth_mm, delta_each_mm=0):
    return max(0.0, distance_mm - voxel_radius_mm - 2 * delta_each_mm - depth_mm)

def depth_budget(distance_mm, voxel_radius_mm, buffer_mm, delta_each_mm=0):
    return distance_mm - voxel_radius_mm - 2 * delta_each_mm - buffer_mm

def uniform_shell_decision(distance_mm, radius_mm, material_wall_mm, buffer_mm, delta_each_mm=0):
    values = [distance_mm, radius_mm, material_wall_mm, buffer_mm, delta_each_mm]
    if not all((np.isfinite(v) and v >= 0 for v in values)):
        raise ValueError('Distances and assumed error bounds must be finite and nonnegative')
    low = max(0.0, distance_mm - radius_mm - 2 * delta_each_mm)
    high = distance_mm + radius_mm + 2 * delta_each_mm
    required = material_wall_mm + buffer_mm
    if low >= required:
        decision = 'LOCAL_CONDITION_ALL_ENCLOSURES'
    elif high < required:
        decision = 'LOCAL_CONDITION_NO_ENCLOSURE'
    else:
        decision = 'UNKNOWN_INTERVAL'
    return {'distance_lower_mm': low, 'distance_upper_mm': high, 'max_uniform_offset_lower_mm': low - buffer_mm, 'max_uniform_offset_upper_mm': high - buffer_mm, 'requested_wall_mm': material_wall_mm, 'margin_lower_mm': low - required, 'margin_upper_mm': high - required, 'decision': decision, 'clinical_feasibility': 'UNKNOWN', 'scope': 'Top2mm pulp region only; original outer contour; no cement spacer, manufacturing error, retention or fatigue included'}

def heat_transfer(rdt_mm, alpha_mm2_s, t_s):
    if rdt_mm <= 0 or alpha_mm2_s <= 0 or t_s <= 0:
        return None
    return float(erfc(rdt_mm / (2 * np.sqrt(alpha_mm2_s * t_s))))

def compatible_anatomy(quantity, tooth_group, reference):
    return quantity == reference['compared_quantity'] and tooth_group == reference['tooth_group']

def material_table(rows, cards):
    out = []
    for r in rows:
        for card in cards:
            wall = card['posterior_w_mm'] if int(r['fdi']) % 10 >= 4 else card['anterior_w_mm']
            for x in [0, 0.25, 0.5, 1.0]:
                for delta in [0, 0.15, 0.3]:
                    q = uniform_shell_decision(r['horn_boundary_min_mm'], r['digital_two_surface_radius_mm'], wall, x, delta)
                    out.append({'case': r['case'], 'fdi': r['fdi'], 'tooth_type': r['tooth_type'], 'product': card['id'], 'buffer_mm': x, 'annotation_delta_each_mm': delta, **q})
    return out

def heat_table(rows):
    out = []
    for r in rows:
        for d in [0.5, 1.0, 1.5]:
            nominal = r['horn_boundary_min_mm'] - d
            lower = max(0, nominal - r['digital_two_surface_radius_mm'])
            upper = nominal + r['digital_two_surface_radius_mm']
            for alpha in [0.1, 0.2, 0.3]:
                for t in [1, 5, 10]:
                    out.append({'case': r['case'], 'fdi': r['fdi'], 'depth_mm': d, 'alpha_mm2_s': alpha, 'time_s': t, 'R_nominal_mm': nominal, 'R_digital_lower_mm': lower, 'R_digital_upper_mm': upper, 'T_fraction_nominal': heat_transfer(nominal, alpha, t), 'T_fraction_upper': 1.0 if lower == 0 else heat_transfer(lower, alpha, t), 'T_fraction_lower': heat_transfer(upper, alpha, t), 'scope': 'Homogeneous half-space surface-temperature step; source amplitude and actual pulp temperature UNKNOWN'})
    return out
