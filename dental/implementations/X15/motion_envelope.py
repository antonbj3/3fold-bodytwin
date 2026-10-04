"""Finite paired-section constraint port. Every physical answer is UNDETERMINED."""
import numpy as np

def interval(model, angles_deg, epsilon_mm=0.0):
    theta = np.deg2rad(np.asarray(angles_deg, dtype=float))
    x = np.array(model['x_root_mm'])
    z = np.array(model['height_mm']) - model['pivot_height_mm']
    rotated = np.cos(theta)[..., None] * x + np.sin(theta)[..., None] * z
    lower = np.max(np.array(model['wall_lingual_mm']) + epsilon_mm - rotated, axis=-1)
    upper = np.min(np.array(model['wall_buccal_mm']) - epsilon_mm - rotated, axis=-1)
    return (lower, upper)

def query(model, t_mm, theta_deg, epsilon_mm=0.608):
    (lo, hi) = interval(model, theta_deg, 0.0)
    margin = float(min(t_mm - lo, hi - t_mm))
    theta = np.deg2rad(theta_deg)
    x = np.array(model['x_root_mm'])
    z = np.array(model['height_mm']) - model['pivot_height_mm']
    rot = x * np.cos(theta) + z * np.sin(theta) + t_mm
    m = np.r_[rot - np.array(model['wall_lingual_mm']), np.array(model['wall_buccal_mm']) - rot]
    active = int(np.argmin(m))
    n = len(x)
    endpoint = active % n
    return {'case': model['case'], 'tooth': model['tooth'], 'source': model['source'], 'translation_mm': t_mm, 'angle_deg': theta_deg, 'conditional_margin_mm': margin, 'max_directional_edge_budget_mm': margin if margin > 0 else None, 'epsilon_scenario_mm': epsilon_mm, 'conditional_geometric_status': 'WITHIN_SECTION_CONSTRAINTS' if margin >= epsilon_mm else 'VIOLATES_SECTION_CONSTRAINTS' if margin < -epsilon_mm else 'UNDETERMINED_EDGE_SCENARIO', 'active_height_mm': model['height_mm'][endpoint], 'active_wall': 'lingual' if active < n else 'buccal', 'active_endpoint': endpoint, 'anatomical_status': 'UNDETERMINED', 'biological_status': 'UNDETERMINED', 'unresolved': ['CEJ', 'full root surface', 'wall requery after vertical displacement', 'center of resistance', 'scanner-specific local edge bound', 'bone remodeling'], 'contract': 'Finite fixed paired section planes; not 3D anatomical containment or safe treatment'}
