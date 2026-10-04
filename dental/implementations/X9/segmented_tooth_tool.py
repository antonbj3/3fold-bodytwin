"""Finite planar IPR research envelope for labelled tissue voxels (mm).

Input NPZ: labels (0 outside, 1 enamel, 2 dentin/other/unknown, 3 pulp),
spacing_mm (3,), origin_mm (3,). Coordinates are mesiodistal, buccolingual,
and crown height. Missing anatomical calibration always yields null safety.
"""
import argparse
import json
import numpy as np
from scipy.ndimage import distance_transform_edt
from ipr import write_json

def build_protection(labels, spacing, residual=0.5, boundary=0.3):
    if labels.ndim != 3 or not np.any(labels == 1) or (not np.any(labels >= 2)):
        raise ValueError('Both enamel and protected non-enamel labels are required')
    if np.any((labels < 0) | (labels > 3)):
        raise ValueError('Unknown label identifier')
    (tooth, enamel, non_enamel) = (labels > 0, labels == 1, labels >= 2)
    spacing = np.asarray(spacing, float)
    if spacing.shape != (3,) or np.any(spacing <= 0):
        raise ValueError('Three positive mm spacings are required')
    d_non = distance_transform_edt(~non_enamel, sampling=spacing)
    d_out = distance_transform_edt(tooth, sampling=spacing)
    diagonal = float(np.linalg.norm(spacing))
    protected = non_enamel | enamel & ((d_non <= residual + boundary + diagonal) | (d_non <= d_out + boundary))
    return (protected, {'residual_mm': residual, 'boundary_budget_mm': boundary, 'cell_diagonal_mm': diagonal, 'protected_cells': int(protected.sum()), 'thickness_fraction_closure': 'd_non <= d_out + boundary protects inner enamel'})

def envelope(labels, protected, spacing, origin, side, y_center, z_center, half_width=0.5, half_height=0.4, residual=0.5, boundary=0.3, overcut=0.05):
    (spacing, origin) = (np.asarray(spacing), np.asarray(origin))
    axes = [origin[i] + spacing[i] * np.arange(labels.shape[i]) for i in range(3)]
    sx = side * axes[0]
    diagonal = float(np.linalg.norm(spacing))
    footprint = (np.abs(axes[1][:, None] - y_center) <= half_width) & (np.abs(axes[2][None, :] - z_center) <= half_height)
    reach = residual + boundary + diagonal
    expanded = (np.abs(axes[1][:, None] - y_center) <= half_width + reach) & (np.abs(axes[2][None, :] - z_center) <= half_height + reach)
    occupied_x = np.any((labels > 0) & footprint[None, :, :], axis=(1, 2))
    px = np.any(protected & expanded[None, :, :], axis=(1, 2))
    if not occupied_x.any() or not px.any():
        return ({'status': 'REJECT_MISSING_FOOTPRINT_OR_PROTECTED_TISSUE', 'scenario_capacity_mm': None, 'calibrated_safe_ipr_mm': None}, None)
    support = float(np.max(sx[occupied_x]))
    obstacle = float(np.max(sx[px]))
    guard = overcut + diagonal / 2
    cap = max(0.0, support - obstacle - guard)

    def feasible(r):
        if r == 0:
            return True
        return not np.any(protected & expanded[None, :, :] & (sx[:, None, None] > support - r - guard))
    (lo, hi) = (0.0, max(5.0, support - float(np.min(sx))))
    for _ in range(60):
        mid = (lo + hi) / 2
        if feasible(mid):
            lo = mid
        else:
            hi = mid
    actual_removed = (labels > 0) & footprint[None, :, :] & (sx[:, None, None] > support - cap)
    protected_overlap = int(np.count_nonzero(actual_removed & protected)) if cap > 0 else 0
    output = {'status': 'CONDITIONAL_VOXEL_MODEL_ONLY', 'scenario_capacity_mm': cap, 'calibrated_safe_ipr_mm': None, 'nominal_plane_coordinate_mm': side * (support - cap), 'support_coordinate_mm': side * support, 'control_capacity_mm': lo, 'parity_error_mm': abs(cap - lo), 'protected_overlap_cells': protected_overlap, 'generic_0_5_pass': feasible(0.5), 'injected_plus_2_mm_rejected': not feasible(cap + 2), 'removed_enamel_mm3': float(np.count_nonzero(actual_removed & (labels == 1)) * np.prod(spacing)), 'removed_non_enamel_mm3': float(np.count_nonzero(actual_removed & (labels >= 2)) * np.prod(spacing)), 'cell_diagonal_mm': diagonal, 'tool_overcut_mm': overcut, 'research_only': True, 'anatomical_calibration': 'UNKNOWN'}
    return (output, actual_removed)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True)
    p.add_argument('--side', choices=['mesial', 'distal'], required=True)
    p.add_argument('--y-mm', type=float, required=True)
    p.add_argument('--z-mm', type=float, required=True)
    p.add_argument('--residual-mm', type=float, default=0.5)
    p.add_argument('--boundary-mm', type=float, default=0.3)
    p.add_argument('--tool-overcut-mm', type=float, default=0.05)
    p.add_argument('--output', default='ipr_output.json')
    args = p.parse_args()
    if min(args.residual_mm, args.boundary_mm, args.tool_overcut_mm) < 0:
        p.error('Budgets must be nonnegative')
    with np.load(args.input, allow_pickle=False) as d:
        (labels, spacing, origin) = (d['labels'], d['spacing_mm'], d['origin_mm'])
    (protection, diag) = build_protection(labels, spacing, args.residual_mm, args.boundary_mm)
    (result, _) = envelope(labels, protection, spacing, origin, 1 if args.side == 'mesial' else -1, args.y_mm, args.z_mm, residual=args.residual_mm, boundary=args.boundary_mm, overcut=args.tool_overcut_mm)
    result['diagnostics'] = diag
    write_json(args.output, result)
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
