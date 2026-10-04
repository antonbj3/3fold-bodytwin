"""Fail-closed PS5 knee-local FDK entry point and input audit.

The trial adapter requires calibrated, full-cycle five-DOF inputs. A stance-only
L1 operator and an unregistered STL mesh cannot establish those inputs.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MESH = Path('external_media')
EVENTS = {'ngait_og_ss1': (5.116667, 6.258333),
          'rightturn6': (5.150000, 6.466667)}
DOF = ('anterior_mm', 'proximal_mm', 'lateral_mm',
       'internal_rotation_deg', 'varus_deg')


class FDKInputError(RuntimeError):
    """Raised before outcome access when the physical input contract is unmet."""
    def __init__(self, audit):
        self.audit = audit
        super().__init__('PS5 five-DOF FDK input gate failed: ' + ', '.join(audit['missing']))


def _sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def input_audit(trial):
    """Inspect model coverage and registration; never reads eTibia outcome."""
    trial = str(trial).removeprefix('PS5__')
    if trial not in EVENTS:
        raise ValueError(f'unsupported PS5 trial: {trial}')
    key = f'PS5__{trial}'
    start, end = EVENTS[trial]
    grid = np.linspace(start, end, 101)
    prep = ROOT / 'results/L1/prep' / f'{key}.npz'
    n12 = ROOT / 'results/N12b/cloud_out/out' / f'{key}.npz'
    z = np.load(prep, allow_pickle=False)
    k = np.load(n12, allow_pickle=False)
    l1time = np.asarray(z['t'], float)
    iktime = np.asarray(k['t'], float)
    nearest = np.argmin(abs(grid[:, None] - l1time[None, :]), axis=1)
    halfstep = float(np.median(np.diff(l1time)) / 2)
    l1covered = abs(grid - l1time[nearest]) <= halfstep + 1e-6
    ikcovered = (grid >= iktime[0] - 1e-6) & (grid <= iktime[-1] + 1e-6)
    mesh_paths = [MESH / 'Femoral Component.stl', MESH / 'Tibial Insert.stl']
    source_paths = [prep, n12, *mesh_paths]
    # No accepted calibration or 5-DOF muscle operator has been supplied for
    # these trials. Explicit booleans keep them distinct from file existence.
    flags = {
        'full_cycle_flexion_IK': bool(ikcovered.all()),
        'full_cycle_six_axis_external_wrench': False,
        'mesh_to_bone_segment_registration': False,
        'five_secondary_DOF_muscle_operator': False,
        'PS5_insert_material_contact_calibration': False,
        'PS5_trial_matched_fluoroscopy': False,
    }
    return {
        'trial': trial, 'events_s': [start, end], 'grid_points': 101,
        'dof': list(DOF), 'flags': flags,
        'missing': [name for name, ok in flags.items() if not ok],
        'l1_rows': int(len(l1time)),
        'l1_time_s': [float(l1time[0]), float(l1time[-1])],
        'l1_near_grid_points': int(l1covered.sum()),
        'l1_uncovered_grid_percent': np.flatnonzero(~l1covered).tolist(),
        'ik_rows': int(len(iktime)),
        'ik_time_s': [float(iktime[0]), float(iktime[-1])],
        'ik_grid_points_in_range': int(ikcovered.sum()),
        'ik_missing_grid_percent': np.flatnonzero(~ikcovered).tolist(),
        'mesh_present': {p.name: p.is_file() for p in mesh_paths},
        'source_sha256': {str(p): _sha(p) for p in source_paths if p.is_file()},
        'outcome_read': False,
    }


def fdk_knee(trial):
    """Return a trial's five-DOF FDK result once its physical input gate passes.

    For the current PS5 trials the required calibration and full-cycle operator
    do not exist; fail with machine-readable audit rather than emit a force
    trace whose coordinate frame and secondary equilibrium are undefined.
    """
    audit = input_audit(trial)
    if audit['missing']:
        raise FDKInputError(audit)
    raise NotImplementedError('validated FDK frame solver is not yet supplied')


if __name__ == '__main__':
    print(json.dumps({trial: input_audit(trial) for trial in EVENTS}, indent=2))
