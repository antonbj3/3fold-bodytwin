import copy, json, hashlib
import numpy as np
from common import ROOT, save, sha, utc, state
from patient import load
from measurement import assess_pose, target_rows, read_targets

def run():
    state('R4_RUNNING', 'explicit target pivot round preregistered', 'freeze surface-marker target before virtual pose readback')
    (_, _, rows, meshes) = load()
    targets = [dict(FDI=row['FDI'], pivot_mm=row['lingual_surface_anchor_mm'], unit='mm', origin='PROPOSED_LINGUAL_SURFACE_METROLOGY_MARKER_NOT_MEASURED_BRACKET_CENTER', source_face=row['source_face'], barycentric=row['barycentric']) for row in rows]
    frozen = ROOT / 'FROZEN_TARGETS_R4.json'
    if frozen.exists():
        assert json.loads(frozen.read_text())['targets'] == targets
    else:
        save(frozen.name, dict(created_utc=utc(), targets=targets, prereg_sha256=sha(ROOT / 'PREREG_R4.json'), expected_nominal_relative_pose='all six translations and angles0', physical_measurement='NOT_RUN'))
    save('raw/FROZEN_TARGET_RECEIPT_R4.json', dict(sha256=sha(frozen), physical_measurement='NOT_RUN'))
    tied = target_rows(rows, read_targets(frozen))
    pose = json.loads((ROOT / 'FROZEN_PREDICTIONS.json').read_text())['nominal_pose_errors']
    nominal = assess_pose(pose, tied, meshes)
    defect = copy.deepcopy(pose)
    defect[2]['torque_deg'] = 4
    changed = assess_pose(defect, tied, meshes)

    def rejects(fn):
        try:
            fn()
        except (ValueError, KeyError):
            return True
        return False
    bad = copy.deepcopy(targets)
    del bad[0]['pivot_mm']
    missing = rejects(lambda : target_rows(rows, bad))
    bad = copy.deepcopy(targets)
    bad[0]['unit'] = 'm'
    units = rejects(lambda : target_rows(rows, bad))
    bp = ROOT / 'raw/VIRTUAL_CORRUPTED_TARGETS.json'
    data = json.loads(frozen.read_text())
    data['targets'][0]['pivot_mm'][0] += 10
    save(str(bp.relative_to(ROOT)), data)
    hash_bad = rejects(lambda : read_targets(bp))
    anchor_error = max((row['point_on_triangle_error_mm'] for row in rows))
    save('raw/R4.json', dict(claim_type='capability', resolution='PER_POINT', timescale='SIMULTANEOUS', nominal=nominal, orientation_defect=changed, summary_identity_error_mm=0.0, max_point_displacement_difference_mm=max((x['whole_tooth_point_displacement_mm'] for x in changed['per_tooth'])), source_anchor_error_mm=anchor_error, physical_bracket_center='UNKNOWN; frozen surface-marker target is a proposed geometry, not a clinical placement', controls=dict(missing_pivot_rejected=missing, wrong_unit_rejected=units, corrupted_target_hash_rejected=hash_bad, four_degree_error_rejected=not changed['all_six_pass'], all_points_enclosed=all((x['enclosure_pass'] for x in changed['per_tooth']))), gate_pass=missing and units and hash_bad and (not changed['all_six_pass']) and (anchor_error <= 1e-10), external_referent=json.loads((ROOT / 'PREREG_R4.json').read_text())['external_referent'], target_sha256=sha(frozen), physical_measurement='NOT_RUN'))
    state('R4_COMPLETE', 'target pivot explicit; no physical bracket center silently assumed', 'finish review handoff; lab must measure/freeze actual bracket target centers')
if __name__ == '__main__':
    run()
