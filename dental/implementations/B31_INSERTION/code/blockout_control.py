"""Open straight-down shadow rule and a more permissive inside-edit relaxation.
No licensed exocad output is consumed or reproduced.
"""
import pathlib, json, time
import numpy as np
from sweep import *
R = pathlib.Path(__file__).resolve().parents[1]

def straight_down_shadow(radius):
    return np.maximum.accumulate(np.asarray(radius)[::-1])[::-1]

def run(root=R):
    root = pathlib.Path(root)
    tick = time.perf_counter()
    profile = np.array([2.0, 3.0, 1.0, 0.0])
    shadow = straight_down_shadow(profile)
    protected = shadow.copy()
    protected[0] = profile[0]
    assert np.array_equal(shadow, [3, 3, 1, 0]) and protected[0] < shadow[0]
    co = json.load(open(root / 'FROZEN_COHORT.json'))
    rows = []
    for r in json.load(open(root / 'rounds/R1/COHORT.json'))['rows']:
        if r['full_path_status'] != 'COLLISION':
            continue
        rec = next((q for q in co['rows'] if q['key'] == r['key']))
        with np.load(rec['mesh_path']) as m, np.load(rec['public_path']) as p:
            keep = np.flatnonzero(m['roles'] != 1)
            tri = m['vertices'][m['faces'][keep]]
            sides = {}
            for (side, q) in r['neighbours'].items():
                if q['path']['status'] == 'COLLISION':
                    z = check_surface(tri, Obstacle(p[side]), np.array([0.0, 0.0, 40.0]), face_ids=keep, roles=m['roles'][keep])
                    assert z['status'] == 'COLLISION'
                    sides[side] = z
            rows.append(dict(key=r['key'], status='COLLISION_AFTER_REMOVING_ALL_INTAGLIO_FACES', paths=sides, scope='Surface-only relaxed control permits even nonclosed crown. Any realistic inner-only block-out preserving exterior/annulus is a subset of this permissive edit class and cannot erase the witness. This is not a deliverable repaired crown.'))
    out = dict(claim_type='capability', source='https://wiki.exocad.com/wiki/index.php/Designing_the_inside_of_the_crown#Advanced', open_rule='Downward cavity shadow; protected marginal zone overrides it locally', exact_profile_fixture=dict(radius_mm=profile.tolist(), blockout_radius_mm=shadow.tolist(), protected_radius_mm=protected.tolist(), protected_zone_incompatibility=True, kind='our_own_fixture', resolution='PHENOMENOLOGICAL'), cohort_relaxed_inside_edit_control=rows, separate_preparation_mesh_control='UNKNOWN; no separate physical/prepared solid in R4public_b. Rule-level and invariant exterior tests run; no fictional CAD result or certified spacer fit.', seconds=time.perf_counter() - tick)
    (root / 'raw/BLOCKOUT_CONTROL.json').write_text(json.dumps(out, indent=2) + '\n')
    print('OPEN_BLOCKOUT_CONTROL: exact shadow fixture; all5cases still collide after maximal inside-surface removal')
if __name__ == '__main__':
    run()
