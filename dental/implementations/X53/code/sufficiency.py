import sys
from common import *
sys.path.insert(0, str(ROOT / 'tests'))
from test_collision import cup, tool
from collision import search

def run():
    points = np.array([[x, y, 0.0] for x in [-0.2, 0.0, 0.2] for y in [-0.2, 0.0, 0.2]])
    t = tool(0.6, 3)
    normal = np.array([0.0, 0.0, 1.0])
    states = []
    for h in [1.0, 6.0]:
        scene = cup(height=h)
        found = [search(scene, p, normal, t, 4)['status'] == 'FOUND' for p in points]
        states.append(dict(height_mm=h, summary_tip_diameter_mm=t['diameter_mm'], summary_local_target_curvature_per_mm=0.0, found_fraction=float(np.mean(found)), point_queries=found))
    error = abs(states[0]['summary_tip_diameter_mm'] - states[1]['summary_tip_diameter_mm'])
    diff = abs(states[0]['found_fraction'] - states[1]['found_fraction'])
    assert error == 0.0 and diff >= read(ROOT / 'PREREG_R1.json')['metrics']['sufficiency_downstream_fraction_min']
    out = dict(round='R1', claim_type='capability', resolution='PER_POINT', states=states, summary_identity_error=error, downstream_access_fraction_difference=diff, rejected_summary='tip diameter + local planar curvature', minimum_sufficient_extension_for_this_pair='neck-reach versus obstacle clearance along the allowed approach; full pose-dependent assembly clearance for general scenes', external_referent=dict(kind='our_own_fixture', locator='tests/test_collision.py:cup; raw/SOURCE_CHECKS.json gives nominal published L3', compared_quantity='constructed pocket straight-access witness; fixture is NOT physical reference', refutes_us=True), numerical_clearance_rigorous_interval_enclosure='MISSING; independent dense Lipschitz control tested to 1e-6 mm')
    dump(ROOT / 'raw/SUFFICIENCY_R1.json', out)
    return out
if __name__ == '__main__':
    print(run())
