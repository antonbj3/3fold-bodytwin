from common import *

def run():
    p = np.array([0.25, 0.0, 0.0])
    a = p.copy()
    b = -p
    s1 = [float(np.linalg.norm(a)), 0.0]
    s2 = [float(np.linalg.norm(b)), 0.0]
    downstream = [dict(tracking_error_mm=float(np.linalg.norm(v - p)), completion=float(v @ p / (p @ p))) for v in [a, b]]
    assert np.array_equal(s1, s2)
    RA = Rotation.from_euler('z', 10, degrees=True).as_matrix()
    RB = Rotation.from_euler('z', -10, degrees=True).as_matrix()
    (ta, tb) = [float(np.degrees(Rotation.from_matrix(R).magnitude())) for R in [RA, RB]]
    assert ta == tb
    rotational_gap = float(np.degrees(Rotation.from_matrix(RB @ RA.T).magnitude()))
    f1 = np.array([1.0, 0.0, 0.0])
    f2 = np.array([0.0, 1.0, 0.0])
    M = np.diag([0.2, 0.05, 0.01])
    u1 = M @ f1
    u2 = M @ f2
    assert np.linalg.norm(f1) == np.linalg.norm(f2)
    measured = np.eye(4)
    world_A = np.eye(4)
    world_B = np.eye(4)
    world_B[1, 3] = 0.25
    scanner_A = np.eye(4)
    scanner_B = np.linalg.inv(world_B)
    e = float(np.max(np.abs(scanner_A @ world_A - scanner_B @ world_B)))
    assert e == 0
    out = dict(claim_type='capability', resolution='PER_TOOTH', external_referent=dict(kind='closed_form', locator='SE(3) invariance G^-1(GH)=H; https://doi.org/10.1107/S0567739476001873', compared_quantity='Signed displacement, scan-gauge invariance and exact downstream residual', refutes_us=True), scalar_motion=dict(summary_A=s1, summary_B=s2, identity_error=0.0, downstream=downstream, downstream_error_difference_mm=0.5, downstream_completion_difference=2.0, sufficient=False, minimum_extension='Signed projection onto the planned vector plus perpendicular residual for this tracking question; full SE3 for arbitrary future point/direction queries'), scalar_rotation=dict(summary_A_deg=ta, summary_B_deg=tb, identity_error=abs(ta - tb), tracking_gap_A_deg=0.0, tracking_gap_B_deg=rotational_gap, sufficient=False, minimum_extension='Oriented relative SO3 pose for arbitrary plan tracking; signed rotation-axis projection plus residual for this fixed question'), scalar_force=dict(summary_A_N=1.0, summary_B_N=1.0, identity_error=0.0, predicted_planned_axis_motion_mm=[float(u1[0]), float(u2[0])], downstream_difference_mm=float(abs(u1[0] - u2[0])), sufficient=False, model_status='PHENOMENOLOGICAL illustrative mobility only; no patient coefficient', physical_rigorous_enclosure='MISSING; u=MF is exact only inside this specified illustrative linear fixture', minimum_extension='Force direction and exposure; crown moment for rotations; measured mobility/force scale'), absolute_gauge=dict(scan_A=measured.tolist(), scan_B=measured.tolist(), identity_error=e, skull_intrusion_mm=[0.0, 0.25], downstream_difference_mm=0.25, sufficient=False, minimum_extension='Externally fixed oriented reference that constrains all six rigid modes; three noncollinear stable points suffice in the rigid idealization'), limitations='Constructed states provide insufficiency under stated mechanics; they are not the external patient outcome reference.')
    write(ROOT / 'raw/SUFFICIENCY.json', out)
    return out
if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
