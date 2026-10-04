import json, hashlib, datetime, pathlib, shutil
R = pathlib.Path(__file__).resolve().parents[1]
P = R.parent / 'PROOF_LANE_FULL_CROWN_R4'
O = R.parent / 'LANE_X95_PERSONA_10X'

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def freeze(n, j):
    p = R / n
    if p.exists():
        raise RuntimeError('immutable already exists ' + str(p))
    p.write_text(json.dumps(j, indent=2, ensure_ascii=False) + '\n')
    p.with_suffix('.sha256').write_text(sha(p) + '\n')

def run():
    t = datetime.datetime.now(datetime.timezone.utc).isoformat()
    inputs = json.loads((P / 'FROZEN_INPUTS_B.json').read_text())
    meshes = json.loads((P / 'FROZEN_PREDICTIONS_C.json').read_text())
    exports = json.loads((P / 'FROZEN_EXPORTS.json').read_text())
    rows = []
    for r in inputs['records']:
        m = next((x for x in meshes['rows'] if x['key'] == r['key'] and x['method'] == 'exact_margin'))
        rows.append(dict(**r, mesh_path=m['mesh_path'], mesh_sha256=m['mesh_sha256']))
    assert len(rows) == 18 and len(exports['exports']) >= 3
    for r in rows:
        assert sha(r['public_path']) == r['public_sha256']
        assert sha(r['mesh_path']) == r['mesh_sha256']
    freeze('FROZEN_COHORT.json', dict(frozen_utc=t, rows=rows, exports=exports['exports'], sources=inputs['sources'], source_manifests={str(P / n): sha(P / n) for n in ['FROZEN_INPUTS_B.json', 'FROZEN_PREDICTIONS_C.json', 'FROZEN_TEST_SET.json', 'FROZEN_EXPORTS.json']}, selection='All eighteen exact_margin C meshes, no outcome selection. All three existing selected exports retained separately.', frame='local mm: source OBJ @ frozen frame.T minus public_b origin_world; no lower arch used', path={'orientation': 'identity', 'translation_mm': [0, 0, 40], 's_interval': [0, 1], 'pose': 'final crown + s*translation; insertion from 1 to 0'}, preparation='R4 cavity cone is an intaglio-defining virtual die, not a separate measured prepared solid. Coincident seated surface is intentional; no new spacer invented. Prep-volume qualification remains UNKNOWN.', licence='Teeth3DS dataset-object licence UNKNOWN. Private read-only processing; no redistribution claim.'))
    shutil.copyfile(O / 'PREREG_B02.json', R / 'PREREG_B02_ORIGINAL.json')
    freeze('PREREG_R1.json', dict(round='R1', frozen_utc=t, claim_type='capability', parent_plan_sha256=sha(R / 'PREREG_B02_ORIGINAL.json'), capability='Locate full-crown straight-path neighbour collisions and certify clearance only with whole interval coverage; assess whether inside-only repair can preserve fixed margin.', obstacle='Endpoint minima and local intaglio normals do not encode whole-sweep crown-neighbour intersection.', changed_operation='Full triangle swept prisms with outward-rounded interval separating projections and exact rational barycentric intersection witnesses proposed by SciPy HiGHS.', consumer=['lab', 'cad', 'gp'], cohort_sha256=sha(R / 'FROZEN_COHORT.json'), metrics={'coverage': 'All18+3exports retained, every missing/incomplete path UNKNOWN unless a certified collision already refutes it', 'boolean_reference_agreement': 'exact, no tolerance', 'margin_max_mm': 0.01, 'falsifier_bulge_mm': 0.2, 'summary_identity_error_required': 0.0}, decision='Digital collision requires exact validated witness. PASS requires every swept pair rigorously separated AND solid/start coverage. Open obstacle surfaces yield surface-only clearance, total path UNKNOWN. Physical insertion always UNKNOWN.', strongest_equally_informed_control='Independent barycentric triangle pair LP over s in [0,1], exact primal reconstruction; open exocad straight-axis shadow/block-out rule, not licensed CAD.', falsifiers=['Identical endpoint minimum with0.2mm intermediate bulge must change path decision', 'Invalid barycentric weight, wrong point or omitted obstacle must be rejected', 'Any unsupported rotation/missing scene information must abstain'], full_cost={'preparation_budget_h': 0.5, 'sweep_and_repair_budget_h': 2.5, 'independent_geometry_budget_h': 1, 'export_budget_h': 0.5, 'fit': 0, 'discovery_query_validation_fallback': 'charged elapsed CPU/RSS; historical data acquisition and lab costs UNKNOWN', 'threads_max': 4, 'intermediate_bytes_max': 3000000000, 'GPU': False}, axis='First round fixed local+z; successor may search declared finite straight axes; rotations UNKNOWN', repair_class='Inside only; exterior and complete exterior marginal curve fixed. An exterior/annulus collision is an impossibility witness for this class. No clinical feasible repair claim.'))
    leaves = [('source triangles', 'EXTERNALLY_MEASURED', 'Teeth3DS frozen OBJ/labels; published dataset surface, scanner calibration UNKNOWN; immutable coordinates used as exact digital dyadics'), ('virtual cone', 'CONSTITUTIVE_CLOSURE', 'R4 public_b .65 XY scaling plus apex; this is the existing construction, not a measured preparation'), ('rigid translation', 'CONSTITUTIVE_CLOSURE', '40mm local+z path, no compliance/load; start safety checked'), ('swept triangle convexity', 'DERIVED_UNDER_ASSUMPTIONS', 'conv(A,A+D) exactly equals union(A+sD), s in[0,1]; affineness sufficient, rotations outside scope'), ('separation bound', 'DERIVED_UNDER_ASSUMPTIONS', 'Any axis with outward-rounded min(B)-max(sweep A)>0 proves empty intersection; float axis need not be true facet normal'), ('collision witness', 'DERIVED_UNDER_ASSUMPTIONS', 'Nonnegative exact rational weights alpha,beta,s,t, sum alpha=sum beta=s+t=1 and A alpha+sD=B beta prove an actual intersection'), ('solid and scanner truth', 'UNKNOWN', 'Open segmented neighbour patches are measured surfaces, not certified closed solids; independent scene/crown metrology needed'), ('physical insertion', 'UNKNOWN', 'Matched prep/crown/neighbor scans, registration band and physical insertion trial absent; no force, remake or10x claim')]
    freeze('DECOMPOSITION_R1.json', dict(chain='desired usable insertion -> entire path free + unchanged margin + smallest inside change -> nonlocal swept scene obstruction -> triangle prism + repair-invariant surface -> exact certificates', equations=['T(s)A=sum_i alpha_i A_i+sD', 'sum alpha=sum beta=s+t=1; alpha,beta,s,t>=0', 'max_{A+sD} n.x < min_B n.x => disjoint for all s', 'unchanged exterior collision => every inside-only repair has same collision'], representation='original case, face id, region, barycentric point, exact rational pose, mm, PER_POINT, SIMULTANEOUS', leaves=[dict(name=n, status=s, stop=a) for (n, s, a) in leaves], linear_sensitivity='None reported. Affine kinematics have rigorous rational and interval enclosures; physical sensitivities UNKNOWN.'))
    freeze('FROZEN_PREDICTIONS.json', dict(frozen_utc=t, prereg_sha256=sha(R / 'PREREG_R1.json'), cohort_sha256=sha(R / 'FROZEN_COHORT.json'), predictions=[dict(key=r['key'], full_path='UNKNOWN until checked; no positivity forecast', physical_path='UNKNOWN') for r in rows], falsifier_prediction='Both endpoint summaries identical exactly; bulge fixture collision only for bulged path', future_measurement='None executed or available; this is prospective digital audit contract.'))
    (R / 'CURRENT_WORK_STATE.json').write_text(json.dumps(dict(lane='X95-whole-insertion-path', updated_utc=t, status='PREREGISTERED_R1', latest_gate='Cohort18 and3exports frozen before geometry', next_operation='First own exact-summary fixture calculation, then primary-source countercheck'), indent=2) + '\n')
    print('FROZEN18cases3exports', sha(R / 'FROZEN_COHORT.json'))
if __name__ == '__main__':
    run()
