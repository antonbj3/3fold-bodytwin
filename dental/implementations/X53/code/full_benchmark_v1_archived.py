"""R4: exact halfspace support for ALL original GenCAD V2 task intaglios.

This does not turn a roof into a full crown. It removes sampling from the
limited planar-intaglio whole-tool query, in the same no-fixture scenario.
"""
from dental_release.paths import expand as _release_expand
import time, collections, trimesh
from fractions import Fraction as Q
from common import *
from collision import Scene, search
from run_r1 import lib

def certificate(design, tool):
    inner = design['inner_vertices']
    outer = design['outer_vertices']
    z = [Q(str(v[2])) for v in inner]
    if not z or min(z) != max(z):
        return dict(status='UNKNOWN_NONPLANAR_INTAGLIO')
    plane = z[0]
    if any((Q(str(v[2])) < plane for v in outer)):
        return dict(status='FAIL_RETAINED_VERTEX_BELOW_SUPPORTING_PLANE')
    r = Q(str(tool['diameter_mm'])) / 2
    ell = Q(str(tool['neck_reach_mm']))
    gauge = Q(str(tool['gauge_mm']))
    rs = Q(str(tool['shank_mm'])) / 2
    rh = Q(str(tool['holder_diameter_mm'])) / 2
    margins = dict(ball='0', shank=str(r + ell - rs), holder=str(r + gauge - rh))
    good = r > 0 and ell > 0 and (gauge >= ell) and (r + ell >= rs) and (r + gauge >= rh)
    return dict(status='EXACT_SCENARIO_FULL_PLANAR_INTAGLIO' if good else 'FAIL_ASSEMBLY_HALFSPACE_SUPPORT', plane_z_mm=str(plane), support_margins_mm=margins, direction=[0, 0, -1], minimum_green_tool_mm=tool['diameter_mm'] / 0.8, scope='All points on the serialized planar inner triangles; all retained vertices above supporting plane; no fixture/workspace; declared capsule envelope. Externally measured assembly UNKNOWN.')

def run():
    path = ROOT / 'PREREG_R4.json'
    if not path.exists():
        pre = dict(round='R4', claim_type='capability', capability='Answer the whole-tool question on every original GenCAD V2 task and remove point sampling from its planar-intaglio subcase', obstacle='R1 full-crown study covers only31 full caps; original576 roof tasks include unknown and abstained cases and must remain visible', changed_operation='Exact rational halfspace support of all tool capsules and serialized mesh vertices: below-plane approach certifies every point of planar intaglio', consumer='PROOF_LANE_GENCAD_V2 task milling check; DENT-MFG-ASBUILT-DEVIATION', resolution='PER_POINT', timescale='SIMULTANEOUS', metrics={'task_count': 576, 'planarity': 'all serialized inner z coordinates machine-identical and exact decimal rational equal', 'retained_support': 'every outer vertex z >= inner plane, exact rational', 'tool_support': 'r+neck >= shank_radius; r+gauge >= holder_radius', 'injected_bad_gauge_green_mm': 1.0}, decision_criteria='Every requested task gets exact scoped certificate, retained ABSTAIN, or UNKNOWN; compare numerical whole-tool pose at four source inner vertices for each certifiable design. Injected1mm gauge must reject support.', strongest_equal_information_control='Existing Scene continuous swept-capsule predicate with identical mesh/tool/direction, no additional data', external_referent=dict(kind='closed_form', locator=_release_expand('https://mathworld.wolfram.com/Sphere.html; @DENTAL_IMPLEMENTATIONS@/GENCAD_V2/gencad_bench_v2/checks.py'), compared_quantity='sphere support sup z=c_z+r and piecewise-linear triangle support; not physical CAM', refutes_us=True), full_cost={'preparation': 'read576 serialized task/output hashes, no source copies', 'fit': 0, 'discovery': 'exact rational supporting-plane and capsule inequalities', 'validation': 'four continuous whole-tool queries per certified design plus an injected early holder', 'questions': 0, 'fallback': 'UNKNOWN on nonplanar/invalid shape; original missing/abstentions retained', 'physical': 'none; same unknown fixture/holder/mount'}, limitations=['Only original roof-model intaglio. Exterior/axial margins/retention remain uncaptured by this proof.', 'No material-compatibility, actual installed neck/holder, workspace or CAM-path generation certified.', 'This is a capability coverage/proof continuation, not novelty or superiority.'])
        dump(path, pre)
        path.with_suffix('.sha256').write_text(sha(path) + '\n')
        dump(ROOT / 'DECOMPOSITION_R4.json', dict(leaves=[dict(leaf='support plane', status='DERIVED_UNDER_ASSUMPTIONS', equation='sup_triangle z = max_vertex z; sup_capsule z=max(endpoint z)+radius', stop_argument='Exact convex support; rational arithmetic eliminates rounding from this narrowly vertical containment test.'), dict(leaf='actual retained material', status='UNKNOWN', equation='Roof mesh above plane, real axial/fixture material not supplied', stop_argument='Cannot extrapolate to a full crown or a mounted blank.')]))
    assert sha(path) == path.with_suffix('.sha256').read_text().strip()
    state('R4_FULL_BENCHMARK_RUNNING', sha(path), 'Exact planar whole-tool support and four numerical poses on all576 original tasks')
    start = time.perf_counter()
    base = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V2'))
    rows = []
    tools = lib('vhf', 5)
    tool = tools[0]
    fnd = DENTAL / 'results/PROOF_LANE_GENCAD_V2'
    import sys
    sys.path.insert(0, str(fnd))
    from gencad_bench_v2.geometry import shell
    public = sorted((p for p in (base / 'public').glob('*.json') if any((p.stem.endswith('_' + level) for level in ['easy', 'normal', 'hard', 'boundary']))))
    for taskfile in public:
        t = read(taskfile)
        designfile = base / 'predictions/constraint_optimizer' / taskfile.name
        if not designfile.exists():
            rows.append(dict(task_id=t['task_id'], status='UNKNOWN_MISSING_DESIGN'))
            continue
        d = read(designfile)
        row = dict(task_id=t['task_id'], family=t['family'], level=t['level'] if 'level' in t else taskfile.stem.rsplit('_', 1)[1], source_task_sha256=sha(taskfile), source_design_sha256=sha(designfile), source_design_file=str(designfile), resolution='PER_POINT', external_and_full_restoration_access='UNKNOWN_UNMODELED_AXIAL_AND_FIXTURE')
        if d.get('status') != 'DESIGN':
            row.update(status=d.get('status', 'UNKNOWN'), reason=d.get('reason'))
            rows.append(row)
            continue
        if d.get('units') != 'mm' or d.get('frame') != t['frame']:
            row.update(status='INVALID_FRAME_OR_UNITS')
            rows.append(row)
            continue
        cert = certificate(d, tool)
        row.update(status=cert['status'], certificate=cert)
        if cert['status'] == 'EXACT_SCENARIO_FULL_PLANAR_INTAGLIO':
            inner = np.asarray(d['inner_vertices'])
            outer = np.asarray(d['outer_vertices'])
            xy = outer[:, :2]
            (v, f) = shell(xy, outer[:, 2], inner[:, 2], np.asarray(d['faces'], int))
            m = trimesh.Trimesh(v, f, process=True)
            m.fix_normals()
            scene = Scene(m.vertices, m.faces)
            normal = np.array([0, 0, -1.0])
            selected = sorted(set([0, len(inner) // 3, 2 * len(inner) // 3, len(inner) - 1]))
            controls = []
            for i in selected:
                w = search(scene, inner[i], normal, tool, 4)
                controls.append(w['status'] == 'FOUND')
            assert all(controls)
            corrupt = dict(tool, gauge_mm=0.8)
            bad = certificate(d, corrupt)
            assert bad['status'] == 'FAIL_ASSEMBLY_HALFSPACE_SUPPORT'
            row.update(numerical_wholetool_controls=controls, injected_holder_rejected=True, per_inner_vertex_minimum_green_mm=[0.6] * len(inner), inner_vertex_count=len(inner), whole_inner_triangle_surface_certified=True)
        rows.append(row)
    assert len(rows) == 576
    counts = dict(collections.Counter((r['status'] for r in rows)))
    result = dict(round='R4', claim_type='capability', rows=rows, requested=576, counts=counts, seconds=time.perf_counter() - start, physical_gate='UNKNOWN', external_referent=read(path)['external_referent'], all_numerical_controls_pass=True, all_injections_rejected=True, prereg_sha256=sha(path), resolution='PER_POINT', summary_cannot_certify='This exact plane containment cannot certify missing axial/fixture geometry; no global restoration claim')
    dump(ROOT / 'RESULTS_R4.json', result)
    freeze(ROOT / 'FROZEN_PREDICTIONS_R4.json', dict(result_sha256=sha(ROOT / 'RESULTS_R4.json'), prereg_sha256=sha(path), code_sha256=sha(Path(__file__)), physical_measurement_performed=False))
    print('fullbenchmark', counts, result['seconds'])
if __name__ == '__main__':
    import collections
    run()
