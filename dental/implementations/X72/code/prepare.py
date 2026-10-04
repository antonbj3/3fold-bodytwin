from dental_release.paths import expand as _release_expand
from common import *
import subprocess
F = _release_expand('@DENTAL_EXTERNAL_ROOT@/engines/3fold-field-engine')
M = _release_expand('@DENTAL_EXTERNAL_ROOT@/engines/3fold-motion-engine')
groups = [(F, '1216032', ['src/field_engine/contact_port_v1.py', 'src/field_engine/contact_port_adapter.py', 'docs/contact_port_v1.md', 'docs/contact_port_v1.schema.json', 'tests/test_contact_port_v1.py']), (F, 'e590aa3', ['src/field_engine/experimental/farkas_witness.py', 'src/field_engine/experimental/_farkas_check.py', 'docs/FARKAS_WITNESS_API.md', 'tests/test_farkas_witness.py']), (F, '7abf765', ['src/field_engine/experimental/normal_cone.py', 'docs/normal_cone.md', 'tests/test_normal_cone.py']), (M, '6e2defd', ['src/motion_engine/ncp/contact_port_adapter.py', 'docs/contact_port_adapter.md', 'tests/test_contact_port_adapter.py']), (M, '95ec803', ['src/motion_engine/ncp/family_sweep.py', 'src/motion_engine/ncp/latent_contact.py', 'tests/test_family_sweep.py', 'tests/test_latent_contact.py'])]

def main():
    DATA.mkdir(parents=True, exist_ok=True)
    rows = []
    for (repo, commit, files) in groups:
        full = subprocess.check_output(['git', '-C', repo, 'rev-parse', commit], text=True).strip()
        for file in files:
            cmd = ['git', '-C', repo, 'show', full + ':' + file]
            b = subprocess.check_output(cmd)
            path = ROOT / 'vendor' / (file[4:] if file.startswith('src/') else commit + '/' + file)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b)
            rows.append(dict(repo=repo, commit=full, file=file, export=str(path.relative_to(ROOT)), sha256=sha(path), command=cmd))
    for p in ['field_engine/__init__.py', 'field_engine/experimental/__init__.py', 'motion_engine/__init__.py', 'motion_engine/ncp/__init__.py']:
        (ROOT / 'vendor' / p).write_text('')
    freeze(ROOT / 'VENDOR_MANIFEST.json', rows)
    support = []
    for (repo, commit, files) in [(F, '1216032', ['src/field_engine/univariate_equilibrium.py']), (M, '6e2defd', ['src/motion_engine/ncp/ncp_ref.py', 'src/motion_engine/ncp/branch_set.py'])]:
        full = subprocess.check_output(['git', '-C', repo, 'rev-parse', commit], text=True).strip()
        for file in files:
            cmd = ['git', '-C', repo, 'show', full + ':' + file]
            b = subprocess.check_output(cmd)
            path = ROOT / 'vendor' / file[4:]
            path.write_bytes(b)
            support.append(dict(repo=repo, commit=full, file=file, export=str(path.relative_to(ROOT)), sha256=sha(path), command=cmd))
    freeze(ROOT / 'VENDOR_SUPPORT_MANIFEST.json', support)
    repo = M
    commit = '6e2defd'
    file = 'src/motion_engine/ncp/contact_rank.py'
    full = subprocess.check_output(['git', '-C', repo, 'rev-parse', commit], text=True).strip()
    cmd = ['git', '-C', repo, 'show', full + ':' + file]
    path = ROOT / 'vendor' / file[4:]
    path.write_bytes(subprocess.check_output(cmd))
    freeze(ROOT / 'VENDOR_SUPPORT_MANIFEST_V2.json', [*support, dict(repo=repo, commit=full, file=file, export=str(path.relative_to(ROOT)), sha256=sha(path), command=cmd)])
    for lane in [_release_expand('X54'), _release_expand('X2'), _release_expand('X10'), _release_expand('X49'), _release_expand('X18')]:
        for p in (BASE / lane / 'code').rglob('*.py'):
            record(p)
        for p in (BASE / lane).glob('PREREG*.json'):
            record(p)
        record(BASE / lane / 'results.json')
        record(BASE / lane / 'RESULTS.md')
    originals = [('LANE_X54_DEFORMABLE_CONTACT/code/run_r3.py', 'run_r3_count_guard.py'), ('LANE_X2_OCCLUSION_B2B/code/summarize.py', 'summarize_18.py')]
    for (src, dst) in originals:
        s = record(BASE / src).read_text()
        if 'run_r3' in src:
            s = s.replace("result=dict(round='R3'", "assert len(rows)>0 and all(len(r['checks']['controls'])>0 for r in rows), 'EMPTY_CONTROL_SELECTION'\n    result=dict(round='R3'")
        else:
            s = s.replace('.05/12', '.05/18')
        p = ROOT / 'corrected_sources' / dst
        p.parent.mkdir(exist_ok=True)
        p.write_text(s)
    specs = {'X54': dict(capability='Exact lower-gap and infeasibility witnesses for stored winning triangle pairs', obstacle='Vacuous all and floating LP cannot certify geometry', changed_operation='Nonempty guard; exact rational vertex LP and Farkas dual check', consumer='regional gap/contact force port', metric='Every named pair gets primal optimum and strict dual for 0.001 mm too-low gap', tolerance=0, gate='all pairs verified; empty list and tampered dual rejected', control='scipy linprog same stored triangle pair', falsifier='exact optimum differs >1e-8 mm or injected bound passes'), 'X2': dict(capability='Recompute all 18 Angle contrasts with corrected multiplicity', obstacle='family frozen for12 but evaluated18', changed_operation='same RNG, 2000 replicates, alpha=.05/18', consumer='population Angle contrast', metric='changed significance and original five-pp effect gate per contrast', tolerance=0, gate='18 contrasts per round; exact replay of old intervals; injection wrong family rejected', control='original alpha .05/12 from identical raw cohorts', falsifier='input or seed drift; family guard accepts12'), 'CORNERS': dict(capability='Distinguish licensed box guarantees from sampled checks', obstacle='corner sampling alone has no interior guarantee', changed_operation='derive coordinate monotonicity per quantity; exact family_sweep for affine fixed-XY PT path; unsigned norms labelled sampling', consumer='X10 seated margin and X49 design gate', metric='exact analytic bound containment / sweep replay', tolerance=0, gate='every reported guaranteed quantity has license; unlicensed quantities downgraded; corrupt cert rejected', control='original corners and sampled perturbations', falsifier='interior exceeds claimed envelope or unsigned norm called multiaffine'), 'PORT': dict(capability='Transport X10 same seating cases using Field and Motion contact port v1', obstacle='dental bespoke semantic port', changed_operation='canonical Gap in mm, conditional event, provenance, branch coverage replayer', consumer='motion contact port', metric='bounds roundtrip exactly; 120um threshold classifications and semantic refusals', tolerance=0, gate='all6 cases roundtrip; source, unit and replay mutations rejected', control='X10 existing signed_gap.certify_port same cases', falsifier='encoding loses condition or sigma UNKNOWN or physical promotion'), 'X18': dict(capability='Conditional stress range with local body force budgets, friction and inclined normals', obstacle='100N independently on every tooth and vertical-only simplex', changed_operation='global100N scene budget, per-body20N caps; Coulomb inner/outer polytope stress extrema; exact normal cone escape direction', consumer='crown stress/force-measurement decision', metric='old, locally budgeted, tilted and frictional upper stress and uncertainty ratio', tolerance=1e-08, gate='FE basis parity <=1e-8 relative; direct mixed solve and conservation; inner<=outer; cap/cone/tensor corruptions rejected', control='original vertical 100N force simplex identical roof; conventional direct linear FE', falsifier='budget/cone violation survives; physics accuracy cannot be claimed without stress facit')}
    for (name, s) in specs.items():
        s.update(claim_type='capability', frozen_utc=now(), lane='X72-falt-fixes', external_referent=dict(kind='published_code', locator='Field/Motion commits in VENDOR_MANIFEST.json; source datasets/papers retain own locators', compared_quantity=s['metric'], refutes_us=True), cost=dict(preparation='measured run time plus existing source-read time unmeasured', fit='none', discovery='unmeasured human/agent reasoning', validation='measured per correction', questions='zero human questions', fallback='UNKNOWN is explicit; no empirical substitute'), physical_scope='CONDITIONAL_MODEL, not clinical or newly measured', resolution={'X2': 'POPULATION', 'X54': 'PER_SURFACE_REGION', 'CORNERS': 'PER_SURFACE_REGION', 'PORT': 'PER_SURFACE_REGION', 'X18': 'PER_SURFACE_REGION'}[name])
        if name == 'X18':
            s.update(force_budget_closure=dict(global_vertical_N=100, per_body_vertical_N=20, friction_mu='1/5', friction_polygon_sides=8), unchanged_decision_ratio=1.5, rigorous_FE_roundoff_enclosure='MISSING; interval arithmetic for FE and eigenvalues not implemented')
        freeze(ROOT / ('PREREG_' + name + '.json'), s)
    leaves = [dict(idea='triangle projected gap', mechanism='equal XY barycentric weights', equation='A x=b,x>=0; minimize c.x; y.A<=c gives lower bound', operation='rational basis enumeration and exact Farkas sum', representation='decimal JSON coordinates interpreted exactly', status='DERIVED_UNDER_ASSUMPTIONS', stop='all6 variables retained; rounded geometry/calibration and global omitted pairs UNKNOWN'), dict(idea='simultaneous contrasts', mechanism='union bound', equation='P(any error)<=sum_i alpha_i=18*(.05/18)', operation='same seeded bootstrap quantiles', representation='population area and conditional mechanical shares', status='DERIVED_UNDER_ASSUMPTIONS', stop='bootstrap marginal coverage not finite-sample exact; annotation/independence UNKNOWN'), dict(idea='seating bounds', mechanism='maximum of contact lifts', equation='h=max(0,-b_i/a_i); g_j=max(0,b_j,max_i!=j[b_j-a_j*b_i/a_i])', operation='own coordinate increases; other coordinates decrease; extrema licensed', representation='independent signed clearance box and positive fixed projections', status='DERIVED_UNDER_ASSUMPTIONS', stop='taper, full-wall bounds, shared datum and manufacturing calibration are required'), dict(idea='geometry error', mechanism='triangle inequality or barycentric convexity', equation='|d(A,B)-d(A0,B0)|<=eA+eB; |min g-min g0|<=eU+eL', operation='set-error or fixed XY axial enclosure', representation='whole-surface deterministic Hausdorff bound or per-vertex signed axial bound', status='DERIVED_UNDER_ASSUMPTIONS', stop='published medians and SD are not such bounds; FP rounding enclosure MISSING'), dict(idea='contact port', mechanism='preserve units/provenance/conditional completeness', equation='1um=1/1000mm; physical UNKNOWN stays UNKNOWN', operation='encode/decode and full consumer replay', representation='v1 Gap+Witness+BranchSet', status='DERIVED_UNDER_ASSUMPTIONS', stop='codec cannot derive seating/normalization/calibrated sigma'), dict(idea='linear roof stress', mechanism='elastic equilibrium and Coulomb forces', equation='K u=f; sigma=D B u; ||t||<=mu n; lambda_max is convex', operation='stress basis, product of local capped cones, vertex maximization', representation='tetrahedra and registered proximity patches', status='CONSTITUTIVE_CLOSURE', stop='E,nu, thickness, support, uniform patch traction,mu and20N cap PHENOMENOLOGICAL; need measured local load-pressure and compliance; FE continuum error UNKNOWN')]
    freeze(ROOT / 'DECOMPOSITION.json', dict(chain='reliable mesh→gap→force→stress→design decision', leaves=leaves))
    state('PREREG_FROZEN', 'Five correction contracts frozen before numerical evaluation', 'Run exact X54/X2/corner/port corrections and locally budgeted X18 FE')
if __name__ == '__main__':
    main()
