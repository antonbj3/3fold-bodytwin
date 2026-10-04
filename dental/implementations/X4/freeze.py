from dental_release.paths import expand as _release_expand
import datetime, hashlib, json, os
from pathlib import Path
ROOT = Path(__file__).resolve().parent
OLD = Path(__import__('os').environ.get('DENTAL_PROJECT_ROOT', str(ROOT.parents[1]))) / 'results/LANE_NEXT_G_MANDIBLE_POSTOP'
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X4_fibula_planner'))

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f'.tmp.{os.getpid()}')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def state(milestone, next_operation, gate=None):
    if gate is None:
        current = ROOT / 'CURRENT_WORK_STATE.json'
        gate = json.loads(current.read_text()).get('latest_gate', 'NOT_RUN') if current.exists() else 'NOT_RUN'
    write(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X4-fibula-planner', updated_utc=now(), milestone=milestone, latest_gate=gate, next_operation=next_operation, review_state='PENDING_INDEPENDENT_REVIEW', data_dir=str(DATA), intermediate_limit_bytes=3000000000))

def freeze(round_id, information, change, gates):
    out = ROOT / f'PREREG_{round_id}.json'
    if out.exists():
        return
    man = json.loads((OLD / 'CASE_MANIFEST.json').read_text())
    cases = sorted((c for (c, r) in man['roles'].items() if all((k in r and r[k]['exists'] for k in ('Pre', 'Post', 'Original')))))
    write(ROOT / 'SOURCE_MANIFEST.json', dict(dataset=man['source_doi'], descriptor='https://doi.org/10.1038/s41597-025-06048-8', license='CC BY 4.0 per local _metadata.json', case_manifest=str(OLD / 'CASE_MANIFEST.json'), case_manifest_sha256=sha(OLD / 'CASE_MANIFEST.json'), population=cases, n_complete=len(cases), source_file_count=410, n_designated_cases=147, Original='original defective anatomy', Pre='prepared defective anatomy', Post='virtual expert annotation', reuse='Read-only R2 component cache, anatomy proxies independently unreviewed', noncomplete_cases='excluded prospectively; keep 118/147 denominator', predecessor=str(OLD / 'RESULTS.md')))
    reg = dict(id=f'X4_{round_id}', frozen_utc=now(), capability='1-3 straight circular graft surrogate segments, lengths, osteotomy planes, STL and implant-scenario export', obstacle='a completed surface does not specify limited straight grafts; absent bilateral curvature is not identifiable from Pre', changed_operation=change, consumer='CAD/CAM researcher testing printed geometric prototypes', population=cases, information=information, parameters=dict(angle_bins=73, angle_range_degrees=[-100, 100], inferior_quantile=0.08, curve_band_quantile=0.22, gap_fraction=0.5, gap_min_samples=20, radius_mm=6.0, roi_distance_mm=3.0, minimum_segment_length_mm=10.0, max_segments=3, max_curve_error_mm=2.0, cost_tolerance_mm=1e-08, sample_count_per_segment=3000, register_trim=0.9, register_iterations=35), metrics=dict(primary='max directional p95 sampled surface distance: expert-added surface to Pre+graft; graft outside Pre to full Post', secondary=['inferior-rail chord maximum', 'graft length vs expert ROI-rail arc length', 'bend and half-bend miter angle', 'Original to preserved Pre sampled distance', 'condyle landmark error UNKNOWN without labels', 'occlusal plane UNKNOWN without supplied plane'], validity='ROI>=200 points; 1-3 segments each>=10mm; no dropped cases; UNKNOWN fails complete cohort'), gates=gates, strongest_equal_information_control='independent exhaustive breakpoint enumeration with same curve, endpoints, radii and segment count', practice_reference='Pre skull-plane sagittal reflection proxy, not full clinical workflow; uniform arclength split is only a weaker reference', falsifiers=['20mm translation fails distance gate', '10x reported length fails length gate', 'broken cut plane fails plane residual', 'exhaustive contradicts DP', 'missing cases fail completeness'], full_cost=dict(preparation='hash and cache reads timed; predecessor preparation separately charged', fit='all curve extraction and fits timed', discovery='all implementation/debug costs logged, wall time reported', validation='registration plus all reference queries', queries='per case', fallback='UNKNOWN retained; no manual correction; no physical or clinical time saving claimed'), external_referent=dict(kind='published_dataset', locator='https://doi.org/10.6084/m9.figshare.28052240.v2', compared_quantity='virtual expert completed surface; Original defective preservation only', refutes_us=False), scientific_admission=False)
    write(out, reg)
    (ROOT / f'PREREG_{round_id}.sha256').write_text(sha(out) + '  ' + out.name + '\n')
    leaves = [('surface', 'EXTERNALLY_MEASURED', 'CT segmentation and expert virtual CAD annotation with fixed source hashes', 'no intact truth or actual outcome'), ('component mask', 'CONSTITUTIVE_CLOSURE', 'inherited inferior-to-skull connected component proxy', 'no anatomy labels; may include teeth'), ('angle ordering', 'CONSTITUTIVE_CLOSURE', 'theta=atan2(x-sagittal_x,-(y-anchor_y)); superior=z assumed', 'not landmark calibration'), ('inferior rail', 'CONSTITUTIVE_CLOSURE', 'median xy of lowest22%; q08(z)+radius per angular bin', 'proxy contour, not medial axis proof'), ('minimax DP', 'DERIVED_UNDER_ASSUMPTIONS', 'D[k,j]=min_i max(D[k-1,i],E[i,j]), E=max chord distance', 'finite sampled rail only'), ('miter', 'DERIVED_UNDER_ASSUMPTIONS', 'n=(u_left+u_right)/norm(...); ellipse major=r/abs(u.n)', 'equal-radius straight circular surrogate, no cortex or pedicle'), ('donor geometry', 'UNKNOWN', 'r=6mm research scenario, no matched donor CT', 'cannot declare harvest or fixation feasibility'), ('occlusal/condyle', 'UNKNOWN', 'no labelled landmarks; superior-plane scenario only', 'surface extrema are proxies'), ('mechanics/union', 'UNKNOWN', 'no load, fixation or physical union reference', 'uncalibrated FE is not validation')]
    write(ROOT / f'DECOMPOSITION_{round_id}.json', dict(idea='make a target shape actionable with straight grafts', mechanism='concentrate curvature into joints with common cut planes', operation=change, equations=['point-segment projection distance in mm', 'Bellman minimax recurrence', 'cylinder-plane intersection', 'min union of implicit fields'], representation='3D rail plus endpoint/cut-plane records and analytically clipped circular cylinders', assumptions='coordinate and component proxies; no clinical manufacturing claim', leaves=[dict(name=n, status=s, basis=b, stopping_argument=t) for (n, s, b, t) in leaves]))
if __name__ == '__main__':
    DATA.mkdir(parents=True, exist_ok=True)
    freeze('R1', dict(prediction=['Pre only'], evaluation=['Post', 'Original'], target_Post_forbidden_until_freeze=True), 'Pre mirror -> angular inferior rail -> minimax DP -> straight mitered circular graft -> STL -> Post comparison', dict(median_primary_p95_mm_max=3.0, paired_gain_vs_uniform_mirror_mm_min=0.5, all_cases_available=True, solver_max_absolute_disagreement_mm=1e-08, every_gate_fault_injection=True))
    state('R1_preregistered', 'implement and freeze all Pre-only predictions before evaluation')
