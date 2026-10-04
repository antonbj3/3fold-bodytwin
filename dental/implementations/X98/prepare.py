"""Prepare only the bounded, preregistered source snapshot; no numerical exploration."""
from dental_release.paths import expand as _release_expand
import datetime
import hashlib
import json
from pathlib import Path
import shutil
ROOT = Path(__file__).resolve().parent
DENT = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
R = DENT / 'results'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, x):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def load(p):
    return json.loads(Path(p).read_text())

def lines(p):
    return [json.loads(l) for l in Path(p).read_text().splitlines() if l.strip()]
if (ROOT / 'PREREG_R1.json').exists():
    raise SystemExit('R1 is already frozen. Use ./run_all.sh; preparation must not rewrite frozen sources or preregistration.')
manifest = []

def copy(p, local, transform=None):
    p = Path(p)
    dst = ROOT / local
    dst.parent.mkdir(parents=True, exist_ok=True)
    if transform:
        dst.write_text(transform(p.read_text()))
    else:
        shutil.copyfile(p, dst)
    manifest.append(dict(source=str(p), source_sha256=sha(p), local=local, sha256=sha(dst), bytes=dst.stat().st_size, transformation='none' if not transform else 'replace from common import * with import numpy as np'))
copy(R / 'LANE_X8_GUIDE_NERVE_RISK/full_geometry.py', 'implant_safety/vendor/x8_geometry.py')
copy(R / 'LANE_X87_NERVE_MARGIN_TABLE/science.py', 'implant_safety/vendor/x87_science.py')
copy(R / 'LANE_X96_DRILL_OVERSHOOT/code/guide.py', 'implant_safety/vendor/x96_guide.py', lambda s: s.replace('from common import *', 'import numpy as np'))
copy(R / 'PROOF_LANE_COMBINE_IMPLANT_SITE/code/control_geometry.py', 'implant_safety/vendor/proof_lane_control.py')
copy(R / 'LANE_X8_GUIDE_NERVE_RISK/GUIDE_PROFILES.json', 'data/guide_profiles.json')
copy(R / 'LANE_X96_DRILL_OVERSHOOT/PROTOCOLS.json', 'data/protocols.json')
copy(R / 'LANE_XREVIEW_BATCH32/REVIEW_LANE_X96_DRILL_OVERSHOOT.json', 'evidence/UPSTREAM_X96_REVIEW.json')
inter = load(R / 'PROOF_LANE_COMBINE_IMPLANT_SITE/raw/DENSE_BONE_INTERSECTION.json')
plans = {(s['case'], s['fdi']): s for s in lines(R / 'PROOF_LANE_COMBINE_IMPLANT_SITE/raw/SITE_PLANS_R1.jsonl')}
lineage = {(s['case'], s['fdi']): s for s in load(R / 'LANE_X73_CANAL_WALL_FULL/raw/LINEAGE_SITES.json')}
bone = {(s['case'], int(s['site'].split('_')[-1])): s for s in lines(R / 'LANE_X75_CBCT_BONE_SITE/SITE_INPUT_PORTS.jsonl') if s['site'].startswith('X8_FDI_')}
sites = []
for k in inter:
    key = (k['case'], k['fdi'])
    (p, c, b) = (plans[key], lineage[key], bone[key])
    assert p['pose'] == c['pose'], 'different source poses'
    dest = 'data/' + c['patient'] + '_points.npz'
    if not (ROOT / dest).exists():
        assert sha(c['point_artifact']['path']) == c['point_artifact']['sha256']
        copy(c['point_artifact']['path'], dest)
    sites.append(dict(site_id=k['case'] + '/FDI' + str(k['fdi']), case=k['case'], fdi=k['fdi'], image_group=k['image_group'], pose=c['pose'], frame=dict(order='zyx', units='mm', spacing_mm=[0.3, 0.3, 0.3], origin_mm=[0.0, 0.0, 0.0], direction=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), points_local=dest, points_sha256=sha(ROOT / dest), annotation_kind='RELEASE_REVISION_UNION', segment=c['segment'], external_source_bindings=c['source_bindings'], published_nominal_interval=[c['distances']['tf2']['lower_mm'], c['distances']['tf2']['upper_mm']], published_union_interval=[c['union']['lower_mm'], c['union']['upper_mm']], bone_binding=dict(pose=c['pose'], image_member_sha256=p['bone_image_member_sha256'], gray=p['gray'], regional_geometry=p['regional_geometry'], density_g_cm3=b['density_g_cm3'], cortex_mm=b['cortex_mm'], X68=b['X68'], X67=b['X67'], debts=b['debts'], locator=str(R / 'LANE_X75_CBCT_BONE_SITE/SITE_INPUT_PORTS.jsonl'))))
dump(ROOT / 'data/sites.json', sites)
for p in [R / 'PROOF_LANE_COMBINE_IMPLANT_SITE/raw/DENSE_BONE_INTERSECTION.json', R / 'PROOF_LANE_COMBINE_IMPLANT_SITE/raw/SITE_PLANS_R1.jsonl', R / 'LANE_X73_CANAL_WALL_FULL/raw/LINEAGE_SITES.json', R / 'LANE_X75_CBCT_BONE_SITE/SITE_INPUT_PORTS.jsonl']:
    manifest.append(dict(source=str(p), source_sha256=sha(p), bytes=p.stat().st_size, local=None, transformation='selected records in data/sites.json'))
manifest.append(dict(source='derived bounded source bundle', local='data/sites.json', sha256=sha(ROOT / 'data/sites.json'), bytes=(ROOT / 'data/sites.json').stat().st_size))
dump(ROOT / 'SOURCE_MANIFEST.json', manifest)
now = datetime.datetime.now(datetime.timezone.utc).isoformat()
prereg = dict(schema='dental-prereg-v1', round='X98-R1', frozen_utc=now, claim_type='capability', capability='One pose/system/length query with source-resolved annotation, tool and guide terms; digital interval vs 2 mm; physical unknowns remain explicit', obstacle='Reviewed parts use separate frames, pose caches, depth datums, confidence contracts and unknown physical ports', changed_operation='Compose reviewed geometry and bounds at exact query pose; telescope observed annotation/tool losses; never add channel sigma to revision union', consumer='Implantology researcher and release repository JSON API', selection=dict(sites='all 22 same-pose dense-revision + bone intersections frozen by PROOF_LANE', guides=['freehand', 'pilot_guided', 'fully_guided', 'dynamic_navigation'], systems=['Straumann_BLX_2024_VeloDrill', 'ProofLane_EV_Guided_2017'], lengths='original X8 physical cylinder lengths', confidence=0.95, restricted_system='ZimVie datum example requires actual length12.6mm, label13mm and platform offset1mm; never transferred to arbitrary lengths'), metrics=dict(nominal_replay_max_error_mm=1e-06, geometry_bracket_max_width_mm=1e-06, control_max_error_mm=1e-06, guide_control_max_error_mm=1e-10, source_identity_error=0, mutation_rejection_fraction=1.0, physical_claims_without_calibration=0), decision_criterion='All numerical and contract gates pass; any unsupported term produces UNKNOWN, never patient approval. Capability remains partial absent independent anatomy/pose/outcome.', strongest_equally_informed_control='Reviewed PROOF_LANE independent 6D SLSQP cylinder-box on every query; X87 Cantelli numerical inversion, X96 Markov/Cantelli scalar guide bound', practice='Fixed nominal annotated cylinder clearance >=2 mm, without guide/tool/revision information', falsifiers=['source/pose/datum mixup', 'sum entry + apex instead of max', 'recomputed answer differs from independent geometry >1e-6mm', 'unknown anatomy reported as physical certificate', 'injected wrong distance accepted'], external_referent=dict(kind='published_dataset', locator='https://doi.org/10.1016/j.media.2026.104095; X73 source member hashes; X96 primary protocol PDFs', compared_quantity='Cylinder-to-published-canal-voxel gap and source-reported extra tool depth; not true nerve safety', refutes_us=False), sufficiency_test='Same exact nominal gap 2.5mm to axial versus lateral voxel; same1mm tool extension; compare swept gap; retain direction/location as minimum geometry', full_cost=dict(preparation='timed source bundling; reasoning/source-read elapsed uninstrumented', fit=0, discovery='no new search/fit; previous full-cost debts inherited', validation='all 22 x2 tool geometries plus four guides; measured', queries_to_human=0, fallback='physical uncertainty UNKNOWN for every site; retain all rejected system/pose ports'), inputs_sha256=sha(ROOT / 'SOURCE_MANIFEST.json'), resources=dict(threads=1, gpu=False, intermediate_limit_bytes=3000000000), numerical_scope='Convex bracket under real arithmetic; rigorous IEEE directed rounding absent. No affine sensitivity reported.')
dump(ROOT / 'PREREG_R1.json', prereg)
(ROOT / 'PREREG_R1.json.sha256').write_text(sha(ROOT / 'PREREG_R1.json') + '\n')
dump(ROOT / 'DECOMPOSITION.json', dict(idea='Source-resolved minimum gap over planned cylinder, retained walls and possible guided displacement', mechanism='Changing occupied sets changes their minimum Euclidean separation', equation='d0=d(C_L,A_TF2); d1=d(C_L,A_union); d2=d(C_(L+h),A_union); conditional interval [max(0,d2_lower-G),d2_upper+G]; G=max(E,A)+2(r+h)sin(theta/2)', operation='Call X8 cylinder/voxel bracket, X87 guide component bounds and X96 protocol datums; carry X75 ports at matched pose', representation='PER_POINT voxel boxes + rigid PER_TOOTH pose; PHENOMENOLOGICAL guide moments; SIMULTANEOUS geometry', leaves=[dict(name=n, status=s, provenance=p, stop_argument=t) for (n, s, p, t) in [('published occupancy', 'EXTERNALLY_MEASURED', 'X73 source mask bytes', 'digital labels observed; deeper anatomical correctness UNKNOWN'), ('cylinder distance', 'DERIVED_UNDER_ASSUMPTIONS', 'X8 full_geometry.py', 'convex set minimum defines query; real arithmetic bracket only'), ('entry/apex/angle sample moments', 'EXTERNALLY_MEASURED', 'Varga2020 Table4 and Wu2020 Figure4', 'source sample endpoints; do not infer signed depth'), ('population transport and drill/implant error match', 'CONSTITUTIVE_CLOSURE', 'X87/X96', 'replace by paired registered drill/implant measurements; no data here'), ('tool extra depth', 'EXTERNALLY_MEASURED', 'X96 primary manuals', 'source specifies maximum or exact datum example; real tool shape UNKNOWN'), ('true canal wall and physical joint coverage', 'UNKNOWN', 'no independent wall/achieved drill pose', 'requires registered independent wall and signed achieved gap loss'), ('bone density/torque/temperature', 'UNKNOWN', 'X75/X68/X67 ports', 'gray is not density; no calibration or regional heat sensors')]]))
dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X98-implant-safety-module', updated_utc=now, stage='R1_FROZEN', latest_gate='PREREG_BEFORE_NUMERICS', next_operation='Implement and run composed module; own calculation precedes web countercheck'))
print(json.dumps(dict(selected_sites=len(sites), source_bundle_bytes=sum(((ROOT / r['local']).stat().st_size for r in manifest if r.get('local'))), prereg_sha256=sha(ROOT / 'PREREG_R1.json'))))
