from dental_release.paths import expand as _release_expand
from common import *
from stats import summarize
from importlib.util import spec_from_file_location, module_from_spec
from scipy.optimize import minimize
start = time.perf_counter()
check()
rows = [json.loads(l) for l in (DATA / 'R1_GEOMETRY.jsonl').read_text().splitlines()]
ss = summarize(1.5, [r for r in rows if r['extra_mm'] == 1.5])
spec = spec_from_file_location('independent_control', DENT / 'results/PROOF_LANE_COMBINE_IMPLANT_SITE/code/control_geometry.py')
mod = module_from_spec(spec)
spec.loader.exec_module(mod)
selected = sorted({r['image_group'] for r in rows}, key=lambda g: hashlib.sha256(('X96-control' + g).encode()).hexdigest())[:8]
controls = []
import zipfile
with zipfile.ZipFile(tf2_io.ZIP) as z:
    for group in selected:
        rr = sorted([r for r in rows if r['image_group'] == group and r['extra_mm'] == 1.5], key=lambda r: (r['case'], r['fdi']))[0]
        arrays = []
        for src in rr['point_sources']:
            if 'path' in src:
                with np.load(src['path']) as a:
                    arrays.extend((a[k].astype(float) * 0.3 for k in a.files if k.endswith('_voxels_zyx')))
            else:
                (lab, sp, _) = tf2_io.read_mha_bytes(z.read(src['member']))
                arrays = [np.argwhere(np.isin(lab, [3, 4])) * 0.3]
        pt = np.unique(np.vstack(arrays), axis=0)
        e = np.array(rr['entry_zyx_mm'])
        a = np.array(rr['axis_zyx'])
        L = rr['length_mm'] + 1.5
        r = rr['radius_mm']
        from full_geometry import project_cylinder
        dc = np.linalg.norm(pt - project_cylinder(pt, e, a, L, r), axis=1)
        ix = np.where(dc - 0.3 * np.sqrt(3) / 2 <= rr['union_upper_mm'] + 1e-06)[0]
        vals = [mod.control_box(pt[i], e, a, L, r, np.full(3, 0.3)) for i in ix]
        d = min((v[0] for v in vals))
        err = max(0, rr['union_lower_mm'] - d, d - rr['union_upper_mm'])
        controls.append(dict(image_group=group, case=rr['case'], fdi=rr['fdi'], candidate=[rr['union_lower_mm'], rr['union_upper_mm']], control_mm=d, error_mm=err, all_optimizer_success=all((v[1] for v in vals)), injected_plus1mm_rejected=abs(rr['union_upper_mm'] + 1 - d) > 1e-06))
e = np.zeros(3)
a = np.array([1.0, 0, 0])
points = np.array([[12.0, 0, 0], [5.0, 4.0, 0]])
from full_geometry import project_cylinder
d0 = np.linalg.norm(points - project_cylinder(points, e, a, 10.0, 2.0), axis=1)
dh = np.linalg.norm(points - project_cylinder(points, e, a, 11.5, 2.0), axis=1)
suff = dict(summary='nominal minimum clearance', summary_values_mm=d0.tolist(), identity_error_mm=float(abs(d0[0] - d0[1])), downstream_extra1p5_clearance_mm=dh.tolist(), downstream_difference_mm=float(abs(dh[0] - dh[1])), smallest_extension='signed distance loss for frozen pose; wall position/direction for new queries', external_referent=dict(kind='our_own_fixture', locator='code/summarize_r1.py', compared_quantity='sufficiency witness only', refutes_us=True))
dump(ROOT / 'raw/R1_CONTROLS.json', controls)
dump(ROOT / 'raw/SUFFICIENCY_R1.json', suff)
res = dict(round='R1', claim_type='information_link', outcome='DIGITAL_DRILL_SWEEP_CONSEQUENCE_MEASURED_PHYSICAL_FREQUENCY_UNKNOWN', summary=ss, gates=dict(coverage=len(rows) == 2 * len(certs()), practice_change=ss['all']['fraction'] >= 0.05, control=all((r['error_mm'] <= 1e-06 and r['injected_plus1mm_rejected'] for r in controls)), summary_sufficiency=suff['downstream_difference_mm'] == 0), external_referent=load(ROOT / 'PREREG_R1.json')['external_referent'], cost=dict(geometry=load(ROOT / 'raw/R1_COST.json'), validation_wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss), formal_floating_point_enclosure='MISSING')
dump(ROOT / 'rounds/R1.json', res)
(ROOT / 'HANDOFF_R1.md').write_text('# R1 digitalt borrsvep\n\n' + json.dumps(res, indent=2) + _release_expand("\n\nFirst local pre-web search measurement: @DENTAL_CASE_ID@ / 44 , nominally 2.962229729532198 mm , extra 1.5 mm fullradiesvep 1.9092399240467286 mm . R1 is a prescribed housing, not a manufacturer-identified drilling length. The next design carries protocol version and reference point for length.\n"))
state('R1_COMPLETE', res['gates'], 'Freeze manufacturer/protocol datum and guide partial identification')
print(json.dumps(res['summary'], indent=2))
