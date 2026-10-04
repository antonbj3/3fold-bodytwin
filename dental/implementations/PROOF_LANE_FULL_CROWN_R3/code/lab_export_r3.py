from common_r3 import *
import sys, time, shutil, zipfile, xml.etree.ElementTree as ET
sys.path.insert(0, str(OLD / 'code'))
sys.path.insert(0, str(OLD / 'code/vendor'))
import export_format as fmt
import trimesh

def check3mf(path, v, f):
    (vv, ff) = fmt.parse(path)
    if not np.array_equal(vv, v) or not np.array_equal(ff, f):
        raise ValueError('coordinate/index mismatch')
    m = trimesh.Trimesh(vv, ff, process=False)
    if not m.is_watertight or not m.is_winding_consistent or m.volume <= 0:
        raise ValueError('topology mismatch')
    return True

def run():
    start = time.perf_counter()
    candidates = []
    for tag in ['A', 'B', 'C']:
        for r in read(ROOT / 'rounds' / f'{tag}.json')['rows']:
            if r['status'] == 'SCORED' and r['digital_closed_shell']:
                candidates.append(dict(round=tag, **r))
    out = DATA / 'lab_exports'
    out.mkdir(exist_ok=True)
    selected = []
    controls = []
    for family in ['molar_crown', 'premolar_crown', 'anterior_crown']:
        rr = [r for r in candidates if r['family'] == family]
        if not rr:
            selected.append(dict(family=family, status='NO_CLOSED_CANDIDATE'))
            continue
        r = min(rr, key=lambda r: (-int(r['anatomy_pass'] and r['functional_pass_sampled']), -int(r['positive_native_band_restored']), -int(r['functional_pass_sampled']), -sum(r['gates'].values()), r['reconstruction_p95_mm'], r['round'], r['participant'], r['key']))
        source = DATA / r['round'] / r['participant'] / r['key']
        m = npz(source / 'mesh.npz')
        p = npz(DATA / 'inputs' / (r['key'] + '.npz'))
        v = m['vertices'] @ p['source_R'].T + p['source_base']
        f = m['faces']
        dest = out / family
        dest.mkdir(exist_ok=True)
        shutil.copyfile(source / 'crown.stl', dest / 'crown.stl')
        fmt.three_mf(dest / 'crown.3mf', v, f, 'RESEARCH SPECIMEN; full anatomy and physical function unresolved')
        check3mf(dest / 'crown.3mf', v, f)
        loaded = trimesh.load(dest / 'crown.stl', process=False)
        actual = np.asarray(loaded.triangles)
        expected = v[f]
        err = float(abs(actual - expected).max())
        assert err <= 2e-05
        welded = trimesh.load(dest / 'crown.stl', process=True)
        stl_closed = bool(welded.is_watertight and welded.is_winding_consistent and (welded.volume > 0))
        selected.append(dict(family=family, round=r['round'], key=r['key'], participant=r['participant'], status='DIGITAL_JOINT_PASS' if r['anatomy_pass'] and r['functional_pass_sampled'] else 'REJECTED_RESEARCH_SPECIMEN', anatomy_p95_mm=r['reconstruction_p95_mm'], anatomy_pass=r['anatomy_pass'], functional_pass_sampled=r['functional_pass_sampled'], positive_native_band_restored=r['positive_native_band_restored'], gates=r['gates'], stl_roundtrip_max_error_mm=err, stl_welded_closed=stl_closed, three_mf_bitidentical=True, files={name: dict(path=str(dest / name), sha256=sha(dest / name), bytes=(dest / name).stat().st_size) for name in ['crown.stl', 'crown.3mf']}, scope='Virtual joint preparation; physical manufacture/fit/load not tested', resolution='PER_TOOTH'))
        with zipfile.ZipFile(dest / 'crown.3mf') as z:
            contents = {n: z.read(n) for n in z.namelist()}
        for fault in ['unit', 'triangle']:
            root = ET.fromstring(contents['3D/3dmodel.model'])
            if fault == 'unit':
                root.set('unit', 'inch')
            else:
                tr = root.find('.//{' + fmt.NS + '}triangles')
                tr.remove(list(tr)[0])
            temp = out / (family + '__bad_' + fault + '.3mf')
            with zipfile.ZipFile(temp, 'w', zipfile.ZIP_DEFLATED) as z:
                for (name, data) in contents.items():
                    z.writestr(name, ET.tostring(root) if name == '3D/3dmodel.model' else data)
            rejected = False
            try:
                check3mf(temp, v, f)
            except (AssertionError, ValueError):
                rejected = True
            assert rejected
            controls.append(dict(family=family, fault=fault, rejected=rejected))
            temp.unlink()
    result = dict(selected=selected, controls=controls, seconds=time.perf_counter() - start, selection_prereg_sha256=sha(ROOT / 'PREREG_EXPORT_SELECTION.json'), physical_test='NOT_RUN')
    dump(ROOT / 'LAB_EXPORTS.json', result)
    if not (ROOT / 'FROZEN_LAB_PREDICTIONS.json').exists():
        freeze(ROOT / 'FROZEN_LAB_PREDICTIONS.json', result)
    return result
if __name__ == '__main__':
    run()
