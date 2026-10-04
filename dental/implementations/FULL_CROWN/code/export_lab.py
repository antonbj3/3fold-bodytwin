from fc_common import *
import zipfile, xml.etree.ElementTree as ET, shutil
import trimesh
NS = 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'

def three_mf(path, v, f, title):
    model = ET.Element('model', unit='millimeter', xmlns=NS)
    metadata = ET.SubElement(model, 'metadata', name='Title')
    metadata.text = title
    res = ET.SubElement(model, 'resources')
    obj = ET.SubElement(res, 'object', id='1', type='model')
    mesh = ET.SubElement(obj, 'mesh')
    verts = ET.SubElement(mesh, 'vertices')
    for p in v:
        ET.SubElement(verts, 'vertex', **{k: format(float(a), '.17g') for (k, a) in zip(['x', 'y', 'z'], p)})
    faces = ET.SubElement(mesh, 'triangles')
    for (a, b, c) in f:
        ET.SubElement(faces, 'triangle', v1=str(a), v2=str(b), v3=str(c))
    build = ET.SubElement(model, 'build')
    ET.SubElement(build, 'item', objectid='1')
    payload = {'3D/3dmodel.model': ET.tostring(model, encoding='utf-8', xml_declaration=True), '[Content_Types].xml': b'<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>', '_rels/.rels': b'<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>'}
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for (name, data) in payload.items():
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 3, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, data)

def parse(p):
    with zipfile.ZipFile(p) as z:
        root = ET.fromstring(z.read('3D/3dmodel.model'))
    assert root.attrib['unit'] == 'millimeter'
    ns = {'m': NS}
    vv = root.findall('.//m:vertex', ns)
    ff = root.findall('.//m:triangle', ns)
    return (np.array([[float(x.attrib[k]) for k in ['x', 'y', 'z']] for x in vv]), np.array([[int(x.attrib[k]) for k in ['v1', 'v2', 'v3']] for x in ff]))

def run():
    rows = []
    for tag in ['R6', 'R7']:
        for r in read(ROOT / 'rounds' / f'{tag}.json')['rows']:
            if r.get('status') == 'SCORED' and r.get('certificate', {}).get('shell_components') == 1:
                rows.append(dict(r, round=tag))
    dest = DATA / 'lab_exports'
    dest.mkdir(exist_ok=True)
    selected = []
    checks = []
    for family in ['molar_crown', 'premolar_crown', 'anterior_crown']:
        rr = [r for r in rows if r['family'] == family]
        if not rr:
            selected.append(dict(family=family, status='NO_VALID_COMPLETE_PAIR'))
            continue
        r = min(rr, key=lambda x: (x['reconstruction_p95_mm'], x['key'], x['round']))
        base = DATA / (r['round'] + '_predictions') / r['participant'] / r['key']
        m = npz(base / 'mesh.npz')
        p = npz(DATA / 'public' / (r['key'] + '.npz'))
        out = dest / family
        out.mkdir(exist_ok=True)
        v = m['vertices'] @ p['source_R'].T + p['source_base']
        f = m['faces']
        three_mf(out / 'crown.3mf', v, f, 'Research crown; full-prescan joint design; see CASE.json for failed gates')
        (cv, cf) = parse(out / 'crown.3mf')
        tm = trimesh.Trimesh(cv, cf, process=False)
        checks.append(dict(family=family, coordinate_identity=bool(np.array_equal(cv, v)), faces_identical=bool(np.array_equal(cf, f)), closed=bool(tm.is_watertight and tm.is_winding_consistent and (tm.volume > 0)), removed_face_rejected=not trimesh.Trimesh(cv, cf[:-1], process=False).is_watertight))
        die = trimesh.load(base / 'die.stl', process=True)
        three_mf(out / 'die.3mf', die.vertices, die.faces, 'Research matching die; not a clinical preparation')
        (dv, df) = parse(out / 'die.3mf')
        checks[-1].update(die_vertices_identical=bool(np.array_equal(dv, die.vertices)), die_faces_identical=bool(np.array_equal(df, die.faces)), die_closed=bool(trimesh.Trimesh(dv, df, process=False).is_watertight))
        for name in ['crown.stl', 'die.stl', 'CERTIFICATE.json']:
            shutil.copy2(base / name, out / name)
        card = dict(key=r['key'], family=family, round=r['round'], information_track='FULL_PREOPERATIVE_SCAN_JOINT_PREPARATION', source_fdi=r['source_fdi'], p95_mm=r['reconstruction_p95_mm'], anatomy_gate_pass=r['reconstruction_p95_mm'] <= 0.35, certificates=r['certificate'], function=r['function'], clinical_and_physical_status='UNVALIDATED_RESEARCH_SPECIMEN', licence='Local source-specific research; Bits2Bites CC BY-NC-SA per brief, Bite2Text exact licence UNKNOWN; no redistribution assumption')
        dump(out / 'CASE.json', card)
        selected.append(dict(card, path=str(out)))
    good = next(dest.rglob('crown.3mf'), None)
    if good is not None:
        bad = DATA / 'negative_controls/bad_unit.3mf'
        bad.parent.mkdir(exist_ok=True)
        with zipfile.ZipFile(good) as z:
            content = {n: z.read(n) for n in z.namelist()}
        content['3D/3dmodel.model'] = content['3D/3dmodel.model'].replace(b'unit="millimeter"', b'unit="inch"')
        with zipfile.ZipFile(bad, 'w') as z:
            for (name, value) in content.items():
                info = zipfile.ZipInfo(name, date_time=(2026, 10, 3, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                z.writestr(info, value)
        try:
            parse(bad)
            reject = False
        except AssertionError:
            reject = True
        checks.append(dict(family='fault_unit', bad_unit_rejected=reject))
    manifest = {str(p.relative_to(dest)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in dest.rglob('*') if p.is_file()}
    record = dict(selected=selected, checks=checks, all_pass=all((all((v for (k, v) in c.items() if k != 'family')) for c in checks)), files=manifest)
    if not (ROOT / 'FROZEN_LAB_EXPORTS.json').exists():
        freeze(ROOT / 'FROZEN_LAB_EXPORTS.json', record)
    else:
        old = read(ROOT / 'FROZEN_LAB_EXPORTS.json')
        assert old['files'] == manifest, 'frozen lab exports drift'
    dump(ROOT / 'raw/EXPORT_VALIDATION.json', record)
    print('3MF export checks', record['all_pass'])
    return record
if __name__ == '__main__':
    run()
