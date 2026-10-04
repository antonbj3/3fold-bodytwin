from geometry import *
import zipfile, xml.etree.ElementTree as ET
NS = 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'

def three_mf(path, v, f, title):
    model = ET.Element('model', unit='millimeter', xmlns=NS)
    ET.SubElement(model, 'metadata', name='Title').text = title
    res = ET.SubElement(model, 'resources')
    obj = ET.SubElement(res, 'object', id='1', type='model')
    mm = ET.SubElement(obj, 'mesh')
    vv = ET.SubElement(mm, 'vertices')
    for p in v:
        ET.SubElement(vv, 'vertex', **{k: format(float(a), '.17g') for (k, a) in zip(['x', 'y', 'z'], p)})
    ff = ET.SubElement(mm, 'triangles')
    for (a, b, c) in f:
        ET.SubElement(ff, 'triangle', v1=str(a), v2=str(b), v3=str(c))
    ET.SubElement(ET.SubElement(model, 'build'), 'item', objectid='1')
    payload = {'3D/3dmodel.model': ET.tostring(model, encoding='utf-8', xml_declaration=True), '[Content_Types].xml': b'<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>', '_rels/.rels': b'<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>'}
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for (name, b) in payload.items():
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 4, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, b)

def load3mf(p):
    with zipfile.ZipFile(p) as z:
        root = ET.fromstring(z.read('3D/3dmodel.model'))
    assert root.attrib['unit'] == 'millimeter'
    ns = {'m': NS}
    return (np.array([[float(x.attrib[k]) for k in ['x', 'y', 'z']] for x in root.findall('.//m:vertex', ns)]), np.array([[int(x.attrib[k]) for k in ['v1', 'v2', 'v3']] for x in root.findall('.//m:triangle', ns)]))

def run():
    st = time.perf_counter()
    rows = read(ROOT / 'rounds/C.json')['rows']
    out = []
    (ROOT / 'exports').mkdir(exist_ok=True)
    for family in ['molar', 'premolar', 'anterior']:
        choices = [r for r in rows if r['method'] == 'exact_margin' and r['family'] == family and (r['status'] == 'SCORED') and r['gates']['closed']]
        r = min(choices, key=lambda r: r['in_situ']['p95_mm'])
        m = npz(r['mesh_path'])
        (v, f) = (m['vertices'], m['faces'])
        mm = trimesh.Trimesh(v, f, process=False)
        stem = ROOT / 'exports' / family
        mm.export(str(stem.with_suffix('.stl')))
        three_mf(stem.with_suffix('.3mf'), v, f, 'RESEARCH GEOMETRY ONLY: wall gate failed; no verified preparation or bite')
        raw = stem.with_suffix('.stl').read_bytes()
        dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
        tt = np.frombuffer(raw, dtype=dt, count=int.from_bytes(raw[80:84], 'little'), offset=84)['v']
        err = float(np.linalg.norm(tt - v[f], axis=2).max())
        (v3, f3) = load3mf(stem.with_suffix('.3mf'))
        er3 = float(np.linalg.norm(v3 - v, axis=1).max())
        assert err <= 1e-05 and er3 <= 1e-12 and np.array_equal(f3, f)
        out.append(dict(family=family, key=r['key'], mesh_path=r['mesh_path'], mesh_sha256=r['mesh_sha256'], stl_path=stem.with_suffix('.stl'), stl_sha256=sha(stem.with_suffix('.stl')), three_mf_path=stem.with_suffix('.3mf'), three_mf_sha256=sha(stem.with_suffix('.3mf')), stl_roundtrip_max_mm=err, three_mf_roundtrip_max_mm=er3, predictions={k: r[k] for k in ['in_situ', 'natural_floor', 'margin_curve_sampled_max_mm', 'wall', 'contact', 'proximal', 'gates']}, status='NOT_FABRICATION_QUALIFIED', selection='minimum in-situ p95 among frozen C exact_margin closed candidates of this type; evaluated-cohort selection'))
    freeze(ROOT / 'FROZEN_EXPORTS.json', dict(prereg_sha256=sha(ROOT / 'PREREG_EXPORTS.json'), exports=out, seconds=time.perf_counter() - st, physical_measurement='NOT_PERFORMED'))
    freeze(ROOT / 'FROZEN_PREDICTIONS.json', dict(claim_type='capability', external_referent=REFERENT, export_manifest_sha256=sha(ROOT / 'FROZEN_EXPORTS.json'), predictions=out, physical_measurement='NOT_PERFORMED; these files precede any future measurement', licence='Dataset licence binding UNKNOWN; private local research only, no redistribution claim'))
if __name__ == '__main__':
    run()
