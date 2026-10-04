from cadlib import *
import zipfile, xml.etree.ElementTree as ET
NS = 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
ET.register_namespace('', NS)

def tag(x):
    return '{' + NS + '}' + x

def write(path, vertices, faces, roles, metadata):
    vertices = np.asarray(vertices)
    faces = np.asarray(faces)
    roles = np.asarray(roles)
    if roles.shape != (len(faces),) or not set(roles).issubset({0, 1, 2}):
        raise ValueError('Invalid region assignment')
    root = ET.Element(tag('model'), unit='millimeter', attrib={'{http://www.w3.org/XML/1998/namespace}lang': 'en-US'})
    ET.SubElement(root, tag('metadata'), name='Title').text = 'Research crown; physical qualification UNKNOWN'
    ET.SubElement(root, tag('metadata'), name='Description').text = json.dumps(clean(metadata), separators=(',', ':'))
    res = ET.SubElement(root, tag('resources'))
    bm = ET.SubElement(res, tag('basematerials'), id='1')
    for (name, color) in [('exterior', '#C8DEF2FF'), ('intaglio', '#ED9368FF'), ('cervical_rim', '#65B68FFF')]:
        ET.SubElement(bm, tag('base'), name=name, displaycolor=color)
    obj = ET.SubElement(res, tag('object'), id='2', type='model', pid='1', pindex='0')
    mesh = ET.SubElement(obj, tag('mesh'))
    vs = ET.SubElement(mesh, tag('vertices'))
    for v in vertices:
        ET.SubElement(vs, tag('vertex'), **{k: format(float(x), '.17g') for (k, x) in zip('xyz', v)})
    ts = ET.SubElement(mesh, tag('triangles'))
    for (f, r) in zip(faces, roles):
        ET.SubElement(ts, tag('triangle'), v1=str(f[0]), v2=str(f[1]), v3=str(f[2]), pid='1', p1=str(r), p2=str(r), p3=str(r))
    build = ET.SubElement(root, tag('build'))
    ET.SubElement(build, tag('item'), objectid='2')
    content = b'<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>'
    rel = b'<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>'
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for (name, data) in [('[Content_Types].xml', content), ('_rels/.rels', rel), ('3D/3dmodel.model', ET.tostring(root, encoding='utf-8', xml_declaration=True))]:
            zi = zipfile.ZipInfo(name, date_time=(2026, 10, 4, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(zi, data)

def readback(path):
    with zipfile.ZipFile(path) as z:
        rt = ET.fromstring(z.read('3D/3dmodel.model'))
    if rt.attrib.get('unit') != 'millimeter':
        raise ValueError('Explicit millimeter unit required')
    v = np.array([[float(x.attrib[k]) for k in 'xyz'] for x in rt.iter(tag('vertex'))])
    f = np.array([[int(x.attrib[k]) for k in ['v1', 'v2', 'v3']] for x in rt.iter(tag('triangle'))])
    r = []
    for x in rt.iter(tag('triangle')):
        a = x.attrib
        if a.get('pid') != '1' or len({a.get(k) for k in ['p1', 'p2', 'p3']}) != 1:
            raise ValueError('Ambiguous face role')
        r.append(int(a['p1']))
    if not set(r).issubset({0, 1, 2}):
        raise ValueError('Unknown face role')
    if f.min() < 0 or f.max() >= len(v):
        raise ValueError('Bad facet index')
    return (v, f, np.array(r))
