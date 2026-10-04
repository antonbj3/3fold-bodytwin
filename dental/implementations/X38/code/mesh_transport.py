"""Bounded, explicit STL / single-mesh 3MF transport. No mesh repair."""
from pathlib import Path
import hashlib, io, json, re, zipfile
import xml.etree.ElementTree as ET
import numpy as np
import trimesh
from scipy.spatial import cKDTree
CORE = 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
SETS = 'http://schemas.microsoft.com/3dmanufacturing/trianglesets/2021/07'
X38 = 'https://3fold.local/x38/transport/v1'
UNITS = {'mm': 1.0, 'millimeter': 1.0, 'micron': 0.001, 'cm': 10.0, 'centimeter': 10.0, 'm': 1000.0, 'meter': 1000.0, 'inch': 25.4, 'foot': 304.8}
MAX_BYTES = 150000000
MAX_FACES = 300000
TOL = 1e-05
ET.register_namespace('', CORE)
ET.register_namespace('t', SETS)

class Refused(ValueError):
    pass

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def digest(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def dump(p, x):
    Path(p).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def physical_mesh(tri):
    t = np.asarray(tri, float)
    if t.ndim != 3 or t.shape[1:] != (3, 3) or (not 0 < len(t) <= MAX_FACES) or (not np.isfinite(t).all()):
        raise Refused('empty/nonfinite/oversized triangle array')
    (v, inv) = np.unique(t.reshape(-1, 3), axis=0, return_inverse=True)
    return trimesh.Trimesh(v, inv.reshape(-1, 3), process=False)

def validate_triangles(t):
    m = physical_mesh(t)
    if np.any(np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) <= 1e-12):
        raise Refused('degenerate facet')
    if len(np.unique(np.sort(m.faces, axis=1), axis=0)) != len(t):
        raise Refused('duplicate facet')
    if not m.is_watertight:
        raise Refused('mesh is not watertight')
    if not m.is_winding_consistent:
        raise Refused('inconsistent winding')
    if m.volume <= 0:
        raise Refused('nonpositive oriented volume')
    return dict(faces=len(t), vertices=len(m.vertices), watertight=True, winding_consistent=True, volume_mm3=float(m.volume), extent_mm=np.ptp(m.vertices, axis=0).tolist(), geometry_resolution='PER_POINT', volume_resolution='PER_TOOTH', self_intersections='UNKNOWN', repair='none; exact vertex dedup only')

def _xml(blob):
    if b'<!DOCTYPE' in blob.upper() or b'<!ENTITY' in blob.upper():
        raise Refused('XML entity declarations unsupported')
    return ET.fromstring(blob)

def load_mesh(path, stl_unit=None):
    p = Path(path)
    if p.stat().st_size > MAX_BYTES:
        raise Refused('input exceeds150MB')
    if p.suffix.lower() == '.stl':
        if stl_unit not in UNITS:
            raise Refused('STL unit must be asserted explicitly')
        blob = p.read_bytes()
        n = int.from_bytes(blob[80:84], 'little') if len(blob) >= 84 else 0
        if len(blob) == 84 + 50 * n:
            if not 0 < n <= MAX_FACES:
                raise Refused('invalid STL face count')
            dt = np.dtype([('normal', '<f4', 3), ('v', '<f4', (3, 3)), ('attr', '<u2')])
            tri = np.frombuffer(blob, dtype=dt, count=n, offset=84)['v'].astype(float)
        else:
            try:
                s = blob.decode('ascii')
            except UnicodeDecodeError as e:
                raise Refused('truncated binary/non-ASCII STL') from e
            lines = [x.strip().split() for x in s.splitlines() if x.strip()]
            if not lines or lines[0][0] != 'solid' or lines[-1][0] != 'endsolid':
                raise Refused('invalid ASCII STL envelope')
            blocks = lines[1:-1]
            if not blocks or len(blocks) % 7:
                raise Refused('malformed ASCII facets')
            faces = []
            for i in range(0, len(blocks), 7):
                b = blocks[i:i + 7]
                if b[0][:2] != ['facet', 'normal'] or len(b[0]) != 5 or b[1] != ['outer', 'loop'] or (b[5] != ['endloop']) or (b[6] != ['endfacet']) or any((len(x) != 4 or x[0] != 'vertex' for x in b[2:5])):
                    raise Refused('malformed ASCII facet')
                faces.append([[float(x) for x in line[1:]] for line in b[2:5]])
                if len(faces) > MAX_FACES:
                    raise Refused('too many facets')
            tri = np.asarray(faces, float)
        return (tri * UNITS[stl_unit], dict(format='STL', native_unit=None, asserted_unit=stl_unit, labels=None, contract=None))
    if p.suffix.lower() != '.3mf':
        raise Refused('supported formats: STL and3MF')
    with zipfile.ZipFile(p) as z:
        infos = z.infolist()
        if sum((a.file_size for a in infos)) > MAX_BYTES or len(infos) > 64:
            raise Refused('3MF expanded input exceeds budget')
        if len({a.filename for a in infos}) != len(infos):
            raise Refused('duplicate ZIP member')
        rel = _xml(z.read('_rels/.rels'))
        models = [x.get('Target') for x in rel if x.get('Type') == 'http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel']
        if len(models) != 1:
            raise Refused('requires one OPC model relationship')
        name = models[0].lstrip('/')
        if '..' in Path(name).parts:
            raise Refused('unsupported model target')
        r = _xml(z.read(name))
    if r.tag != '{' + CORE + '}model':
        raise Refused('unsupported3MF namespace')
    if r.get('requiredextensions'):
        raise Refused('required3MF extensions unsupported')
    unit = r.get('unit', 'millimeter')
    if unit not in UNITS:
        raise Refused('unsupported3MF unit')
    ns = {'c': CORE, 't': SETS}
    obs = r.findall('c:resources/c:object', ns)
    items = r.findall('c:build/c:item', ns)
    if len(obs) != 1 or len(items) != 1 or items[0].get('objectid') != obs[0].get('id') or (obs[0].find('c:components', ns) is not None):
        raise Refused('only one built mesh object supported; components/multi-object UNKNOWN')
    obj = obs[0]
    mesh = obj.find('c:mesh', ns)
    if mesh is None or obj.get('type', 'model') != 'model':
        raise Refused('requires solid model mesh')
    v = np.array([[float(e.get(k)) for k in ['x', 'y', 'z']] for e in mesh.findall('c:vertices/c:vertex', ns)])
    fs = mesh.findall('c:triangles/c:triangle', ns)
    f = np.array([[int(e.get(k)) for k in ['v1', 'v2', 'v3']] for e in fs])
    if not 0 < len(f) <= MAX_FACES or v.ndim != 2 or v.shape[1] != 3 or (f.min() < 0) or (f.max() >= len(v)):
        raise Refused('invalid3MF mesh indices')
    transform = items[0].get('transform')
    if transform:
        a = np.array([float(x) for x in transform.split()])
        if a.shape != (12,) or not np.isfinite(a).all():
            raise Refused('invalid build transform')
        a = a.reshape(4, 3)
        if np.linalg.det(a[:3]) <= 0:
            raise Refused('singular/reflected transform unsupported')
        v = v @ a[:3] + a[3]
    labels = None
    groups = mesh.findall('t:trianglesets/t:triangleset', ns)
    if groups:
        identifiers = [g.get('identifier') for g in groups]
        if any((not x for x in identifiers)) or len(set(identifiers)) != len(identifiers):
            raise Refused('invalid triangle-set identifiers')
        labels = [''] * len(f)
        for g in groups:
            label = g.get('name')
            if not label:
                raise Refused('unnamed triangle set')
            ids = [int(x.get('index')) for x in g.findall('t:ref', ns)]
            for b in g.findall('t:refrange', ns):
                (lo, hi) = (int(b.get('startindex')), int(b.get('endindex')))
                if not 0 <= lo <= hi < len(f):
                    raise Refused('invalid triangle-set range')
                ids.extend(range(lo, hi + 1))
            for i in ids:
                if not 0 <= i < len(f) or labels[i]:
                    raise Refused('overlapping/invalid triangle sets')
                labels[i] = label
        if any((not x for x in labels)):
            raise Refused('incomplete triangle-set partition')
    metadata = r.findall('c:metadata', ns)
    if len({x.get('name') for x in metadata}) != len(metadata):
        raise Refused('duplicate metadata name')
    md = {x.get('name'): x.text for x in metadata}
    contract = md.get('x38:contract_sha256')
    materials = r.findall('c:resources/c:basematerials/c:base', ns)
    if materials:
        mg = r.findall('c:resources/c:basematerials', ns)
        if len(mg) != 1 or len(materials) != 1 or obj.get('pid') != mg[0].get('id') or (obj.get('pindex') != '0'):
            raise Refused('unsupported/invalid material assignment')
        for e in fs:
            if e.get('pid', obj.get('pid')) != mg[0].get('id') or any((e.get(k, '0') != '0' for k in ['p1', 'p2', 'p3'])):
                raise Refused('triangle material assignment changed')
    return (v[f] * UNITS[unit], dict(format='3MF', native_unit=unit, asserted_unit=None, labels=labels, contract=contract, material_names=[x.get('name') for x in materials], build_transform_applied=bool(transform)))

def write_stl(path, tri, ascii=False):
    t = np.asarray(tri, float)
    normal = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    if ascii:
        with Path(path).open('w') as f:
            f.write('solid X38\n')
            for (row, n) in zip(t, normal):
                f.write('  facet normal ' + ' '.join((format(x, '.17g') for x in n)) + '\n    outer loop\n')
                for v in row:
                    f.write('      vertex ' + ' '.join((format(x, '.17g') for x in v)) + '\n')
                f.write('    endloop\n  endfacet\n')
            f.write('endsolid X38\n')
    else:
        dt = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
        rows = np.zeros(len(t), dtype=dt)
        rows['n'] = normal
        rows['v'] = t
        Path(path).write_bytes(b'X38; coordinates require sidecar unit'.ljust(80, b' ') + len(t).to_bytes(4, 'little') + rows.tobytes())

def write_3mf(path, tri, labels, material, contract):
    m = physical_mesh(tri)
    model = ET.Element('{' + CORE + '}model', {'unit': 'millimeter', '{http://www.w3.org/XML/1998/namespace}lang': 'en-US', 'xmlns:x38': X38})
    ET.SubElement(model, '{' + CORE + '}metadata', {'name': 'x38:contract_sha256', 'preserve': 'true'}).text = contract
    resources = ET.SubElement(model, '{' + CORE + '}resources')
    mats = ET.SubElement(resources, '{' + CORE + '}basematerials', {'id': '1'})
    ET.SubElement(mats, '{' + CORE + '}base', {'name': material, 'displaycolor': '#CCCCCCFF'})
    obj = ET.SubElement(resources, '{' + CORE + '}object', {'id': '2', 'type': 'model', 'pid': '1', 'pindex': '0'})
    mesh = ET.SubElement(obj, '{' + CORE + '}mesh')
    verts = ET.SubElement(mesh, '{' + CORE + '}vertices')
    for v in m.vertices:
        ET.SubElement(verts, '{' + CORE + '}vertex', {k: format(x, '.17g') for (k, x) in zip(['x', 'y', 'z'], v)})
    faces = ET.SubElement(mesh, '{' + CORE + '}triangles')
    for row in m.faces:
        ET.SubElement(faces, '{' + CORE + '}triangle', {k: str(x) for (k, x) in zip(['v1', 'v2', 'v3'], row)})
    sets = ET.SubElement(mesh, '{' + SETS + '}trianglesets')
    for (j, label) in enumerate(sorted(set(labels))):
        g = ET.SubElement(sets, '{' + SETS + '}triangleset', {'name': label, 'identifier': 'x38:region' + str(j)})
        ids = np.flatnonzero(np.asarray(labels) == label)
        start = prev = int(ids[0])
        for k in list(ids[1:]) + [None]:
            if k is not None and k == prev + 1:
                prev = int(k)
                continue
            ET.SubElement(g, '{' + SETS + '}refrange', {'startindex': str(start), 'endindex': str(prev)})
            if k is not None:
                start = prev = int(k)
    ET.SubElement(ET.SubElement(model, '{' + CORE + '}build'), '{' + CORE + '}item', {'objectid': '2'})
    content = b'<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>'
    rel = b'<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>'
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for (name, data) in [('[Content_Types].xml', content), ('_rels/.rels', rel), ('3D/3dmodel.model', ET.tostring(model, encoding='utf-8', xml_declaration=True))]:
            info = zipfile.ZipInfo(name, (2026, 10, 3, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, data)

def correspondence(reference, current, tol=TOL):
    """Unique oriented per-facet bijection; no nearest-face relabel on remeshing."""
    a = np.asarray(reference, float)
    b = np.asarray(current, float)
    if a.shape != b.shape or not np.isfinite(b).all():
        raise Refused('facet count/nonfinite geometry changed')
    if np.array_equal(a, b):
        return (np.arange(len(a)), 0.0)
    tree = cKDTree(a.mean(1))
    near = tree.query_ball_point(b.mean(1), tol, workers=1)
    match = np.empty(len(b), int)
    errors = np.empty(len(b))
    for (i, ids) in enumerate(near):
        if len(ids) > 64:
            raise Refused('ambiguous local facet correspondence')
        good = []
        for j in ids:
            e = min((float(np.linalg.norm(np.roll(a[j], k, axis=0) - b[i], axis=1).max()) for k in range(3)))
            if e <= tol:
                good.append((j, e))
        if len(good) != 1:
            raise Refused('changed/reversed/ambiguous triangle at index ' + str(i))
        (match[i], errors[i]) = good[0]
    if len(np.unique(match)) != len(a):
        raise Refused('correspondence is not bijective')
    return (match, float(errors.max()))
