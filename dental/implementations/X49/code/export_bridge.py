"""Reuse X38 transport; native metadata binds source, report and sinter value."""
import json
import math
import tempfile
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from parents import HERE, transport as T, exportgate as X, sha
NS = 'https://3fold.local/x49/design-gate/v1'

def check_native(asset, sidecar, expected_sha, source_sha, report_sha, factor):
    x = X.check(asset, sidecar, expected_sha)
    with zipfile.ZipFile(asset) as z:
        root = T._xml(z.read('3D/3dmodel.model'))
    meta = {e.get('name'): e.text for e in root.findall('{' + T.CORE + '}metadata')}
    required = {'x49:source_sha256': source_sha, 'x49:design_report_sha256': report_sha, 'x49:sinter_factor': format(factor, '.17g'), 'x49:sinter_convention': 'manufacturing_length / final_sintered_length'}
    if any((meta.get(k) != v for (k, v) in required.items())):
        raise ValueError('native X49 source/report/sinter metadata changed')
    x['native_x49_metadata'] = required
    return x

def export(paths, report, out, contract, factor):
    if report.get('inputs', {}).get('crown', {}).get('sha256') != sha(paths['crown']):
        raise ValueError('Design report crown SHA256 does not match export source')
    diag = report['inputs']['crown']['diagnostics']
    crown_faults = any((w.get('mesh') == 'crown' for w in report['rules']['mesh_health'].get('witnesses', [])))
    if crown_faults or not diag['watertight'] or (not diag['winding_consistent']) or (report['rules']['scale']['status'] != 'PASS'):
        return dict(status='FAIL', reason='Export requires healthy crown and explicit plausible scale', resolution='PER_TOOTH')
    if factor is None or not math.isfinite(factor) or factor <= 0:
        return dict(status='UNKNOWN', reason='Explicit positive sinter factor required; no inferred process compensation', resolution='PER_TOOTH')
    meta = json.loads(Path(contract).read_text()) if contract else {}
    unit = report.get('declared_units')
    if unit not in ('mm', 'um'):
        return dict(status='UNKNOWN', reason='No authoritative units', resolution='PER_TOOTH')
    unit = 'micron' if unit == 'um' else unit
    out = Path(out)
    if out.exists():
        raise ValueError('export output exists; preserve previous artifacts and use a new directory')
    with tempfile.TemporaryDirectory(prefix='x49-export-') as tmp:
        tmp = Path(tmp)
        rpath = tmp / 'DESIGN_GATE.json'
        T.dump(rpath, report)
        active = json.loads((HERE / 'ACTIVE_ROUND.json').read_text())
        pred = HERE / active['predictions_file']
        (tri, _) = T.load_mesh(paths['crown'], unit)
        labels = ['unclassified_surface'] * len(tri)
        names = [[] for _ in tri]
        for label in ('intaglio', 'exterior', 'occlusal', 'margin'):
            for i in meta.get('regions', {}).get(label, []):
                names[i].append(label)
        labels = ['+'.join(sorted(a)) if a else 'unclassified_surface' for a in names]
        rp = tmp / 'REGIONS.json'
        T.dump(rp, dict(labels=labels, semantics=meta.get('region_semantics', 'UNKNOWN; caller labels, no anatomical validation')))
        x = X.export(paths['crown'], out, unit, rp, report['material'], factor, 'final_sintered', pred, rpath)
    path = out / 'model.3mf'
    with zipfile.ZipFile(path) as z:
        members = {n: z.read(n) for n in z.namelist()}
    root = T._xml(members['3D/3dmodel.model'])
    root.set('xmlns:x49', NS)
    source_sha = sha(paths['crown'])
    report_sha = sha(out / 'DESIGN_GATE.json')
    values = {'x49:source_sha256': source_sha, 'x49:design_report_sha256': report_sha, 'x49:sinter_factor': format(factor, '.17g'), 'x49:sinter_convention': 'manufacturing_length / final_sintered_length'}
    for (k, v) in values.items():
        ET.SubElement(root, '{' + T.CORE + '}metadata', {'name': k, 'preserve': 'true'}).text = v
    members['3D/3dmodel.model'] = ET.tostring(root, encoding='utf-8', xml_declaration=True)
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for (n, b) in members.items():
            z.writestr(n, b)
    side = out / 'model.3mf.json'
    c = json.loads(side.read_text())
    c['asset_sha256'] = sha(path)
    T.dump(side, c)
    h = sha(side)
    (out / 'model.3mf.json.sha256').write_text(h + '\n')
    checked = check_native(path, side, h, source_sha, report_sha, factor)
    checked.pop('regions_in_imported_face_order')
    receipt = dict(status='PASS', scope='transport only', asset=str(path), asset_sha256=sha(path), sidecar_sha256=h, design_report_sha256=report_sha, source_sha256=source_sha, native_metadata=values, roundtrip=checked, upstream_verdict=report['verdict'], commercial_CAM_import='UNKNOWN', full_XSD_validation='UNKNOWN', sinter_calibration='UNKNOWN')
    T.dump(out / 'EXPORT_RECEIPT.json', receipt)
    return dict(status='PASS', reason='3MF geometry, unit, per-facet label, sinter and source/report hash round-trip preserved', resolution='PER_POINT', receipt=receipt, manufacturing_release='UNKNOWN; upstream FAIL/UNKNOWN retained')
