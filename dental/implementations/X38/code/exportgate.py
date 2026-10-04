"""CLI for geometry-bound laboratory export. PASS concerns transport only."""
import argparse, json, math, os, shutil, sys
from pathlib import Path
from mesh_transport import *

def read_contract(sidecar, expected):
    if sha(sidecar) != expected:
        raise Refused('sidecar hash differs from frozen expected SHA256')
    c = json.loads(Path(sidecar).read_text())
    payload = c['payload']
    if digest(payload) != c['contract_sha256']:
        raise Refused('contract digest differs')
    if payload['unit'] != 'mm':
        raise Refused('contract unit must be mm')
    if payload['geometry_tolerance_mm'] != TOL:
        raise Refused('unsupported tolerance')
    pred = Path(sidecar).parent / payload['frozen_predictions']['path']
    if sha(pred) != payload['frozen_predictions']['sha256']:
        raise Refused('frozen predictions changed')
    if payload.get('design_gate'):
        g = Path(sidecar).parent / payload['design_gate']['path']
        if sha(g) != payload['design_gate']['sha256']:
            raise Refused('design-gate report changed')
    return c

def check(asset, sidecar, expected, stl_unit=None, converted=False):
    c = read_contract(sidecar, expected)
    p = c['payload']
    matching_hash = sha(asset) == c['asset_sha256']
    if not matching_hash and (not converted):
        raise Refused('asset hash differs; use --converted for explicit geometry revalidation')
    (tri, meta) = load_mesh(asset, stl_unit)
    diagnostics = validate_triangles(tri)
    (match, err) = correspondence(p['triangles_mm'], tri)
    labels = [p['regions'][int(i)] for i in match]
    if meta['format'] == '3MF':
        if meta['contract'] is not None and meta['contract'] != c['contract_sha256']:
            raise Refused('native3MF contract changed')
        if meta['labels'] is not None and meta['labels'] != labels:
            raise Refused('native3MF regions disagree with geometric sidecar')
        if meta['material_names'] and meta['material_names'] != [p['material']['name']]:
            raise Refused('native3MF material changed')
        if not converted and (meta['contract'] is None or meta['labels'] is None):
            raise Refused('native3MF metadata missing')
    return dict(status='PASS', claim_type='capability', scope='tessellation-preserving transport only', asset_hash_exact=matching_hash, converted=converted, geometry_max_correspondence_mm=err, geometry_tolerance_mm=TOL, unit_mm_verified=True, region_source='native3MF+sidecar' if meta.get('labels') is not None else 'geometry-bound sidecar', region_counts={s: labels.count(s) for s in sorted(set(labels))}, regions_in_imported_face_order=labels, anatomical_region_validity=p['region_semantics'], material=p['material'], sinter=p['sinter'], frozen_predictions=p['frozen_predictions'], design_gate=p.get('design_gate'), manufacturing_release='UNKNOWN; transport does not assess design/physical calibration', diagnostics=diagnostics)

def export(source, out, stl_unit, regions, material, factor, geometry_state, predictions, design_gate=None):
    if not material.strip() or not math.isfinite(factor) or factor <= 0:
        raise Refused('explicit material and positive finite sinter factor required')
    (tri, meta) = load_mesh(source, stl_unit)
    validate_triangles(tri)
    if regions:
        r = json.loads(Path(regions).read_text())
        labels = r['labels']
        semantics = r['semantics']
    elif meta.get('labels') is not None:
        labels = meta['labels']
        semantics = 'caller must independently validate anatomical meaning'
    else:
        labels = ['unclassified_surface'] * len(tri)
        semantics = 'UNKNOWN; no anatomical labels supplied'
    if len(labels) != len(tri) or any((not isinstance(x, str) or not x.strip() or len(x) > 100 for x in labels)):
        raise Refused('regions must name every input face exactly once')
    json.loads(Path(predictions).read_text())
    report = None
    if design_gate:
        report = json.loads(Path(design_gate).read_text())
        entries = report.get('inputs', {})
        if not isinstance(entries, dict) or not any((isinstance(v, dict) and v.get('sha256') == sha(source) for v in entries.values())):
            raise Refused('design-gate report belongs to a different source mesh or lacks input SHA256')
    out = Path(out)
    if out.exists():
        raise Refused('output exists; choose new directory to preserve prior frozen export')
    out.mkdir(parents=True)
    shutil.copyfile(predictions, out / 'FROZEN_PREDICTIONS.json')
    payload = dict(schema='x38-physical-triangle-contract-v1', unit='mm', geometry_tolerance_mm=TOL, triangles_mm=tri.tolist(), regions=labels, region_semantics=semantics, source_asset_sha256=sha(source), material=dict(name=material, physical_properties='UNKNOWN; caller label only'), sinter=dict(factor=factor, convention='manufacturing_length / final_sintered_length', geometry_state=geometry_state, operation='metadata_only; no implicit scaling', calibration='UNKNOWN; verify same batch'), frozen_predictions=dict(path='FROZEN_PREDICTIONS.json', sha256=sha(predictions)), design_gate=None)
    if design_gate:
        shutil.copyfile(design_gate, out / 'DESIGN_GATE.json')
        payload['design_gate'] = dict(path='DESIGN_GATE.json', sha256=sha(design_gate), source_asset_sha256=sha(source), verdict=report.get('verdict', 'UNKNOWN'), scope='source hash matched; upstream report retained; does not grant transport or manufacture approval')
    cd = digest(payload)
    reports = []
    manifest = []
    for (name, fmt) in [('model.stl', 'binary'), ('model_ascii.stl', 'ascii'), ('model.3mf', '3mf')]:
        path = out / name
        if fmt == '3mf':
            write_3mf(path, tri, labels, material, cd)
        else:
            write_stl(path, tri, fmt == 'ascii')
        side = out / (name + '.json')
        dump(side, dict(payload=payload, contract_sha256=cd, asset_sha256=sha(path)))
        h = sha(side)
        side.with_suffix(side.suffix + '.sha256').write_text(h + '\n')
        reports.append({k: v for (k, v) in check(path, side, h, 'mm' if fmt != '3mf' else None).items() if k != 'regions_in_imported_face_order'})
        manifest.append(dict(asset=name, asset_sha256=sha(path), sidecar=side.name, sidecar_sha256=h))
    dump(out / 'EXPORT_RECEIPT.json', dict(status='PASS', scope='transport only', files=manifest, roundtrip=reports))
    return dict(status='PASS', output=str(out), files=manifest, scope='transport only; manufacturer/CAM validation remains UNKNOWN')

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='command', required=True)
    ex = sub.add_parser('export')
    ex.add_argument('source')
    ex.add_argument('--out', required=True)
    ex.add_argument('--stl-unit', choices=sorted(UNITS))
    ex.add_argument('--regions')
    ex.add_argument('--material', required=True)
    ex.add_argument('--sinter-factor', required=True, type=float)
    ex.add_argument('--geometry-state', required=True, choices=['final_sintered', 'compensated_green'])
    ex.add_argument('--predictions', required=True)
    ex.add_argument('--design-gate')
    ck = sub.add_parser('check')
    ck.add_argument('asset')
    ck.add_argument('--sidecar', required=True)
    ck.add_argument('--expected-sidecar-sha256', required=True)
    ck.add_argument('--stl-unit', choices=sorted(UNITS))
    ck.add_argument('--converted', action='store_true')
    ck.add_argument('--output')
    a = ap.parse_args()
    try:
        if a.command == 'export':
            x = export(a.source, a.out, a.stl_unit, a.regions, a.material, a.sinter_factor, a.geometry_state, a.predictions, a.design_gate)
        else:
            x = check(a.asset, a.sidecar, a.expected_sidecar_sha256, a.stl_unit, a.converted)
            x.pop('regions_in_imported_face_order')
        if getattr(a, 'output', None):
            dump(a.output, x)
        print(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False))
        return 0
    except (ValueError, KeyError, OSError, zipfile.BadZipFile, ET.ParseError) as e:
        x = dict(status='REFUSED', reason=str(e), manufacturing_release='UNKNOWN')
        if getattr(a, 'output', None):
            dump(a.output, x)
        print(json.dumps(x, ensure_ascii=False))
        return 2
if __name__ == '__main__':
    sys.exit(main())
