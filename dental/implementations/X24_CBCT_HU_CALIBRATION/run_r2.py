import os
import json, hashlib, time, gzip, struct, zipfile, datetime
from pathlib import Path
import numpy as np
from calibration import fit_affine, apply
from frozen_coefficients_match import models_match
from measure_phantom import P, SOURCE, DATA, sample, prepare_rois, sha

def parse_nii(blob):
    b = gzip.decompress(blob)
    dim = struct.unpack('<8h', b[40:56])
    dt = struct.unpack('<h', b[70:72])[0]
    vo = int(struct.unpack('<f', b[108:112])[0])
    (sl, ic) = struct.unpack('<2f', b[112:120])
    A = np.array([struct.unpack('<4f', b[i:i + 16]) for i in [280, 296, 312]] + [(0.0, 0.0, 0.0, 1.0)])
    m = {2: np.uint8, 4: np.int16, 8: np.int32, 16: np.float32, 64: np.float64, 512: np.uint16, 256: np.int8}
    (nx, ny, nz) = dim[1:4]
    img = np.frombuffer(b, dtype=m[dt], count=nx * ny * nz, offset=vo).reshape(nz, ny, nx).astype(np.float32)
    return (img * (sl if sl else 1) + ic, A)

def weighted_median(v, w):
    ii = np.argsort(v)
    cw = np.cumsum(w[ii])
    n = int(cw[-1])
    return float((v[ii][np.searchsorted(cw, (n - 1) // 2 + 1)] + v[ii][np.searchsorted(cw, n // 2 + 1)]) / 2)

def cached_sample(img, cache):
    rows = []
    raw = {}
    for m in ['AIR', 'PMP', 'LDPE', 'Polystyrene', 'Delrin', 'Teflon']:
        ids = cache[m + '_image_flat_index']
        w = cache[m + '_legacy_multiplicity']
        v = img.ravel()[ids].astype(float)
        mean = float(v @ w / w.sum())
        sd = float(np.sqrt((v - mean) ** 2 @ w / w.sum()))
        rows.append({'material': m, 'gray_median': weighted_median(v, w), 'gray_SD_voxels': sd, 'n_finite': int(w.sum()), 'n_unique_image_points': len(ids), 'resolution': 'PER_SURFACE_REGION'})
        for (k, a) in [('image_flat_index', ids), ('gray', v), ('legacy_multiplicity', w)]:
            raw[m + '_' + k] = a
    return (rows, raw)

def lagrange(x, y, g):
    g = np.asarray(g, float)
    out = np.zeros_like(g)
    for i in range(3):
        fac = np.ones_like(g)
        for j in range(3):
            if j != i:
                fac *= (g - x[j]) / (x[i] - x[j])
        out += y[i] * fac
    return out

def calc(model, g, variant):
    if variant == 'quadratic3_candidate':
        return apply(model['quadratic'], g)
    if variant == 'piecewise_linear3_diagnostic':
        x = np.array(model['x'])
        y = np.array(model['y'])
        g = np.asarray(g, float)
        return np.where(g <= x[1], y[0] + (g - x[0]) * (y[1] - y[0]) / (x[1] - x[0]), y[1] + (g - x[1]) * (y[2] - y[1]) / (x[2] - x[1]))
    return apply(model[variant], g)
OUT = Path(os.environ.get('X24_OUTPUT_DIR', str(P)))
OUT.mkdir(parents=True, exist_ok=True)

def run():
    t0 = time.perf_counter()
    pre = json.loads((P / 'PREREG_R2.json').read_text())
    ph = sha(P / 'PREREG_R2.json')
    assert ph == (P / 'PREREG_R2.sha256').read_text().strip()
    bas = json.loads((P / 'PHANTOM_BASELINE.json').read_text())
    ref = next((b for b in bas if b['device'] == 'Varian'))
    V = {r['material']: r['gray_median'] for r in ref['inserts']}
    models = {}
    for b in bas:
        G = {r['material']: r['gray_median'] for r in b['inserts']}
        x = [G[m] for m in pre['anchors']]
        y = [V[m] for m in pre['anchors']]
        c = np.polynomial.polynomial.polyfit(x, y, 2)
        models[b['device']] = {'x': x, 'y': y, 'quadratic': c.tolist(), 'air_delrin2_diagnostic': fit_affine([G['AIR'], G['Delrin']], [V['AIR'], V['Delrin']]).tolist(), 'air_polystyrene2_R1': fit_affine([G['AIR'], G['Polystyrene']], [V['AIR'], V['Polystyrene']]).tolist()}
    frozen = {'frozen_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': ph, 'reference': V, 'models': models, 'new_measurements_read_at_freeze': 0}
    fp = P / 'FROZEN_RESPONSE_R2.json'
    if not fp.exists():
        fp.write_text(json.dumps(frozen, indent=2) + '\n')
        (P / 'FROZEN_RESPONSE_R2.sha256').write_text(sha(fp) + '\n')
        fp.chmod(292)
    else:
        old = json.loads(fp.read_text())
        assert models_match(old['models'], models) and old['reference'] == V
        models = old['models']
    assert sha(fp) == (P / 'FROZEN_RESPONSE_R2.sha256').read_text().strip()
    new = []
    cachepath = P / 'PHANTOM_NEW_R2.json'
    if cachepath.exists():
        new = json.loads(cachepath.read_text())
        assert [r['zip_member'] for r in new] == pre['new_scans']
    else:
        with zipfile.ZipFile(SOURCE / 'cbct_study.zip') as z:
            for name in pre['new_scans']:
                dev = name.split('/')[0]
                b = next((r for r in bas if r['device'] == dev))
                blob = z.read(name)
                (img, A) = parse_nii(blob)
                same = img.shape == tuple(b['image_shape_zyx']) and np.array_equal(A, np.array(b['image_affine']))
                if same:
                    with np.load(b['raw_samples']['path']) as cache:
                        (rows, raw) = cached_sample(img, cache)
                else:
                    (rois, _) = prepare_rois(SOURCE / 'extracted' / dev)
                    (rows, raw) = sample(img, A, rois)
                arr = DATA / ('R2_' + dev + '_' + hashlib.sha256(name.encode()).hexdigest()[:12] + '.npz')
                np.savez_compressed(arr, **raw)
                new.append({'device': dev, 'scan': 'R2_new', 'zip_member': name, 'compressed_member_sha256': hashlib.sha256(blob).hexdigest(), 'inserts': rows, 'affine_matches_baseline': same, 'image_affine': A.tolist(), 'image_shape_zyx': list(img.shape), 'raw_samples': {'path': str(arr), 'bytes': arr.stat().st_size, 'sha256': sha(arr)}})
                (P / 'PHANTOM_NEW_R2_partial.json').write_text(json.dumps(new, indent=2) + '\n')
                print(len(new), dev, Path(name).parent.name, flush=True)
                del img, blob, raw
        cachepath.write_text(json.dumps(new, indent=2) + '\n')
    rows = []
    drift = []
    checkmax = 0.0
    mono = []
    for (dev, model) in models.items():
        xx = np.linspace(model['x'][0], model['x'][-1], 1001)
        mono.append({'device': dev, 'min_derivative': float(np.min(model['quadratic'][1] + 2 * model['quadratic'][2] * xx))})
    for b in bas + new:
        G = {r['material']: r['gray_median'] for r in b['inserts']}
        mod = models[b['device']]
        for m in pre['held_out']:
            for variant in pre['variants']:
                h = float(calc(mod, G[m], variant))
                row = {'device': b['device'], 'scan': b['scan'], 'source': b.get('zip_member', b.get('source')), 'material': m, 'variant': variant, 'H_reference': V[m], 'gray': G[m], 'H_calibrated': h, 'error_HU': h - V[m], 'practice_error_HU': G[m] - V[m], 'pass_40HU': abs(h - V[m]) <= 40, 'resolution': 'PER_SURFACE_REGION'}
                if variant == 'quadratic3_candidate':
                    checkmax = max(checkmax, abs(h - float(lagrange(mod['x'], mod['y'], G[m]))))
                    row['injected_target_plus100_rejected'] = abs(h - (V[m] + 100)) > 40
                rows.append(row)
        if b['scan'] == 'R2_new':
            for m in pre['anchors']:
                error = float(calc(mod, G[m], 'quadratic3_candidate') - V[m])
                drift.append({'device': b['device'], 'source': b['zip_member'], 'material': m, 'error_HU': error, 'pass_40HU': abs(error) <= 40})
    ele = [r for r in rows if r['device'] != 'Varian']
    cand = [r for r in ele if r['variant'] == 'quadratic3_candidate']
    variants = {}
    for var in pre['variants']:
        rr = [r for r in ele if r['variant'] == var]
        variants[var] = {'n': len(rr), 'max_abs_HU': max((abs(r['error_HU']) for r in rr)), 'median_abs_HU': float(np.median([abs(r['error_HU']) for r in rr])), 'pass_fraction': sum((r['pass_40HU'] for r in rr)) / len(rr), 'resolution': 'PER_SURFACE_REGION'}
    imp = float(np.median([abs(r['practice_error_HU']) for r in cand]) / np.median([abs(r['error_HU']) for r in cand]))
    gate = {'heldout40HU': all((r['pass_40HU'] for r in cand)), 'new_scan_anchor_drift40HU': all((r['pass_40HU'] for r in drift if r['device'] != 'Varian')), 'improvement4x': imp >= 4, 'monotonic': all((r['min_derivative'] > 0 for r in mono)), 'equal_information_check': checkmax <= 1e-08, 'error_injection': all((r['injected_target_plus100_rejected'] for r in cand))}
    out = {'round': 'R2', 'claim_type': 'information_link', 'outcome': 'PASS_SCANNER_REFERENCE_TRANSFER' if all(gate.values()) else 'NEGATIVE_THREE_ANCHOR_SUFFICIENCY', 'external_referent': pre['external_referent'], 'prereg_sha256': ph, 'frozen_response_sha256': sha(fp), 'gates': gate, 'variants': variants, 'practice_median_abs_HU': float(np.median([abs(r['practice_error_HU']) for r in cand])), 'median_improvement_factor': imp, 'equal_information_max_HU': checkmax, 'rows': rows, 'anchor_drift': drift, 'monotonicity': mono, 'cost': {'wall_s': time.perf_counter() - t0, 'new_scans': len(new), 'new_physical_acquisitions': 0, 'raw_sample_bytes': sum((r['raw_samples']['bytes'] for r in new))}, 'dropout': {'selected_new_scans': len(pre['new_scans']), 'rejected_new_scans': 0, 'new_ROIs': 6 * len(new), 'rejected_ROIs': 0, 'denominator_note': 'Varian rows kept as reference stability diagnostics, excluded from Elekta transfer efficacy'}, 'status_absolute_HA_or_bone_E': 'UNKNOWN_NOT_MEASURED'}
    (OUT / 'results_R2.json').write_text(json.dumps(out, indent=2) + '\n')
    (OUT / 'HANDOFF_R2.md').write_text('# R2 : extended scanner reference\n\n' + out['outcome'] + '\n' + json.dumps(variants, indent=2) + '\n\nAll material canned. Next operation: consumer sensitivity on existing implant profiles, explicit conditional scenarios and K42 refrain without paired mechanically reference. No new physical water/HA - measurement .\n')
    s = json.loads((P / 'CURRENT_WORK_STATE.json').read_text())
    s.update(phase='R2_COMPLETE', latest_gate=gate, next_operation='freeze consumer sensitivity R3')
    (OUT / 'CURRENT_WORK_STATE.json').write_text(json.dumps(s, indent=2) + '\n')
    print(out['outcome'], gate, variants, flush=True)
if __name__ == '__main__':
    run()
