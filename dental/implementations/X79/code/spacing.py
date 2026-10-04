from common import *
from data import reports
import re, zipfile, numpy as np, collections
from scipy.spatial import cKDTree, distance

def observations(rr):
    out = []
    excluded = []
    for (c, reps) in rr.items():
        if not reps:
            continue
        r = reps[0]
        for sentence in re.split('(?<=[.!?])\\s+', r['source_text']):
            if not re.search('\\bspacing\\b|\\bspaces\\b|\\bdiastema\\w*', sentence, re.I):
                continue
            for clause in re.split(';|\\bwhereas\\b|\\bwhile\\b|\\band (?=(?:mild|moderate|severe) crowding)', sentence, flags=re.I):
                if not re.search('\\bspacing\\b|\\bspaces\\b|\\bdiastema\\w*', clause, re.I):
                    continue
                if re.search('extraction space|space reserve|space closure|eruption space|space for.*erupt|space for eruption', clause, re.I):
                    excluded.append(dict(case_id=c, clause=clause, reason='WRONG_SPACE_QUANTITY'))
                    continue
                if re.search('\\bmay\\b|\\bmight\\b|\\bpossible\\b|\\bappear\\b|\\bseem\\b', clause, re.I):
                    excluded.append(dict(case_id=c, clause=clause, reason='HEDGED_SPACING'))
                    continue
                neg = bool(re.search('no (?:significant )?(?:spacing|spaces|diastema)|(?:spacing|spaces|diastema\\w*) (?:is |are )?absent|absence of (?:spacing|spaces|diastema)', clause, re.I))
                jaws = []
                if re.search('upper|maxillary', clause, re.I):
                    jaws.append('upper')
                if re.search('lower|mandibular', clause, re.I):
                    jaws.append('lower')
                if re.search('both (?:dental )?arch', clause, re.I):
                    jaws = ['upper', 'lower']
                teeth = sorted(set((int(f) for f in re.findall('\\b([1-4][1-8])\\b', clause))))
                if not jaws and teeth:
                    jaws = sorted(set(('upper' if f // 10 <= 2 else 'lower' for f in teeth)))
                if not jaws:
                    excluded.append(dict(case_id=c, clause=clause, reason='ARCH_SCOPE_UNKNOWN'))
                    continue
                for jaw in jaws:
                    out.append(dict(case_id=c, jaw=jaw, polarity=int(not neg), named_teeth=teeth, source_member=r['member'], source_sha256=r['sha256'], clause=clause, resolution='PER_ARCH' if len(teeth) != 2 else 'PER_TOOTH', timescale='SIMULTANEOUS'))
    return (out, excluded)

def decode_points(blob):
    n = int.from_bytes(blob[80:84], 'little') if len(blob) >= 84 else -1
    if n < 0 or len(blob) != 84 + 50 * n:
        raise ValueError('Requires original binary STL for exact X11 vertex-order mapping')
    dt = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    return np.unique(np.frombuffer(blob, dtype=dt, count=n, offset=84)['v'].reshape(-1, 3), axis=0)

def main():
    st = time.perf_counter()
    cpu = time.process_time()
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    rr = reports(m['test'], 'ios')
    (obs, ex) = observations(rr)
    write('raw/R4_SPACING_SOURCE_OBSERVATIONS.json', obs)
    write('raw/R4_SPACING_SOURCE_EXCLUSIONS.json', ex)
    positive = [r for r in obs if r['polarity']]
    keys = sorted({(r['case_id'], r['jaw']) for r in positive})
    results = {}
    checks = []
    errors = []
    with zipfile.ZipFile(ZIP) as z:
        for (i, (c, jaw)) in enumerate(keys):
            p = X11 / c / 'labels+landmarks.json'
            meta = read(p)
            a = meta['arches'][jaw]
            try:
                member = c + '/ios/ios_' + jaw + '.stl'
                blob = z.read(member)
                assert digest(blob) == a['input_sha256']
                pts = decode_points(blob).astype(np.float64)
                lp = Path(a['labels']['path'])
                assert sha(lp) == a['labels']['sha256']
                labels = np.load(lp)['labels']
                assert len(pts) == len(labels)
                height = pts @ np.array(a['frame']['occlusal_unit'])
                points = {}
                for (f, n) in a['label_counts'].items():
                    f = int(f)
                    if f == 0 or n < 100 or f % 10 == 8:
                        continue
                    ix = np.flatnonzero(labels == f)
                    q = np.median(height[ix])
                    points[f] = pts[ix[height[ix] >= q]]
                pairs = []
                qs = [1, 2] if jaw == 'upper' else [3, 4]
                adj = [(q * 10 + k, q * 10 + k + 1) for q in qs for k in range(1, 7)] + [(qs[0] * 10 + 1, qs[1] * 10 + 1)]
                for (f1, f2) in adj:
                    if f1 not in points or f2 not in points:
                        continue
                    (p1, p2) = (points[f1], points[f2])
                    (ds, ix) = cKDTree(p2).query(p1, k=1, workers=1)
                    k = int(np.argmin(ds))
                    gap = float(ds[k])
                    if i < 5:
                        s1 = p1[np.linspace(0, len(p1) - 1, min(64, len(p1)), dtype=int)]
                        s2 = p2[np.linspace(0, len(p2) - 1, min(64, len(p2)), dtype=int)]
                        kd = float(np.min(cKDTree(s2).query(s1, workers=1)[0]))
                        br = float(np.min(distance.cdist(s1, s2)))
                        checks.append(dict(case_id=c, jaw=jaw, pair=[f1, f2], identity_error_mm=abs(kd - br), wrong_gap_injection_rejected=abs(kd - (br + 1.0)) > 1e-09))
                    pairs.append(dict(FDI_pair=[f1, f2], pointset_minimum_gap_mm=gap, point_witness_A_mm=p1[k].tolist(), point_witness_B_mm=p2[int(ix[k])].tolist(), pointset_sizes=[len(p1), len(p2)], surface_enclosure='MISSING; not a continuous/anatomical surface minimum'))
                results[c + ':' + jaw] = dict(case_id=c, jaw=jaw, pairs=pairs, maximum_adjacent_pointset_minimum_mm=max((p['pointset_minimum_gap_mm'] for p in pairs), default=None), source_member=member, source_sha256=digest(blob), labels_sha256=a['labels']['sha256'], resolution='PER_TOOTH point-set minima; PER_ARCH maximum', timescale='SIMULTANEOUS')
            except (AssertionError, KeyError, ValueError) as e:
                errors.append(dict(case_id=c, jaw=jaw, error=str(e)))
            if i % 10 == 0:
                print('spacing arch', i + 1, '/', len(keys), flush=True)
    write('raw/R4_SPACING_GEOMETRY.json', results)
    joined = []
    for r in positive:
        g = results.get(r['case_id'] + ':' + r['jaw'])
        if not g or g['maximum_adjacent_pointset_minimum_mm'] is None:
            continue
        pairs = g['pairs']
        exact = [p for p in pairs if set(p['FDI_pair']) == set(r['named_teeth'])] if len(r['named_teeth']) == 2 else []
        q = exact[0]['pointset_minimum_gap_mm'] if exact else g['maximum_adjacent_pointset_minimum_mm']
        joined.append(dict(**r, candidate_gap_mm=q, geometry_scope='NAMED_PAIR' if exact else 'ARCH_ANY_ADJACENT', candidate_spacing=q > 0.5, geometry_claim='Pointset scenario only; scoped report is not a numeric gap measurement'))
    unique = {(r['case_id'], r['jaw']): r for r in joined}
    v = list(unique.values())
    maxerr = max((r['identity_error_mm'] for r in checks), default=None)
    valid = bool(checks) and maxerr <= 1e-09 and all((r['wrong_gap_injection_rejected'] for r in checks))
    result = dict(claim_type='capability', outcome='UNKNOWN_BINARY_ACCURACY_NO_NEGATIVE_FACIT', source_positive_arches=len(keys), compared_arches=len(v), source_negative_observations=sum((r['polarity'] == 0 for r in obs)), positive_only_agreement=sum((r['candidate_spacing'] for r in v)) / len(v) if v else None, sensitivity='Not population sensitivity: selectively reported positive spacing only', specificity='UNKNOWN_NO_EXPLICIT_NO_SPACING_SOURCE', geometry_errors=errors, pointset_comparator=dict(n=len(checks), max_residual_mm=maxerr, pass_gate=valid), physical_uncertainty_enclosure='MISSING', external_referent=read('PREREG_R4_SPACING_POINTSETS.json')['external_referent'], cost=cost(st, cpu))
    write('raw/R4_SPACING_JOINED.json', joined)
    write('raw/R4_POINTSET_CONTROLS.json', checks)
    write('raw/R4_RESULTS.json', result)
    (ROOT / 'HANDOFF_R4.md').write_text('R4: positive-only spacing point-set agreement ' + str(result['positive_only_agreement']) + '; specificity UNKNOWN because no explicit no-spacing reference exists. Point-set comparator ' + str(valid) + '. Next physical operation: blinded per-tooth restoration/presence and pair-spacing labels with a complete negative assessment and repeated named raters; no threshold retuning.\n')
    state('R4_DECIDED', result['outcome'], 'Package rerunnable demo, figure, matched annotation specification and missing graph coverage')
    print(json.dumps(result), flush=True)
if __name__ == '__main__':
    main()
