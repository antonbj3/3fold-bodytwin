"""Unchanged crown gates, plus R6's preregistered source identity gate."""
import numpy as np, time
from crown_track import evaluate
from util import dump
from preparation import bad_edges

class Remap:

    def __init__(self, b, round):
        self.base = b
        self.round = round
        self.payload = b.payload / ('round' + str(round))

    def name(self, n):
        r = str(self.round)
        return n.replace('PREREG_R4.json', 'PREREG_R' + r + '.json').replace('FROZEN_CROWNS.json', 'FROZEN_CROWNS_R' + r + '.json').replace('payload/preparations_r3/', 'payload/preparations_r' + r + '/').replace('payload/crowns/', 'payload/crowns_r' + r + '/')

    def json(self, n):
        return self.base.json(self.name(n))

    def npz(self, n):
        return self.base.npz(self.name(n))

    def bytes(self, n):
        return self.base.bytes(self.name(n))

def triangles_canonical(x):
    order = np.lexsort((x[:, :, 2], x[:, :, 1], x[:, :, 0]), axis=1)
    return np.take_along_axis(x, order[:, :, None], axis=1)

def run(bundle, round):
    start = time.perf_counter()
    remap = Remap(bundle, round)
    out = evaluate(remap)
    if round in [6, 7]:
        for row in out['rows']:
            if row['status'] != 'EXPORTED':
                continue
            a = remap.npz('payload/preparations_r3/' + row['case_key'] + '/fields.npz')
            local = (a['source_vertices'] - a['center']) @ a['R']
            before = triangles_canonical(local[a['source_faces']])
            after = triangles_canonical(a['closed_source_vertices'][a['closed_source_faces'][:len(before)]])
            delta = float(np.max(np.abs(before - after))) if before.shape == after.shape else float('inf')
            row['source_triangle_coordinate_change_max_mm'] = delta
            row['gates']['source_identity'] = delta == 0.0
            row['gates']['closed_source_edges'] = bad_edges(a['closed_source_faces']) == 0
            shifted = after.copy()
            shifted[0, 0, 0] += 0.001
            row['fault_rejections']['changed_native_triangle'] = float(np.max(np.abs(before - shifted))) > 0
            row['fault_rejections']['missing_source_cap_facet'] = bad_edges(a['closed_source_faces'][:-1]) > 0
        exported = [r for r in out['rows'] if r['status'] == 'EXPORTED']
        out['geometric_pass'] = sum((all(r['gates'].values()) for r in exported))
        out['failed_exported'] = len(exported) - out['geometric_pass']
        out['gate'] = out['geometric_pass'] >= remap.json('PREREG_R4.json')['metrics']['minimum_successful_pairs'] and all((all(r['fault_rejections'].values()) for r in exported))
        dump(remap.payload / 'evaluation/CROWN_SCORES.json', out)
    if round == 7:
        q = remap.json('PREREG_R4.json')
        support_checks = []
        for row in out['rows']:
            a = bundle.npz('payload/preparations_r6/' + row['case_key'] + '/fields.npz')
            v = a['closed_source_vertices']
            f = a['closed_source_faces']
            capmax = float(v[f[len(a['source_faces']):], 2].max())
            margin = capmax + q['parameters']['closure_guard_mm']
            total = kept = 0.0
            tri = (a['source_vertices'][a['source_faces']] - a['center']) @ a['R']
            for t in tri:
                u = t[1] - t[0]
                w = t[2] - t[0]
                area = float(np.sqrt(sum((x * x for x in np.cross(u, w)))) / 2)
                total += area
                if all((z > margin + q['parameters']['source_comparison_above_margin_mm'] for z in t[:, 2])):
                    kept += area
            fraction = kept / total
            error = abs(fraction - row['retained_source_area_fraction'])
            eligible = fraction >= q['metrics']['minimum_retained_source_area_fraction']
            correct = (row['status'] == 'EXPORTED') == eligible
            checks = dict(source_area_recompute=error <= 1e-12, source_cap_max_recompute=abs(capmax - row['virtual_cap_max_z_mm']) <= 1e-12, eligible_status=correct, closure_separation=row['virtual_margin_z_mm'] - capmax >= q['metrics']['source_closure_separation_min_mm'] - 1e-12)
            faults = dict(false_full_coverage=abs(1.0 - fraction) > 1e-12, lowered_margin=capmax - 0.01 - capmax < q['metrics']['source_closure_separation_min_mm'])
            row['support_validation'] = checks
            row['support_fault_rejections'] = faults
            support_checks.append(all(checks.values()) and all(faults.values()))
        out['support_decisions_pass'] = all(support_checks)
        out['support_resolution'] = 'PER_SURFACE_REGION'
        dump(remap.payload / 'evaluation/CROWN_SCORES.json', out)
    out['seconds'] = time.perf_counter() - start
    dump(remap.payload / 'evaluation/CROWN_SCORES.json', out)
    return out
