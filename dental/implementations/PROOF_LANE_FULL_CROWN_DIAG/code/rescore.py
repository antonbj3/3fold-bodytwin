from common import *
from local_pairs import fit, sample, metrics
import resource
resource.setrlimit(resource.RLIMIT_AS, (3584 * 1024 ** 2, 3584 * 1024 ** 2))
sys.path.insert(0, str(R3 / 'code'))
from memory_height import height
sys.path.insert(0, str(BASE / 'PROOF_LANE_GENCAD_V6/code'))
import importlib.util

def loadmodule(name, p):
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m
saved_common = sys.modules['common']
try:
    sys.modules['common'] = loadmodule('v6common_diagnostic', BASE / 'PROOF_LANE_GENCAD_V6/code/common.py')
    from reuse import grid
    from contact import compare, gates, enclosure
finally:
    sys.modules['common'] = saved_common

def run():
    pr = read(ROOT / 'PREREG_RESCORE.json')
    start = time.perf_counter()
    rows = []
    for (p, h) in pr['input_hashes'].items():
        assert sha(p) == h, p
    thresholds = pr['thresholds_mm']
    spec = read(R3 / 'PREREG_A.json')['metrics']['v6']
    for rec in pr['records']:
        t0 = time.perf_counter()
        key = rec['key']
        old = rec['source_score']
        r = {k: v for (k, v) in rec.items() if k != 'source_score'}
        r.update(resolution='PER_TOOTH', clinical_status='NOT_ESTABLISHED', actual_finish_line='UNKNOWN', loaded_occlusion='UNKNOWN', signed_proximal_contact='UNKNOWN')
        if rec['status'] == 'SCORED':
            p = npz(PUBLIC / (key + '.npz'))
            a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
            m = npz(rec['mesh_path'])
            src = a['source_triangles']
            src = src[src.mean(1)[:, 2] >= float(p['margin_z'])]
            ext = m['vertices'][m['faces'][m['face_roles'] == 0]]
            fixed = metrics(ext, src, 2048)
            (R, u, starts) = fit(sample(ext, 1536), sample(src, 1536))
            align = ext @ R.T + u
            aligned = {str(n): metrics(align, src, n) for n in [2048, 8192]}
            (xy, ff) = grid(src)
            ceiling = height(a['antagonist_triangles'], xy, True)
            rg = ceiling - height(src, xy)
            pg = ceiling - height(ext, xy)
            contact = compare(xy, ff, pg, rg)
            cg = gates(contact, spec)
            bound = enclosure(xy, ff, pg, rg, 0.05)
            assert cg == old['contact_gates'], (key, cg, old['contact_gates'])
            wall = old['wall']
            prox = old['proximal_patches']
            functional = {**old['gates'], 'spatial_contact': all(cg.values())}
            nominal = all((v for (k, v) in functional.items() if k != 'wall_cover_bound'))
            cover = all(functional.values())
            threshold = thresholds[rec['family']]
            shape_pass = aligned['8192']['p95_mm'] <= threshold if threshold is not None else None
            r.update(status='RESCORED', in_situ=fixed, old_p95_mm=old['reconstruction_p95_mm'], oracle_aligned=aligned, shape_threshold_mm=threshold, shape_ratio=aligned['8192']['p95_mm'] / threshold if threshold else None, shape_envelope_pass=shape_pass, aligned_R=R, aligned_t=u, alignment_start_rms_mm=starts, contact=contact, contact_gates=cg, pattern_enclosure=bound, wall=wall, proximal_patches=prox, technical_planarity_mm=old['margin_planarity_mm'], nominal_function_pass=nominal, cover_function_pass=cover, positive_native_band=contact['reference']['area_mm2'] > 0, shape_and_nominal_function=bool(shape_pass and nominal), shape_and_cover_function=bool(shape_pass and cover), metric_provenance=dict(shape='freshly computed from frozen meshes', contact='freshly computed; exact gate parity asserted', wall_proximal_planarity='hash-bound inherited evidence; no replay claim'), continuous_numeric_enclosure='MISSING; probe density comparison only')
        r['seconds'] = time.perf_counter() - t0
        rows.append(r)
        dump(ROOT / 'raw/RESCORE.json', dict(rows=rows, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
        print(rec['tag'], key, r['status'], r.get('shape_ratio'), r.get('nominal_function_pass'), flush=True)
        state('RESCORING', str(len(rows)) + '/' + str(len(pr['records'])), 'Finish all frozen designs and scalar diagnostic')
if __name__ == '__main__':
    run()
