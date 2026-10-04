from common import *
from geometry import load_case, wall_and_removal, export_stl, gap
from optimize import solve
import time

def run():
    check_lock()
    start = time.perf_counter()
    pr = read(ROOT / 'PREREG_R1.json')
    (q, _) = parents()
    rows = []
    pred = []
    for r in pr['cohort']:
        t = load_case(r)
        row = dict(key=r['key'], outputs={})
        prediction = dict(key=r['key'], outputs={})
        for (name, regional) in [('regional', True), ('ordinary_mesh', False)]:
            (d, opt) = solve(t, r, regional)
            (v, gates) = wall_and_removal(t, d, r)
            p = DATA / 'R1' / r['key'] / (name + '.npz')
            p.parent.mkdir(exist_ok=True, parents=True)
            np.savez_compressed(p, vertices=v, faces=t['faces'], face_roles=t['face_roles'], relief_mm=d, predicted_gap=gap(t, d))
            e = ROOT / 'exports' / r['key'] / (name + '.stl')
            export_stl(e, v, t['faces'])
            row['outputs'][name] = dict(path=str(p), sha256=sha(p), STL_path=str(e), STL_sha256=sha(e), optimization=opt, wall_and_removal=gates)
            contact = q.compare(t['xy'], t['grid_faces'], gap(t, d), t['reference_gap'])
            prediction['outputs'][name] = dict(artifact_sha256=sha(p), area_integral_abs_displacement_mm3=gates['area_integral_abs_displacement_mm3'], wall_lower_mm=gates['wall_lower_mm'], nominal_contact_symdiff_mm2=contact.get('symdiff_mm2'), negative_gap_area_mm2=contact.get('negative_gap_area_mm2'), force_N=None, chair_adjustment_min=None)
            print('constructed', r['key'], name, 'J', opt['phase2_area_integral_mm3'], 'diff', contact.get('symdiff_mm2'), 'subset', gates['removal_status'], flush=True)
        rows.append(row)
        pred.append(prediction)
        dump(ROOT / 'raw/CONSTRUCTION_PROGRESS.json', dict(rows=rows, phase='NOT_YET_VALIDATED'))
    dump(ROOT / 'raw/CONSTRUCTIONS_R1.json', dict(rows=rows, seconds=time.perf_counter() - start))
    if not (ROOT / 'FROZEN_PREDICTIONS.json').exists():
        freeze(ROOT / 'FROZEN_PREDICTIONS.json', dict(claim_type='capability', prereg_sha256=sha(ROOT / 'PREREG_R1.json'), input_lock_sha256=sha(ROOT / 'INPUT_LOCK_R1.json'), predictions=pred, physical_measurement='NOT_RUN', source_information='Source gaps were already consumed by optimizer; retrospective construction test, not held-out prediction', next_validation='Independent actual-mesh ray recomputation and same-information unchanged v4 area check'))
    else:
        old = read(ROOT / 'FROZEN_PREDICTIONS.json')['predictions']
        for (a, b) in zip(old, pred):
            for n in a['outputs']:
                if a['outputs'][n]['artifact_sha256'] != b['outputs'][n]['artifact_sha256']:
                    raise ValueError('Output drift after freeze: ' + a['key'] + ' ' + n)
    state('R1_CONSTRUCTIONS_FROZEN', 'Ten regional and10 ordinary-mesh LP outputs hash-frozen before independent mesh validation', 'Optimize full uniform parameter interval; validate all actual outputs, faults and continuous common-height boxes')
if __name__ == '__main__':
    run()
