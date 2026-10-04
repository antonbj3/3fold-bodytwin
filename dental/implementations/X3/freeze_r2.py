from crownbench import *
import copy

def freeze():
    if (P / 'PREREG_R2.json').exists():
        return
    parent = json.load(open(P / 'RESULTS_R1.json'))
    pr = copy.deepcopy(json.load(open(P / 'PREREG_R1.json')))
    pr.update(round='R2_VISIBLE_GINGIVAL_RING', frozen_utc=now(), capability=pr['capability'], parent={'round': 'R1', 'decision': parent['decision'], 'results_sha256': sha(P / 'RESULTS_R1.json'), 'raw_sha256': sha(P / 'RAW_METRICS_R1.jsonl')}, obstacle='Bilateral and neighbor-gap placement has no patient-local cervical anchor; lower average gap can conceal cusp/pose mistakes. Change the positional observations, retain independent held-out original anatomy.', changed_operation='Extract a local boundary component from the surviving masked arch: >=20 vertices,>=50% gingiva labels,planar centroid distance<=6mm from mirrored hint,span<=25mm,q95 radius<=12mm. Fit the visible homolog cervical rim to that ring by two-sided robust six-parameter ICP, then reconstruct SDF. Mesh control uses the identical ring observation and fitted parameters.', methods=['mirror', 'population_mean_sdf', 'neighbor_mesh', 'neighbor_sdf', 'ring_mesh', 'ring_sdf'])
    pr['data']['selection'] = 'Eight reserve patients fixed in SPLIT.json before R1 outcomes; no reserve geometry read in R1. Training prior and calibration constants unchanged. Anatomical eligibility identical; ring extraction failure remains a failure, no silent fallback.'
    pr['data']['counts']['scheduled_R2_teeth'] = 64
    pr['data']['mask_site'] = 'R2 now uses observed surviving gingival hole rim. It is not a removed-tooth border; the synthetic deletion footprint still encodes cervical site information. Deployment on healed extraction sites is outside this claim.'
    pr['ring_fit'] = {'observations': 'max256 observed masked-rim and256 homolog-boundary points, largest source boundary component', 'parameters': 'local xyz translation, rotation about visible arch normal,xy scales', 'bounds': 'translation+-3mm, rotation+-0.40rad,scales0.70..1.30', 'objective': 'bidirectional nearest-rim xyz residual /0.30mm; regularizer translation0.35/mm,angle0.35/0.20rad,scales0.35/0.20; scipy least_squares soft_l1 f_scale1,max80eval; initialrim centroidoffset clipped+-2.5mm,rotation0,scales1', 'consumer': 'patient-local crown placement prior plus benchmark methods using actual remaining-mask input', 'control': 'conventional robust mesh ICP with exactly same observations, fit and degrees of freedom; SDF reconstruction is compared to this direct mesh'}
    pr['gates']['data'] = 'At least6 of8 sealed test patients,3 tooth classes,>=80% of64 scheduled teeth completed. Coverage change is frozen for this smaller reserve before its numerical work; same error thresholds and decision criteria as R1.'
    pr['gates']['strong_control'] = 'ring_mesh is the identical-information conventional ICP control;0.05mm occlusal mean difference tolerance. No algorithmic novelty claim for a TIE.'
    pr['falsifiers'] += ['observed mask ring cannot be identified in>=80% of scheduled cases', 'synthetic hole-rim fit improves cervical boundary but degrades held-out upper surface', 'ring fit fails to beat mirror by frozen10% and0.50mm targets']
    pr['full_cost']['fit'] = 'All R1 baseline fits plus observed rim extraction and robust bidirectional registration; mesh and SDF share the fit but standalone fit cost charged to both'
    pr['full_cost']['validation'] = 'Sealed8-patient reserve untouched in R1; frozen predictions prior to hidden-reference measurement; same injection controls and numerical bound. Adaptive construction chosen after R1, evaluated on a distinct reserve.'
    pr['successor'] = 'If ring registration fails, next information port is independently labeled cervical/preparation margin or registered bite, not another parameterization of same fit. If ring helps, test failure of that assumption on partial fracture masks and actual preparation pairs.'
    pr['implementation_manifest'] = [{'path': str(p), 'sha256': sha(p)} for p in [pathlib.Path(__file__).resolve().parent / n for n in ['ring_method.py', 'crownbench.py', 'generate.py', 'evaluate.py']]]
    put(P / 'PREREG_R2.json', pr)
    decomp = copy.deepcopy(json.load(open(P / 'DECOMPOSITION.json')))
    decomp['round'] = 'R2'
    decomp['leaves'] += [{'name': 'visible_ring_observation', 'status': 'EXTERNALLY_MEASURED', 'equation': 'B=boundary(F_visible) selected by local centroid/extent/gingival-label filters; target vertices absent', 'source': 'Same published real IOS, masked mesh hash-bound in public reserve contexts', 'stop': 'Mask boundary is a deterministic transformation of independent acquired geometry, but the cervical footprint is inherited from synthetic deletion; it is not a measured prepared margin or healed socket.'}, {'name': 'ring_registration', 'status': 'CONSTITUTIVE_CLOSURE', 'equation': 'min_p sum rho(||T_p(b_source)-nearest_B||/.30mm)^2 + reverse residual + regularization', 'assumptions': 'Cervical boundary of contralateral tooth can register to observed hole rim. Translation,axial rotation andxy scales are sufficient; no held-out-surface fit.', 'stop': 'Natural asymmetry/crowding/missing gingiva can invalidate correspondence. Strong direct-mesh ICP receives the identical observations and can refute the SDF claim.'}]
    put(P / 'DECOMPOSITION_R2.json', decomp)
    handoff = '# R1 closed; R2 preregistered\n\n' + json.dumps(parent['decision'], indent=2) + '\n\nR1 metric controls: ' + str(parent['metric_controls_all_pass']) + '. Raw metrics and frozen output hashes preserved. No R1 metric or threshold changed.\n\nNext construction changes the information used: surviving gingival hole rim, robust cervical registration, and identical-information mesh control. Eight separately sealed patients are frozen in SPLIT.json. PREREG_R2.json precedes their mask preparation and generation.\n\nGraph coverage proposal remains PENDING_INDEPENDENT_REVIEW; no source graph modified.\n'
    put(P / 'HANDOFF_R1.json', {'decision': parent['decision'], 'result_sha256': sha(P / 'RESULTS_R1.json'), 'next': 'R2 observed ring on sealed reserve'})
    (P / 'HANDOFF_R1.md').write_text(handoff)
    (P / 'HANDOFF.md').write_text(handoff)
    state('R2_PREREG_FROZEN', parent['decision'], 'prepare sealed reserve contexts then ring fit without hidden tooth data')
    print('R2_FROZEN', sha(P / 'PREREG_R2.json'), flush=True)
if __name__ == '__main__':
    freeze()
