from common import *
from model_fit import dataset, fit_arrays

def run():
    prereg = {'claim_type': 'algorithm', 'frozen_utc': now(), 'kind': 'stronger_control_before_hidden_query', 'changed_operation': 'Direct equally informed standard multioutput kernel ridge without the candidate12-mode output truncation, paired with the same full contact field and beta1 branch. Same gamma1/ridge0.1 chosen before test and same208 dev donor set.', 'candidate_changes': 'NONE; original FROZEN_GENERATOR.json remains immutable', 'metric_and_gate_changes': 'NONE; strongest controls must refute same frozen primary goal', 'name': 'control_full_kernel', 'external_referent': json.load(open(ROOT / 'PREREG_R1.json'))['external_referent'], 'full_cost': {'fit': '9 full-output model fits counted', 'discovery': 'no new candidate tuning; direct counterpart included from equations', 'validation': 'same exact inequalities and native metrics; one joint query', 'queries': '0 before extension; no second test query', 'fallback': 'same LP/ABSTAIN'}, 'assumption': 'Kernel matrix and full scene features identical to candidate; only output rank truncation removed', 'stop_argument': 'This is the direct conventional computation underlying candidate, not a novel algorithm claim.'}
    dump(ROOT / 'PREREG_CONTROL_EXTENSION.json', prereg)
    (ROOT / 'PREREG_CONTROL_EXTENSION.sha256').write_text(sha(ROOT / 'PREREG_CONTROL_EXTENSION.json') + '  PREREG_CONTROL_EXTENSION.json\n')
    (ds, reject, access) = dataset(split()['dev_all'])
    base = DATA / 'FULL_KERNEL_CONTROL'
    base.mkdir(exist_ok=True)
    cost = []
    for (fam, rr) in ds.items():
        start = time.perf_counter()
        m = fit_arrays(rr, 'kernel', (1.0, 0.1), output_modes=0)
        p = base / (fam + '.npz')
        np.savez_compressed(p, **m)
        cost.append(dict(family=fam, seconds=time.perf_counter() - start, bytes=p.stat().st_size, sha256=sha(p)))
    dump(ROOT / 'raw/FULL_KERNEL_FIT.json', dict(fits=cost, reject=reject, reference_access_keys=access, test_references='NOT_ACCESSED'))
    cfg = json.load(open(ROOT / 'FINAL_CONFIG.json'))
    cfg['participants'].append('control_full_kernel')
    cfg['control_extension'] = 'standard full output kernel; same contact field beta1'
    dump(ROOT / 'EVAL_CONFIG.json', cfg)
    print('Direct full-output kernel control frozen; candidate unchanged')
if __name__ == '__main__':
    run()
