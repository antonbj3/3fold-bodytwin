"""Fault injection at the prospective lab interfaces; examples are explicitly fixtures."""
import thread_budget
import copy, csv, json, tempfile
from pathlib import Path
import numpy as np
from common import ROOT, save, sha
from lot_cli import freeze as lot_freeze, evaluate
from qc_cli import freeze as qc_freeze, score

def main():
    check = []
    rng = np.random.default_rng(5104)
    examples = ROOT / 'examples'
    examples.mkdir(exist_ok=True)
    proto_sha = 'a' * 64
    context = {'lot_id': 'FIXTURE_NOT_LAB', 'configuration': 'fixed100N example', 'protocol_locator': 'our_own_fixture', 'protocol_sha256': proto_sha, 'quantity': 'cycles', 'unit': 'cycles', 'proof_threshold': 5000000, 'planned_ids': ['test_%02d' % i for i in range(1, 47)], 'p_bad': 0.1, 'p_good': 0.01, 'power': 0.8, 'alpha': 0.05}
    records = [dict(specimen_id=i, lot_id=context['lot_id'], configuration=context['configuration'], protocol_sha256=proto_sha, value=5000000, failure=0, unit='cycles', locator='our_own_fixture:' + i) for i in context['planned_ids']]
    records[0].update(value=1000000, failure=1)
    protocol = {'protocol_sha256': proto_sha, 'frame_id': 'external_fiducial_fixture', 'regions': ['MG'], 'weeks': 52, 'weekly_n': 5, 'alpha': 0.05, 'source_kind': 'OUR_OWN_FIXTURE', 'slope_bound_um_per_day': None, 'hard_error_d_um': None, 'gap0_um': {'MG': 50}, 'gap_limit_um': 120, 'calibration_locator': None, 'model_scope': 'COMMON_ADDITIVE_SIGNED_REGION'}

    def row(i, r, c):
        return dict(measurement_id=str(i), region='MG', reference_um=str(r), coupon_um=str(c), unit='um', frame_id=protocol['frame_id'], protocol_sha256=proto_sha, locator='our_own_fixture:' + str(i))
    baseline = [row(i, float(rng.normal(0, 10)), float(rng.normal(0, 10))) for i in range(100)]
    future = [row(i, float(rng.normal(0, 1)), float(rng.normal(0, 1))) for i in range(5)]
    for (name, o) in [('lot_context.json', context), ('lot_measurements.json', records), ('qc_protocol.json', protocol)]:
        save('examples/' + name, o)
    for (name, rs) in [('qc_baseline.csv', baseline), ('qc_week1.csv', future)]:
        with (examples / name).open('w') as f:
            w = csv.DictWriter(f, fieldnames=list(rs[0]))
            w.writeheader()
            w.writerows(rs)
    with tempfile.TemporaryDirectory(prefix='X51-cli-', dir=ROOT / 'raw') as td:
        p = Path(td)
        lf = p / 'lot_plan.json'
        qf = p / 'qc_prediction.json'
        state = p / 'qc_state.json'
        lot_freeze(context, lf)
        out = evaluate(lf, records)
        wrong = copy.deepcopy(records)
        wrong[1].update(value=500, failure=1)
        check.append(dict(name='LOT_ACCEPTANCE_BOUNDARY', nominal_pass=out['verdict'] == 'ACCEPT_RISK_CRITERION_ONLY', injected_error_rejected=evaluate(lf, wrong)['verdict'] == 'NOT_CERTIFIED'))
        pending = copy.deepcopy(records)
        pending[1].update(value=500, failure=0)
        check.append(dict(name='CENSOR_NO_FALSE_SURVIVAL', nominal_pass=evaluate(lf, pending)['unknown_endpoints'] == 1, injected_error_rejected=evaluate(lf, pending)['verdict'] == 'NOT_CERTIFIED'))
        check.append(dict(name='NO_OPTIONAL_ACCEPTANCE', nominal_pass=evaluate(lf, records[:5])['verdict'] == 'INCOMPLETE_FIXED_PLAN', injected_error_rejected=evaluate(lf, records[:5])['verdict'] != 'ACCEPT_RISK_CRITERION_ONLY'))

        def rejects(fn):
            try:
                fn()
                return False
            except ValueError:
                return True
        dup = copy.deepcopy(records)
        dup[1]['specimen_id'] = dup[0]['specimen_id']
        badunits = copy.deepcopy(records)
        badunits[0]['unit'] = 'N'
        check.append(dict(name='LOT_UNITS_AND_DUPLICATES', nominal_pass=True, injected_error_rejected=rejects(lambda : evaluate(lf, dup)) and rejects(lambda : evaluate(lf, badunits))))
        qc_freeze(baseline, protocol, qf)
        q = score(qf, future, 1, state)
        altered = copy.deepcopy(future)
        altered[0]['coupon_um'] = '100'
        check.append(dict(name='QC_WEEK_LOCK', nominal_pass=score(qf, future, 1, state) == q, injected_error_rejected=rejects(lambda : score(qf, altered, 1, state))))
        uncalibrated = copy.deepcopy(protocol)
        uncalibrated.update(source_kind='LAB_MEASUREMENT', slope_bound_um_per_day=0.0, hard_error_d_um=0.0)
        qbad = p / 'uncalibrated_prediction.json'
        qc_freeze(baseline, uncalibrated, qbad)
        injected = score(qbad, future, 1, p / 'uncalibrated_state.json')
        check.append(dict(name='PHYSICAL_UNKNOWN_HOLD', nominal_pass=q['verdict'] == 'HOLD_UNKNOWN_OR_HORIZON', injected_error_rejected=injected['regions']['MG']['horizon_guard']['verdict'] == 'HOLD_UNKNOWN', injected_protocol='LAB_MEASUREMENT and zero numerical bounds without calibration locator'))
        badframe = copy.deepcopy(future)
        badframe[0]['frame_id'] = 'bestfit_crown'
        check.append(dict(name='REGISTRATION_FRAME', nominal_pass=True, injected_error_rejected=rejects(lambda : score(qf, badframe, 2, state))))
        qf.write_text(qf.read_text().replace('OUR_OWN_FIXTURE', 'LAB_MEASUREMENT'))
        check.append(dict(name='FROZEN_QC_HASH', nominal_pass=True, injected_error_rejected=rejects(lambda : score(qf, future, 2, state))))
        save('raw/CONSUMER_OUTPUTS.json', {'lot': out, 'qc': q, 'examples_kind': 'our_own_fixture, no external physical facit'})
    save('raw/CONSUMER_CHECKS.json', {'checks': check, 'all_pass': all((c['nominal_pass'] and c['injected_error_rejected'] for c in check)), 'external_referent': {'kind': 'our_own_fixture', 'locator': 'examples/*.json,examples/*.csv', 'compared_quantity': 'input contract and CLI output', 'refutes_us': True}})
    if not all((c['nominal_pass'] and c['injected_error_rejected'] for c in check)):
        raise RuntimeError('Consumer fault gate failed')
    print('Prospective consumers:', len(check), 'rejecting controls PASS')
if __name__ == '__main__':
    main()
