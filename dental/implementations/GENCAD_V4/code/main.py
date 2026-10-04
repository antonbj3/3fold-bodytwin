from common import *
from integrity import verify

def run():
    start = time.perf_counter()
    receipt = verify()
    state('REPLAYING_RELEASE', 'Input/code/runtime hashes verified', 'Execute six frozen rounds and full-preparation replay')
    import discrimination, pareto_probe, replay_preparation, score_panel, score_whole, literature, report, support_probe
    discrimination.run()
    pareto_probe.run()
    support_probe.run()
    replay_preparation.run()
    score_panel.run()
    score_whole.run()
    score_whole.run('R5')
    literature.evaluate()
    verify()
    report.run()
    verify()
    dump(ROOT / 'raw/DEMO_RECEIPT.json', dict(input_integrity=receipt, seconds=time.perf_counter() - start, semantic_digest=read(ROOT / 'raw/SEMANTIC_DIGEST.json')['sha256'], physical_measurement='NOT_RUN'))
    state('DEMO_REPRODUCED', 'Seven rounds re-evaluated; full-crown anatomical gates remain failed', 'Independent review and source-supported boundary/registered prepared specimen')
    print('DEMO_REPRODUCED')
if __name__ == '__main__':
    run()
