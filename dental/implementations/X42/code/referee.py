"""Lane-owned trusted one-query evaluation; never writes upstream checkout.

The official V3 Bundle/evaluator check original release bytes. The lane owns
only the submission manifest, its ledger and aggregate extensions. The
participant process cannot access this file, ledger or the source references.
"""
from common import *
from scorer import score_submission
from release_anchor import BENCHMARK_SHA256
from analysis_tools import paired, decide
import gzip, csv, io, collections, fcntl

def rows_from_official(path):
    raw = list(csv.DictReader(io.StringIO(gzip.decompress(Path(path).read_bytes()).decode())))
    rows = []
    for r in raw:
        a = {'case_key': r['case_key'], 'family': r['family'], 'level': r['level'], 'method': r['participant'], 'verdict': r['verdict'], 'dataset': r['dataset'], 'split': r['split'], 'feasibility': r['feasibility']}
        for (k, old) in [('rmse_mm', 'height_rmse_mm'), ('contact_error_mm2', 'contact_area_error_mm2')]:
            a[k] = float(r[old]) if r[old] else None
        rows.append(a)
    return rows

def eval_once(predroot, receipt, replay=False):
    ledger = ROOT / 'raw/QUERY_LEDGER.json'
    lock = ROOT / 'raw/QUERY_LEDGER.lock'
    out = DATA / ('replay_evaluation' if replay else 'evaluation')
    with lock.open('a+') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        if ledger.exists():
            old = json.load(open(ledger))
            if not replay or old['prediction_freeze_sha256'] != receipt:
                raise RuntimeError('one frozen hidden-query budget exhausted; only same receipt replay allowed')
        elif replay:
            raise RuntimeError('cannot replay an unmeasured freeze')
        else:
            dump(ledger, dict(query_count=1, charged_utc=now(), benchmark_sha256=BENCHMARK_SHA256, prediction_freeze_sha256=receipt, policy='Joint selected construction + frozen matched controls; no generation or selection after test feedback', upstream_mutated=False, official_api='score_submission; lane-owned ledger instead of upstream intake mutation'))
    result = score_submission(BENCH, predroot, receipt, out)
    rows = rows_from_official(out / 'SCORED_ROWS.csv.gz')
    sub = [r for r in rows if r['split'] == 'test']
    cfg = json.load(open(ROOT / 'FINAL_CONFIG.json'))
    candidate = cfg['selected_candidate']
    contrasts = []
    for control in cfg['participants']:
        if control == candidate:
            continue
        contrasts.extend(paired(sub, candidate, control))
    perlevel = []
    for level in ['easy', 'normal', 'hard', 'boundary']:
        lr = [r for r in sub if r['level'] == level]
        for control in cfg['participants']:
            if control != candidate:
                perlevel.extend([dict(d, level=level) for d in paired(lr, candidate, control)])
    family_level = []
    for ds in ['Bite2Text', 'Bits2Bites']:
        ss = [r for r in rows if r['dataset'] == ds and (r['split'] == 'test' if ds == 'Bite2Text' else True)]
        for fam in sorted({r['family'] for r in ss}):
            for level in ['easy', 'normal', 'hard', 'boundary']:
                for name in cfg['participants']:
                    rr = [r for r in ss if r['family'] == fam and r['level'] == level and (r['method'] == name)]
                    group = collections.defaultdict(list)
                    for r in rr:
                        group[r['case_key']].append(int(r['verdict'] == 'PASS'))
                    y = np.array([np.mean(v) for v in group.values()])
                    rng = np.random.default_rng(6142)
                    boot = np.array([y[rng.integers(len(y), size=len(y))].mean() for _ in range(2000)])
                    family_level.append(dict(dataset=ds, family=fam, level=level, method=name, counts=dict(collections.Counter((r['verdict'] for r in rr))), case_clusters=len(y), strict_pass_fraction=float(y.mean()), case_bootstrap_95=np.quantile(boot, [0.025, 0.975]).tolist(), resolution='POPULATION'))
    gates = {d['control']: decide(d) for d in contrasts if d['family'] == 'all'}
    scope = dict(primary_test='588 Bite2Text patient-case clusters; no independent-person claim beyond upstream case documentation', difficulty_levels='easy/normal/hard/boundary are requirement scenarios, not levels1/2/3 of evidence', level1='Roof wall, nominal planar film, static continuous nonpenetration, fixed insertion; ideal planar intaglio only. Full milling/sidewalls/margins NOT_SCORED.', level3='Native source-height and static occlusal proximity area only; approximal contact and substance removal NOT_SCORED.', measured_mirror='NOT_AVAILABLE in equal-information public interface; cannot claim defeat', clinical='No physical or clinical validation', multiple_comparisons='Primary all-family paired gate preregistered; per-family/level intervals descriptive without simultaneous multiplicity adjustment')
    agg = dict(claim_type='algorithm', selected_candidate=candidate, external_referent=json.load(open(ROOT / 'PREREG_R1.json'))['external_referent'], benchmark_sha256=BENCHMARK_SHA256, prediction_freeze_sha256=receipt, official_evaluation=result, primary_test_summary=summarize(sub), paired_contrasts=contrasts, paired_by_level=perlevel, family_level=family_level, gates=gates, scope=scope, hidden_queries_charged=1, replay_of_same_frozen_predictions=replay)
    dump(ROOT / 'raw' / ('REPLAY_AGGREGATES.json' if replay else 'FINAL_AGGREGATES.json'), agg)
    return agg
