"""Explain the frozen PASS count against the official feasibility oracle.

This is a post-measurement diagnostic of the already preregistered obstruction,
not another estimator, threshold, query or physical certificate. No native
reference coordinates are read. The scope is the release's accepted scoring
contract: the official scorer rejects PASS on its INFEASIBLE rows.
"""
from common import *
import csv, gzip, collections

def certificate():
    counts = collections.Counter()
    reasons = collections.Counter()
    missed = collections.Counter()
    extra = collections.Counter()
    with gzip.open(DATA / 'evaluation/SCORED_ROWS.csv.gz', 'rt') as f:
        rows = [r for r in csv.DictReader(f) if r['split'] == 'test' and r['participant'] == 'constraint_optimizer']
    grouped = collections.defaultdict(list)
    for r in rows:
        grouped[r['case_key']].append(r)
    for (key, rr) in grouped.items():
        tt = {t['task_id']: t for t in tasks(key)}
        scenes = {}
        for r in rr:
            orig = next((t for t in tt.values() if t['family'] == r['family'] and t['level'] == r['level']))
            reason = None
            if r['feasibility'] != 'FEASIBLE':
                reason = 'official_' + r['feasibility']
            else:
                gf = orig['geometry_file']
                if gf not in scenes:
                    scenes[gf] = scene(orig)
                t = scenes[gf]
                if not t['A'].shape[0] and t['family'] != 'veneer':
                    reason = 'measured_antagonist_overlap_absent'
            attainable = reason is None
            counts['tasks'] += 1
            counts['operational_upper_bound'] += int(attainable)
            counts['LP_PASS'] += int(r['verdict'] == 'PASS')
            if reason:
                reasons[reason] += 1
            if attainable and r['verdict'] != 'PASS':
                missed[r['verdict']] += 1
            if not attainable and r['verdict'] == 'PASS':
                extra[reason] += 1
    out = dict(claim_type='capability', diagnostic_only=True, derived_from_preregistered_obstruction='DECOMPOSITION.json optimal pass ceiling leaf', counts=dict(counts), excluded_by_reason=dict(reasons), unattained_by_LP=dict(missed), unexpected_PASS=dict(extra), operational_ceiling_saturated=not missed and (not extra), external_referent={'kind': 'published_code', 'locator': str(BENCH / 'code/scorer.py'), 'compared_quantity': 'PASS on official fixed roof feasibility and availability contract', 'refutes_us': True}, assumptions=['Official evaluator rejects PASS on INFEASIBLE rows', 'Flat intaglio can give wall and ideal milling eligibility in this contract', 'No overlapping antagonist gives UNKNOWN except veneer by the original rule', 'No tolerance or clinical accuracy guarantee beyond the scorer'], does_not_prove='Full3D physical infeasibility; mathematical exactness of floating-point LP; film/wall feasibility outside fixed interface', resolution='PER_TOOTH; counts aggregate to POPULATION', new_hidden_queries=0, source_rows_sha256=sha(DATA / 'evaluation/SCORED_ROWS.csv.gz'))
    dump(ROOT / 'raw/PASS_CEILING_CERTIFICATE.json', out)
    return out
if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
