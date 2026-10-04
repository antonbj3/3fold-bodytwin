from common import *
from stats import collapse, summarize
from guide import bound, displacement95, scalar_control
from scipy.stats import beta
start = time.perf_counter()
check()
rows = [json.loads(l) for l in (DATA / 'R2_GEOMETRY.jsonl').read_text().splitlines()]
r0 = [json.loads(l) for l in (DATA / 'R1_GEOMETRY.jsonl').read_text().splitlines() if json.loads(l)['extra_mm'] == 0]
assert len(rows) == 4 * len(certs()), 'R2 geometry incomplete; run geometry first'
profiles = load(X8 / 'GUIDE_PROFILES.json')['profiles']
protocols = load(ROOT / 'PROTOCOLS.json')['protocols']

def cal(g):
    return int(hashlib.sha256(('X96-cal:' + g).encode()).hexdigest(), 16) % 100 < 60
scores = {}
for r in r0:
    if r['paired_revisions']:
        scores[r['image_group']] = max(scores.get(r['image_group'], 0), max(0, r['nominal_tf2_upper_mm'] - r['union_lower_mm']))
cal_scores = [v for (g, v) in scores.items() if cal(g)]
test_scores = [v for (g, v) in scores.items() if not cal(g)]
rank = int(np.ceil((len(cal_scores) + 1) * 0.95))
q = sorted(cal_scores)[rank - 1] if rank <= len(cal_scores) else None
ok = sum((v <= q for v in test_scores)) if q is not None else 0
n = len(test_scores)
ci = [float(beta.ppf(0.025, ok, n - ok + 1)) if ok else 0, float(beta.ppf(0.975, ok + 1, n - ok)) if ok < n else 1.0] if n else [None, None]
calibration = dict(split='sha256(X96-cal:<image_group>) mod100 <60 calibration; remaining test', calibration_groups=len(cal_scores), test_groups=n, quantile_rank=rank, revision_budget95_mm=q, heldout_success=ok, heldout_coverage=ok / n if n else None, heldout_binomial_ci95=ci, empirical95_gate=bool(n and ok / n >= 0.95), calibration_scores=[dict(image_group=g, score_mm=v, split='calibration' if cal(g) else 'test') for (g, v) in sorted(scores.items())], scope='known-image-group maximum revision loss only at frozen nominal poses; not independent wall error', resolution='POPULATION')
dump(ROOT / 'raw/R2_REVISION_CALIBRATION.json', calibration)
tables = []
guide_rows = []
per_site = []
controls = []
for pr in protocols:
    h = pr['extra_mm']
    rr = [r for r in rows if r['extra_mm'] == h]
    tab = summarize(h, rr)
    tab['protocol'] = pr['id']
    tables.append(tab)
    cc = collapse(rr)
    for p in profiles:
        probs = bound([r['union_lower_mm'] for r in rr], p, h)
        for (r, u) in zip(rr, probs):
            per_site.append(dict(case=r['case'], fdi=r['fdi'], image_group=r['image_group'], protocol=pr['id'], guide=p['id'], clearance_mm=r['union_lower_mm'], moment_upper=float(u), true_geometric_event_probability='UNKNOWN', physical_nerve_injury_probability='UNKNOWN', paired=r['paired_revisions'], resolution='PER_TOOTH', time_scale='SIMULTANEOUS'))
        bypose = {(r['case'], r['fdi']): float(v) for (r, v) in zip(rr, probs)}
        for cohort in ['all', 'paired']:
            selected = [c for c in cc if c['accepted'] and (cohort == 'all' or c['paired'])]
            groups = {}
            for c in selected:
                u = min(1, sum((bypose[r['case'], r['fdi']] for r in c['rows'] if r['nominal_tf2_lower_mm'] >= 2)))
                groups.setdefault(c['image_group'], []).append(u)
            nums = np.array([sum(v) for v in groups.values()])
            dens = np.array([len(v) for v in groups.values()])
            rng = np.random.default_rng(9601)
            ix = rng.integers(0, len(nums), size=(2000, len(nums)))
            boot = nums[ix].sum(1) / dens[ix].sum(1)
            gd = displacement95(p, h)
            m = 2 + h + q + gd if q is not None else None
            guide_rows.append(dict(protocol=pr['id'], guide=p['id'], cohort=cohort, sites=len(selected), image_groups=len(groups), physical_frequency='UNKNOWN', conditional_geometric_frequency_interval=[0.0, float(nums.sum() / dens.sum())], upper_bound_bootstrap95=np.quantile(boot, [0.025, 0.975]).tolist(), guide_displacement95_mm=gd, required_nominal_margin_joint90_mm=m, revision_budget95_mm=q, formula='2+h+A95+G95(h); union bound .05+.05; fixed source moments exact and revision-score exchangeability assumed', scope='PHENOMENOLOGICAL sufficient margin, not minimal or clinically required; all cohort lacks revision info for unpaired sites', physical_validation=False, resolution='PHENOMENOLOGICAL'))
        for d in [2.0, 3.0, 5.0, 10.0, 20.0]:
            v = float(bound([d], p, h)[0])
            c = float(scalar_control(d, p, h))
            controls.append(dict(protocol=pr['id'], guide=p['id'], clearance_mm=d, error=abs(v - c), injected_plus0p1_rejected=bool(abs(v + 0.1 - c) > 1e-06)))
import csv
with (ROOT / 'raw/PER_SITE_GUIDE_R2.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=list(per_site[0]))
    w.writeheader()
    w.writerows(per_site)
dump(ROOT / 'raw/R2_GUIDE_CONTROLS.json', controls)
st = (ROOT / 'sources/straumann_blx.txt').read_text()
at = (ROOT / 'sources/proof_lane_ev_guided.txt').read_text()
zt = (ROOT / 'sources/zimvie_t3.txt').read_text()
import re
source_checks = dict(blx=re.search('tip is up to 0\\.5\\s*mm longer', st) != None, proof_lane=re.search('tip can be up to 1\\s*mm longer', at) != None, zimvie_tip=re.search('3\\.85 mm\\s+N/A\\s+1\\.2 mm', zt) != None, zimvie_datum=re.search('13 mm\\s+12\\.6 mm\\s+1 mm\\s+0\\.4 mm\\s+13\\.7 mm', zt) != None)
from decimal import Decimal as D
aligned = D('13.7') + D('1.2') - D('12.6') - D('1.0')
source_checks['datum_arithmetic'] = aligned == D('1.3')
source_checks['injected_tip_only_rejected'] = aligned != D('1.2')
source_checks['injected_blx1p5_rejected'] = protocols[0]['extra_mm'] != 1.5
dump(ROOT / 'raw/R2_SOURCE_CONTROLS.json', source_checks)
res = dict(round='R2', claim_type='information_link', outcome='PROTOCOL_DATUM_AND_GUIDE_CONDITIONAL_BOUNDS_COMPUTED_TRUE_POPULATION_FREQUENCY_UNKNOWN', protocol_tables=tables, guide_tables=guide_rows, calibration={k: v for (k, v) in calibration.items() if k != 'calibration_scores'}, gates=dict(coverage=len(rows) == 4 * len(certs()), source_controls=all(source_checks.values()), guide_controls=all((r['error'] <= 1e-06 and r['injected_plus0p1_rejected'] for r in controls)), practice_consequence_all_protocols=all((t['all']['fraction'] >= 0.05 for t in tables)), revision_transfer=calibration['empirical95_gate']), cost=dict(geometry=load(ROOT / 'raw/R2_COST.json'), analysis_validation_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss), external_referent=load(ROOT / 'PREREG_R2.json')['external_referent'], zimvie_tip_only_ablation=summarize(1.2, [r for r in rows if r['extra_mm'] == 1.2]), formal_floating_point_enclosure='MISSING')
dump(ROOT / 'rounds/R2.json', res)
(ROOT / 'HANDOFF_R2.md').write_text('# R2 protocol and direction\n\n' + json.dumps({k: v for (k, v) in res.items() if k != 'guide_tables'}, indent=2) + "\n\nNew sources changed direction: Universal 1.5 mm is not an identified manufacturer's protocol; version-bound dimensions and actual apex reference are needed. ZimVie tip-only 1.2 mm is separated from reference-adjusted 1.3 mm. Physical guide frequency UNKNOWN ; the entire tool's targeted bugs and independent canal wall are missing. The next design will limit the possible shape of the drill point instead of calling the full radius housing real drill.\n")
state('R2_COMPLETE', res['gates'], 'Freeze tool-shape identified interval and axis witness')
print(json.dumps(dict(tables=[dict(protocol=t['protocol'], all=t['all'], paired=t['paired_revisions']) for t in tables], calibration=res['calibration'], gates=res['gates']), indent=2))
