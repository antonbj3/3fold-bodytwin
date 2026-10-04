import json, hashlib, datetime, sys
from pathlib import Path
r = sys.argv[1]
p = json.loads(Path('PREREG_R1v2.json').read_text())
p['round'] = r
p['frozen_at_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
p['parent_artifact'] = 'raw/R1_RESULTS.json'
p['claim_type'] = 'information_link'
if r == 'R2':
    p['capability'] = 'Recover an unambiguous reference-status decision by observing only the joint distribution of the two non-invasive tests, without additional reference observations.'
    p['obstacle'] = 'Marginal tables omit dependence; an unlabeled joint frequency may remove it without opening a tooth.'
    p['changed_operation'] = 'Add exact observed joint event counts C1P1=264 and C1P0=116, C1Pmissing=0 to the 18-atom polytope. Do NOT use Table2 PPV/SN/SP.'
    p['input'] = 'R1 Figure1 reference/marginal counts plus Table2 joint event totals only; joint reference counts excluded.'
    p['metrics']['decision_gate'] = 'all four posterior widths <=0.10; fixed threshold grid remains 0.1,...,0.9; t=0.8 highlighted as declared benchmark, not clinical.'
    p['strongest_equally_informed_control'] = 'Independent 18-atom LP with same unlabeled joint totals. Compare acquired-information result against R1 without these totals, not LP superiority.'
    p['falsifier'] = 'Unlabeled joint observation leaves width>0.10 for any pattern, or endpoint control disagreement>1e-8.'
    p['full_cost']['external_acquisition'] = 'Two published event counts reused; hypothetical prospective joint tally is NOT_RUN and its time UNKNOWN. No destructive verification supplied.'
else:
    raise ValueError(r)
Path(f'PREREG_{r}.json').write_text(json.dumps(p, indent=2) + '\n')
Path(f'PREREG_{r}.sha256').write_text(hashlib.sha256(Path(f'PREREG_{r}.json').read_bytes()).hexdigest() + '\n')
d = json.loads(Path('DECOMPOSITION_R1v2.json').read_text())
d['round'] = r
d['changed_operation'] = p['changed_operation']
d['leaves'].append({'leaf': 'observed_joint_event_totals', 'status': 'EXTERNALLY_MEASURED', 'basis': 'PMC4884138 Table2 joint positive counts 264 and116.', 'stop_argument': 'This gives event denominator, not reference-conditioned numerator; no PPV is smuggled into an unlabeled measurement.'})
Path(f'DECOMPOSITION_{r}.json').write_text(json.dumps(d, indent=2) + '\n')
