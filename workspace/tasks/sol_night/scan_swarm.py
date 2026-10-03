"""Compact summary of new The swarm/reserve_worker results since the previous scan, for the coordinator's combinatorial view.
Writes tasks/build_night/swarm/digest_<time>.md and updates the marker. Usage: python3 scan_swarm.py [--since-min M] [--max N]"""
import json, sys, time, itertools, collections
from pathlib import Path

STORE = Path('/mnt/games-240/research/bunny48_20260926/bodytwin')
D = Path('./tasks/build_night'); OUT = D / 'swarm'; OUT.mkdir(exist_ok=True)
MARK = OUT / 'last_scan'

def rd(p, n=None):
    try: t = p.read_text(errors='replace'); return t[:n] if n else t
    except Exception: return ''

def js(p):
    try: return json.loads(p.read_text())
    except Exception: return {}

args = sys.argv[1:]
since = float(MARK.read_text()) if MARK.exists() else time.time() - 3 * 3600
if '--since-min' in args: since = time.time() - 60 * float(args[args.index('--since-min') + 1])
cap = int(args[args.index('--max') + 1]) if '--max' in args else 60
now = time.time()

rows = []
for r in list(STORE.glob('*/RESULTS.md')) + list(Path('./results').glob('BT-XSEED-*/RESULTS.md')):
    try: m = r.stat().st_mtime
    except Exception: continue
    if m <= since: continue
    d = r.parent; job = js(d / 'JOB.json'); route = js(d / 'MODEL_ROUTE.json')
    model = str(route.get('model') or route.get('selected_model') or '?').split('/')[-1]
    res = rd(r)
    head = ' '.join(l.strip() for l in res.splitlines()[1:12] if l.strip() and not l.startswith('#'))[:420]
    rows.append(dict(id=d.name, t=time.strftime('%H:%M', time.localtime(m)), model=model, cat=job.get('category', '?'),
                     target=job.get('target_id', '?'), keys=job.get('source_keys', []), decision=(job.get('decision') or '')[:200],
                     head=head, followups=(d / 'FOLLOWUPS.json').exists() or (d / 'SWARM_QUEUE_ADD.json').exists()))
rows.sort(key=lambda x: x['t'])

L = [f"# Warm scan {time.strftime('%d/%m %H:%M')} — {len(rows)} nya resultat sedan {time.strftime('%H:%M', time.localtime(since))}", '']
L.append('Modeller: ' + ', '.join(f'{k} {v}' for k, v in collections.Counter(r['model'] for r in rows).most_common()))
L.append('Kategorier: ' + ', '.join(f'{k} {v}' for k, v in collections.Counter(r['cat'] for r in rows).most_common()))
# Share of new jobs (last 2 h) with an external facit according to research_value/external_referent (Anton 1/10).
try:
    import sys as _sys; _sys.path.insert(0,'~/research/AGENT_DASHBOARD_20260930'); import research_value as _rv
    _jobs=[p for p in STORE.glob('*/JOB.json') if time.time()-p.stat().st_mtime<7200]
    # 2/10 19:05 (anton-5f): this line counted ONLY the structured external_referent field, so the
    # requirement that free_controller._normalize injects into the job's BODY read as 'saknas'. At
    # 18:57 it printed 'saknas 62 av 72' and I nearly booked a regression; measured directly, all
    # 116 jobs from the last two hours carry the requirement - 47 in the field and 69 in the body.
    # Counting a requirement as absent because it arrived by the other of its two routes is the same
    # miscount class as the four we logged today, in an instrument I own.
    def _req(path):
        raw = json.loads(path.read_text())
        kind = (_rv.assess(raw).get('external_referent') or {}).get('kind')
        if kind:
            return kind
        return 'i brieftext' if 'external_referent' in str(raw.get('body', '')) else 'saknas helt'
    _k=collections.Counter(_req(p) for p in _jobs)
    L.append("External reference REQUIRED in new jobs (2 h), field or brief text:"
             + ', '.join(f'{k} {v}' for k,v in _k.most_common()) + f' av {len(_jobs)}')
except Exception as _e:
    L.append(f'External reference observation: could not be measured ({type(_e).__name__})')
# 2/10 (anton-5f): the line above reads JOB.json, i.e. what the PARENT proposed - that is never filled in afterwards
# and therefore always looks bad. What is decisive is what the job DELIVERS in results.json. Measure both, otherwise
# every tick reads as a regression (I made that mistake at 03:55 and 05:10).
try:
    # 3/10 04:45: BT-DATG-jobb is DATAINSAMLING with its own file format and lacks JOB.json completely.
    # They legitimately have no external_referent — they SAMLAR external data. To count them in the same
    # denominators such as research jobs made the reference statistics misleading: 56 of 121 s went out to
    # missing reference, while every research job with job file supplied it (69 of 69).
    _done=[p for p in STORE.glob('*/results.json')
           if time.time()-p.stat().st_mtime<7200 and (p.parent/'JOB.json').exists()]
    _datg=[p for p in STORE.glob('*/results.json')
           if time.time()-p.stat().st_mtime<7200 and not (p.parent/'JOB.json').exists()]
    _d=collections.Counter()
    for _p in _done:
        try: _e2=json.loads(_p.read_text()).get('external_referent')
        except Exception: _d["unreadable"]+=1; continue
        _d[(_e2.get('kind') if isinstance(_e2,dict) else "inget field") or 'tomt']+=1
    _ok=sum(v for k,v in _d.items() if k not in ("inget field",'tomt',"unreadable"))
    L.append(f'Yttre facit LEVERERAT (results.json, 2 h): {_ok} av {len(_done)} — ' + ', '.join(f'{k} {v}' for k,v in _d.most_common()))
except Exception as _e:
    L.append(f'Delivered reference observations: could not be measured ({type(_e).__name__})')
# 2/10 21:10 (anton-5f): etiketten omdopt. 'Okorsade par' raknar par som aldrig NAMNTS ihop i nagot
# job's key list, which is not the same as them being able to JAMFORAS. An integration job in math
# pain that 900 of 903 couple of years UNDECIDED therefore that 40 of 43 families do not get their artifacts in
# inputs/, medan den har raden samtidigt sa '0 par aldrig namnda ihop, av 903'. Jag laste det som full tackning
# and quoted it further. The line must say what it matters.
L += ['', '## Nyckelpar som aldrig NAMNTS ihop (ej samma sak som jamforbara)', '']
# Keys from new constructive results; pairs that have never occurred together in any job in STATE.json.
try: J = json.loads(Path('./tasks/free48/STATE.json').read_text())['jobs']
except Exception: J = {}
co = set()
for r in J.values():
    ks = sorted(set(r.get('source_keys') or []))
    co.update(itertools.combinations(ks, 2))
cons = [r for r in rows if r['cat'] in ('mechanism', 'integration', 'calibration', 'exploratory')]
freq = collections.Counter(k for r in cons for k in set(r['keys']))
ex = {k: next(r['id'] for r in cons if k in r['keys']) for k in freq}
cand = [(freq[a] + freq[b], a, b) for a, b in itertools.combinations(sorted(freq), 2) if (a, b) not in co]
cand.sort(reverse=True)
L.append(f'{len(freq)} nycklar i {len(cons)} konstruktiva resultat; {len(cand)} okorsade par av {len(freq)*(len(freq)-1)//2}.')
for n, a, b in cand[:20]: L.append(f'- {a} × {b}  (nya resultat {freq[a]}+{freq[b]}; ex {ex[a]} / {ex[b]})')
L += ['', '## Resultat (nyast sist)', '']
for r in rows[-cap:]:
    L.append(f"### {r['id']} · {r['t']} · {r['model']} · {r['cat']} · {r['target']} · {'+'.join(r['keys'])}")
    L.append(f"Beslut: {r['decision']}")
    L.append(f"Resultat: {r['head']}")
    L.append('')
p = OUT / f'digest_{time.strftime("%m%d_%H%M")}.md'; p.write_text('\n'.join(L))
MARK.write_text(str(now)); print(p, len(rows))
