"Quick refill (Anton 24/9 16:3x: maximize The_swarm NOW). Reproduktions-/robusthetspaket BT-Q-<src> from existing results.\nDeterministic, no LLM. Copies the source folder's small files (≤30 MB total, no >10 MB-filer) as inputs/."
import os, shutil, sys, json
W=''; R=W+'/results'
N=int(sys.argv[1]) if len(sys.argv)>1 else 200
queued={l.split()[2] for l in open(W+'/tasks/lanes/bt_queue.txt') if len(l.split())==3}
srcs=sorted([d for d in os.listdir(R) if d.startswith(('BT-','N','CX-')) and not d.startswith(('BT-Q-','BT-R-','BT-P-','BT-X-','BT-AG-'))
       and os.path.isfile(f'{R}/{d}/RESULTS.md') and os.path.isfile(f'{R}/{d}/results.json')], key=lambda d:-os.path.getmtime(f'{R}/{d}/RESULTS.md'))
out=[]; prof='ABC'; k=0
for s in srcs:
    if len(out)>=N: break
    for fam,q in (('Q',"REPRODUCTION: recalculate maximum 3 key figures from inputs/src/results.json with OWN code from raw data/code in inputs/src/. Deviation >1 % relative = flag. Can't it be recalculated from what's there: UNKNOWN with exactly what's missing."),
                  ('S',"ROBUSTHET: select ONE analysis choice on which the conclusion in inputs/src/RESULTS.md rests (threshold ±20 %, other normalization, omit one person/sample, other seed). Rerun with your own code. Does the conclusion hold? countertest: permutation/placebo.")):
        j=f'BT-{fam}-{s}'
        if j in queued or os.path.exists(f'{R}/{j}'): continue
        files=[]; tot=0
        for root,_,fs in os.walk(f'{R}/{s}'):
            if '__pycache__' in root or '/.' in root: continue
            for f in fs:
                p=os.path.join(root,f); sz=os.path.getsize(p)
                if sz>10e6 or f.startswith('.'): continue
                files.append((p,sz)); tot+=sz
        if tot>30e6 or not files: break
        d=f'{R}/{j}'; os.makedirs(d+'/inputs/src')
        for p,_ in files:
            t=d+'/inputs/src/'+os.path.relpath(p,f'{R}/{s}'); os.makedirs(os.path.dirname(t),exist_ok=True); shutil.copy2(p,t)
        shutil.copy2(W+'/tasks/NIGHT_PREAMBLE.md',d+'/inputs/NIGHT_PREAMBLE.md')
        open(d+'/BRIEF.md','w').write(f'''# {j} — {('reproduktion' if fam == 'Q' else 'robusthet')} av {s}\nRead inputs/NIGHT_PREAMBLE.md. Arbeta bara i denna katalog; ≤ 60 s per run; 1 wire, < 1 GB.\nBuilds on:{s} (kopia i inputs/src/: RESULTS.md, results.json, code and data).\n{q}\nPREREG.md + PREREG.sha256 before calculation. Deliver results.json (source/deduction/UNKNOWN separated) and RESULTS.md starting with the line "{j}"; first paragraph: HOLDS / AVVIKER / UNKNOWN med tal.\n''')
        out.append(f'{prof[k%3]} swarm {j}'); k+=1
with open(W+'/tasks/lanes/bt_queue.txt','a') as f: f.write('\n'.join(out)+'\n')
print(len(out),"queued")
