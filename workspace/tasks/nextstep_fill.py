"N-package: performs a result EGET named next step (ur dess RESULTS.md) with source files as input. reserve_worker C.\nValue controlled (the result itself points out the step), deterministic. Skipping meta-sources and already lined up."
import os, re, shutil, sys
W=''; R=W+'/results'; N=int(sys.argv[1]) if len(sys.argv)>1 else 100
queued={l.split()[2] for l in open(W+'/tasks/lanes/bt_queue.txt') if len(l.split())==3}
meta=re.compile(r'^(CX-(BOOKKEEP|SWARMGEN|PLANNER|CLOUDLANES|PACKETFACTORY|PARAM|FIELDSHARE)|BT-(AN-G|Q-|S-|N-|R-|P-|X-|AG-)|FX-)')
pat=re.compile("(?is)(?:next steg|next step)[^:\\n]*[:：]?\\s*(.{40,600}?)(?:\\n\\n|\\Z)")
out=[]; k=0
srcs=sorted([d for d in os.listdir(R) if not meta.match(d) and os.path.isfile(f'{R}/{d}/RESULTS.md')], key=lambda d:-os.path.getmtime(f'{R}/{d}/RESULTS.md'))
for s in srcs:
    if len(out)>=N: break
    j=f'BT-N-{s}'
    if j in queued or os.path.exists(f'{R}/{j}'): continue
    m=pat.search(open(f'{R}/{s}/RESULTS.md',errors='ignore').read())
    if not m: continue
    step=' '.join(m.group(1).split())
    files=[]; tot=0
    for root,_,fs in os.walk(f'{R}/{s}'):
        if '__pycache__' in root or '/.' in root: continue
        for f in fs:
            p=os.path.join(root,f); sz=os.path.getsize(p)
            if sz<=10e6 and not f.startswith('.'): files.append(p); tot+=sz
    if tot>40e6 or not files: continue
    d=f'{R}/{j}'; os.makedirs(d+'/inputs/src')
    for p in files:
        t=d+'/inputs/src/'+os.path.relpath(p,f'{R}/{s}'); os.makedirs(os.path.dirname(t),exist_ok=True); shutil.copy2(p,t)
    shutil.copy2(W+'/tasks/NIGHT_PREAMBLE.md',d+'/inputs/NIGHT_PREAMBLE.md')
    open(d+'/BRIEF.md','w').write(f'''# {j} — performs the next step as {s} yourself pointed out\nRead inputs/NIGHT_PREAMBLE.md (§0 riktning, §2 nollmodeller). Arbeta bara i denna katalog; ≤ 60 s per run; 1–2 threads, < 1 GB.\nBuilds on:{s} (kopia i inputs/src/: RESULTS.md, results.json, kod, data).\nNext step according to the source: "{step}"\nTask: implements this step as far as possible with the files in inputs/src/. If the step requires data that is not available here: do the largest delimited part that goes, and write exactly which data is missing (UNKNOWN) — Don't make it up.\nPREREG.md + PREREG.sha256 before calculation (hypotes, kriterium med tal, motprov)Deliver results.json and RESULTS.md starting with the line "{j}"; first paragraph: what the step showed, with numbers.\n''')
    out.append(f'C reserve_worker {j}')
with open(W+'/tasks/lanes/bt_queue.txt','a') as f: f.write('\n'.join(out)+'\n')
print(len(out),"N-package queued (reserve_worker C)")
