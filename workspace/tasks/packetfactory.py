#!/usr/bin/env python3
"""Deterministic, bounded BodyTwin packet generator and validator."""
import argparse
import csv
import hashlib
import json
import math
import re
import shutil
from pathlib import Path

def _private_input_pattern(public_pattern):
    """Combine public exclusions with required operator-maintained private terms."""
    import json
    from pathlib import Path
    import re
    config = Path.home() / ".bodytwin" / "private_source_terms.json"
    terms = json.loads(config.read_text())
    if not isinstance(terms, list) or not terms or any(not isinstance(t, str) or not t.strip() for t in terms):
        raise ValueError("A nonempty private-source exclusion list is required")
    return "(?:" + public_pattern + ")|(?:" + "|".join(re.escape(t) for t in terms) + ")"


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
PREAMBLE = ROOT / 'tasks/NIGHT_PREAMBLE.md'
MAX_BYTES = 50_000_000

# Only metrics directly recoverable from the named frozen input are eligible.
R = [
 ('BT-B50','inputs/tibia_params.csv','validation.b17_input_row_count','rows'),
 ('BT-B137','inputs/BT-B110/tibia_params_v2.csv','reference.subject_count','subjects'),
 ('BT-B178','inputs/BT-B137/inputs/BT-B110/tibia_params_v2.csv','reference.finite_pooled_count','finite:tibial_torsion_abs_deg'),
 ('BT-B126','inputs/BT-B66/params.csv','scope.n_input_rows','rows'),
 ('BT-B48','inputs/B25_landmarks.csv','cohort.side_counts.L','side:L'),
]
P = [
 ('BT-B92','inputs/B42_params.csv','ccd_deg'),
 ('BT-B121','inputs/BT-B42/params.csv','ccd_deg'),
 ('BT-B174','inputs/BT-B126/inputs/BT-B66/params.csv','ccd_deg'),
 ('BT-B64','inputs/B25_landmarks.csv','head_radius_mm'),
 ('BT-B142','inputs/BT-B85/landmarks.csv','x_mm'),
]
X = [
 ('BT-B50','inputs/tibia_params.csv','BT-B137','inputs/BT-B110/tibia_params_v2.csv','tibial_length_mm'),
 ('BT-B50','inputs/tibia_params.csv','BT-B178','inputs/BT-B137/inputs/BT-B110/tibia_params_v2.csv','tibial_torsion_abs_deg'),
 ('BT-B92','inputs/B42_params.csv','BT-B121','inputs/BT-B42/params.csv','ccd_deg'),
 ('BT-B92','inputs/B66_params.csv','BT-B174','inputs/BT-B126/inputs/BT-B66/params.csv','ccd_deg'),
 ('BT-B48','inputs/B25_landmarks.csv','BT-B64','inputs/B25_landmarks.csv','head_radius_mm'),
]

def field(data, path):
    for key in path.split('.'):
        data = data[key]
    return data

def rows(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def metric(data, kind):
    if kind == 'rows': return len(data)
    if kind == 'subjects': return len({r['subject'] for r in data})
    if kind.startswith('side:'): return sum(r.get('side') == kind[5:] for r in data)
    if kind.startswith('finite:'):
        key = kind[7:]
        return sum(math.isfinite(float(r[key])) for r in data if r.get(key))
    raise ValueError(kind)

def safe_source(src, paths):
    d = RESULTS / src
    if not (d/'RESULTS.md').is_file() or not (d/'results.json').is_file(): return False
    for rel in paths:
        p = (d/rel).resolve()
        if not p.is_relative_to(d.resolve()) or not p.is_file() or p.is_symlink(): return False
    # No credentials/private patient records; allow only the explicitly selected small tables.
    if any(re.search(_private_input_pattern('credential|secret|token|password|patient|collaborator|health'), rel, re.I) for rel in paths): return False
    inp=d/'inputs'
    if inp.exists() and any(re.search(r'credential|secret|token|password|patient|collaborator|\.env$|\.dcm$',str(p.relative_to(inp)),re.I) for p in inp.rglob('*') if p.is_file()):return False
    return True

def used_audits():
    out=set()
    for p in RESULTS.glob('BT-AN-*/inputs/REQUIREMENT.json'):
        try:
            value=json.loads(p.read_text()).get('audit_rows',[])
            out.update(value if isinstance(value,list) else [value])
        except (ValueError,OSError): pass
    for p in RESULTS.glob('BT-AG-*/inputs/REQUIREMENT.json'):
        try:
            value=json.loads(p.read_text()).get('audit_rows',[])
            out.update(value if isinstance(value,list) else [value])
        except (ValueError,OSError): pass
    return out

def audit_candidates():
    used=used_audits(); out=[]
    for line in (ROOT/'notes/RESULTS_INDEX.md').read_text().splitlines():
        m=re.match(r'\| A(\d+) \|',line)
        if not m or 'ej granskad' not in line or int(m[1]) in used: continue
        ids=re.findall(r'results/([^/| ]+)/',line)
        if not ids: continue
        src=ids[0]; d=RESULTS/src
        needed=['RESULTS.md','results.json','PREREG.md','PREREG.sha256']
        if not safe_source(src,needed):continue
        extra=[p for p in (d/'inputs').rglob('*') if p.is_file()] if (d/'inputs').exists() else []
        if not extra or any(p.is_symlink() or re.search(r'credential|secret|token|password|patient|collaborator|health|\.h5$|\.dcm$',str(p.relative_to(d)),re.I) for p in extra):continue
        if sum(p.stat().st_size for p in extra)+sum((d/n).stat().st_size for n in needed)>MAX_BYTES//3:continue
        try:
            actual=hashlib.sha256((d/'PREREG.md').read_bytes()).hexdigest()
            if actual != (d/'PREREG.sha256').read_text().split()[0]:continue
            data=json.loads((d/'results.json').read_text())
            def numeric_count(value):
                if isinstance(value,dict):return sum(map(numeric_count,value.values()))
                if isinstance(value,list):return sum(map(numeric_count,value))
                return int(isinstance(value,(float,int)) and not isinstance(value,bool) and math.isfinite(value))
            if numeric_count(data)<3:continue
        except (ValueError,OSError,IndexError):continue
        out.append((int(m[1]),src,line))
    return out

def specs():
    for src,rel,f,k in R:
        if safe_source(src,[rel]) and field(json.loads((RESULTS/src/'results.json').read_text()),f)==metric(rows(RESULTS/src/rel),k):
            yield ('R',f'BT-R-{src}',[(src,rel),(src,'results.json')],dict(field=f,metric=k,source=src))
    for src,rel,col in P:
        if safe_source(src,[rel]):
            rr=rows(RESULTS/src/rel)
            if len(rr)>=5 and col in rr[0] and ('subject' in rr[0] or 'mesh_id' in rr[0]):
                yield ('P',f'BT-P-{src}',[(src,rel),(src,'results.json')],dict(column=col,source=src,choice='leave_one_subject_out'))
    for a,pa,b,pb,col in X:
        if safe_source(a,[pa]) and safe_source(b,[pb]):
            ar,br=rows(RESULTS/a/pa),rows(RESULTS/b/pb)
            if ar and br and col in ar[0] and col in br[0]:
                yield ('X',f'BT-X-{a}-{b}',[(a,pa),(a,'results.json'),(b,pb),(b,'results.json')],dict(column=col,sources=[a,b]))
    aud=audit_candidates()
    for i in range(0,len(aud)-2,3):
        group=aud[i:i+3]
        ident=f'BT-AG-{group[0][0]}-{group[-1][0]}'
        files=[]
        for _,s,_ in group:
            files.extend((s,n) for n in ('results.json','RESULTS.md','PREREG.md','PREREG.sha256'))
            files.extend((s,str(p.relative_to(RESULTS/s))) for p in (RESULTS/s/'inputs').rglob('*') if p.is_file())
        yield ('AG',ident,files,dict(audit_rows=[n for n,_,_ in group],sources=[s for _,s,_ in group],register_lines=[l for _,_,l in group]))

def generate(spec):
    family,ident,files,req=spec; d=RESULTS/ident
    if d.exists(): return False
    inp=d/'inputs'; inp.mkdir(parents=True)
    try:
        shutil.copy2(PREAMBLE,inp/'NIGHT_PREAMBLE.md')
        for src,rel in files:
            target=inp/src/rel
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(RESULTS/src/rel,target)
        req.update(family=family,id=ident,files=[f'{src}/{rel}' for src,rel in files])
        (inp/'REQUIREMENT.json').write_text(json.dumps(req,ensure_ascii=False,indent=2)+'\n')
        if family=='AG':
            (inp/'RESULTS_INDEX_ROWS.md').write_text('\n'.join(req['register_lines'])+'\n')
            question=f"Granska oberoende registerraderna {req['audit_rows']} against copied source files."
            criterion="3/3 rows: recalculate at least one numerical claim per row; deviation at most 1 % relative. Otherwise UNKNOWN."
            counter="Also compare against the source's own conclusion without using it as independent ground truth."
        elif family=='R':
            question=f"Recalculate {req['field']} from the raw table with its own code."
            criterion="Relative deviation at most 1 % against the frozen result field; flag larger deviations."
            counter="Recalculate after omitting one row; the number should then respond according to the metric definition."
        elif family=='P':
            question=f"Is the median of {req['column']} stable when a person is left out?"
            criterion="All personwise omissions change the median by at most 20 % relative to the full table."
            counter="Compare with a 20 % higher threshold and report the largest change and denominator."
        else:
            question=f"Transfer the same median analysis of {req['column']} between the two frozen tables."
            criterion="Both tables have at least 5 finite values; report median ratio and overlap without assuming a common cohort."
            counter="Compare min/max and person counts; report UNKNOWN for external generalization."
        refs=', '.join('inputs/'+x for x in req['files'])
        (inp/'DATA_SUFFICIENCY.md').write_text(f'# Data adequacy — {ident}\nLadda {refs}. Require finite values and fields indicated in: inputs/REQUIREMENT.json. Run the python3 tasks/packetfactory.py check {ident}`. Saknat data ger UNKNOWN.\n')
        (d/'BRIEF.md').write_text(f"# {ident}\nRead inputs/NIGHT_PREAMBLE.md. Endast paketet, 1 thread.\nBuilds on:{', '.join(dict.fromkeys((src for src, _ in files)))}.\nQuestion: {question}\nKriterium: {criterion}\nMotprov: {counter}\nUnderlag: {refs}, inputs/REQUIREMENT.json, inputs/DATA_SUFFICIENCY.md.\nSkriv PREREG.md and PREREG.sha256 before analysis. Separate source, derivation and hypothesis. Deliver results.json and RESULTS.md with first row`# {ident}`.\n")
        errors=check(ident)
        if errors: raise ValueError('; '.join(errors))
        return True
    except Exception:
        shutil.rmtree(d)
        raise

def check(ident):
    d=RESULTS/ident; errors=[]
    try:
        req=json.loads((d/'inputs/REQUIREMENT.json').read_text())
        brief=(d/'BRIEF.md').read_text()
        if not brief.startswith('# '+ident+'\n'):errors.append('brief-ID')
        if "Builds on:" not in brief or 'PREREG' not in brief:errors.append('brief-kontrakt')
        if not (d/'inputs/NIGHT_PREAMBLE.md').is_file() or not (d/'inputs/DATA_SUFFICIENCY.md').is_file():errors.append('preamble/dataprov')
        if any(x.is_symlink() for x in d.rglob('*')):errors.append('symlink')
        for rel in req['files']:
            p=(d/'inputs'/rel).resolve()
            if not p.is_relative_to((d/'inputs').resolve()) or not p.is_file():errors.append('saknad '+rel)
        if req['family']=='R':
            src=req['source']; rel=req['files'][0].removeprefix(src+'/')
            observed=metric(rows(d/'inputs'/src/rel),req['metric'])
            claimed=field(json.loads((d/'inputs'/src/'results.json').read_text()),req['field'])
            if abs(observed-claimed)>max(1e-9,abs(claimed)*.01):errors.append('R-avvikelse')
        elif req['family'] in ('P','X'):
            for rel in (x for x in req['files'] if x.endswith('.csv')):
                rr=rows(d/'inputs'/rel); col=req['column']
                vals=[float(r[col]) for r in rr if r.get(col) and math.isfinite(float(r[col]))]
                if len(vals)<5:errors.append("insufficient numeric data "+rel)
        elif req['family']=='AG':
            if len(req['audit_rows'])!=3 or len(set(req['audit_rows']))!=3:errors.append("audit-coverage")
            for src in req['sources']:
                sd=d/'inputs'/src
                if hashlib.sha256((sd/'PREREG.md').read_bytes()).hexdigest() != (sd/'PREREG.sha256').read_text().split()[0]:errors.append('audit-hash '+src)
                data=json.loads((sd/'results.json').read_text())
                def nums(value):
                    if isinstance(value,dict):
                        for v in value.values():yield from nums(v)
                    elif isinstance(value,list):
                        for v in value:yield from nums(v)
                    elif isinstance(value,(float,int)) and not isinstance(value,bool):yield value
                if sum(math.isfinite(x) for x in nums(data))<3:errors.append('audit-numerik '+src)
            if len((d/'inputs/RESULTS_INDEX_ROWS.md').read_text().splitlines())!=3:errors.append('registerrader')
        else:errors.append('familj')
        size=sum(p.stat().st_size for p in d.rglob('*') if p.is_file())
        if size>MAX_BYTES:errors.append("over 50 MB")
    except Exception as e:errors.append(f'{type(e).__name__}: {e}')
    return errors

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('command',choices=['inventory','generate','check']); ap.add_argument('id',nargs='?');ap.add_argument('--limit',type=int,default=150);ap.add_argument('--family',choices=['AG','R','P','X']);a=ap.parse_args()
    if a.command=='check':
        ids=[a.id] if a.id else [p.name for p in RESULTS.iterdir() if p.is_dir() and re.match(r'BT-(R|P|X|AG)-',p.name)]
        out={i:check(i) for i in ids};print(json.dumps(out,ensure_ascii=False,indent=2));raise SystemExit(bool(any(out.values())))
    pool=list(specs())
    priority={'AG':0,'R':1,'P':2,'X':3}
    pool.sort(key=lambda s:(priority[s[0]], -(RESULTS/s[3].get('source','BT-B50')/'results.json').stat().st_mtime if s[0]=='R' else 0, s[1]))
    counts={f:sum(s[0]==f for s in pool) for f in ('AG','R','P','X')}
    if a.command=='inventory':print(json.dumps(dict(available=counts,total=sum(counts.values())),indent=2));return
    selected=[s for s in pool if (not a.family or s[0]==a.family) and not (RESULTS/s[1]).exists()][:a.limit]
    made=[]
    for s in selected:
        if generate(s):made.append(s[1])
    print(json.dumps(dict(generated=made,available=counts),ensure_ascii=False))

if __name__=='__main__':main()
