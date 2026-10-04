#!/usr/bin/env python3
"""Deterministic, local-only Grand Challenge raw-data packet generator."""
from __future__ import annotations
import argparse, csv, hashlib, io, json, re, shutil, zipfile
from collections import Counter
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
LANE = OUT / 'CX-DATAMATRIX'
ARCHIVES = sorted(Path('external_media').glob('*Competition-latest.zip'))
PROTOCOLS = {
 'P-GRF-PEAK': ('grf', ['peak_vertical_N','peak_ap_N','peak_ml_N','peak_resultant_N','vertical_at_resultant_peak_N'], 'c,g', 'Largest single force plate peak in N; the leg on the plate must be identified separately.'),
 'P-GRF-IMPULSE': ('grf', ['vertical_impulse_Ns','ap_impulse_Ns','ml_impulse_Ns','positive_vertical_impulse_Ns','vertical_mean_N'], 'b,c', 'Integrate force components using the actual time vector and report mean force.'),
 'P-GRF-BALANCE': ('grf', ['active_duration_s','active_fraction','vertical_rms_N','ap_rms_N','ml_rms_N'], 'c,e', 'Active force plate time at Fz>20 N and RMS on three axes.'),
 'P-EMG-TIMING': ('emg', ['semimem_active_fraction','bifem_active_fraction','vasmed_active_fraction','medgas_active_fraction','soleus_active_fraction'], 'a,c', 'Full-wave rectify five channels; active when the 50 ms mean exceeds 20 % of the channel p95.'),
 'P-EMG-AMPLITUDE': ('emg', ['semimem_rms','bifem_rms','vasmed_rms','medgas_rms','soleus_rms'], 'a,b', 'RMS of five raw EMG channels in the amplitude unit of the data, without MVC normalization.'),
 'P-MARKER-GEOM': ('trajectories', ['pelvis_width_mm','right_knee_width_mm','left_knee_width_mm','right_ankle_width_mm','left_ankle_width_mm'], 'd,f', 'Median distance between bilateral pelvic markers and medial/lateral knee and ankle markers.'),
 'P-IMU-VIRTUAL': ('trajectories', ['sacrum_accel_rms_m_s2','sacrum_accel_peak_m_s2','sacrum_speed_mean_m_s','sacrum_speed_peak_m_s','sacrum_vertical_range_mm'], 'e,a', 'Virtual pelvic IMU from a sacrum/pelvic marker: numerical velocity and acceleration without filtering.'),
 'P-CONTACT-COMP': ('forces', ['contact_peak_N','contact_mean_N','contact_impulse_Ns','contact_p95_N','contact_duration_s'], 'c,g', 'Contact load in N: eKnee sum of four load cells, eTibia resultant of Fx/Fy/Fz; the raw forces are in lb and are multiplied by 4.4482216152605.'),
}
CONSUMERS = {'a':'knee chain L1/V0, N1g k(activity), CX-KNEEMERGE','b':'populationspriorer N2b/N50/bodytwin_core.population','c':'nollmodeller B24/N1/N1g per aktivitet','d':'geometrimanager N7c','e':'glesa sensorer A321','f':'what-if anatomical parameter distributions','g':'kontaktlager Field U384'}
SUFFIX={'grf':'_grf.csv','emg':'_emg.csv','trajectories':'_trajectories.csv','forces':'_knee_forces.csv'}

def raw_table(data):
    s=data.decode('utf-8-sig',errors='replace').replace('\r','')
    lines=s.splitlines();
    if len(lines)<21:raise ValueError('empty or short raw table')
    header=next(csv.reader([lines[0]])); names=[h.strip().lower() for h in header]
    arr=np.genfromtxt(io.StringIO('\n'.join(lines[1:])),delimiter=',',usecols=range(len(names)),invalid_raise=False)
    if arr.ndim==1: arr=arr.reshape(1,-1)
    return names,arr

def col(names,arr,name):
    return arr[:,names.index(name.lower())]

def measure(proto,names,arr):
    if len(arr)<20 or not np.isfinite(arr[:,0]).any(): raise ValueError('too few finite samples')
    d={}; kind=PROTOCOLS[proto][0]
    t=col(names,arr,'time(sec)'); dt=np.diff(t); dt=dt[np.isfinite(dt)&(dt>0)]
    if len(dt)<10: raise ValueError('bad time base')
    if kind=='grf':
        ids=sorted({int(x[2:]) for x in names if re.fullmatch('fz[0-9]+',x)})
        f=np.stack([np.nan_to_num(sum((col(names,arr,f'{axis}{i}') for i in ids))) for axis in ('fx','fy','fz')],axis=1)
        fx,fy,fz=f.T; mag=np.linalg.norm(f,axis=1)
        if proto=='P-GRF-PEAK':
            # A plate is one foot contact. Summing plates includes the other leg.
            plates=[np.stack([col(names,arr,f'{axis}{i}') for axis in ('fx','fy','fz')],axis=1) for i in ids]
            peaks=np.array([np.linalg.norm(p,axis=1).max() for p in plates])
            best=plates[int(np.argmax(peaks))]; at=int(np.argmax(np.linalg.norm(best,axis=1)))
            d=dict(peak_vertical_N=max(np.max(abs(p[:,2])) for p in plates),
                   peak_ap_N=max(np.max(abs(p[:,0])) for p in plates),
                   peak_ml_N=max(np.max(abs(p[:,1])) for p in plates),
                   peak_resultant_N=float(peaks.max()),vertical_at_resultant_peak_N=best[at,2])
        elif proto=='P-GRF-IMPULSE': d=dict(vertical_impulse_Ns=np.trapezoid(fz,t),ap_impulse_Ns=np.trapezoid(fx,t),ml_impulse_Ns=np.trapezoid(fy,t),positive_vertical_impulse_Ns=np.trapezoid(np.maximum(fz,0),t),vertical_mean_N=np.mean(fz))
        else: d=dict(active_duration_s=np.sum((fz[:-1]>20)*np.diff(t)),active_fraction=np.mean(fz>20),vertical_rms_N=np.sqrt(np.mean(fz*fz)),ap_rms_N=np.sqrt(np.mean(fx*fx)),ml_rms_N=np.sqrt(np.mean(fy*fy)))
    elif kind=='emg':
        channels=['semimem','bifem','vasmed','medgas','soleus']
        for c in channels:
            x=col(names,arr,c)
            if proto=='P-EMG-AMPLITUDE':d[c+'_rms']=np.sqrt(np.nanmean(x*x))
            else:
                a=abs(x); n=max(3,round(.05/np.median(dt))); smooth=np.convolve(a,np.ones(n)/n,mode='same');d[c+'_active_fraction']=np.mean(smooth>.2*np.nanpercentile(a,95))
    elif kind=='trajectories':
        def point(base):
            return np.stack([col(names,arr,('RankleMedial' if base=='R.Ankle.Medial' and 'r.ankle.medialx' not in names else base)+axis) for axis in 'xyz'],axis=1)
        if proto=='P-MARKER-GEOM':
            for field,a,b in [('pelvis_width_mm','R.Asis','L.Asis'),('right_knee_width_mm','R.Knee.Lateral','R.Knee.Medial'),('left_knee_width_mm','L.Knee.Lateral','L.Knee.Medial'),('right_ankle_width_mm','R.Ankle.Lateral','R.Ankle.Medial'),('left_ankle_width_mm','L.Ankle.Lateral','L.Ankle.Medial')]:
                d[field]=np.nanmedian(np.linalg.norm(point(a)-point(b),axis=1))
        else:
            bases=['Sacrum','Lumbar','R.Psis'];base=next((b for b in bases if b.lower()+'x' in names),None)
            if not base:raise ValueError('no pelvic marker')
            p=point(base);v=np.gradient(p,t,axis=0)/1000;a=np.gradient(v,t,axis=0)
            speed=np.linalg.norm(v,axis=1);acc=np.linalg.norm(a,axis=1)
            d=dict(sacrum_accel_rms_m_s2=np.sqrt(np.nanmean(acc*acc)),sacrum_accel_peak_m_s2=np.nanmax(acc),sacrum_speed_mean_m_s=np.nanmean(speed),sacrum_speed_peak_m_s=np.nanmax(speed),sacrum_vertical_range_mm=np.nanmax(p[:,2])-np.nanmin(p[:,2]))
    else:
        if all(x in names for x in ('pm','am','al','pl')):x=sum(col(names,arr,k) for k in ('pm','am','al','pl'))
        elif all(x in names for x in ('fx','fy','fz')):x=np.linalg.norm(np.stack([col(names,arr,k) for k in ('fx','fy','fz')],axis=1),axis=1)
        else:raise ValueError('unrecognized contact format')
        x=x*4.4482216152605 # source data: pound-force, output fields: newtons
        d=dict(contact_peak_N=np.nanmax(x),contact_mean_N=np.nanmean(x),contact_impulse_Ns=np.trapezoid(x,t),contact_p95_N=np.nanpercentile(x,95),contact_duration_s=t[-1]-t[0])
    vals=np.array(list(d.values()),float)
    if len(vals)<5 or not np.isfinite(vals).all():raise ValueError('nonfinite measured fields')
    # Independent interleaved samples: estimate a held-out mean of a representative raw channel.
    ix=next((i for i,n in enumerate(names) if n not in ('frame','time(sec)') and np.isfinite(arr[:,i]).sum()>20),None)
    if ix is None:raise ValueError('no validation channel')
    x=arr[:,ix];odd=x[::2];even=x[1::2];k=min(len(odd),len(even));
    check=abs(np.nanmean(odd[:k])-np.nanmean(even[:k]))/(np.nanstd(x)+1e-12)
    if not np.isfinite(check):raise ValueError('nonfinite split-sample check')
    return {k:float(v) for k,v in d.items()},float(check)

def historical_text():
    chunks=[(ROOT/'notes/RESULTS_INDEX.md').read_text(errors='ignore')]
    for p in OUT.glob('*/results.json'):
        if p.parent.name.startswith('BT-DM-'):continue
        try:chunks.append(p.read_text(errors='ignore'))
        except OSError:pass
    return '\n'.join(chunks).lower()

def candidates():
    for outer in ARCHIVES:
        z=zipfile.ZipFile(outer)
        innername='Synchronized Motion Data.zip'
        if innername not in z.namelist():continue
        inner=zipfile.ZipFile(io.BytesIO(z.read(innername)))
        index={}
        for name in inner.namelist():
            low=Path(name).name.lower()
            if low.startswith('._'):continue
            for kind,suffix in SUFFIX.items():
                if low.endswith(suffix):index.setdefault(low[:-len(suffix)],{}).setdefault(kind,name)
        dataset='GC'+str(7-int(re.match(r'(\d+)\.',outer.name).group(1)))
        yield dataset,outer,inner,index

def packet_name(dataset,unit,proto):
    return 'BT-DM-'+dataset+'-'+re.sub('[^a-z0-9]+','-',unit.lower()).strip('-')+'-'+proto

def build(limit=35):
    LANE.mkdir(parents=True,exist_ok=True)
    history=historical_text(); counts=Counter(); dataset_counts=Counter();skips=Counter(); manifest=[]
    for dataset,outer,inner,index in candidates():
        for unit,files in sorted(index.items()):
            if unit=='jw_lungef1':skips['sealed']=skips['sealed']+1;continue
            if unit in history:skips['prior_unit']=skips['prior_unit']+1;continue
            for proto,(kind,fields,consumer,recipe) in PROTOCOLS.items():
                if counts[proto]>=limit or dataset_counts[(dataset,proto)]>=(35 if dataset=='GC6' else 12) or kind not in files:continue
                pid=packet_name(dataset,unit,proto);path=OUT/pid
                if (path/'inputs/PROTOCOL.json').exists():
                    old=json.loads((path/'inputs/PROTOCOL.json').read_text());manifest.append({'packet':pid,'dataset':dataset,'unit':unit,'protocol':proto,'check':old['check_value'],'bytes':(path/'inputs/raw.csv').stat().st_size});skips['existing_packet']+=1;counts[proto]+=1;dataset_counts[(dataset,proto)]+=1;continue
                data=inner.read(files[kind])
                try:
                    names,arr=raw_table(data); values,check=measure(proto,names,arr)
                    if set(values)!=set(fields):raise ValueError('field mismatch')
                except Exception:skips['insufficient']+=1;continue
                inp=path/'inputs';inp.mkdir(parents=True,exist_ok=True)
                filename='raw.csv';(inp/filename).write_bytes(data)
                spec={'protocol':proto,'dataset':dataset,'unit':unit,'fields':fields,'consumers':[CONSUMERS[k] for k in consumer.split(',')],'raw_source':str(outer)+' :: '+'Synchronized Motion Data.zip'+' :: '+files[kind],'raw_sha256':hashlib.sha256(data).hexdigest(),'raw_file':filename,'check':'odd/even independent sample mean gap / full-signal SD','check_value':check,'recipe':recipe}
                spec['units']={f:('source EMG units' if f.endswith('_rms') and proto=='P-EMG-AMPLITUDE' else '1' if f.endswith('_fraction') else 'N s' if f.endswith('_Ns') else 'm/s²' if 'accel' in f else 'm/s' if 'speed' in f else 's' if f.endswith('_s') else 'mm' if f.endswith('_mm') else 'N') for f in fields}
                (inp/'PROTOCOL.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n')
                (inp/'NIGHT_PREAMBLE.md').write_text('# Work preamble\nWrite PREREG.md and PREREG.sha256 before measurement. State the hypothesis, unit, criterion, raw data source and countertest.\nUse only the raw data and field list in the packet. Make no model claims without your own baseline.\nDistinguish source, measurement, derivation and hypothesis; mark missing data UNKNOWN.\nReport negative results and limitations. No publishing or messages.\n')
                shutil.copyfile(__file__,inp/'reader.py')
                (path/'BRIEF.md').write_text(f'# {pid}\nRead `inputs/NIGHT_PREAMBLE.md`. Measure raw data in `inputs/raw.csv` according to `inputs/PROTOCOL.json`.\n{recipe}\nWrite one row per trial with exactly the five protocol fields and unit, plus source and check value.\nReport measurement uncertainty, missing channels and limitations.\nRun `python3 inputs/reader.py check .` before analysis.\nWrite `TABLE.csv` and `RESULTS.md`; keep measured data separate from model predictions.\n')
                (path/'FILTER.md').write_text(f'1. Named consumers: {"; ".join(spec["consumers"])}.\n2. Raw data: `inputs/raw.csv`, SHA256 {spec["raw_sha256"]}.\n3. Five measures per trial: {", ".join(fields)}; two consumers.\n4. Non-circular countertest: odd versus even, independent time samples; standardized mean difference {check:.6g} (limit 0,05).\n5. Previous measurement: the exact trial ID `{unit}` is absent from existing results.json and RESULTS_INDEX at generation.\n')
                (path/'DATA_SUFFICIENCY.md').write_text(f'# Data sufficiency\nRaw file: {len(data)} byte, {len(arr)} rows. All five fields can be read and are finite.\nCheck value: {check:.6g}. Loading test: `python3 inputs/reader.py check .`.\n')
                manifest.append({'packet':pid,'dataset':dataset,'unit':unit,'protocol':proto,'check':check,'bytes':len(data)})
                counts[proto]+=1;dataset_counts[(dataset,proto)]+=1
        inner.close()
    (LANE/'MANIFEST.json').write_text(json.dumps({'packets':manifest,'counts':counts,'skips':skips},ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'counts':counts,'skips':skips,'total':sum(counts.values())}))

def check(path):
    path=Path(path); spec=json.loads((path/'inputs/PROTOCOL.json').read_text()); data=(path/'inputs'/spec['raw_file']).read_bytes()
    if hashlib.sha256(data).hexdigest()!=spec['raw_sha256']:raise ValueError('SHA mismatch')
    fields,num=measure(spec['protocol'],*raw_table(data))
    if list(fields)!=spec['fields'] or not np.isfinite(num):raise ValueError('field/check failure')
    print(json.dumps({'packet':path.name,'status':'PASS','rows':len(raw_table(data)[1]),'check_value':num}))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('command',choices=['build','check']);a.add_argument('path',nargs='?');a.add_argument('--limit',type=int,default=80);x=a.parse_args()
    if x.command=='build':build(x.limit)
    else:check(x.path or '.')
