import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
R=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(R/'.mplconfig')
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();dest=a.out.resolve()
if not dest.is_relative_to(R):raise ValueError('Lane-local figures only')
for ext in ['png','pdf']:
 if dest.with_suffix('.'+ext).exists():raise FileExistsError(dest)
nom=np.load(R/'r1/nominal_complete_v4/fields.npz');abl=np.load(R/'r1/no_seal_oxygen_v1/fields.npz')
fig,ax=plt.subplots(1,3,figsize=(14,4.2));col=['#126b87','#cc813a']
ax[0].plot(nom['blood_t_min'],nom['seal'],color=col[0],label='Persistent seal')
ax[0].plot(nom['blood_t_min'],nom['blood_m3']*1e9/120,color=col[1],label='Blood /120 µL')
ax[0].set(xlabel='Minutes',ylabel='Fraction / scaled volume',title='Radius + shear + seal memory',xlim=(0,120),ylim=(0,1.1));ax[0].legend(fontsize=8)
for data,label,c in [(nom,'Seal → atmosphere',col[0]),(abl,'Atmospheric coupling off',col[1])]:
 ax[1].plot(data['x']*1e3,data['early_pressure'],label=label,color=c)
 ax[2].plot(data['days'],data['strength_pct'],label=label,color=c)
ax[1].axhline(10,color='gray',ls=':',label='Synthetic threshold')
ax[1].set(xlabel='Distance from incision (mm)',ylabel='O2 pressure (Torr)',title='120 min; FV oxygen inventory',xlim=(0,1));ax[1].legend(fontsize=8)
rows=json.loads((R/'sources/LANE_SURGICAL_RESPONSE/r2/world/summary.json').read_text())['curves'][0]['rows']
ax[2].scatter([x['day'] for x in rows],[x['breaking_pct'] for x in rows],marker='x',s=22,c='#a43d4c',label='Rat source (unmatched assay)')
ax[2].set(xlabel='Days',ylabel='Synthetic % intact strength',title='Inventory proxy; native rupture unknown',xlim=(0,90));ax[2].legend(fontsize=7)
for x in ax:x.spines[['top','right']].set_visible(False);x.grid(alpha=.16)
fig.suptitle('Surgical synthesis R1 · conditional closures · PENDING_INDEPENDENT_REVIEW',fontsize=13,y=1.02)
fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(dest.with_suffix('.'+ext),dpi=160,bbox_inches='tight')
print(dest)
