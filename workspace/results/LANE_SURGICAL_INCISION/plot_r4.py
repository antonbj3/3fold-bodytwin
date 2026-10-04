import os
from pathlib import Path
P=Path(__file__).resolve().parent
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='2'
os.environ['MPLCONFIGDIR']=str(P/'r4/mpl_config')
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=json.loads((P/'SOURCE_TARGETS_R4_v2.json').read_text());F=json.loads((P/'r4/fracture_transfer.json').read_text());T=json.loads((P/'r4/total_force_transfer.json').read_text());B=D['Barnett2016']
fig,ax=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
gauges=['16','18','21','25'];xs=range(4);js=[B['first_minus_repeat_effective_J_means_J_m2'][g][1]/1000 for g in gauges]
for x,g,j in zip(xs,gauges,js):
 band=F['per_tool_plotted_bar_plus_reading_bands_J_m2'][g]
 ax[0].plot([x,x],[band[0]/1000,band[1]/1000],color='.70',lw=5)
ax[0].errorbar(list(xs),js,yerr=.05,fmt='o',color='#215ca0',capsize=4,label='Figure means ± reading bound')
ax[0].axhline(F['training_J_J_m2']/1000,ls='--',color='#bc5324',label='16G calibration, no refit')
ax[0].set(xticks=list(xs),xticklabels=gauges,xlabel='Gauge (bevel also changes)',ylabel='Effective paired-work J (kJ/m²)',title='A. New-tool transfer fails, 20 mm/s')
ax[0].legend(fontsize=8);ax[0].text(.02,.02,'Grey: reported plot bars + reading\nNot a confidence interval',transform=ax[0].transAxes,fontsize=8)
rows=T['heldout'];xs=list(range(len(rows)))
ax[1].errorbar(xs,[r['observed_peak_N'] for r in rows],yerr=.03,fmt='o',capsize=4,color='#215ca0',label='Observed peak ± reading')
ax[1].plot(xs,[r['models']['linear_D']['prediction_N'] for r in rows],'s--',color='#bc5324',label='16G one-anchor F ∝ D')
ax[1].set(xticks=xs,xticklabels=[str(r['gauge']) for r in rows],xlabel='Heldout gauge',ylabel='Total peak force (N)',title='B. One-gauge force prediction')
ax[1].legend(fontsize=8)
O=D['Owen2022']['experiments']['4'];keys=['21_new','21_reuse12','21_reuse36','21_reuse100'];xs=[1,12,36,100];ys=[O[k]['median_N'] for k in keys]
ax[2].plot(xs,ys,'o-',color='#215ca0',label='Published medians, new sites')
ax[2].axhline(ys[0],ls='--',color='#bc5324',label='Same diameter, no wear state')
ax[2].set(xlabel='Tool puncture count',ylabel='Piglet 21G peak force (N)',title='C. Diameter cannot encode wear',ylim=(.4,1.5))
ax[2].legend(fontsize=8);ax[2].text(.02,.02,'Owen2022, n=7 per group\nNot same-hole repeat experiment',transform=ax[2].transAxes,fontsize=8)
fig.suptitle('INCISION R4 — primary skin data; approximate figure readings; pending independent review',fontsize=11)
for ext in ['pdf','png']:
 file=P/f'r4/needle_series_overview.{ext}'
 if file.exists():raise FileExistsError(file)
 fig.savefig(file,dpi=180)
plt.close(fig)
