import os
from pathlib import Path
ROOT=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(ROOT/'r3/mplconfig')
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from diagnose_r2 import load_reference,targets
for suffix in ['pdf','png']:
 if (ROOT/f'SHARED_INVENTORY_R3.{suffix}').exists():raise FileExistsError(suffix)
with np.load(ROOT/'r3/experiment_v2/all_R3/fields.npz') as f:t=f['days'];strength=f['strength_pct']
with np.load(ROOT/'r3/experiment_v2/all_R3/shared_reference.npz') as f:C=f['C'];UIM=f['UIM'];chem=f['chemical_strength_pct'];ch=f['chemical_marks']
ref=load_reference()
fig,axes=plt.subplots(1,3,figsize=(13.5,4),layout='constrained')
a=axes[0];a.plot(t,C,label='FV C',linewidth=3);a.plot(t,UIM.sum(1),'--',label='sum U/I/M');a.plot(t,UIM[:,2],label='mature M');a.set(xlabel='Day',ylabel='Collagen per fixed reference',title='Shared inventory: trace error <6e-15');a.legend(fontsize=8)
a=axes[1];a.plot(t,strength,label='Autonomous U/I/M (RMSE 21.81)');a.plot(t,chem,label='Autonomous chemistry (44.34)');a.plot(ref['days'],ref['strength_pct'],':',label='External HP diagnostic (8.31)')
rows=targets();a.scatter([q['day'] for q in rows],[q['target_pct'] for q in rows],c='black',s=18,label='Separate rat source targets')
a.set(xlabel='Day',ylabel='Strength proxy (% intact)',title='Mass consistency does not identify maturity');a.legend(fontsize=7)
a=axes[2];hp=np.divide(ch[:,1],C,out=np.zeros_like(C),where=C>0);a.plot(t,hp,label='Shared chemistry H/C');a.axhline(.01445,ls='--',c='orange',label='Conditional site ceiling .01445');a.plot(ref['days'],ref['HP'],':',label='External mouse HP input');a.set(xlabel='Day',ylabel='HP mol/mol collagen',title='Fixed sites/yield cannot reach HP .043');a.legend(fontsize=8)
fig.suptitle('R3 synthetic closures — pending independent review; native prediction UNKNOWN',fontsize=11)
for suffix in ['pdf','png']:fig.savefig(ROOT/f'SHARED_INVENTORY_R3.{suffix}',dpi=160)
