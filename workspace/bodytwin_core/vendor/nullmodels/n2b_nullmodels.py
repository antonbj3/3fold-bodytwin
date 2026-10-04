"""N2b: obligatory null models for the collaborator's box lift (no in-vivo facit exists for the collaborator).
BT-B24: BW x OrthoLoad median peak %BW, activity 'Lifting' (results/BT-B24/null_models.json); BW = sum of N14b segment masses.
N1 (N12): F = k |GRF| per leg, k = 2.2 and 2.4 (N12 README: fitted LOSO on Grand Challenge knee, gait); |GRF| per foot = |sum of
GRF_Prediction contact forces| from the collaborator's h5 (Fout of the contact actuators times their direction is not needed: we use the
magnitude of GRF_Prediction_<side>/FGlobal).  -> n2b_nullmodels.json"""
import json, numpy as np, h5py
an=json.load(open('n2b_anthro.json')); BW=an['body_mass_kg_sum_segments']*9.81
d=json.load(open('../BT-B24/null_models.json')); it=d if isinstance(d,list) else next(v for v in d.values() if isinstance(v,list))
b24={v['facet_id']:dict(pctBW=v['population_null']['coefficient_pctBW'],peak_N=v['population_null']['coefficient_pctBW']/100*BW,
      loso_rmse_N=v['error']['rmse_N']) for v in it if 'Lifting' in v['facet_id']}
f=h5py.File('external_mount','r')
grf={}
for side in ('Right','Left'):
    base='Output/_Main/EnvironmentModel/ForcePlates/GRF_Prediction_%s/Contacts'%side
    tot=f['Output/_Main/EnvironmentModel/ForcePlates/GRF_Prediction_%s/FGlobal'%side][:].reshape(-1,3)
    grf[side]=np.linalg.norm(tot,axis=1)
out=dict(BW_N=BW,BT_B24_lifting=b24,
         GRF_peak_N={s:float(v.max()) for s,v in grf.items()},
         N1_peak_N={s:{'k2.2':float(2.2*v.max()),'k2.4':float(2.4*v.max())} for s,v in grf.items()},
         note='knee: BT-B24 has no Lifting facet; N1 k from knee gait (N12) applied to hip/knee as a scale baseline only')
json.dump(out,open('n2b_nullmodels.json','w'),indent=1); print(json.dumps(out,indent=1))
