"""All inherited defaults remain visible, with domain transfers labelled synthetic."""
from dataclasses import asdict
from .ports import synth
from .hemostasis import CO,PL
from .vendor import response_r1 as r

def manifest():
 q={k:synth(v['value'],v['unit'],'Frozen Q036/model.py parameter_table; '+v['source'],v['rationale']+'; cardiac ischemia proxy transferred to wound: SYNTETISKT').dict() for k,v in r.q036.parameter_table(r.QP).items()}
 c={}
 for k,v in dict(CO['LOCKED'],TF=1000.,XIIa=0.,S_amp=5.).items():
  u='nM' if k.startswith('Km_') or k in ['Ki_tfpi','TF','XIIa'] else '1' if k.startswith('S_') else '1/min' if k in ['k6_base','k9_base','k10','kdiss_VIIIa','kclear_XIa'] else '1/(nM min)'
  c[k]=synth(v,u,'Frozen coagulation_hemostasis.py rhs','Illustrative/tuned 15-species closure, not original full34-species kinetics').dict()
 for k in ['kAT_IIa','kAT_Xa','kAT_IXa']:c[k]=synth(CO[k],'1/(nM min)','Frozen anticoagulation parameter conversion','Source-proxy rates; wound-transfer assumption').dict()
 return {'Q036_parameters':q,'Q036_state_units':r.q036.STATE_UNITS,'coagulation_rates':c,
 'coagulation_initial_nM':{k:synth(v,'nM','Frozen coagulation y0_vec','Specified uniform-concentration reservoir; no whole-blood advection').dict() for k,v in zip(CO['IDX'],CO['y0_vec']())},
 'platelet_parameters':{k:synth(PL[k],'1/min' if k.startswith('K') else '1','Frozen platelet_hemostasis.py','Illustrative adhesion/aggregation/coverage closure').dict() for k in ['K_ON0','K_AGG0','KOFF_A','KOFF_P','W_A','W_P','THETA_C']},
 'FV_default_parameters':asdict(r.Params()),'FV_scope':'All unoverridden FV Params are synthetic/proxy; units in frozen response_r1.py. pcap60Torr,healthy51Torr,km3Torr,flowtransition15um,Dfibro2e-8m²/day; long_rhs literal coefficients have normalized-state/day units except .3 and oxygen/gap gates dimensionless.',
 'long_rhs_laws':'Ndot=2.4 pulse footprint+1.5dead exp(-t/4)-1.8N; Mdot=.8N-.35M; Fdot=.45M h(1-F)/(1+gap/.0005)-.06F; Cdot=kdep Fh(1-C)-.008MC; Xdot=kmature Ch(1-X); Vdot=(.20Fh+.005)(1-V)-.025MV',
 'additional_closures':'Two-pathway shear sigmoid low600/900 and high6000/20000 s^-1 with0.05/0.95 proxy anchors; radius_floor.05; Poiseuille; shape0.3/0.7plateletcoverage; irreversible seal and atmospheric permeability link; 10Torrhypoxia threshold;75% isotropic reference and affine e^TQe. All synthetic transfer, no measured joint posterior.'}
