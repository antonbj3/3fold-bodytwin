"""Explicit default scenario; all numbers are assumptions unless stated otherwise."""
from .ports import synth, unknown

def scenario_config():
 def p(v,u,scope):return synth(v,u,'PREREG_R1.json / documented scenario',scope).dict()
 d={'tool':{'kind':'scalpel','edge_radius':p(5e-6,'m','Local edge radius; no calibrated radius→work law'),'speed':p(.02,'m/s','Tangential path speed; not used to extrapolate toughness')},
 'path':p([[0.,0.,0.],[.01,0.,0.]],'m','Straight planar incision path'),
 'cut_depth':p(.002,'m','Penetration depth, single geometric crack plane'),
 'layers':[{'name':'dermis','thickness':p(.0015,'m','Synthetic dermal thickness'),'gamma_cut':p(150.,'J/m²','Synthetic BINDINGS severance scenario; not needle J')},
           {'name':'subcutis','thickness':p(.0005,'m','Synthetic subcutis thickness'),'gamma_cut':p(150.,'J/m²','Explicit assumed subcutis cutting law; not borrowed empirical dermis')}],
 'parameters':{}}
 values={'gap':(.0002,'m'),'biological_width':(.00025,'m'),'perfusion_width':(.00025,'m'),
 'radius_bins':([8e-6,25e-6],'m'),'density_bins':([1.18e6,2e4],'1/m²'),
 'driving_pressure':(4000.,'Pa'),'blood_viscosity':(.0035,'Pa s'),'vascular_path':(.001,'m'),
 'domain':(.008,'m'),'grid_dx':(5e-6,'m'),'early_dt':(.1,'min'),'early_until':(120.,'min'),
 'late_dt':(.1,'day'),'late_until':(90.,'day'),'oxygen_storage':(30.,'mlO2/(m³ Torr)'),
 'oxygen_permeability':(1.3e-6,'mlO2/(m min Torr)'), 'oxygen_consumption':(1470.,'mlO2/(m³ min)'),
 'face_conductance':(.003,'mlO2/(m² min Torr)'), 'face_pressure':(160.,'Torr'),
 'seal_oxygen_coupling':(1.,'1'),'aggregation_factor':(1.,'1'),
 'k_deposit':(.1,'1/day'),'k_mature_legacy':(.04,'1/day'),
 'k_U_I':(.1,'1/day'),'k_I_M':(.05,'1/day'),'turnover_U_I':(.03,'1/day'),
 'birth_scale':(.1,'1/day'),'bridge_fraction':(1.,'1'),'load_angle':(0.,'degree')}
 for k,(v,u) in values.items():d['parameters'][k]=p(v,u,'Synthetic '+k+'; joint scenario, no identified native posterior')
 d['parameters']['oxygen_permeability']['source']='RESPONSE R1 Stücker viable-skin proxy; transfer to incision SYNTHETIC'
 d['parameters']['oxygen_consumption']['source']='RESPONSE R1 dermal consumption proxy; dry-volume/lineage transfer SYNTHETIC'
 d['review_state']='PENDING_INDEPENDENT_REVIEW';d['scientific_admission']=False
 return d

def evidence_config():
 d=scenario_config()
 for layer in d['layers']:layer['gamma_cut']=unknown('J/m²','INCISION/BINDINGS final ports','Native layer scalpel work not acquired').dict()
 for name in ['gap','biological_width','perfusion_width','radius_bins','density_bins','k_U_I','k_I_M','turnover_U_I','bridge_fraction']:
  u=d['parameters'][name]['unit'];d['parameters'][name]=unknown(u,'Final source ports','Native '+name+' remains unmeasured').dict()
 return d
