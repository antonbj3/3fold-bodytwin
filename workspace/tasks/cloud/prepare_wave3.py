#!/usr/bin/env python3
"""Prepare reserve public-data/own-code lanes for paced launches after wave 2."""
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parent/'lanes'
ITEMS={
'CLOUD-W3-NHANES-NONLINEAR':('NHANES adjusted nonlinear DXA response','Use bundled public CDC XPT files. Compare preregistered linear adjusted baseline with splines for age, leg lean and height in 1999–2000 training, 2001–02 holdout. Tune only within training. Quantify calibration and sex residuals.','PASS if adjusted+DXA nonlinear fit improves holdout RMSE >=5% against tuned adjusted no-DXA model with paired PSU-bootstrap lower bound >0.'),
'CLOUD-W3-NHANES-CALIB':('NHANES sex and age calibration of strength models','Use bundled CDC XPT files. Audit held-out residuals by sex and age cells for adjusted with/without DXA, avoiding changing cohort or target. Build prespecified calibration slopes and uncertainty.','PASS if adding DXA reduces absolute mean error in every sex-by-age cell and overall RMSE; else FAIL.'),
'CLOUD-W3-KKT-DEGEN':('Degenerate active-set derivative certificates','Original synthetic parametric QP only. Analyze strict complementarity failure, nearly collinear columns, changing peak frame, and numerical tolerance. Compare one-sided directional derivatives and naive frozen-active KKT.','PASS if all seeded degenerate cases are flagged or bounded and interior derivative relative error <0.1%; no claim about A369 actual matrices.'),
'CLOUD-W3-BATCH-CERT':('Batched NNLS stopping certificates','Original synthetic NNLS only. Compare stopping rules based on projected gradient, duality gap and objective change on ill-conditioned ensembles. Independent SciPy reference and timing.','PASS if a proposed stopping rule attains <=1e-6 scaled KKT and <=1e-7 objective gap in >=99% of >=100 cases.'),
'CLOUD-W3-PATELLA-CONTACT':('Patella contact switching and moment arm','Original planar pulley geometry. Derive continuous tendon path and moment arm across contact onset/exit, test analytic versus finite-difference virtual work and define angle convention.','PASS if energy continuity and derivative convergence hold away from transition and all transitions are detected; public measured peak remains UNKNOWN absent data.'),
'CLOUD-W3-DISC-IDENT':('Disc pressure parameter identifiability','Original postural equilibrium model with explicit area, muscle arm and co-contraction. Use synthetic postures to determine which combinations pressure-only observations identify; add one independent imaging measure.','PASS mathematical test if profile likelihood and Fisher rank agree on all prespecified designs; empirical pressure validation UNKNOWN absent same-level public measurements.'),
'CLOUD-W3-OED-CORREL':('Correlated measurement error in minimal protocol','Build synthetic design with systematic landmark and imaging biases. Derive posterior covariance and show how naive independence changes selected measurement; calibrate held-out intervals.','PASS if covariance-aware 90% intervals cover 85–95% and naive intervals fail in high-correlation case across >=1000 trials.'),
'CLOUD-W3-BONE-HEAL':('Fatigue accumulation with healing counterexample','Derive and simulate Miner versus damage-healing ODE under repeated bouts. Use source-backed ranges if available; report unknown if unavailable. Compare sequence dependence and ranking.','PASS mathematical test if unit checks and analytic limiting cases agree <1e-5 and a ranking reversal is reproducible; empirical injury risk UNKNOWN.'),
'CLOUD-W3-CELL-MESH':('Cell geometry volume-to-field reduction','Use a public segmented volume if verified, else synthetic labeled geometry. Compare voxel, tetrahedral or finite-volume geometric quantities, boundary flux and 1-voxel uncertainty.','PASS public-data criterion only if version, labels and voxel scale verified and volume discrepancy <=1%; otherwise synthetic PARTIAL.'),
'CLOUD-W3-CELL-COUPLE':('Cell metabolism and diffusion coupling','Construct a minimal original coupled ATP-consumption and diffusion toy cell with organelle geometry; test conservation and nondimensional sensitivity. Seek matched public structure/function reference; distinguish unmatched species/cell lines.','PASS numerical criterion if coupled mass/energy balance <1e-6 and mesh convergence <2%; empirical calibration UNKNOWN without matched data.'),
'CLOUD-W3-VASCULAR-SAMPLE':('Minimal sampling for vascular transient','Original two-compartment synthetic model with volume and solute conservation. Optimize plasma-only versus added interstitial time points under known noise; compute held-out parameter coverage.','PASS if added interstitial points raise Fisher rank and 90% interval coverage is 85–95% in >=1000 simulations; real-world optimal schedule UNKNOWN.'),
'CLOUD-W3-NATO-ERROR':('NATO chain error-budget certification','Build original activity-load-strain-fatigue chain with interval arithmetic or Monte Carlo, clearly separate source-backed from assumed ranges. Compare analytic and sampled bounds, find dominant link and reversal conditions.','PASS mathematical criterion if dimensions and limiting cases pass and interval bound contains all >=100000 seeded samples; public risk mapping UNKNOWN.'),
}
COMMON='Use only public primary data and original code. No restricted model data, TLEM, the collaborator data, LHDL, private matrices, other repositories, paid extra usage, subagents, push, email or publishing. Cite sources actually checked. Distinguish synthetic proof from empirical validity. Never invent missing data. Deliver RESULTS.md, results.json and runnable code; final message print each small file verbatim between =====FILE <path>===== and =====END FILE=====, total <=200 KB.'
for lane,(title,work,criterion) in ITEMS.items():
 d=ROOT/lane;d.mkdir(parents=True,exist_ok=True)
 source_note = 'Read SOURCE_NOTES.md for one checked public primary reference.\n\n' if lane in ('CLOUD-W3-DISC-IDENT','CLOUD-W3-BONE-HEAL','CLOUD-W3-PATELLA-CONTACT','CLOUD-W3-CELL-MESH') else ''
 brief=f'# {lane}: {title}\n\n{COMMON}\n\n{source_note}## Bounded task\n{work}\n\nReport baseline, uncertainty, counterexamples, exact commands, source/version and limitations.\n'
 prereg=f'# PREREG {lane}\n\nQuestion: {title}.\n\nFrozen criterion: {criterion}\n\nIf required data are unavailable, report UNKNOWN/PARTIAL.\n'
 for name,content in [('BRIEF.md',brief),('PREREG.md',prereg)]:
  p=d/name
  if p.exists() and p.read_text()!=content: raise RuntimeError(str(p))
  p.write_text(content)
 (d/'PREREG.sha256').write_text(hashlib.sha256(prereg.encode()).hexdigest()+'  PREREG.md\n')
 if lane.startswith('CLOUD-W3-NHANES'):
  import shutil
  for name in ('nhanes_public.zip','MANIFEST.sha256','a359.py'):
   p=d/name
   if not p.exists():shutil.copy2(ROOT/'CLOUD-A359DATA'/name,p)
 print(lane)
