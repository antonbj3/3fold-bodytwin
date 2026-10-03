#!/usr/bin/env python3
"""Prepare the eight public-data BodyTwin cloud briefs; never copy source files."""
from pathlib import Path
import hashlib
import subprocess

ROOT = Path(__file__).resolve().parent / 'lanes'
COMMON = '''# Rules
Use only public data that you download yourself and code you write here. No restricted model data, TLEM, the collaborator data, LHDL, internal A matrices, private repository or private file. The register statements below are unverified claims to challenge, not data or ground truth. Cite primary source URLs and DOI where applicable; state what was actually downloaded. Use a reproducible, small Python implementation with fixed random seeds and no hidden inputs. Distinguish derivation, synthetic test, and empirical validation. If a source or required data cannot be obtained, report UNKNOWN rather than inventing a result. Do not push, publish, email, or call subagents.

Final message format: print exact UTF-8 contents of RESULTS.md, results.json, and small code files with `=====FILE <relative path>=====` then `=====END FILE=====`. Combined files ≤200 KB. Results must be independently checkable from cited sources and code. The cloud VM receives this git checkout only.
'''
LANES = {
'CLOUD-A359': (
 'Independent NHANES strength versus DXA reproduction',
 'Register A359 reports NHANES 1999–2002 MSX + DXA, age ≥50, held-out 2001–02 n≈1485: lean-leg×height RMSE 28.7 Nm; mass-only 36.5; DXA leg lean improves 19.1% in holdout. These numbers are a claim to audit, not input data.',
 'Strongest baseline: prespecified mass/height/sex/age model and landmark-based LengthMass and LengthMassFat formulas described only as formulas, without restricted model data files. Counterexamples: BMI fat correction reportedly worsened CV; sex residuals +18/−17 Nm and MSDPF-to-torque lever-arm inference may bias comparison.',
 'Download public CDC NHANES strength and DXA tables yourself; identify exact files, variables, merge keys, exclusions and units. Fit on 1999–2000, test once on 2001–02; add grouped CV in training. Compare mass-only, lean-leg, lean-leg×height with same folds and complete-case sample. Report n, missingness, RMSE, uncertainty by bootstrap, sex residuals and sensitivity to lever-arm conversion.',
 'Primary: a CDC dataset/data dictionary URL for each table. Pass if held-out DXA model improves mass-only RMSE by ≥10% on the same persons and no leakage; otherwise fail/UNKNOWN. Never claim the register number reproduced without matching cohort and target.'),
'CLOUD-A367': (
 'Public patellar tendon moment arm versus knee angle',
 'Register A367: BodyTwin simplified arm at 0/45/90° = 51.4/42.7/35.6 mm, maximum at 0°. It cites Krevolin et al. 2004 DOI 10.1016/j.jbiomech.2003.09.010 as a peak around 45°, but literature RMSE is UNKNOWN.',
 'Strongest baseline: constant arm and simple monotonic three-point interpolation. Counterexample: a public mechanism can have a maximum at 45° only under a particular angle convention, tendon routing and measurement definition; check all three.',
 'Find open primary knee model or digitizable public table with DOI. Build an explicit smooth low-parameter r(theta) with maximum near 45° and positive arms across 0–90°, fit with held-out angles or a second source. Show parameters, curve/table, uncertainty, units, angle convention and a runnable plotting/evaluation script.',
 'Pass only if a public quantitative source supports a peak 30–60° and the model beats constant and monotonic baselines on held-out points. Otherwise UNKNOWN/fail; do not fabricate digitized values.'),
'CLOUD-B5DISC': (
 'L5/S1 lever-arm and intradiscal pressure audit',
 'BodyTwin plan B5 proposes a thin enclosed disc layer K=(M/t)A and Wilke in vivo intradiscal pressure as an external target. No validated local numerical result is supplied.',
 'Strongest baseline: published measured pressure by posture/activity and a simple axial-force divided by effective disc area model. Counterexamples: nucleus pressure is not uniform stress; muscle co-contraction, orientation and lever arms can dominate body-weight-only estimates.',
 'Find public Wilke primary paper DOI and numerical pressure values and another primary source for disc geometry/load. Derive a transparent moment equilibrium for at least neutral standing and flexion, propagate area, lever-arm and muscle-force uncertainty, compare predicted MPa to measured MPa. Code the model and unit tests.',
 'Pass if model improves an axial body-weight-only baseline for at least two prespecified postures without fitted posture-specific parameters; otherwise fail/UNKNOWN. Treat clinical interpretation as out of scope.'),
'CLOUD-B1OED': (
 'Minimal individual measurement protocol by optimal design',
 'Register A72: 5 radiographic femur measures cut surface RMS 2.684→2.094 mm in 35 LOSO folds; two external measures did not beat affine. A358: shared knee-physics terms had N_eff≈1 and lacked individual geometry. A275 prioritizes individual HJC and imaged muscle attachments; values are unreviewed.',
 'Strongest baseline: equal-cost greedy selection by independent Gaussian Fisher information. Counterexamples: correlated systematic errors destroy assumed independence; DXA lean-mass prediction A359 may not identify moment arms; X-ray cannot directly measure strength.',
 'Derive a small transparent OED over candidate DXA, X-ray, landmarks, strength and knee angle measurements. Use only the register facts above as qualitative constraints, not fabricated patient matrices. State assumed Jacobian/noise and show sensitivity across plausible ranges, correlations and costs. Provide Python script and a ranked protocol under at least two budgets.',
 'Pass only if top recommendation survives ≥80% of stated sensitivity scenarios and beats equal-cost random and one-modality baselines on posterior target variance. Synthetic evidence must be labeled synthetic.'),
'CLOUD-NATO': (
 'Activity to tissue load to epsilon-N fatigue to injury risk',
 'Register A126: internal unreviewed chain predicts huge uncertainty; density–modulus accounts for reported 0.676 of log-damage variance, and a QCT relation shrinkage was claimed 22.8×. Do not reuse its internal patient inputs.',
 'Strongest baseline: exposure count times constant strain amplitude with Miner accumulation, compared with a public bone ε-N reference. Counterexamples: healing/remodeling, sequence dependence, muscle force and density–modulus uncertainty can reverse risk ranking.',
 'Use a public ε-N/bone-fatigue primary source with DOI and units. Implement an explicitly synthetic activity→load→strain→cycles→damage chain, uncertainties per link and variance/sensitivity decomposition; compare constant-load baseline. Show at least one counterexample in which ranking changes under plausible source-backed uncertainty.',
 'Pass if all links have dimensional checks, source-backed parameter ranges, uncertainty budget, and a reproducible ranking-sensitivity test. Do not call damage probability observed injury risk without epidemiology; label risk mapping UNKNOWN.'),
'CLOUD-N1GAUDIT': (
 'Independent constraint-net N_eff audit',
 'Register A324: N1g activity×GRF baseline RMSE 0.378 BW; five physics additions failed; JW transplanted geometry was a counterexample. A358: certified clipping 0.448/0.447, shared physics terms N_eff≈1, placebo ≈ real. These are claims, not raw data.',
 'Strongest baseline: N1g alone; compare mathematical estimators using the same latent shared input. Counterexamples: perfectly correlated constraints yield only one source, but independent measurements or distinct person-specific Jacobians may raise effective rank.',
 'Prove under specified assumptions why the effective information rank can be 1, distinguish rank from sample-size N_eff, and construct numerical counterexamples. Simulate a simple shared-geometry versus individual-geometry model and exact placebo; give necessary individual measurements for a useful certificate. Include code that computes eigenvalues, variance and calibration.',
 'Pass if derivation states the covariance/Jacobian assumptions, simulation reproduces rank-one collapse, and at least one individual-information counterexample increases rank and passes an independent test. Do not claim the actual BodyTwin N_eff is proven without matrices.'),
'CLOUD-VASCULAR': (
 'Vascular–interstitial Starling/lymph transient',
 'BodyTwin plan BT-B suggests a coupled vascular–interstitial transient. No private physiological code or measurement series is provided in this bundle.',
 'Strongest baseline: one-compartment constant filtration and fixed lymph drainage. Counterexamples: revised Starling principle/glycocalyx, nonlinear lymph response and conservation failure can invalidate a simple linear law.',
 'Find a public primary source with DOI and a quantitative transient or steady target. Build a two-volume mass-conserving model with pressure/oncotic filtration, lymph return and at least one transported solute; simulate a bounded perturbation. Compare baseline and revised model, check nonnegative volumes, mass conservation and parameter identifiability. Include runnable code.',
 'Pass if conservation residual <1e-8 relative in numerical tests and a public quantitative target is matched better than the baseline without extra fitted free parameters; otherwise report the partial result and UNKNOWN.'),
'CLOUD-B4CURVE': (
 'Strength–angle curve as identifiable individual signature',
 'Register A359 reports individual strength poorly predicted by mass laws relative to DXA lean mass; A367 reports wrong patellar-arm peak in a simplified model. Neither establishes individual curve identifiability.',
 'Strongest baseline: population mean torque–angle curve with individual amplitude only. Counterexamples: tendon moment arm, optimal fiber length and activation can produce nearly collinear curve changes; measurement noise and velocity confound identity.',
 'Locate primary public multi-angle knee-strength curves with DOI or an open dataset. Derive a low-dimensional torque(theta) model factoring capacity, force-length and moment arm. Compute Fisher rank/profile likelihood for measurements at 2, 3 and ≥5 angles; cross-validate held-out angles or clearly label synthetic simulations. Recommend minimum angles and report unresolved parameter combinations.',
 'Pass only if a real public multi-person dataset supports an individual signature with held-out improvement over amplitude-only and the FIM remains full rank under source-backed noise; else UNKNOWN/fail, with mathematical findings separate.'),
}

def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    for lane, (title, register, baseline, work, criterion) in LANES.items():
        folder = ROOT / lane
        folder.mkdir(exist_ok=True)
        brief = (f'# {lane}: {title}\n\n{COMMON}\n## Register context\n{register}\n\n'
                 f'## Strongest baseline and prior counterexamples\n{baseline}\n\n'
                 f'## Work and deliverables\n{work}\n\n'
                 'Deliver RESULTS.md, results.json, and small runnable Python code. Report provenance, '
                 'method, sample size or synthetic status, baseline, uncertainty, limitations and exact commands.\n')
        prereg = (f'# PREREG {lane}\n\nQuestion: {title}.\n\n'
                  f'Criterion: {criterion}\n\n'
                  'Freeze all thresholds before data retrieval. Record deviations in RESULTS.md; do not edit this file. '
                  'An inaccessible source yields UNKNOWN, never an invented pass.\n')
        for name, content in [('BRIEF.md', brief), ('PREREG.md', prereg)]:
            path = folder / name
            if not path.exists():
                path.write_text(content)
        digest = hashlib.sha256((folder / 'PREREG.md').read_bytes()).hexdigest()
        (folder / 'PREREG.sha256').write_text(f'{digest}  PREREG.md\n')
        if not (folder / '.git').exists():
            subprocess.run(['git', 'init', '-q', str(folder)], check=True)
        subprocess.run(['git', 'add', 'BRIEF.md', 'PREREG.md', 'PREREG.sha256'], cwd=folder, check=True)
        if subprocess.run(['git', 'rev-parse', '--verify', 'HEAD'], cwd=folder,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode:
            subprocess.run(['git', '-c', 'user.name=BodyTwin Cloud',
                            '-c', 'user.email=bodytwin-cloud@localhost',
                            'commit', '-qm', 'Preregister public cloud lane'], cwd=folder, check=True)
        print(lane)

if __name__ == '__main__':
    main()
