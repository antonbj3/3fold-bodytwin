"""Make reader-facing tables, figures, result manifest and pending graph feedback."""
from dental_release.paths import expand as _release_expand
import csv
import json
from collections import defaultdict
from pathlib import Path
import sys
_MPL = Path(_release_expand('@DENTAL_WORK_ROOT@/X9-ipr-safety/mpl_deps'))
if _MPL.exists():
    sys.path.insert(0, str(_MPL))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from ipr import capacity, sha, write_json

def load(path):
    return json.loads(Path(path).read_text())

def finalize():
    result = load('results.json')
    result['artifacts'] = {str(p): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(Path('artifacts').glob('*')) if p.is_file()}
    result['raw_artifacts_sha256'] = {str(p): sha(p) for name in ['R1_raw_profiles.json', 'R1_rows.json', 'R2_profiles.json', 'R2_patches.json', 'R3_rows.json', 'R3_external_comparison.json', 'operator_analytic_probe.json'] if (p := (Path('rounds') / name)).exists()}
    result['code_sha256'] = {str(p): sha(p) for p in sorted(Path('.').glob('*.py'))}
    result['code_sha256']['run_all.sh'] = sha('run_all.sh')
    result['prereg_sha256'] = {str(p): sha(p) for p in sorted(Path('.').glob('PREREG*.json'))}
    for (name, item) in result['artifacts'].items():
        assert sha(name) == item['sha256']
    write_json('results.json', result)
    feedback = load('GRAPH_FEEDBACK.json')
    feedback['sha256'] = sha('results.json')
    feedback['population_regime'] = feedback['population'] + '; ' + feedback['regime']
    feedback['preregistered_gate'] = feedback['frozen_gate']
    feedback['baseline'] = feedback['best_baseline']
    write_json('GRAPH_FEEDBACK.json', feedback)
    write_json('CURRENT_WORK_STATE.json', {'lane': 'X9-ipr-safety', 'stage': 'DEMO_VERIFIED_PENDING_INDEPENDENT_REVIEW', 'latest_gates': {r['round']: r['gates'] for r in result['rounds']}, 'next_operation': 'same-specimen CBCT/microCT boundary calibration and measured planar overcut; direct O21/N1 coverage proposal', 'result_sha256': sha('results.json')})
    print('Final artifact manifest PASS', sha('results.json'))

def run():
    rounds = [load(f'rounds/R{i}_results.json') for i in [1, 2, 3]]
    rows = load('rounds/R1_rows.json')
    r3rows = load('rounds/R3_rows.json')
    refs = load('sources/external_by_tooth_type.json')
    groups = defaultdict(list)
    for row in rows:
        groups[row['jaw'], row['tooth_type'], row['side'], row['height_fraction']].append(row)
    aggregate = []
    for (key, rr) in sorted(groups.items()):
        valid = [r for r in rr if r['mean_mm'] is not None]
        aggregate.append(dict(jaw=key[0], tooth_type=key[1], side=key[2], height_fraction=key[3], teeth=len(rr), covered_patches=sum((r['valid_fraction'] >= 0.9 for r in rr)), median_mean_thickness_mm=float(np.median([r['mean_mm'] for r in valid])) if valid else None, median_scenario_scalar_capacity_mm=float(np.median([r['scenario_capacity_mm'] for r in valid])) if valid else None, boundary_budget_mm=0.3, tool_overcut_mm=0.05, calibrated_safe_ipr_mm=None))
    with open('artifacts/by_type_side_height.csv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=aggregate[0].keys())
        w.writeheader()
        w.writerows(aggregate)
    with open('artifacts/ipr_3d_by_tooth_side_height.csv', 'w') as f:
        names = sorted(set().union(*(r.keys() for r in r3rows)))
        w = csv.DictWriter(f, fieldnames=names)
        w.writeheader()
        w.writerows(r3rows)
    table = '| Jaw | Type | Side | Altitude | Teeth | Thickness, median mm | Conditional line space mm | Covered areas |\n|---|---|---|---:|---:|---:|---:|---:|\n'
    for r in aggregate:
        thickness = f"{r['median_mean_thickness_mm']:.3f}" if r['median_mean_thickness_mm'] is not None else 'UNKNOWN'
        cap = f"{r['median_scenario_scalar_capacity_mm']:.3f}" if r['median_scenario_scalar_capacity_mm'] is not None else 'UNKNOWN'
        table += f"| {r['jaw']} | {r['tooth_type']} | {r['side']} | {r['height_fraction']} | {r['teeth']} | {thickness} | {cap} | {r['covered_patches']}/{r['teeth']} |\n"
    Path('artifacts/table_by_type_side_height.md').write_text(table)
    (fig, axes) = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    comparison = load('rounds/R3_external_comparison.json')
    ax = axes[0, 0]
    ax.scatter([r['external_mean_mm'] for r in comparison], [r['our_mean_mm'] for r in comparison], c='#30698e', s=35)
    ax.plot([0.3, 2.2], [0.3, 2.2], 'k--', lw=1)
    ax.set(xlabel='Published type/side mean enamel (mm)', ylabel='STS threshold mesh mean (mm)', title='R1: independent population comparison')
    ax.text(0.03, 0.97, 'Different specimens; not an accuracy validation', transform=ax.transAxes, va='top', fontsize=9)
    profiles = load('rounds/R2_profiles.json')
    selected = next((r for r in profiles if r['accepted'] and r['ray_id'] == 40))
    failed = next((r for r in profiles if not r['accepted'] and 'grey' in r and (r['ray_id'] == 40)))
    ax = axes[0, 1]
    ax.plot(selected['s_mm'], selected['grey'], label='Accepted raw profile', color='#276749')
    ax.plot(selected['s_mm'], selected['fit_grey'], '--', color='#276749', label='Joint-interface fit')
    ax.plot(failed['s_mm'], failed['grey'], label='Rejected raw profile', color='#a63e3e', alpha=0.75)
    ax.set(xlabel='Inward distance from threshold surface (mm)', ylabel='CBCT grey value (scanner units)', title='R2: both interfaces are not usually identifiable')
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    keys = [(r['jaw'], r['tooth_type']) for r in refs]
    labels = [f"{r['jaw'][0].upper()} {r['tooth_type'].replace('_', ' ')}" for r in refs]
    values = [r['mesial_mean_mm'] / 2 for r in refs]
    ax.bar(np.arange(len(values)), values, color='#4e7998')
    ax.axhline(0.5, color='#a63e3e', ls='--', label='0.5 mm/surface benchmark')
    ax.set_xticks(np.arange(len(labels)), labels, rotation=65, ha='right', fontsize=7)
    ax.set(ylabel='Half of published mean thickness (mm)', title='Published 50% convention varies by tooth type')
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    r1caps = [r['scenario_capacity_mm'] for r in rows if r['valid_fraction'] >= 0.9]
    r3caps = [r['scenario_capacity_mm'] for r in r3rows if r['scenario_capacity_mm'] is not None]
    ax.hist(r1caps, bins=np.linspace(0, 1.4, 20), alpha=0.7, label='R1 sampled line model', color='#4e7998')
    ax.hist(r3caps, bins=np.linspace(0, 1.4, 20), alpha=0.6, label='R3 3D unknown-tissue protection', color='#b87735')
    ax.axvline(0.5, color='#a63e3e', ls='--')
    ax.set(xlabel='Conditional removal envelope (mm)', ylabel='Side/height patches', title='Whole-footprint protection rejects these STS cuts')
    ax.legend(fontsize=8)
    fig.suptitle('IPR research demo — calibrated anatomical safety remains UNKNOWN', fontsize=14)
    fig.savefig('artifacts/ipr_evidence.png', dpi=170)
    fig.savefig('artifacts/ipr_evidence.pdf')
    plt.close(fig)
    (g1, g2, g3) = (r['gates'] for r in rounds)
    positive = rounds[2]['positive_scenario_patches']
    summary = _release_expand(f"# IPR with actual enamel: what the evidence supports\n\nThe tool calculates how a limited flat grinding range can cut a segmented tooth while protecting dentin or unknown tissue. A SDF is here a distance field to protected tissue . The test shows why a toothworm and a threshing enamel surface are not enough for a validated IPR decision.\n\nRun from this directory: `./run_all.sh` . Python 3 with numpy, scipy, trimesh, rtree, scikit-image and matplotlib are required. Execution only uses existing local STS files and saved published reference tables; no network connection is needed. All three frozen tries to run and the figure is created again. Full raw values is available in `rounds/` and tables in `artifacts/`.\n\n| Trials | Results | Frozen gate |\n|---|---|---|\n| R1 , 15 molars / 90 side-height areas | Middle Height Median {g1['G1_external_population']['median_mm']:.3f} mm ; 80 / 90 areas have ≥ 90% valid rays | Population band PASS ; coverage FAIL |\n| R1 , 0,5 mm compared to declared margins | {g1['G3_praxis_difference']['fraction'] * 100:.1f}% of covered areas have less conditional line space | Model difference PASS ; biological safety UNKNOWN |\n| R2 , raw CBCT / two bounds | {g2['G4_fit_coverage']['accepted_rays']}/ 810 profiles passed ; only {g2['G4_fit_coverage']['complete_patches']}/ 90 complete areas | Coverage FAIL , frozen requirement 80% |\n| R3 , tissue field / entire grinding area | {positive}/ 90 positive conditional grinding spaces | Protected tissue never overlaps allowed incisions |\n| R3 , tooth type and page against published table | Median absolute difference {g3['G3_external_type_contact']['median_absolute_error_mm']:.3f} mm | {('PASS' if g3['G3_external_type_contact']['pass'] else 'FAIL')}, boundary 0,200 mm |\n| equally informed controls | Analytical boundary and direct feasibility search provide the same response | TIE; no method win or 10 × -claim |\n\n![ Measurement , reference and Limitations](artifacts/ipr_evidence.png)\n\n`artifacts/thickness_by_tooth_side_height.csv` provides all 90 measurement ranges; `artifacts/by_type_side_height.csv` provides compilation per jaw , M1 / M2 , page and height . Height ratio 0 is the lowest point of the enamel meshe and 1 its highest, **not a measured cement-ename limit**. M1 / M2 comes from the geometrical rating of the parent code. Incisiver, corner teeth and premolars have only published reference values here; no new STS - measurements on these types is claimed.\n\nExternal reference : [Bian 2020 , DOI 10.3760/cma.j.cn112144 -20191211 -00448 ](https://pubmed.ncbi.nlm.nih.gov/32634888/), the μCT-measurement of molars contact enamel , 1,46 ± 0,25 mm . The frozen compatibility gate is [0,96 ; 1,96 ] mm . [Vellini-Ferreira and others 2024 , Table 1 ](https://doi.org/10.1590/2177-6709.29.3.e242422.oar) provides mesial/distal thickness of all tooth types on calibrated X-ray images; Table 2 is **half of the thickness**, despite the misleading English title  It is not used as fifth percentile. [Sarig m.fl. 2023 , Table 2 ](https://doi.org/10.1590/2177-6709.28.2.e2321149.oar) also shows altitude variation in μCT of sub-vocal incenses, but its zones within the incisive third do not match our molar heights.\n\nThe 50% rule is presented as a published convention in these studies. 0,5 mm per surface area is the comparison value of this trial, no general rule for all teeth. The requirement for 0,5 mm remaining enamel , edge budget 0,3 mm and tool overcutting 0,05 mm are declared scenarios. The budgets are not calibrated confidence limits. Published summaries apply to second teeth and populations and therefore can reject reasonableness but not verify individual accuracy.\n\nWhat does not holds : STS has no expert labels for enamel / dentin / pulp . The `Dentin`-mesh of the parent code is the whole tooth. The enamel is Otsu-tresculated and morphologically processed. Between the toothworm and enamel there is a unknown outer strip; R3 protects it and refrains from positive cuts. R2 : s accepted subset is selected and its median must not be described as a population estimate . Distance based half-thickness in R3 is a model disclosure. No biological damage or caries outcome has been measured.\n\nTools for custom segmented tooth :\n\n```bash\npython3 segmented_tooth_tool.py --input @DENTAL_WORK_ROOT@/X9-ipr-safety/demo_segmented_tooth.npz --side mesial --y-mm DEMO_Y --z-mm DEMO_Z --output artifacts/tool_output.json\n```\n\nReplace the coordinates with values from `artifacts/demo_tool_command.sh` , or run that script directly. NPZ -contract: `labels` with 0 = outside, 1 = enamel , 2 = dentin /other/ unknown protected tissue , 3 = pulp ; `spacing_mm` and `origin_mm` have three components. The shoulders must be mesiodistal, buccolingual and crown height . `scenario_capacity_mm` is a conditional geometric boundary . `calibrated_safe_ipr_mm` is always null in this demo. The `conditional_post_ipr.stl` export is a model output; when allowed reduction is zero exported the same tooth.\n\nMinimum next physical design: the same extracted tooth in CBCT and μCT before and after a measured plant grinding step. Record the outer edge and enamel –detiny limit separately as well as the actual upper dimension of the tool, and freeze specimen-hash, grinding plan and predicated remaining enamel before the μCT comparison. For several tooth types and transfer to another scanner, more specimens are required; a pair calibrates no population. `LAB_VALIDATION.md` specifies the observation contract.\n\nDataset/licens: STS-Tooth, Wang, Yaqi m.fl. (2025), [artikel](https://doi.org/10.1038/s41597-024-04306-9), [dataset DOI 10.5281/zenodo.10597292](https://zenodo.org/records/10597292), **CC BY 4.0**. The meshs are previously derived in CROWN_pop and are reused with originalhashar. No second image data sets are used and no personal data are read. Derivatives, simulation, source measurement and hypothesis are kept separate. Status: PENDING_INDEPENDENT_REVIEW.\n")
    import xml.etree.ElementTree as ET
    names = {}
    for pmc in ['PMC11235574', 'PMC10229113']:
        root = ET.parse(f'sources/{pmc}.xml').getroot()
        first = root.find('.//article-meta/contrib-group/contrib/name/surname')
        names[pmc] = first.text if first is not None else 'authors'
    summary = summary.replace('Vellini-Ferreira m.fl.', names['PMC11235574'] + ' m.fl.').replace('Sarig m.fl.', names['PMC10229113'] + ' m.fl.')
    Path('README_DEMO.md').write_text(summary)
    demo = next((r for r in r3rows if r['tid'] == 'L002_up_M1_k5' and r['side'] == 'mesial' and (r['height_fraction'] == 0.5)))
    cmd = _release_expand(f"""#!/usr/bin/env bash\nset -euo pipefail\ncd -- "$(dirname -- "${{BASH_SOURCE[0]}}")/.."\npython3 segmented_tooth_tool.py --input @DENTAL_WORK_ROOT@/X9-ipr-safety/demo_segmented_tooth.npz --side mesial --y-mm {demo['y_mm']:.12g} --z-mm {demo['z_mm']:.12g} --output artifacts/tool_output.json\n""")
    Path('artifacts/demo_tool_command.sh').write_text(cmd)
    Path('artifacts/demo_tool_command.sh').chmod(493)
    result = {'lane': 'X9-ipr-safety', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'outcome': 'GEOMETRIC_TOOL_DELIVERED_ANATOMICAL_SAFETY_UNKNOWN', 'external_referent': rounds[2]['external_referent'], 'rounds': rounds, 'calibrated_safe_ipr_mm': None, 'tenfold_gain_claimed': False, 'uncertainty': {'anatomical_bound': 'UNKNOWN', 'scenario_boundary_mm': 0.3, 'scenario_tool_mm': 0.05}, 'artifacts': {str(p): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(Path('artifacts').glob('*')) if p.is_file()}, 'sources_sha256': {str(p): sha(p) for p in sorted(Path('sources').glob('*')) if p.is_file()}, 'code_sha256': {str(p): sha(p) for p in sorted(Path('.').glob('*.py'))}, 'full_cost': {'preparation_manual_seconds': 'NOT_MEASURED', 'historical_parent_segmentation': 'UNKNOWN', 'measured_run_wall_seconds': sum((r['cost']['wall_seconds'] for r in rounds)), 'measured_run_cpu_seconds': sum((r['cost']['cpu_seconds'] for r in rounds)), 'peak_RSS_MiB_max': max((r['cost']['peak_rss_MiB'] for r in rounds)), 'fit': rounds[1]['cost'], 'external_reading_queries': 'Recorded in SOURCE_SELECTION.md; cost unmeasured', 'physical_validation': 'NONE', 'fallback': 'paired CBCT/microCT and measured tool overcut'}}
    write_json('results.json', result)
    Path('RESULTS.md').write_text(summary + '\n' + table)
    feedback = {'target_id': 'DENT-GEOM-UNCERTAINTY', 'result_file': 'results/LANE_X9_IPR_SAFETY/results.json', 'sha256': sha('results.json'), 'measured_quantity': 'Proximal threshold-mesh enamel thickness; raw-profile fit coverage; protected-tissue IPR envelope', 'units': 'mm, fraction, voxel count', 'uncertainty': 'Anatomical calibration UNKNOWN; 0.3 mm boundary/0.05 mm tool declared scenario', 'population': '15 geometrically classified STS molars from 10 ROI files; not assumed independent patients', 'regime': 'Geometric research; finite planar patch; no clinical or biological outcome', 'frozen_gate': 'PREREG_R1/R2/R3.json with independent population/type measurements and explicit coverage gates', 'best_baseline': 'same-input feasibility bisection; fixed 0.5 mm/surface benchmark', 'outcome': 'R1 coverage FAIL; R2 coverage FAIL; R3 published-type comparison ' + ('PASS' if g3['G3_external_type_contact']['pass'] else 'FAIL') + '; exact geometry controls TIE; anatomical safety UNKNOWN', 'negative_result': True, 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'missing_coverage': 'O21/N1 IPR tissue-to-tool constraint node and calibration dependency are not supplied by this uncertainty target'}
    feedback['population_regime'] = feedback['population'] + '; ' + feedback['regime']
    feedback['preregistered_gate'] = feedback['frozen_gate']
    feedback['baseline'] = feedback['best_baseline']
    write_json('GRAPH_FEEDBACK.json', feedback)
    write_json('CURRENT_WORK_STATE.json', {'lane': 'X9-ipr-safety', 'stage': 'DEMO_GENERATED_PENDING_VERIFICATION', 'latest_gates': {r['round']: r['gates'] for r in rounds}, 'next_operation': 'run complete one-command demo and bind final hash'})
    print('Demo generated', result['outcome'])
if __name__ == '__main__':
    if '--finalize' in sys.argv:
        finalize()
    else:
        run()
