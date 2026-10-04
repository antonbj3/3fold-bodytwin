"""Reader-facing figures, result package and a bounded graph-coverage proposal."""
from dental_release.paths import expand as _release_expand
import sys, datetime, resource
sys.path.insert(0, _release_expand('@DENTAL_IMPLEMENTATIONS@/X7/code/vendor_plotting'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry import *

def figure(rounds):
    final = rounds[-1]
    rn = final['round']
    rows = json.loads((H / 'raw' / ('RESULTS_' + rn + '_ROWS.json')).read_text())
    row = next((r for r in rows if r['status'] == 'SCORED' and r['type'] == 'M1'))
    preds = json.loads((H / 'raw' / ('PREDICTIONS_' + rn + '.json')).read_text())
    p = next((r for r in preds if r['case'] == row['case'] and r['fdi'] == row['fdi']))
    z = np.load(p['file'])
    suffix = '_R2' if rn == 'R2' else ''
    ref = np.load(D / (str(row['case']) + '_' + str(row['fdi']) + suffix + '_reference.npz'))
    (fig, axes) = plt.subplots(2, 3, figsize=(13, 8))
    xy = z['xy']
    for (ax, arm, name) in zip(axes[0], ['original', 'practice', 'informed'], ['Original scan reference', 'Generic donor crown', 'Measured-antagonist crown']):
        zz = ref['original_z'] if arm == 'original' else z['z_' + arm]
        im = ax.scatter(xy[:, 0], xy[:, 1], c=zz, s=14, cmap='viridis')
        mask = ref['0.1_' + arm + '_contact']
        ax.scatter(xy[mask, 0], xy[mask, 1], s=22, facecolors='none', edgecolors='red', label='|gap|≤0.1 mm')
        ax.set_aspect('equal')
        ax.set_title(name + '\n' + str(row['arms']['0.1'][arm]['patch_count']) + ' patches')
        ax.set_xlabel('RAS x [mm]')
        ax.set_ylabel('RAS y [mm]')
        ax.legend(fontsize=8)
        fig.colorbar(im, ax=ax, label='Height [mm]', shrink=0.6)
    types = ['P1', 'P2', 'M1', 'M2']
    xx = np.arange(4)
    summary = final['summary']['by_type']
    axes[1, 0].bar(xx - 0.18, [summary[t]['practice']['IoU'] for t in types], 0.36, label='Generic')
    axes[1, 0].bar(xx + 0.18, [summary[t]['informed']['IoU'] for t in types], 0.36, label='Antagonist')
    axes[1, 0].set_xticks(xx, types)
    axes[1, 0].set_ylabel('Median contact-map IoU')
    axes[1, 0].set_ylim(0, 1)
    axes[1, 0].legend()
    axes[1, 1].bar(xx - 0.18, [summary[t]['practice']['absolute_area_error_mm2'] for t in types], 0.36, label='Generic')
    axes[1, 1].bar(xx + 0.18, [summary[t]['informed']['absolute_area_error_mm2'] for t in types], 0.36, label='Antagonist')
    axes[1, 1].set_xticks(xx, types)
    axes[1, 1].set_ylabel('Median projected area error [mm²]')
    axes[1, 1].legend()
    fe = json.loads((H / 'raw/FE_RESULTS.json').read_text())
    sphere = []
    for t in types:
        sphere.append(np.median([r['answers']['sphere']['relative_peak_tensile_error'] * 100 for r in fe['rows'] if r['type'] == t and 'relative_peak_tensile_error' in r['answers']['sphere']]))
    axes[1, 2].bar(xx, sphere, color='slateblue')
    axes[1, 2].set_xticks(xx, types)
    axes[1, 2].set_ylabel('Sphere-load stress deviation [%]')
    axes[1, 2].set_title('Same roof, normalized100 N\nConditional FE; physical accuracy UNKNOWN', fontsize=10)
    fig.suptitle('X18: registered opposing geometry constrains crown contacts\n' + rn + ' example case' + str(row['case']) + ', predicted FDI' + str(row['fdi']) + '; contact count is threshold/grid dependent', fontsize=12)
    fig.tight_layout()
    fig.savefig(H / 'figures/crown_antagonist.png', dpi=170)
    fig.savefig(H / 'figures/crown_antagonist.pdf')
    plt.close(fig)

def package():
    rounds = []
    for rn in ['R1', 'R2']:
        p = H / 'rounds' / (rn + '.json')
        if p.exists():
            j = json.loads(p.read_text())
            j['round'] = rn
            rounds.append(j)
    pr = json.loads((H / 'PREREG_R1.json').read_text())
    ver = json.loads((H / 'raw/VERIFICATION.json').read_text())
    final = rounds[-1]
    figure(rounds)
    sources = json.loads((H / 'raw/DONOR_FIT.json').read_text())
    fe = json.loads((H / 'raw/FE_RESULTS.json').read_text())
    files = [dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(D.glob('*')) if p.is_file()]
    total = sum((x['bytes'] for x in files))
    assert total < 3000000000
    result = dict(lane='X18-crown-antagonist', claim_type='information_link', review_state='PENDING_INDEPENDENT_REVIEW', capability='Runnable posterior roof design and scan-referenced contact scoring with continuous projected clearance and conditional FE cross-loads', external_referent=pr['external_referent'], rounds=rounds, physical_force_accuracy='UNKNOWN', physical_stress_accuracy='UNKNOWN', whole_patient_crown_fit='UNKNOWN: cervical-site roof proxy, no full preparation/axial crown', target_tooth_identity='UNKNOWN: X11 transfer accuracy not independently measured in Bits2Bites', numerical_checks=dict(passed=ver['passed'], total=ver['total'], failed_initial_checks='raw/failed_checks/VERIFICATION_R1.json,13/14; occlusal boundary leakage fixed only in R2'), conditional_FE=fe, data_manifest_file='raw/DATA_MANIFEST.json', data_bytes=total, full_cost=dict(donor_fit_seconds=sources['wall_seconds'], R1_preparation_design=json.loads((H / 'raw/PREDICT_COST.json').read_text()), R1_evaluation=rounds[0]['cost'], R2_preparation_design=json.loads((H / 'raw/PREDICT_R2_COST.json').read_text()) if len(rounds) > 1 else None, R2_evaluation=final['cost'] if len(rounds) > 1 else None, failed_attempts='raw/failed_checks/R2_THREADS.json; system matplotlib ABI and graph interface dispatch failures', upstream_X11='Existing pretrained target labels; model training/120patient Teeth3DS validation are external prior preparation costs. No new training charged as0; inherited training total UNKNOWN here. See upstream COST files.', physical_acquisition='Dataset already on disk; scanner/operator costs UNKNOWN; no lab force/strain acquisition', discovery='Current R1+R2 walltimes, checks and interrupted attempt retained; LLM reasoning time/tokens UNKNOWN', fallback='Missing-site and empty-contact cases explicitly UNKNOWN; no fabricated pressure', peak_rss_kib=max(json.loads((H / 'raw/PREDICT_COST.json').read_text())['peak_rss_kib'], json.loads((H / 'raw/PREDICT_R2_COST.json').read_text())['peak_rss_kib']) if len(rounds) > 1 else None), negative_results=['R1 contact recovery gate failed', 'R1 initial site leakage control failed', 'FDI domain transfer and actual contact pressure unvalidated', 'Norms do not identify an individual contact map'], next_operation='Full axial/preparation crown plus independently marked contact pressure at registered pose; freeze before lab measurement. If geometry gate fails, acquire target-domain tooth and cervical annotations before further optimization.')
    if (H / 'rounds/R3.json').exists():
        result['force_uncertainty_capability'] = json.loads((H / 'rounds/R3.json').read_text())
    if (H / 'rounds/R4.json').exists():
        result['personal_morphology_information_link'] = json.loads((H / 'rounds/R4.json').read_text())
        result['full_cost']['R4_design'] = json.loads((H / 'raw/PREDICT_R4_COST.json').read_text())
    if (H / 'raw/SDF_EXPORT.json').exists():
        result['sdf_export'] = json.loads((H / 'raw/SDF_EXPORT.json').read_text())
    if (H / 'raw/VERIFICATION_R2.json').exists():
        result['R2_numerical_checks'] = json.loads((H / 'raw/VERIFICATION_R2.json').read_text())
    if (H / 'raw/VERIFICATION_R4.json').exists():
        result['R4_numerical_checks'] = json.loads((H / 'raw/VERIFICATION_R4.json').read_text())
    if (H / 'raw/PRESSURE_PORT_VERIFICATION.json').exists():
        result['future_pressure_port_checks'] = json.loads((H / 'raw/PRESSURE_PORT_VERIFICATION.json').read_text())
    result['uncertainty'] = dict(quantified_sensitivity='Saved0.05/0.1/0.2mm proximity bands and±0.05mm antagonist pose shifts in raw per-tooth rows', physical_contacts='UNKNOWN pressure and scan registration bias; proximity threshold is not physical contact detection', tooth_type='Conditional on unvalidated target-domain X11 labels', FE='Uncalibrated support/material/traction closures and unmeasured discretization error', population='Descriptive small deterministic samples; no population confidence guarantee')
    if (H / 'raw/PREDICTIONS_R2.json').exists():
        predictions = json.loads((H / 'raw/PREDICTIONS_R2.json').read_text())
        designed = [r for r in predictions if r['status'] == 'DESIGNED']
        result['R2_optimizer_diagnostics'] = dict(sites_attempted=len(predictions), sites_designed=len(designed), LP_success=sum((r['control']['success'] for r in designed)), LP_failures=[dict(case=r['case'], fdi=r['fdi'], message=r['control']['message']) for r in designed if not r['control']['success']], QP_optimality_not_established=sum((not a['qp']['success'] for r in designed for a in r['attempts'])), hard_facet_clearance_evaluated_separately=True)
    dump(H / 'raw/DATA_MANIFEST.json', dict(files=files, total_bytes=total))
    result['data_manifest_sha256'] = sha(H / 'raw/DATA_MANIFEST.json')
    dump(H / 'results.json', result)
    text = 'Original surface contact pattern in registered Bits2Bites - bite rejects both generated crown roof constructions . The measured antagonist removes projected overlaps, but no frozen contact reproduction gate is clear. The code can now follow the same scan through crown shape, contact witnesses, conditional FE, STL and a distance field. Physical contact force and stress accuracy are UNKNOWN ; a full patient crown is not yet delivered.\n\n'
    for r in rounds:
        text += r['round'] + ' — frozen contact reproduction gate ' + r['summary']['decision'] + '.\n\n| Tooth type | n | IoU without/against antagonist | Area error without/against [mm²] | Location error without/against [ mm ] | Max overlap without/against [ mm ] |\n|---|---:|---:|---:|---:|---:|\n'
        for kind in ['P1', 'P2', 'M1', 'M2']:
            g = r['summary']['by_type'][kind]
            (a, b) = (g['practice'], g['informed'])
            text += '| ' + kind + ' | ' + str(g['n']) + ' | ' + f"{a['IoU']:.3f}/{b['IoU']:.3f}" + ' | ' + f"{a['absolute_area_error_mm2']:.3f}/{b['absolute_area_error_mm2']:.3f}" + ' | ' + f"{a['centroid_error_mm']:.2f}/{b['centroid_error_mm']:.2f}" + ' | ' + f"{a['max_penetration_mm']:.3f}/{b['max_penetration_mm']:.3g}" + ' |\n'
        text += '\n'
    text += "Location error is symmetrical distance between patch centers. Missing a map used the site's diagonal as frozen error penalty; the number is then no measured physical contact shift. The area is projected raster area, not pressure bearing area. Gapband0,1mm; sample box size varies with cervical site, approximately 0,2 – 0,5 mm. Several adjacent sample components are not clinically verified contact points.\n\n"
    text += '![ Contact maps and model output](figures/crown_antagonist.png)\n\n'
    text += 'Run `./run_all.sh` . It reads local original mesher, reproduces the contact evaluation, solves FE - load cases , tests independent LP /Hooke-controls and writes the figure. `./run_all.sh --design` recreates R2 designs in the same data directory; execution preserves previous freeze files in raw/replays. Numerical reproduction is not a new blinded validation.\n\n'
    text += 'What is not holds: a clinical complete crown is not delivered. The original Cervical 30% is a virtual preparation/site proxy; axial sidewall, cement, preparation and manufacturing access are not optimized. FDI in these scanss is predicted and lacks independent dental facits. The registration error of the antagonist and the real force distribution are unknown. The gap has been tested for entire projected triangle facets, with floating point arithmetic ; this is no arbitrary 3 D collision certification.\n\n'
    text += 'FE - discretization error has not been quantified by refinement. Stress tables includes eight illustrative positions and is not a population estimate per tooth type. Normal forces acts vertically; friction, real contact normals and motion are missing. R2 : s QP -optimality is not certified when the solar flag is false; all geometric differences and all-facet repair are tested anyway. L1 checker timeout is 3 s per call; timeout is UNKNOWN, no method win.\n\n'
    if 'sdf_export' in result:
        text += '`demo_roof_sdf.npz` provides euclidic distance to the exported closed roof surface on a 0,5 mm grid. Shape change was made with the exact height surface / Triangle representation. The export grid interpolation error is UNKNOWN; it must not replace the contact gate.\n\n'
    text += "FE compares load maps on exactly the same original roof with 100N in each non-empty load case. E= 210GPa , ν = 0,3 , 1,2 mm vertical shell thickness and clipped outer support ring are declared Closures. The stress deviation of the sphere is the model's sensitivity to load map error, not measured patient voltage failure. Tomma contact maps provides UNKNOWN_NO_CONTACT; they are not replaced with invented force. Full raw values and equally informed control outcomes are in results.json. R1 : s direct projection and separately expressed semiplane control provide the same heights; R2 : s L1 - LP runs on the same geometric constraints.\n\n"
    if 'force_uncertainty_capability' in result:
        r3 = result['force_uncertainty_capability']
        text += 'R3 makes the missing knowledge of the power measurable in the model. The same contact geometry and the same 100 N provide the following allowed upper stress levels when force distribution is unknown.\n\n| Predicted tooth | Patcher | Area-viktad topp [MPa] | Sharp top [ MPa ] | Kvot |\n|---|---:|---:|---:|---:|\n'
        for row in r3['rows']:
            if 'upper_over_area_peak' in row:
                text += f"| {row['fdi']} | {row['patch_count']} | {row['area_weighted_peak_MPa']:.2f} | {row['sharp_upper_peak_MPa']:.2f} | {row['upper_over_area_peak']:.2f} |\n"
        text += '\nThe upper level is exactly for the declared linear roof model and force simplex : peak tensile stress is convex in the positive force fractions and reaches its maximum at a corner point. A separately mixed load solution matches the tensor base; injected 2 × tensor is rejected. [Boyd/Vandenberghe, chapter 3 ](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf) supports the convexity step. A total bite force gives no shape channels. N − 1 separate force fractions is sufficient for N patcher; this is not a proven minimum for a single peak voltage. A registered pressure image can provide multiple such channels at once. The model lacks material damage, cement and really support ; the peaks are not verified patient voltages.\n\n'
    if 'personal_morphology_information_link' in result:
        rr = result['personal_morphology_information_link']['summary']
        text += "R4 adds the same patient's corresponding tooth on the second side as a measured shape prior, with four separate arms on the six new cases. Frozen outcome: " + rr['decision'] + '. Too large holes in the predicted tooth surface is rejected; they are not filled with an invented personal tooth.\n\n| Tooth type | Scorad sites | IoU generic without/with antagonist | IoU personal without/with antagonist |\n|---|---:|---:|---:|\n'
        for t in ['P1', 'P2', 'M1', 'M2']:
            r = rr['by_type'][t]
            aa = r['arms']
            values = [aa[a]['IoU'] for a in ['generic_without', 'generic_with', 'personal_without', 'personal_with']]
            formatted = ['UNKNOWN' if x is None else f'{x:.3f}' for x in values]
            text += '| ' + t + ' | ' + str(r['n']) + ' | ' + formatted[0] + '/' + formatted[1] + ' | ' + formatted[2] + '/' + formatted[3] + ' |\n'
        text += '\nThe requirement was at least six scored sites per tooth type. Personal form and antagonist are different information channels; R4 : s separate comparison personal to/without antagonist isolates bite information. When coverage gate fails is available, there is no population gain to report. Anatomic FDI - reference and real preparation surface are still missing.\n\n'
    text += 'Published contact data provides limitations, not a universal crown rule: [Qadeer2023 , DOI10.1111/joor.13451 ](https://doi.org/10.1111/joor.13451) compiles 11 – 70 contacts for the entire dentition depending on the registration method. [Kordaß 2023 , DOI10.1016/j.aanat.2023.152112 ](https://doi.org/10.1016/j.aanat.2023.152112) displays shifting ABC patterns and contact-free posteriora teeth ; this does not justify a fixed number per tooth . R2 :s 1 – 4 patcher/ 0,2 – 5 mm² is a visible design resolution that originality gets reject .\n\n'
    text += "Dataset/licence: [ Bits2Bites_v01 ](https://ditto.ing.unimore.it/bits2bites/) CC BY-NC-SA4.0 according to project briefing and upstreamX 2 . Zip reads member for member, no archives unpacked. No personal data or report texts are used. Bite2Text is not used in the numerical trial; its exact license is UNKNOWN in X11 and no additional data was needed. Teeth3DS is upstream training data for X11 ; authors' [shallenge-repo](https://github.com/abenhamadou/3DTeethSeg_MICCAI_Challenges#license) indicates CC BY-NC-ND4.0 for dataset; MIT for code; exact local snapshot license is not verified. No Teeth3DS -mesher exported. ToothFairy2 and mandible defects are not used.\n\n"
    text += 'Next trial: annotate cervical tape and teeth independently, manufacture a full crown on known die, register antagonist and calibrated pressure image at known total force. Then measure fit and load-shift towards frozen predictions. No lab readings have been performed. All result is waiting for independent review.\n'
    text = text.replace('100N', '100 N').replace('100N', '100 N').replace('cervikala30%', 'cervikala 30 %').replace('gapband0,1', 'gapband 0,1').replace('11–70', '11–70').replace('normalized100 N', 'normalized 100 N')
    (H / 'README_DEMO.md').write_text(text)
    (H / 'RESULTS.md').write_text(text)
    feedback = dict(target_id='DENT-CROWN-MEASURED-ANTAGONIST-LOTO', related_existing_id='DENT-REG-IF-OCCLUSAL', status='PENDING_INDEPENDENT_REVIEW', review_state='PENDING_INDEPENDENT_REVIEW', claim_type='information_link', result_file='results/LANE_X18_CROWN_ANTAGONIST/results.json', sha256=sha(H / 'results.json'), outcome=final['summary']['decision'], measured_quantity='scan-derived posterior contact IoU, area and projected penetration', units=['dimensionless', 'mm2', 'mm'], uncertainty='Target FDI and pose error UNKNOWN; 0.05/0.2mm sensitivity reported; force/stress physical validity UNKNOWN', population_regime='Two disjoint12case Bits2Bites samples; lower P1/P2/M1/M2 with cervical site proxy', preregistered_gate='IoU median gain>=.15, area error ratio<=.8, every informed projected penetration<=1e-7mm', baseline='Same generic donor without registered antagonist; declared software practice proxy', negative_result=final['summary']['decision'] != 'PASS', external_referent=pr['external_referent'], coverage_status='MISSING_SCOPED_TARGET: existing interface has no numerical dispatch permission/baseline; attempted dispatch failed, no source/working graph edited', next_operation=result['next_operation'])
    if 'personal_morphology_information_link' in result:
        feedback['R4_outcome'] = result['personal_morphology_information_link']['summary']['decision']
    dump(H / 'GRAPH_FEEDBACK.json', feedback)
    handoff = 'Read README_DEMO.md and results.json. All source members and saved label hashes checked. No target occlusal surface enters design except the preserved R1 percentile-boundary leakage found by injection; exact-band extraction fixes it in R2.\n\n' + ''.join((r['round'] + ': ' + r['summary']['decision'] + ', n=' + str(r['summary']['by_type']['all']['n']) + ', median IoU gain=' + str(r['summary']['by_type']['all']['median_IoU_gain']) + ', area ratio=' + str(r['summary']['by_type']['all']['area_error_ratio']) + '.\n' for r in rounds)) + '\nFE stress numbers are normalized scenario model differences, no external physical stress reference . Contact maps derive from independent measured geometry, no pressure observations. Full preparation/sidewalls missing. The original source scan reference can return contact recovery and did refereR 1 .\n\nNext construction: independently annotated cervical/tooth region and whole crown axial surfaces, plus a registered calibrated pressure-image measurement. Continue the parent actual-antagonist capability, not a renamed clipping rule. Preserve all failed gates and old prediction hashes. Graphcoverage: proposal DENT-CROWN-MEASURED-ANTAGONIST-LOTO, related interface DENT-REG-IF-OCCLUSAL cannot carry this numerical claim. No generated/native graphs edited.\n\nRun ./run_all.sh. Arrays ' + str(round(total / 1000000.0, 1)) + 'MB under ' + str(D) + '; hashes raw/DATA_MANIFEST.json. No GPU/subagents/mail/publication.\n'
    (H / 'HANDOFF.md').write_text(handoff)
    state('REVIEWABLE_DEMO', final['summary']['decision'], 'Independent tooth/preparation and force measurement; construct full crown geometry next', data_bytes=total, rounds_completed=len(rounds))
    if 'personal_morphology_information_link' in result:
        r4 = result['personal_morphology_information_link']['summary']
        r4handoff = '\nR4 actual measured same-patient homolog prior + antagonist factorial: ' + r4['decision'] + '. Coverage per type ' + str({k: r4['by_type'][k]['n'] for k in ['P1', 'P2', 'M1', 'M2']}) + '. Missing/fragmented source roof returned UNKNOWN; do not lower85%/0.8mm gates. Read PREREG_R4.json, FROZEN_PREDICTIONS_R4.json, raw/RESULTS_R4_ROWS.json. Next load-bearing construction requires validated tooth/cervical/preparation surfaces or a representation that preserves source surface support, before another generic envelope optimizer. Then acquire pressure channels;100N alone cannot calibrate contact shares.\n'
        (H / 'HANDOFF.md').write_text(handoff + r4handoff)
        state('REVIEWABLE_FOUR_CONSTRUCTIONS', r4['decision'], 'Source-support-preserving tooth/preparation representation plus independent pressure channels', data_bytes=total, rounds_completed=4)
    print('Packaged', len(rounds), 'rounds,', round(total / 1000000.0, 1), 'MB')
if __name__ == '__main__':
    package()
