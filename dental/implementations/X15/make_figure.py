from dental_release.paths import expand as _release_expand
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent

def main():
    r2 = [json.loads(l) for l in open(HERE / 'RAW_R2.jsonl')]
    r3 = json.load(open(HERE / 'R3_COUPLED_CURVES.json'))
    r4 = json.load(open(HERE / 'SUMMARY_R4.json'))
    (fig, ax) = plt.subplots(2, 2, figsize=(12, 8))
    ray = next((r for r in r2 if r['case'] == 'ToothFairy2F_001' and r['tooth'] == 44 and (r['height_from_apex_mm'] == 5)))
    p = np.load(ray['profile_path'])
    x = p['x_mm']
    gray = p['gray']
    a = ax[0, 0]
    a.plot(x, gray, color='#246a8b', label='Local CBCT raw gray')
    for (side, name, color) in ((1, 'buccal', '#af4d32'), (-1, 'lingual', '#276d50')):
        m = ray[name]
        for k in ('root_x_mm', 'outer_x_mm'):
            a.axvline(m[k], ls='--', color=color, alpha=0.6)
        for e in ('root_image_edge', 'bone_image_edge'):
            if m[e]['status'] == 'MEASURED_IMAGE_EDGE_SURROGATE':
                a.axvline(side * m[e]['HM_mm'], color=color, alpha=0.8)
    a.set_xlim(-7, 7)
    a.set_xlabel('Buccolingual coordinate along ray (mm)')
    a.set_ylabel('Raw gray (not HU)')
    a.set_title(_release_expand('A  Actual TF2 @DENTAL_CASE_ID@ / 44, 5 mm above apex\nDashed label edges; solid image edge surrogates'), fontsize=10)
    a = ax[0, 1]
    offsets = [r[n]['image_minus_label_mm'] for r in r2 for n in ('buccal', 'lingual') if r[n].get('image_status') == 'MEASURED_IMAGE_CLEARANCE_SURROGATE']
    a.hist(offsets, bins=35, color='#246a8b')
    a.axvline(0, color='black', ls='--')
    a.axvline(np.median(offsets), color='#af4d32')
    a.set_xlabel('Image minus annotation envelope separation (mm)')
    a.set_ylabel('Directional rays')
    a.set_title(f'B  New image information, n={len(offsets)}\nMedian {np.median(offsets):.3f} mm; anatomy truth UNKNOWN', fontsize=10)
    a = ax[1, 0]
    curve = next((c for c in r3 if c['case'] == 'ToothFairy2F_001' and c['tooth'] == 44 and (c['source'] == 'label')))
    theta = np.array(curve['angles_deg'])
    lo = np.array(curve['lower_t_mm'])
    hi = np.array(curve['upper_t_mm'])
    a.fill_between(theta, lo, hi, where=lo <= hi, color='#246a8b', alpha=0.15, label='Annotation section contract')
    a.fill_between(theta, lo + 0.608, hi - 0.608, where=lo + 0.608 <= hi - 0.608, color='#246a8b', alpha=0.55, label='0.608 mm edge scenario')
    for (ti, an) in ((0, 0), (0.25, 0), (0.5, 5), (1, 10)):
        a.plot(an, ti, 'o', color='#af4d32')
    a.axhline(0, color='gray', lw=0.5)
    a.axvline(0, color='gray', lw=0.5)
    a.set_xlabel('Declared rigid rotation (degrees)')
    a.set_ylabel('Translation in buccal direction (mm)')
    a.set_title(_release_expand('C  Joint motion set, @DENTAL_CASE_ID@ / 44\nFinite planes and sampled points; not a treatment limit'), fontsize=10)
    a.legend(fontsize=8)
    a = ax[1, 1]
    n = r4['central_R3_admitted_queries']
    surv = r4['admissions_surviving_sampled_R4']
    a.bar(['Central planes\nadmit', 'Survive moved\n3D samples'], [n, surv], color=['#246a8b', '#af4d32'])
    a.set_ylabel('Research queries at epsilon = 0.608 mm')
    a.set_title(f"D  Actual moved 3D surface requery\nSurvival gate: {r4['survival_gate']}; biological accuracy UNKNOWN", fontsize=10)
    for (i, v) in enumerate((n, surv)):
        a.text(i, v + 0.5, str(v), ha='center')
    fig.suptitle('X15: measurable geometry, changed information, and failed validity gates', fontsize=13)
    fig.tight_layout()
    fig.savefig(HERE / 'BONE_ENVELOPE_FIGURE.png', dpi=160)
    fig.savefig(HERE / 'BONE_ENVELOPE_FIGURE.pdf')
    plt.close(fig)
if __name__ == '__main__':
    main()
