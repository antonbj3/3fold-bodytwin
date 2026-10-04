"""Readable Swedish answer-first report, printable tables and standard plots."""
import html, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parent

def mdtable(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |'] + ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])

def htmltable(headers, rows):
    return '<table><thead><tr>' + ''.join(('<th>' + html.escape(str(x)) + '</th>' for x in headers)) + '</tr></thead><tbody>' + ''.join(('<tr>' + ''.join(('<td>' + html.escape(str(x)) + '</td>' for x in r)) + '</tr>' for r in rows)) + '</tbody></table>'

def render(result, gt, tt):
    from analysis import GUIDE_NAMES, GUIDE_ORDER, TYPE_NAMES
    (ROOT / 'figures').mkdir(exist_ok=True)
    allg = [next((t for t in gt if t['guide'] == g and t['region'] == 'ALL')) for g in GUIDE_ORDER]
    gh = ['Guide type', 'Apex scenario¹ mm [B-band]', 'Whole body² mm [B-band]', 'Changed scenario / 100 [ 95% KI ]', 'Clinical requirement']
    gr = []
    for t in allg:
        p = result['guide_parameters'][t['guide']]
        sb = p['scenario_boundary_sensitivity_range_mm']
        bb = p['whole_body_boundary_sensitivity_range_mm']
        gr.append([GUIDE_NAMES[t['guide']], f"{p['scenario_apex_required_gap_mm']:.2f} [{sb[0]:.2f}–{sb[1]:.2f}]", f"{p['whole_body_sufficient_gap_B0p3_mm']:.2f} [{bb[0]:.2f}–{bb[1]:.2f}]", f"{t['scenario_changed_per100']:.1f} [{t['scenario_changed_CI95_low']:.1f}–{t['scenario_changed_CI95_high']:.1f}]", 'UNKNOWN'])
    th = ['Tooth type (n)', 'Reduktion mm', 'Total hard tissue³ median [P5 – P95 ] mm', 'Dentin <0,5⁴ /100', 'Dentin <1⁴ /100']
    tr = []
    for t in tt:
        tr.append([f"{TYPE_NAMES[t['tooth_type']]} ({t['teeth_n']})", f"{t['reduction_mm']:.1f}", f"{t['proxy_remaining_median_mm']:.2f} [{t['proxy_remaining_P5_mm']:.2f}–{t['proxy_remaining_P95_mm']:.2f}]", f"{t['dentin_below_0p5_model_lower_per100']:.1f}–100", f"{t['dentin_below_1_model_lower_per100']:.1f}–100"])
    regional = []
    for t in gt:
        if t['region'] == 'ALL':
            continue
        regional.append([t['region'], GUIDE_NAMES[t['guide']], t['sites_n'], f"{t['two_mm_accept_per100']:.1f}", f"{t['scenario_changed_per100']:.1f} [{t['scenario_changed_CI95_low']:.1f}–{t['scenario_changed_CI95_high']:.1f}]", f"{t['whole_body_changed_per100']:.1f}", 'UNKNOWN'])
    rh = ['Region', 'Guide', 'Sites', '2 mm accepterar /100', 'Changed apex scenario / 100 [ KI ]', 'Changed whole body limit / 100', 'Regional clinical requirement']
    guide_notes = '¹ Research scenario: 95% probability of ≥ 1 mm local apex margin under gamma distributed error magnitude, isotrope direction, plane channel and X8 : s adopted imagejitter. The comparison with 2 mm uses the entire cylinder label spacing as a simplified local proxy. It does not certify body, drill tracks or real anatomy. 0 – 0,7 mm is a deterministic sensitivity to extra edge budget, not a confidence interval.\n\n² Sufficient conditions for the entire implant body: union/Cantelli-boundary with published random moments treated as exact population moments, radius 2 mm and B=0,3 mm . The tape shows B=0 – 0,7 mm . These assumptions are unvalidated. No empirical 95% requirement or statistical interval for the required clinical distance can be identified. The tip/transition and surgical track of the drill require a separate budget.\n\nThe same population profile is transferred to all regions. Regional differences in changed decisions come from the location geometry; regional guide error is Unknown . Dynamic navigation comes from mixed upper/lower jaw. KI for number of changed scenario labels is 512 case-clusters, conditional on the model; model and anatomy errors are not included. 95% within the scenario model is a different claim than 95% confidence in a measured population.'
    tooth_notes = '³ The difference from the CT vertical outer tooth boundary to the highest pulpavoxlar, minus a uniform local depth reduction. This is a total hard tissue proxy (enamel + dentin ), with crown direction selected by cross-sectional area; no expertly marked cusp/horn or actual preparation surface . P5 – P95 is observed tooth variation, not KI. Front teeth does not have the applicable occlusal definition and is only included as marked proxyes in the raw table.\n\n⁴ Mathematically possible fractional range in the column model, with digital two boundary band ± 0,520 mm and unknown enamel limit. If even the maximum residual hard tissue is < the threshold, the dentin cannot reach the threshold. The other teeth can be on either side: the upper limit is therefore 100% . This is not a measured clinical prevalence. With unknown physical segmentation error, the clinical proportion is completely UNKNOWN. Raw data also contains sensitivity ± 0,15/0,30 mm per boundary .\n\n0,5 mm is Murray 2003 : s histological criteria in young teeth after class V-preparation . 1 mm is Camps 2000 : s study boundary between dentin thickness groups, no universal protection guarantee. Bacteria, age, preparation site and material affect the transfer to occlusive crowns. No biological risk ratio is calculated.'
    practice = "Fixed reduction per material is the comparison's information-poor starting point: 1,0 mm for metal in the published preparation overview, 1,5 – 2,0 mm as PFM scenarios and 2,0 mm as conventional posterial ceramic scenario. The sources state product/surfaced rules and variation; no single depth applies to all modern products. They are scenarios, not material choices or a clinical recommendation. The table shows what existing individual geometry can query if in addition to a fixed depth. Number of real teeth with too little dentin at these depths is still UNKNOWN ."
    sources = '- [Varga 2020 , DOI 10.1111/clr.13578 ](https://doi.org/10.1111/clr.13578), published T4 MANDIBLE / accepted manuscript T3 : observed entry/apex/angle, mean and individual SD . Source values are read from saved manuscript; signed margin and regional 95% -quantile missing.\n- [Wu 2020 , DOI 10.1186 /s 40729 -020 -00272 -0 ](https://pmc.ncbi.nlm.nih.gov/articles/PMC7683639/), Results/ Fig4 : dynamic navigation, same quantities ; other population.\n- [Murray 2003 , DOI 10.1046/j.0143-2885.2003.00609.x ]( https://research.birmingham.ac.uk/en/publications/remaining-dentine-thickness-and-human-pulp-responses/ : measured dentin and histology; 0,5 mm criterion in the given regimen.\n- [Camps 2000 , DOI 10.1016 / S0109 -5641 ( 00 ) 00041 -5 ](https://pubmed.ncbi.nlm.nih.gov/10967193/): class V in premolars, thickness groups at 0,5 and 1 mm ; bacterial contamination affects response.\n- [Published preparation overview](https://pmc.ncbi.nlm.nih.gov/articles/PMC4755680/) and [National Dental PBRN ](https://pmc.ncbi.nlm.nih.gov/articles/PMC6283672/): fixed reduction depths and varying practices.\n- [ToothFairy2 ](https://toothfairy2.grand-challenge.org/dataset/) and [Pulpy3D publication](https://doi.org/10.1007/978-3-031-72111-3_2): published labels, no paired post-preparation measurements.\n- [NIST/Wilk tolerance](https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/tolelimi.htm): independent mathematical reference for the proposed signed loss measurement.'
    npost = result['attrition']['pulp']['posterior_teeth']
    nfront = result['attrition']['pulp']['anterior_out_of_occlusal_scope']
    brief = f"""# Clinical Geometric Issues from Existing Data\n\n**Answer:** four guide system and three reduction depth can now be asked with a command, with site and tooth lines behind two printable tables. External reference are published guide errors and separate tooth /pulpaannotations. Reference validates inputs and geometry on the labels; it does not validate real nerve margin or dentin after preparation. Both clinical finalities are **Unknown**.\n\nThe wizard profile can change many scenario labels compared to the 2 mm rule. For reduction we can count remaining total hard tissue and a possible dentin interval. The enamel limit is missing, so hard tissue must not be rewritten to dentin. This is a research basis for a prosthetic researcher, not a clinical recommendation.\n\nRun `./run_all.sh` from this directory. Python 3 with local NumPy/SciPy/Matplotlib is required. No downloads, package installations, GPU steps or full dataset extractions. The input is a small hash-bound copy of already local X8/X12 - result . [Print: Both Tables and Regional Details](PRINT_TABLES.html); PDF tables are available in `tables/`.\n\n## Nerve margin : goal 95% with at least 1 mm\n\n{mdtable(gh, gr)}\n\n{guide_notes}\n\n[Regional table](tables/GUIDE_MARGIN.csv) and [individual sites ](raw/GUIDE_PER_SITE.csv). The figure shows how scenario and conservative condition differ: [guide_margins.pdf ](figures/guide_margins.pdf).\n\n## Hard tissue after reduction; dentin is partially identified\n\n{mdtable(th, tr)}\n\n{tooth_notes}\n\n{practice}\n\n{npost} posterior teeth included, from 3 449 non-cut lower jaw teeth from 282 case. {nfront} front teeth is outside the anatomical application of the occlusive issue. Three posterior teeth lacks column measurement and is rejected for this query (2 premolar 2 , 1 molar 3 ); they were never zero-imputed. Missing/rejected cases are not included in the denominator of the table. The selection does not represent a clinical population. [All tooth lines and bugband](raw/REDUCTION_PER_TOOTH.csv); [median fall-clustered KI ](tables/REDUCTION_TISSUE.csv).\n\n## What does not holds and smallest next measurement\n\nThere is no matched physically reference for our finals. Azim fossa-chamber roof and Khojastepour cusp-horn rejected as error measurement operator; no residual measurement or anatomical approval is calculated. X8 : s label-picture edge-jitter does not give any boom against real channel. X12 : s digital tapes certify no physical segmentation. Guide profiles and dentin thresholds have second populations other than our virtual sites/teeth. All physical conclusions are PENDING_INDEPENDENT_REVIEW .\n\nThe next decisive design adds observation: registered planned/measured minimum gap loss per guide and region; blind expert marked tooth shaft, horns and DEJ on the same specim before and after defined reduction. [MEASUREMENT_PLAN.md ](MEASUREMENT_PLAN.md) indicates port, number and measurement contract. [ FROZEN_PREDICTIONS.json ](FROZEN_PREDICTIONS.json) binds existing proxy predictions before such a measurement ; measurements has not been performed. 0,5 mm threshold is close to the image resolution and error bands: additional decimals give no new anatomy.\n\n## Reference , verification and license\n\n{sources}\n\nEqually informed controls : separate integral of gamma distribution, skalar union/Cantelli, actual HiGHS -extremization of dentin quantities and binomial/Wilks. Numerical answers are consistent; no algorithm superiorness is claimed. Each numerical control displays a fault injection that folds, in `raw/CONTROLS.json` and `raw/MEASUREMENT_CONTROLS.json`. These are verifications within specified models, not physical validation.\n\nLicens: ToothFairy2 CC BY-SA 4.0. The Pulpy3D paper indicates the original ToothFairy as CC BY - SA , but the explicit license for the pulp expansion is UNKNOWN according to X12 : s source control. This prevents external distribution without a complementary license, and does not imply any new right to publish data. No Bits2Bites, Teeth3DS, STS or mandible defect data are used. The source texts are used local for source verification; each copyright follows the source. No direct identification patient fields are used.\n"""
    (ROOT / 'README_DEMO.md').write_text(brief)
    (ROOT / 'README_CLINICIAN.md').write_text(brief)
    (ROOT / 'tables/GUIDE_MARGIN.md').write_text('# Conditional nerve margin issue\n\n' + mdtable(gh, gr) + '\n\n' + guide_notes + '\n\n' + mdtable(rh, regional) + '\n')
    (ROOT / 'tables/REDUCTION_TISSUE.md').write_text('# Reduction and hard tissue\n\n' + mdtable(th, tr) + '\n\n' + tooth_notes + '\n\n' + practice + '\n')
    style = 'body{font:14px Arial;max-width:1150px;margin:24px auto;color:#12242f}table{border-collapse:collapse;width:100%;margin:18px 0}td,th{border:1px solid #999;padding:7px;text-align:left}th{background:#e8edf0}p{line-height:1.45}.page{break-before:page}@media print{body{margin:0;font-size:10pt}table{font-size:9pt}tr{break-inside:avoid}thead{display:table-header-group}}'
    doc = '<!doctype html><html lang="en"><meta charset="UTF-8"><title>X31 printable research tables</title><style>' + style + '</style> <body> <h1> Nerve margin — conditional research response </h1> <p> Clinical 95% requirement: UNKNOWN . All distance in mm . Geometric margin ≠  nerve injury risk . </p> ' + htmltable(gh, gr) + ''.join(('<p>' + html.escape(p) + '</p>' for p in guide_notes.split('\n\n'))) + '<h2> Regional scenario changes </h2>' + htmltable(rh, regional) + '<h1 class="page"> Reduction — total hard tissue, not measured dentin </h1>' + htmltable(th, tr) + ''.join(('<p>' + html.escape(p) + '</p>' for p in (tooth_notes + '\n\n' + practice).split('\n\n'))) + '<p> Sources : Varga DOI10.1111/clr.13578 ; Wu DOI10.1186/s40729-020-00272-0; Murray DOI10.1046/j.0143-2885.2003.00609.x; Camps DOI10.1016/S0109-5641(00)00041-5; TF2/Pulpy3D. Full Contracts/sources : README_DEMO.md . PENDING_INDEPENDENT_REVIEW.</p></body></html>'
    (ROOT / 'PRINT_TABLES.html').write_text(doc)
    for (name, title, headers, data, foot, sz) in [('GUIDE_MARGIN', 'Nerve margin : conditional research response', gh, gr, 'Clinical 95% requirement UNKNOWN. ¹ Apex gamma/isotropy + image jitter. ² Whole body, moment closure, B=0,3.\nB-band 0 – 0,7 is sensitivity, not KI . KI for changed labels = case cluster, model condition.\nSources : Varga DOI10.1111/clr.13578 T4 ; Wu DOI10.1186/s40729-020-00272-0 Fig4.\nRegional profiles and true channel boundary are missing. No clinical recommendation.\nResolution : margin scenarios PHENOMENOLOGICAL ; shares POPULATION ; location geometry PER_TOOTH .', (11.7, 8.3)), ('REDUCTION_TISSUE', 'Reduction: hard tissue proxy and possible dentin interval', th, tr, '³ Total hard tissue along CT - column minus uniform depth; no expert cusp/ DEJ / preparation .\n⁴ Possible fractional range under digital band ± 0,520 mm and unknown enamel . Physical share UNKNOWN .\n0,5 mm: Murray DOI10.1046/j.0143-2885.2003.00609.x; 1 mm: Camps study bin DOI10.1016/S0109-5641(00)00041-5.\nP5 – P95 is spread. Poor jaw, selected selection. No clinical recommendation.\nResolution : column minimum PER_TOOTH ; Table shares/medians POPULATION ; reduction is an adopted scenario.', (11.7, 8.3))]:
        (fig, ax) = plt.subplots(figsize=sz)
        ax.axis('off')
        fig.text(0.07, 0.93, title, fontsize=15, va='top')
        fig.text(0.07, 0.88, 'PENDING_INDEPENDENT_REVIEW', fontsize=10, va='top')
        box = [0, 0.5, 1, 0.3] if len(data) < 6 else [0, 0.17, 1, 0.64]
        tab = ax.table(cellText=data, colLabels=[x.replace(' /100', '\n/100').replace(' [', '\n[') for x in headers], bbox=box, cellLoc='left')
        tab.auto_set_font_size(False)
        tab.set_fontsize(8)
        for ((i, j), cell) in tab.get_celld().items():
            if i == 0:
                cell.set_facecolor('#dfe8ed')
        fig.text(0.07, 0.06, foot, fontsize=9, va='bottom')
        fig.subplots_adjust(top=0.88, bottom=0.1, left=0.07, right=0.93)
        fig.savefig(ROOT / f'tables/{name}.pdf')
        plt.close(fig)
    (fig, axes) = plt.subplots(1, 2, figsize=(11.5, 4.7), layout='constrained')
    y = np.arange(4)
    for (i, g) in enumerate(GUIDE_ORDER):
        p = result['guide_parameters'][g]
        sb = p['scenario_boundary_sensitivity_range_mm']
        bb = p['whole_body_boundary_sensitivity_range_mm']
        axes[0].plot(sb, [i - 0.12] * 2, c='#147d92', lw=5)
        axes[0].scatter(p['scenario_apex_required_gap_mm'], i - 0.12, c='#147d92', s=35)
        axes[0].plot(bb, [i + 0.12] * 2, c='#a05a2c', lw=5)
        axes[0].scatter(p['whole_body_sufficient_gap_B0p3_mm'], i + 0.12, c='#a05a2c', s=35)
        t = allg[i]
        axes[1].errorbar(t['scenario_changed_per100'], i, xerr=[[t['scenario_changed_per100'] - t['scenario_changed_CI95_low']], [t['scenario_changed_CI95_high'] - t['scenario_changed_per100']]], fmt='o', c='#147d92')
    for ax in axes:
        ax.set_yticks(y, [GUIDE_NAMES[g] for g in GUIDE_ORDER])
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.25)
    axes[0].axvline(2, c='black', ls='--', label='2 mm')
    axes[0].set_xlabel('Planerad etikettmarginal (mm)')
    axes[0].set_title('Turquoise: local apex; brown: whole body\nBand: B sensitivity 0 – 0,7 mm')
    axes[0].legend()
    axes[1].set_xlabel('Changed compared to 2 mm, per 100 sites')
    axes[1].set_title('Apex scenario, 95% case-clustered CI\n2 581 virtual sites ; clinical requirement UNKNOWN')
    for ext in ['png', 'pdf']:
        fig.savefig(ROOT / f'figures/guide_margins.{ext}', dpi=170)
    plt.close(fig)
    (fig, ax) = plt.subplots(figsize=(9.5, 4.7), layout='constrained')
    types = ['premolar1', 'premolar2', 'molar1', 'molar2', 'molar3']
    colors = ['#136b85', '#539e8c', '#c1843d']
    for (i, r) in enumerate([1, 1.5, 2]):
        rs = [next((t for t in tt if t['tooth_type'] == kind and t['reduction_mm'] == r)) for kind in types]
        vals = np.array([t['proxy_remaining_median_mm'] for t in rs])
        low = np.array([t['proxy_remaining_P5_mm'] for t in rs])
        high = np.array([t['proxy_remaining_P95_mm'] for t in rs])
        ax.errorbar(np.arange(5) + (i - 1) * 0.2, vals, yerr=[vals - low, high - vals], fmt='o', c=colors[i], label=f'{r:.1f} mm reduktion')
    ax.set_xticks(range(5), [TYPE_NAMES[k] for k in types])
    ax.set_ylabel('Remaining hard tissue : CT -column proxy ( mm )')
    ax.set_title('Total hard tissue after uniform column reduction\nP5 – P95 = observed variation; dentin and physical precision Unknown')
    ax.legend()
    ax.grid(axis='y', alpha=0.25)
    for ext in ['png', 'pdf']:
        fig.savefig(ROOT / f'figures/reduction_tissue.{ext}', dpi=170)
    plt.close(fig)
    g0 = allg[0]
    findings = f'# X31 : queries if margin and reduction becomes executable\n\nFour guide profiles and three occlusive reduction scenarios are linked to existing site and tooth lines. Published guide units were re-read against original source, and X12 : s annotations are re-used with frozen physical scale. Physical 95% -margin per region and real dentin share is still UNKNOWN : signed losses, real channel limit and DEJ missing.\n\n{mdtable(gh, gr)}\n\n{mdtable(th, tr)}\n\nR1 rejected identification: radial magnitudes are compatible with opposite directions, and the same hard tissue mask with opposite dentin decisions. R2 changed the representation to multiple ranges and conditional inverse query. R3 provides a measuring port; no physical measurement has been performed. All PREREG , negative gates and source snapshots are preserved.\n\nTables use {len(gt)} guide/region lines and {len(tt)} Postior tooth type / Reduction Rows; raw queries is available in raw/. Full interpretation limits, practice, license and sources are in README_DEMO.md . B-bands are not statistically KI, and P5 – P95 is not KI. Case-clustered resampling provides only model-conditioned selection uncertainty. Numerical controls provides consistent answers to the same information; no algorithm superior is claimed.\n\nNo missing measurements was replaced by invented distribution parameters. External anatomic reference for fossa/roof and cusp/horn are rejected on the wrong operator; no anatomical validation is claimed. The biggest obstacle is matched measurement, not the solution algorithm. Status PENDING_INDEPENDENT_REVIEW. Cost, loss and resolution per quantity are available in results.json.\n'
    (ROOT / 'RESULTS.md').write_text(findings)
