from common import *
import collections
import matplotlib.pyplot as plt
FAMS = ['molar_crown', 'premolar_crown', 'anterior_crown']
LABELS = {'molar_crown': 'Molar', 'premolar_crown': 'Premolar', 'anterior_crown': 'Anterior'}

def summarize(rows, tag, fam):
    a = [r for r in rows if r['tag'] == tag and r['family'] == fam]
    s = [r for r in a if r['status'] == 'RESCORED']
    return dict(requested=len(a), scored=len(s), old_anatomy_pass=sum((r['old_p95_mm'] <= 0.35 for r in s)), shape_envelope_pass=sum((bool(r['shape_envelope_pass']) for r in s)), function_nominal=sum((r['nominal_function_pass'] for r in s)), function_cover=sum((r['cover_function_pass'] for r in s)), shape_and_nominal=sum((r['shape_and_nominal_function'] for r in s)), shape_and_cover=sum((r['shape_and_cover_function'] for r in s)), positive_native_nominal=sum((r['nominal_function_pass'] and r['positive_native_band'] for r in s)), median_shape_p95_mm=float(np.median([r['oracle_aligned']['8192']['p95_mm'] for r in s])) if s else None, median_shape_ratio=float(np.median([r['shape_ratio'] for r in s])) if s else None, clinical_status='NOT_ESTABLISHED')

def make_figure(pairs, rows):
    (fig, axs) = plt.subplots(1, 3, figsize=(14, 4.5), gridspec_kw={'width_ratios': [1, 1.6, 0.85]})
    for (i, f) in enumerate(FAMS):
        vals = [r['metrics']['8192']['p95_mm'] for r in pairs if r['family'] == f]
        axs[0].scatter(i + np.linspace(-0.12, 0.12, len(vals)), vals, color=['#245c9f', '#c57918', '#498743'][i], s=38)
        axs[0].plot([i - 0.2, i + 0.2], [max(vals)] * 2, c='black', lw=1.5)
    axs[0].axhline(0.35, color='#af3636', ls='--', label='Old 0.35 mm')
    axs[0].set_xticks(range(3), ['Molar', 'Premolar', 'Anterior'])
    axs[0].set_ylabel('Bilateral p95 (mm), PER_TOOTH')
    axs[0].set_title('Annotated natural teeth\n6 pairs/type; not accepted crowns')
    axs[0].legend(fontsize=8)
    tags = ['R2', 'A', 'B', 'C']
    for (i, t) in enumerate(tags):
        rr = [r for r in rows if r['tag'] == t and r['status'] == 'RESCORED']
        for (j, r) in enumerate(rr):
            axs[1].scatter(i + (j - (len(rr) - 1) / 2) * 0.025, r['shape_ratio'], c=['#245c9f', '#c57918', '#498743'][FAMS.index(r['family'])], marker='o' if r['nominal_function_pass'] else 'x', s=35)
    axs[1].axhline(1, color='black', ls='--')
    axs[1].set_yscale('log')
    axs[1].set_xticks(range(4), ['R2 prior', 'R3 A', 'R3 B', 'R3 C'])
    axs[1].set_ylabel('Oracle shape p95 / bilateral envelope')
    axs[1].set_title('Existing designs, PER_TOOTH\nCircle: nominal function; x: failed')
    axs[2].bar([0, 1], [1, 0], color=['#498743', '#af3636'])
    axs[2].set_xticks([0, 1], ['Surface A', 'Surface B'])
    axs[2].set_ylim(0, 1.2)
    axs[2].set_ylabel('Contact area (mm²)')
    axs[2].set_title('Exact equal-summary witness\np95 = 0.25 mm in both')
    for ax in axs:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=0.2)
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/diagnosis.png', dpi=180)
    fig.savefig(ROOT / 'figures/diagnosis.pdf')
    plt.close(fig)

def run():
    raw = read(ROOT / 'raw/RESCORE.json')
    rows = raw['rows']
    assert len(rows) == 90
    ap = read(ROOT / 'raw/ANNOTATED_PAIRS.json')
    lp = read(ROOT / 'raw/LOCAL_PAIRS.json')
    pr = read(ROOT / 'PREREG_RESCORE.json')
    pub = read(ROOT / 'raw/PUBLISHED_MEASUREMENTS.json')
    summary = {t: {f: summarize(rows, t, f) for f in FAMS} for t in ['R2', 'A', 'B', 'C', 'scalar']}
    pair_summary = {}
    for fam in FAMS:
        rr = [r for r in ap['rows'] if r['family'] == fam]
        v = [r['metrics']['8192']['p95_mm'] for r in rr]
        pair_summary[fam] = dict(n=len(v), median_mm=np.median(v), min_mm=min(v), max_mm=max(v), above_old035=sum((x > 0.35 for x in v)), max_probe_sensitivity_mm=max((r['probe_sensitivity_p95_mm'] for r in rr)))
    main = [r for r in rows if r['tag'] != 'scalar']
    scored = [r for r in main if r['status'] == 'RESCORED']
    res = dict(claim_type=['capability', 'information_link'], review_state='PENDING_INDEPENDENT_REVIEW', main_answer='Whole-outer p95 is insufficient for crown function; 0.35mm is not a validated clinical threshold. Independent annotated bilateral shape scale permits a descriptive screen only.', clinical_acceptable_pair_floor='UNKNOWN; no local same-preparation accepted crown pair/repeat scan identified', external_referent=dict(kind='published_dataset', locator=['https://osf.io/xctdy/', str(V4 / 'payload/whole_private')], compared_quantity='Natural bilateral triangle surface distances and native registered static gap, not clinical acceptance', refutes_us=True), published_referents=pub['records'], resolution='PER_POINT fields -> PER_SURFACE_REGION -> PER_TOOTH. Six independent case clusters per local cohort; no population extrapolation.', pair_summary=pair_summary, thresholds_mm=pr['thresholds_mm'], summary=summary, main_counts=dict(attempted=len(main), scored=len(scored), shape_pass=sum((bool(r['shape_envelope_pass']) for r in scored)), shape_and_function=sum((r['shape_and_nominal_function'] for r in scored)), clinical_pass=0, clinical_status='NOT_ESTABLISHED; zero verified is not zero clinically acceptable'), dropout=dict(main_designs=dict(total=72, rejected=72 - len(scored), fraction=(72 - len(scored)) / 72, reasons=dict(collections.Counter((r.get('source_score', {}).get('reason', r['status']) for r in pr['records'] if r['tag'] != 'scalar' and r['status'] != 'SCORED')))), annotated_pairs=dict(total=18, rejected=sum((not r['valid_for_envelope'] for r in ap['rows'])), screened_arches=6 + len(ap['screening_rejections']), rejected_arches=len(ap['screening_rejections']), clinical_calibration_rejected=18, reason='No accepted restoration labels; no actual preparation'), initial_pairs=dict(total=18, measured=18, rejected_for_clinical_floor=18, reason='inferred labels, fragmentation, no clinical paired-restoration labels'), literature=pub['screening']), sufficiency=read(ROOT / 'raw/SUFFICIENCY.json'), validation=read(ROOT / 'raw/VALIDATION.json'), component_audit=read(ROOT / 'raw/COMPONENT_AUDIT.json'), cost=dict(local_pairs_seconds=lp['seconds'], annotated_pairs_seconds=ap['seconds'], rescore_seconds=raw['seconds'], max_observed_peak_rss_MiB=max(lp['peak_rss_MiB'], ap['peak_rss_MiB'], raw['peak_rss_MiB']), preparation='Included within stage times', fit='Included multistart ICP; no learning or generation', validation='included distance/contact work; separate validation script', discovery='Historical R1-R3 compute and X11 labeling cost UNKNOWN, not omitted as zero', questions=90, queries='72 main +18 equally informed scalar diagnostics', fallback='Failed10 main designs retained', human_lab_measurement='NOT_RUN', tokens='UNKNOWN'), edges=[dict(from_quantity='annotated native surface point', to_quantity='point/triangle distance', resolution='PER_POINT', time_scale='SIMULTANEOUS'), dict(from_quantity='static gap field', to_quantity='regional contact pattern', resolution='PER_POINT', time_scale='SIMULTANEOUS'), dict(from_quantity='design outer+intaglio+finishline', to_quantity='manufactured measurement', resolution='PER_POINT', time_scale='HANDOVER', status='PROPOSED_NOT_RUN')], phenomenological_debts=[dict(quantity='contact band0.1mm; wall0.5mm; pattern1mm2; centroid0.5mm', resolution='PHENOMENOLOGICAL', replacement='Registered loaded contact/scan repeatability; material/process-specific wall coupon and crown tests; accepted designs with blinded functional labels'), dict(quantity='family maximum bilateral p95', resolution='PHENOMENOLOGICAL', replacement='Distribution of accepted restoration differences conditional on same preparation and product')], uncertainty='R3 source and homolog labels are inferred; fragmented/mislabelled targets can confound apparent shape error. No rigorous general mesh rounding/sampling enclosure; four-start ICP selects symmetric point RMS, not p95 directly; attained p95 is not a global minimum; scanner/segmentation uncertainty not separated; upper Teeth3DS to lower R3 transfer not validated', raw_files={str(p.relative_to(ROOT)): sha(p) for p in [ROOT / 'raw/LOCAL_PAIRS.json', ROOT / 'raw/ANNOTATED_PAIRS.json', ROOT / 'raw/RESCORE.json', ROOT / 'raw/PUBLISHED_MEASUREMENTS.json']}, data_artifacts={str(p): dict(sha256=sha(p), bytes=p.stat().st_size) for p in DATA.glob('*.npz')})
    dump(ROOT / 'results.json', res)
    make_figure(ap['rows'], rows)
    table = '| Method | Type | Scored/6 | Shape envelope/6 | Nominal function/6 | Positive native contact + nominal/6 | Shape+function/6 | Cover function/6 | Median shape p95 mm |\n|---|---|---:|---:|---:|---:|---:|---:|---:|\n'
    for tag in ['R2', 'A', 'B', 'C']:
        for fam in FAMS:
            s = summary[tag][fam]
            table += f"| {tag} | {LABELS[fam]} | {s['scored']} | {s['shape_envelope_pass']} | {s['function_nominal']} | {s['positive_native_nominal']} | {s['shape_and_nominal']} | {s['function_cover']} | {s['median_shape_p95_mm']:.3f} |\n"
    pt = '| Natural pair type | n | Median p95 mm | Range mm | Above0.35 mm |\n|---|---:|---:|---:|---:|\n'
    for (fam, s) in pair_summary.items():
        pt += f"| {LABELS[fam]} | 6 | {s['median_mm']:.3f} | {s['min_mm']:.3f}–{s['max_mm']:.3f} | {s['above_old035']}/6 |\n"
    density_counts = {str(n): sum((r['metrics'][str(n)]['p95_mm'] > 0.35 for r in ap['rows'])) for n in [2048, 8192]}
    res['sampling_sensitivity'] = dict(annotated_exceedances035=density_counts, max_pair_p95_difference_mm=max((r['probe_sensitivity_p95_mm'] for r in ap['rows'])), max_rescore_shape_difference_mm=max((abs(r['oracle_aligned']['8192']['p95_mm'] - r['oracle_aligned']['2048']['p95_mm']) for r in scored)), interpretation='sampling sensitivity only; not uncertainty interval or certified convergence')
    dump(ROOT / 'results.json', res)
    text = f"# The shape and function of the crown must be assessed separately\n\nIt is now possible to distinguish design deviation from placement and functional requirements of the preserved R3-designerna. **0,35 mm p95 over the entire outer surface is not a sufficient clinical quality measure.** The actual variation floor between two clinically accepted crowns on the same preparation is still UNKNOWNOur local measurement applies to mirrored natural teeth, not such crown pairs.\n\n{pt}\nMeasurement: six Teeth3DS upper jaws, dataset teeth labels, rigid registration without scaling, two-way point–triangel-p95 med 8192 area distributed centroid probes per area. All 18 pairs retained; verification verified that the component selection did not add any surfaces to 36/36 The values are: PER_TOOTHNo population limits. 8192 probes exceeding: 4/18 value 0,35 mm, med 2048 prober 3/18The biggest change is: {res['sampling_sensitivity']['max_pair_p95_difference_mm']:.3f} mm; this is test density sensitivity, no rigorous error limit. Therefore there is support for: 0,35 Sometimes natural bilateral variation is rejected, but not to arbitrarily raise a clinical limit.\n\nThe first measurement of R3:s egna 18 kontralateraler gav 17/18 over 0,35, but several have large disconnected fragments and unsafe labels. This result is preserved and not used for calibration. The difference between cohorts cannot be causally attributed to labels alone. R3-re-evaluation retains its original reference surfaces; incorrect or fragmented reference can thus contribute to the large design distances. Zero form hits do not isolate generator errors from reference errors.\n\n{table}\nThe table applies to: PER_TOOTH. Shape envelope is a **beskrivande** form check after rigid registration against natural tooth: frozen limits {pr['thresholds_mm']['molar_crown']:.3f}/{pr['thresholds_mm']['premolar_crown']:.3f}/{pr['thresholds_mm']['anterior_crown']:.3f} mm for molar/premolar/framtand, set to the largest independently measured bilateral value. Registration selects the lowest symmetric point-RMS av fyra ICP-starts, then measured p95. This is a achieved value, no lower limit possible p95. Contact is calculated in the original mode. 0,35The gate is still in the raw data. {res['main_counts']['shape_pass']}/72 is able to handle the design control and {res['main_counts']['shape_and_function']}/72 klarar form plus nominell funktion. 10/72 design failure is included in the denominator; the mediar and figure use only 62 measurable designer.\n\n![Diagnos](figures/diagnosis.png)\n\nClinical margin, strained occlusion and signed proximal contact are UNKNOWNPlan margin is not fit to a real preparation limit. Nominal function reuses 0,5 mm wall and other frozen technical limits; these are not material specific clinical requirements. The Cover column uses the older conservative facet boundary; rigorous rounding closure is missing. No line involves clinical approval.\n\nThe inherited nominal the validation gate can accept that both reference and construction lack contact. Therefore, positive contact bands are separately reported in the reference. A has 11/18 nominally approved but only 4/18 with positive reference band and nominal approval: 0 molarer, 3 premolars and 1 framtand. R2 and the scaler control has 0/18 with both characteristics; B and C have 1/18 This does not change a frozen gate and is not a new clinical contact assessment.\n\nSatisfaction test: two conditions have exactly p95=RMS=maximum distance=0,25 mm, identitetsfel 0, men geometrisk kontakt 1 respektive 0 mm². En signerad gapstorhet skiljer just det paret. Ett andra par bevarar dessutom kontaktarea 1 mm² and number of components 1, men kontaktcentrum flyttas 3 mm and symmetrical difference becomes 2 mm². There is at least one location coordinate required. The samples show why the location and signs of the deviation must accompany; they do not prove that our vector captures all clinical function.\n\nLika informerad scalar_control: The results can be found separately in results.json; the same acquired gap and the same revaluation. No algorithm advantage is claimed. Read LITERATURE.md for the exact dimensions of the primary sources, and README_DEMO.md for driving, limitations and next measurement.\n"
    (ROOT / 'RESULTS.md').write_text(text)
    state('DIAGNOSTIC_ROUNDS_COMPLETE', res['main_counts'], 'Freeze review bundle, verify one-command demo, bind graph feedback')
    return res
if __name__ == '__main__':
    run()
