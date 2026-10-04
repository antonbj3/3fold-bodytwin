from common import *
import csv, collections

def table(rounds):
    out = ["| Construction | Tooth type | p95 in mode, median mm | Natural pair, median mm | Form /6 | Marginal /6 | Wall / 6 | Complete /6 |", '|---|---|---:|---:|---:|---:|---:|---:|']
    for (tag, methods) in [('B', ['rigid_control', 'boundary_deform', 'missing_donor']), ('C', ['exact_margin', 'fixed_prep_projection']), ('D2', ['measured_neighbours'])]:
        if tag not in rounds:
            continue
        for method in methods:
            for family in ['molar', 'premolar', 'anterior']:
                r = rounds[tag]['summary'][method][family]
                g = r['gate_counts']
                out.append('| %s/%s | %s | %.3f | %.3f | %d | %d | %d | %d |' % (tag, method, family, r['median_p95_mm'], r['median_natural_floor_mm'], r['shape_pass'], g.get('margin', 0), g.get('wall_geometric_bound', 0), r['technical_joint']))
    return '\n'.join(out)

def figure(rounds):
    import matplotlib.pyplot as plt
    (ROOT / 'figures').mkdir(exist_ok=True)
    (fig, axes) = plt.subplots(1, 3, figsize=(13, 4.2))
    families = ['molar', 'premolar', 'anterior']
    x = np.arange(3)
    width = 0.18
    r = rounds['B']['summary']
    axes[0].bar(x - width * 1.5, [r['rigid_control'][f]['median_natural_floor_mm'] for f in families], width, label='Natural pair, oracle fit')
    choices = [('B', 'rigid_control', 'Collar rigid'), ('C', 'exact_margin', 'Exact margin'), ('B', 'missing_donor', 'No contralateral')]
    for (j, (tag, m, label)) in enumerate(choices):
        if tag in rounds:
            axes[0].bar(x + width * (j - 0.5), [rounds[tag]['summary'][m][f]['median_p95_mm'] for f in families], width, label=label)
    axes[0].set_xticks(x, families)
    axes[0].set_ylabel('Median p95 (mm), PER_TOOTH')
    axes[0].set_title('Original annotated exterior')
    axes[0].legend(fontsize=7)
    if 'C' in rounds:
        rr = [r for r in rounds['C']['rows'] if r['method'] == 'exact_margin' and r['status'] == 'SCORED']
        for (f, col) in zip(families, ['#287271', '#d68c45', '#6457a6']):
            ss = [r for r in rr if r['family'] == f]
            axes[1].scatter([r['natural_floor']['p95_mm'] for r in ss], [r['in_situ']['p95_mm'] for r in ss], color=col, label=f)
        axes[1].plot([0, 0.9], [0, 0.9], ':', color='gray')
        axes[1].set_xlabel('Natural pair after oracle fit (mm)')
        axes[1].set_ylabel('Generated in supplied position (mm)')
        axes[1].set_title('Registration adds an error budget')
        axes[1].legend(fontsize=8)
        labels = ['Shape', 'Closed', 'Margin', 'Wall', 'Proximal\nboth', 'Verified\nfunction']
        counts = [sum((r['gates']['shape'] for r in rr)), sum((r['gates']['closed'] for r in rr)), sum((r['gates']['margin'] for r in rr)), sum((r['gates']['wall_geometric_bound'] for r in rr)), sum((r['gates']['mesial'] and r['gates']['distal'] for r in rr)), 0]
        axes[2].bar(np.arange(len(labels)), counts, color=['#287271'] * 3 + ['#ad4f47'] * 3)
        axes[2].set_xticks(np.arange(len(labels)), labels, rotation=30, ha='right')
        axes[2].set_ylim(0, 20)
        axes[2].set_ylabel('Cases /18')
        axes[2].set_title('Exact-margin variant: separate gates')
        axes[2].text(5, 1, 'UNKNOWN', rotation=90, ha='center', fontsize=8)
    fig.suptitle('Research geometry only — virtual preparation; bite pose unverified', fontsize=11)
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/demo.png', dpi=160)
    plt.close(fig)

def run():
    rounds = {p.stem: read(p) for p in sorted((ROOT / 'rounds').glob('*.json')) if p.stem in ['A', 'B', 'C', 'D2']}
    allrows = [r for (tag, q) in rounds.items() if tag != 'A' for r in q['rows']]
    scored = [r for r in allrows if r['status'] == 'SCORED']
    valid = read(ROOT / 'raw/VALIDATION.json') if (ROOT / 'raw/VALIDATION.json').exists() else None
    costs = {tag: {k: v for (k, v) in q.items() if k in ['seconds', 'peak_rss_MiB']} for (tag, q) in rounds.items()}
    gen = {p.stem: read(p) for p in ROOT.glob('FROZEN_PREDICTIONS_*.json')}
    prep = read(ROOT / 'FROZEN_INPUTS_B.json')
    ownfiles = [p for p in DATA.rglob('*') if p.is_file()]
    disk = sum((p.stat().st_size for p in ownfiles))
    assert disk <= 3000000000
    result = dict(claim_type='capability', status='PARTIAL_EXTERIOR_RECONSTRUCTION_FULL_CAPABILITY_NOT_ESTABLISHED', review_state='PENDING_INDEPENDENT_REVIEW', external_referent=REFERENT, rounds={tag: dict(summary=q.get('summary'), dropout=q.get('dropout'), result_path=ROOT / 'rounds' / (tag + '.json'), sha256=sha(ROOT / 'rounds' / (tag + '.json'))) for (tag, q) in rounds.items()}, records=scored, thresholds=read(ROOT / 'PREREG_B.json')['metrics'], resolution=dict(distances='PER_POINT', contacts_wall='PER_SURFACE_REGION', decisions='PER_TOOTH', counts='descriptive frozen cohort, not population estimates'), time_scale='SIMULTANEOUS', limitations=['Preparation and finish-line are virtual; native annotated boundary is not a clinical margin', 'Lower jaw PCA-normalized separately: bite transform/loaded occlusion UNKNOWN', 'No physical manufactured specimen or measured fit', 'No rigorous floating point/source uncertainty enclosure', 'Mesh watertightness does not prove no self-intersection or insertability', 'Adaptive rounds on same frozen18 teeth; no held-out generalization claim for C/D', 'The shape donor control had same inputs; no new algorithm-superiority claim', 'Changing annotated dataset prevents causal attribution of improvement over R1-R3 to algorithm alone'], dropout=dict(A_planar_rejections=18, post_A_requested=len(allrows), post_A_scored=len(scored), post_A_rejected=len(allrows) - len(scored), fraction=(len(allrows) - len(scored)) / len(allrows) if allrows else None), cost=dict(preparation_seconds=prep['seconds'], generation={k: dict(seconds=v['seconds'], peak_rss_MiB=v['peak_rss_MiB']) for (k, v) in gen.items()}, scoring=costs, validation_seconds=valid['seconds'] if valid else None, peak_rss_MiB=max([prep['peak_rss_MiB']] + [v['peak_rss_MiB'] for v in gen.values()] + [q.get('peak_rss_MiB', 0) for q in rounds.values()]), threads=4, address_space_cap_MiB=3500, own_array_bytes=disk, raw_acquisition_annotation_cost='UNKNOWN historical external cost', cpu_time_seconds='NOT_MEASURED; elapsed seconds per phase recorded', fit='included in generation', discovery='A/B/C/D preserved; see ATTEMPTS.json', query_count=len(allrows), fallback='six same-FDI templates per missing-donor query; costs included'), validation=valid, sufficiency=valid['sufficiency'] if valid else None, large_arrays=[dict(path=p, bytes=p.stat().st_size, sha256=sha(p)) for p in ownfiles if p.stat().st_size > 50000000], source_locks={p.name: sha(p) for p in ROOT.glob('FROZEN_*.json')}, graph_coverage=read(ROOT / 'GRAPH_COVERAGE_PROPOSAL.json'))
    if (ROOT / 'FROZEN_EXPORTS.json').exists():
        result['frozen_exports'] = read(ROOT / 'FROZEN_EXPORTS.json')
    if (ROOT / 'raw/DEMO_REPLAY.json').exists():
        result['demo_replay'] = read(ROOT / 'raw/DEMO_REPLAY.json')
    if (ROOT / 'raw/PREPARATION_ORACLE.json').exists():
        result['preparation_oracle'] = read(ROOT / 'raw/PREPARATION_ORACLE.json')
    if (ROOT / 'raw/FINAL_INTEGRITY.json').exists():
        result['final_integrity'] = dict(path=ROOT / 'raw/FINAL_INTEGRITY.json', sha256=sha(ROOT / 'raw/FINAL_INTEGRITY.json'), checked_mesh_records=108)
    result['cost']['failed_interrupted_scorer_seconds'] = 'UNKNOWN additional time before V1 exception; successful row times retained. No total-cost superiority claim.'
    result['cost']['cpu_time_seconds'] = 'UNKNOWN; per-phase elapsed and peak RSS measured. Preparation-oracle and demo elapsed are recorded in their named result objects.'
    for name in ['demo_replay', 'preparation_oracle']:
        if name in result:
            result['cost']['peak_rss_MiB'] = max(result['cost']['peak_rss_MiB'], result[name]['peak_rss_MiB'])
    result['preserved_implementation_failures'] = dict(D_original_generated=16, D_original_rejected=2, D_original_rejection_fraction=2 / 18, reason='Triangle-soup reconstruction welded distinct near-coincident boundary knots. Original outputs/errors retained; D2 preserves indexed geometry.', scorer_v1='Stopped after30 C scores; history/score_c_v1.log retained. V2 preserves original index topology and completes same frozen set.', validation_v1='Exact centroid==3 assertion failed at3.0000000000000027. Same-summary identity still exactly0; controlled validation correction retained.')
    dump(ROOT / 'results.json', result)
    fields = ['round', 'method', 'key', 'family', 'p95_mm', 'registered_p95_mm', 'natural_floor_mm', 'excess_mm', 'margin_mm', 'wall_sample_mm', 'wall_lower_mm', 'proximal_mesial_mm2', 'proximal_distal_mm2', 'contact_symdiff_mm2', 'contact_count_error', 'contact_centroid_mm', 'contact_negative_area_mm2', 'shape_pass', 'technical_joint']
    with (ROOT / 'PER_TOOTH.csv').open('w') as file:
        w = csv.DictWriter(file, fieldnames=fields)
        w.writeheader()
        for (tag, q) in rounds.items():
            if tag == 'A':
                continue
            for r in q['rows']:
                if r['status'] != 'SCORED':
                    continue
                w.writerow(dict(round=tag, method=r['method'], key=r['key'], family=r['family'], p95_mm=r['in_situ']['p95_mm'], registered_p95_mm=r['registered_shape']['p95_mm'], natural_floor_mm=r['natural_floor']['p95_mm'], excess_mm=r['in_situ_excess_over_floor_mm'], margin_mm=r['margin_curve_sampled_max_mm'], wall_sample_mm=r['wall']['sampled_min_mm'], wall_lower_mm=r['wall']['continuous_lower_mm'], proximal_mesial_mm2=r['proximal']['mesial'].get('area_mm2'), proximal_distal_mm2=r['proximal']['distal'].get('area_mm2'), contact_symdiff_mm2=r['contact'].get('symdiff_mm2'), contact_count_error=r['contact'].get('count_error'), contact_centroid_mm=r['contact'].get('centroid_distance_mm'), contact_negative_area_mm2=r['contact'].get('negative_gap_area_mm2'), shape_pass=r['gates']['shape'], technical_joint=r['technical_joint']))
    text = "# Natural outer anatomy can be maintained; a complete working crown is not yet shown\n\nOriginally listed contralaterals provide measurable recovery of six frozen jaws, three tooth types per jaw. The reference is Teeth3DS original triangles; no clinical acceptance measure is established. Each millimetre value below is a PER_TOOTH-mThe median only sums up the six cases. DIAG:s descriptive independent calibration levels (rounded off in spring PREREG), not clinical boundaries.\n\n" + table(rounds) + "\n\nComplete in the table means the technical connection form+closed shell+margin+geometric wall boundary, not clinical function. Bet registration and actual preparation are missing; the entire ability remains UNKNOWN/unestablished. GenCAD v6:s contact area, component number, position, symmetrical difference and negative gap area stand per tooth in PER_TOOTH.csv and results.json, only in the delivered roe pose.\n\n![Resultat](figures/demo.png)\n\nAlla 18 Plana preparations in A fell; no cases were replaced. B retains the entire original surface but approximates the margin curve. C inserts all curve nodes and maintains exactly the linear limit; the extra wall projection is tested as a separate candidate. D examines measured neighboring surfaces and refrains from antagonist change without verified bet pose. p95/RMS/max=0,25 mm and same contact area=1 mm² samt komponentantal=1; nevertheless, the contact center moves about3 mm and symmetrical difference is2 mm². The identity error is exact0. A mode coordinate separates these two states.\n\nEqually informed control: regular rigid neck registration is actually done, with the same data and preparation. No algorithm advantage is claimed. The improved comparability against R1–R3 is not a paired causal test: data cohort and annotation have changed. See LITERATURE.md, COMMANDS.md, all frozen files and HANDOFF.md.\n"
    if 'preparation_oracle' in result:
        text += "\nOracle-control of preparation: also the original measured exterior triangles miss the virtual wall requirement in %d/ 18 case. The wall misses therefore do not isolate generator errors. The control does not change the preparation and does not give a p95-opportunity kit.\n" % result['preparation_oracle']['source_sampled_failure']
    text += "\nImplementation errors preserved: first D run gav16/18 filer och2/18 interruption after an auxiliary function welded together almost congruent margin knots. D2 preserves the original's indexed triangles. The same correction allowed the C assessment to continue after30 rows; frozen dimensions and limits were not changed.\n"
    (ROOT / 'RESULTS.md').write_text(text)
    figure(rounds)
    return result
if __name__ == '__main__':
    run()
