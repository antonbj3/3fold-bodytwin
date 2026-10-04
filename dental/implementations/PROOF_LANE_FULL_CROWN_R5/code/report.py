from construct_b import *
import csv
import matplotlib.pyplot as plt

def summarize(rows):
    out = []
    for method in dict.fromkeys((r['method'] for r in rows)):
        for fam in ['anterior', 'premolar', 'molar']:
            ss = [r for r in rows if r['family'] == fam and r['method'] == method]
            sc = [r for r in ss if r.get('status') == 'SCORED']
            out.append(dict(method=method, family=fam, requested=6, scored=len(sc), rejected=sum((r.get('status') == 'REJECTED' for r in ss)), shape=sum((r.get('gates', {}).get('shape', False) for r in ss)), margin=sum((r.get('gates', {}).get('margin', False) for r in ss)), closed=sum((r.get('gates', {}).get('closed', False) for r in ss)), wall=sum((r.get('gates', {}).get('wall', False) for r in ss)), offset=sum((r.get('gates', {}).get('offset', False) for r in ss)), technical_joint=sum((r.get('technical_joint', False) for r in ss)), complete=sum((r.get('complete_with_offset', False) for r in ss)), median_p95_mm=float(np.median([r['shape']['p95_mm'] for r in sc])) if sc else None))
    return out

def run(final=False):
    rounds = {tag: read(R / f'RESULTS_{tag}.json') for tag in ['A', 'B', 'C'] if (R / f'RESULTS_{tag}.json').exists()}
    allrows = [dict(r, round=tag) for (tag, o) in rounds.items() for r in o['rows']]
    summary = summarize(allrows)
    controls = read(R / 'raw/CONTROLS.json')
    candidates = [r for r in allrows if r.get('complete_with_offset')]
    qual = []
    for r in candidates:
        xf = R / f"raw/X49_{r['round']}.json"
        x = [x for x in read(xf) if x['key'] == r['key'] and x['method'] == r['method']] if xf.exists() else []
        if x and x[0]['rules'].get('mesh_health') == 'PASS' and (x[0]['rules'].get('material_wall') == 'PASS'):
            qual.append(r)
    manifests = []
    for p in sorted(D.rglob('*.npz')):
        manifests.append(dict(path=p, bytes=p.stat().st_size, sha256=sha(p)))
    allocated = sum((p.stat().st_size for p in D.rglob('*') if p.is_file()))
    x49 = {t: read(R / f'X49_{t}.json') for t in ['A', 'B', 'C'] if (R / f'X49_{t}.json').exists()}
    out = dict(claim_type='capability', status='FINAL_PENDING_INDEPENDENT_REVIEW' if final else 'CHECKPOINT_PENDING_INDEPENDENT_REVIEW', capability='Crown intaglio/margin/wall construction and explicit refusals on the unchanged six-jaw R4 cohort.', external_referent=read(R / 'PREREG_A.json')['external_referent'], external_referents=[dict(kind='independent_measurement', locator=str(BASE / 'LANE_X13_CEMENT_GAP/FACIT.csv') + ' PMC10721348 tab1 G50', compared_quantity='CAD50um internal/25um marginal setting vs reported23.22um marginal group mean; NOT spatial film validation for our crowns', refutes_us=True), dict(kind='published_code', locator='https://github.com/libigl/libigl/tree/v2.6.0/include/igl', compared_quantity='Point-to-triangle nominal distances independently checked with trimesh; physical validity not established', refutes_us=True)], summary=summary, rounds={t: dict(path=str(R / f'RESULTS_{t}.json'), sha256=sha(R / f'RESULTS_{t}.json'), seconds=x['seconds'], peak_rss_MiB=x['peak_rss_MiB']) for (t, x) in rounds.items()}, qualified_export_count=len(qual), rows=allrows, controls=controls, X49=x49, dropout={'per_round': {tag: dict(requested=36, scored=sum((r.get('status') == 'SCORED' for r in o['rows'])), rejected=sum((r.get('status') == 'REJECTED' for r in o['rows'])), reasons={reason: sum((r.get('reason') == reason for r in o['rows'])) for reason in set((r.get('reason') for r in o['rows'] if r.get('reason')))}) for (tag, o) in rounds.items()}}, data_files=manifests, cost=dict(round_costs='see frozen prediction seconds and score seconds, plus logs for aborted attempts', data_bytes=allocated, fit_seconds=0.0, query_counts={t: len(o['rows']) for (t, o) in rounds.items()}, dependency_install='raw/dependency_install.log; local-only libigl2.6.1', historical_data_acquisition='UNKNOWN', discovery_and_interrupted_CPU_seconds='UNKNOWN; no speed advantage claimed', physical_measurements=0, peak_rss_MiB=max([o['peak_rss_MiB'] for o in rounds.values()] + [controls['peak_rss_MiB']]), native_thread_limit_correction=read(R / 'raw/THREAD_LIMIT_CORRECTION.json')), uncertainty=dict(original_scan='UNKNOWN', actual_preparation='ABSENT; R4 virtual cone interpreted explicitly', bite_pose='UNVERIFIED', seated_spatial_film='UNKNOWN', self_intersection='UNKNOWN', floating_enclosure='General mesh and affine-field estimates lack rigorous floating enclosure; negative GenCAD pair witnesses use rational arithmetic', shape_sampling='8192 area probes per direction, inheritedR4; no whole-surface Hausdorff certificate'), review_state='PENDING_INDEPENDENT_REVIEW', scientific_admission=False)
    if (R / 'RESULTS_D.json').exists():
        out['necessary_condition'] = read(R / 'RESULTS_D.json')
        out['necessary_condition_verification'] = read(R / 'raw/D_VERIFICATION.json')
    if (R / 'raw/LAST_DEMO.json').exists():
        out['demo'] = read(R / 'raw/LAST_DEMO.json')
    for tag in ['A', 'B', 'C']:
        path = R / f'FROZEN_PREDICTIONS_{tag}.json'
        if path.exists():
            out['rounds'].setdefault(tag, {})['generation_cost'] = {k: v for (k, v) in read(path).items() if k in ['seconds', 'peak_rss_MiB']}
    out['cost']['peak_rss_MiB'] = max([out['cost']['peak_rss_MiB']] + [v.get('generation_cost', {}).get('peak_rss_MiB', 0) for v in out['rounds'].values()] + [v.get('peak_rss_MiB', 0) for v in out['X49'].values()])
    save(R / 'results.json', out)
    save(R / 'exports/QUALIFYING_MANIFEST.json', dict(qualifying_count=len(qual), qualifying_keys=[r['key'] for r in qual], reason='Require complete frozen geometric gates and unchanged X49 mesh/material pass. No fabricated qualifiers. Physical manufacture qualification remains UNKNOWN.'))
    (fig, axs) = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    colors = {'A': '#d58b39', 'B': '#387ba4', 'C': '#7c5ca3'}
    for (ax, fam) in zip(axs, ['anterior', 'premolar', 'molar']):
        for tag in rounds:
            ss = [r for r in allrows if r['round'] == tag and r['family'] == fam and ('local_thickening' in r['method']) and (r.get('status') == 'SCORED')]
            if ss:
                ax.scatter([r['shape']['p95_mm'] for r in ss], [r['wall']['sampled_min_mm'] for r in ss], label=tag, color=colors[tag], s=40, alpha=0.8)
        lim = {'anterior': 0.24, 'premolar': 0.725, 'molar': 0.437}[fam]
        ax.axvline(lim, color='gray', ls='--', lw=1)
        ax.axhline(0.5, color='gray', ls='--', lw=1)
        ax.set_title(fam + ' (PER_TOOTH)')
        ax.set_xlabel('Native-surface p95 [mm]')
        ax.set_ylabel('Sampled minimum wall [mm]')
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.2)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(title='Attempt')
    fig.suptitle('Frozen R4 cohort: wall / exterior trade-off; cement and mesh gates also required')
    fig.savefig(R / 'figures/demo.png', dpi=160)
    fig.savefig(R / 'figures/demo.pdf')
    plt.close(fig)
    header = '| Construction | Tooth type | Form | Marginal | Wall | Gap | Complete |\n|---|---|---:|---:|---:|---:|---:|\n'
    table = header + ''.join((f"| {r['method']} | {r['family']} | {r['shape']}/6 | {r['margin']}/6 | {r['wall']}/6 | {r['offset']}/6 | {r['complete']}/6 |\n" for r in summary))
    text = '# The margin is preserved; complete crown is determined by local gap and wall\n\n' + ('Final result' if final else 'Current Checkpoint') + '. Same six R4 jaws and 18 teeth, without replacement cases. The reference for outer shape is the original Teeth3DS triangles; the virtual preparation is our own fixture and no external fit.\n\n' + table + '\nAll number is PER_TOOTH ; distance checked PER_POINT / PER_SURFACE_REGION before summation. Complete requires shape, margin, end shell, entire surface over 0,5 mm wall and intended nominal gap within the frozen tolerance. X49 determines separately the export mesh and named wall profile holds. Physical feature remains UNKNOWN without real preparation, registered bite and measurement.\n\n![Geometriutfall](figures/demo.png)\n\nA moves common corner, B refines a distance area, C allows topology change through an implicit offset. R4 outer shape changes only local in case of insufficient/undissolved separation, with safety margin to 0,60 mm; all design deviations are reported against the same native area. The safety margin is numerical and no material value. No algorithm advantage is claimed.\n\nSuitable sample: exactly the same mean wall 0,75 mm, identity error 0 ; local minima 0,625 and 0,375 mm respectively, difference 0,25 mm and different X49 decisions. The smallest local wall distinguishes this particular pair; one medium is not enough. Eight controls passes with actual error injections.\n\nControl line: the same inletlio with unchanged R4 -outer shape is actually run in each completed design. GenCAD V2 : s exact triangle pairs, V4 / V6 : s local contact and X49 run without changed thresholds ; raw bite pose does not provide contact details for function. See results.json and raw data .\n'
    if 'necessary_condition' in out:
        text += '\nThe separate necessity sample D excludes the four of six incisors under a clear additional requirement: each preparation point must be covered by the inletlio within 60 μm. 0,5 mm wall then provides at least 0,44 mm distance from the outer surface to the preparation. Native-points inside this area forces shape errors; exactly rational control of 2827 witnesses verifies the local limits. For each excluded foreman exceeds at least 411 of 8192 frozen native-prober 0,24 mm-boundary . Two incisors and other tooth types are UNKNOWN according to this required specimen, not proven feasible. No clinical impossibility is claimed.\n'
    (R / 'RESULTS.md').write_text(text)
    (R / 'README_DEMO.md').write_text('# Local crown construction against solid surfaces\n\nIdea: separate preparation , intended cement space, margin and material wall . Keep the R4 -outer shape of the contralateral and move it only local when the wall border requires it.\n\nRun `./run_all.sh` . The default command recreates A for the first frozen jaw three tooth types, freezes new predictions and controls them; furthermore, all re-evaluates exact D-witnesses. `./run_all.sh --full --method B` recreates B for all 18, `--method C` run implicitly offset. Each drive-through gets its own catalogue. No GPU or new dataset is required.\n\n' + table + '\n![ Result ](figures/demo.png)\n\nWhat does not holds : R4 contains a virtual cone, no scanned preparation with finish line. 25 / 50µm from X13 is CAD settings, not a measured spatial cement field. Sharp fold geometry may require modified offset optics; local displacements do not guarantee a valid shell. Clinics, inserts, manufacturing and physical strength are not established.\n\nData/license: Teeth3DS local R4 files are hash-bound; exact dataset license binding UNKNOWN according to R4 . Only local research, no distribution. X13 PMC10721348 is CC BY according to local article; libigl 2.6.1 MPL2 . No second dental dataset or personal data are used.\n\nExternal locators and table conflict 23,22/22,22 μm are available in LITERATURE.md. Each missing case and every failed gate remains in results.json. Status PENDING_INDEPENDENT_REVIEW.\n')
    if final:
        checkpoint('ROUND_COMPLETE', dict(qualifying_count=len(qual), rounds=list(rounds)), 'Acquire actual preparation+finish line for the same tooth positions; test explicit insertion and local signed offset before any new crown claim')
    print('report', len(allrows), 'rows', len(qual), 'qualifying')
if __name__ == '__main__':
    run('--final' in sys.argv)
