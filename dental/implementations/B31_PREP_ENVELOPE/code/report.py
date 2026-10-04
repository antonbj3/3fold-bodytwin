from common import *
from collections import Counter

def main():
    r1 = read(ROOT / 'raw/R1_RESULTS.json')
    r2 = read(ROOT / 'raw/R2_RESULTS.json') if (ROOT / 'raw/R2_RESULTS.json').exists() else []
    r3 = read(ROOT / 'raw/R3_RESULTS.json') if (ROOT / 'raw/R3_RESULTS.json').exists() else []
    r4 = read(ROOT / 'raw/R4_RESULTS.json') if (ROOT / 'raw/R4_RESULTS.json').exists() else []
    controls = read(ROOT / 'raw/CONTROLS.json')
    suff = read(ROOT / 'raw/SUFFICIENCY.json')
    path = {r['key']: r for r in r3}
    surface = {r['key']: r for r in r4}
    exact_surface_controls = [q for r in r4 if r.get('exact_witness_file') for q in read(r['exact_witness_file'])['pairs'] if 'exact_primal_residual' in q]
    if exact_surface_controls:
        controls['rows'].append({'name': 'exact_original_shoulder_primal', 'valid': all((q['exact_primal_residual'] == '0' for q in exact_surface_controls)), 'injected_value': 'each recorded exact point shifted0.2mm while barycentric weights unchanged', 'injected_rejected': all((q['injected_0.2mm_witness_shift_rejected'] for q in exact_surface_controls)), 'tested_witnesses': len(exact_surface_controls)})
    summaries = {'R1': {'requested': 18, 'completed': len(r1), 'complete_crowns': 0, 'source_closure_unknown': sum((r['status'] == 'UNKNOWN_SOURCE_CLOSURE' for r in r1)), 'reconstruction_unknown': sum((r['status'] == 'UNKNOWN_SURFACE_RECONSTRUCTION' for r in r1)), 'status_counts': dict(Counter((r['status'] for r in r1)))}, 'R2': {'requested': 18, 'completed': len(r2), 'complete_crowns': 0, 'candidate_shells': sum(('candidate_path' in r for r in r2)), 'closed_shells_nominal': sum((r.get('closed_shell', False) for r in r2)), 'exact_carrier_cone_pass': sum((r.get('canonical_intaglio', {}).get('exact_cone_failures') == 0 for r in r2)), 'nominal_wall_lower_pass': sum((r.get('wall', {}).get('geometric_lower_mm', -1) >= 0.5 for r in r2)), 'source_shoulder_counterwitnesses': sum((r.get('added_source_control', {}).get('numeric_intersection_witnesses', 0) for r in r2)), 'source_parity_control_mismatches': sum((r.get('field', {}).get('parity_control', {}).get('interior_control_mismatches', 0) for r in r2)), 'status_counts': dict(Counter((r['status'] for r in r2)))}, 'R3': {'requested': 18, 'completed': len(r3), 'exact_axis_clear': sum((r['whole_path_status'] == 'EXACT_CLEAR_FOR_FINITE_PREP_AXIS_MODEL' for r in r3)), 'exact_axis_collision': sum((r['whole_path_status'] == 'EXACT_COLLISION_COUNTERWITNESS' for r in r3)), 'unknown_no_candidate': sum((r['whole_path_status'] == 'UNKNOWN_NO_CANDIDATE' for r in r3)), 'collision_pairs': sum((r.get('collision_pairs', 0) for r in r3)), 'LP_disagreements': sum((r.get('LP_control_disagreements', 0) for r in r3))}}
    summaries['R4'] = {'requested': 18, 'completed': len(r4), 'exact_crossing_cases': sum((r['status'] == 'EXACT_EXPORTED_SHELL_CROSSING_COUNTERWITNESS' for r in r4)), 'exact_crossing_pairs': sum((r.get('proper_crossing_pairs', 0) for r in r4)), 'tested_positive_pairs': sum((r.get('tested_pairs', 0) for r in r4)), 'LP_disagreements': sum((r.get('LP_disagreements', 0) for r in r4)), 'uninspected_negative_pairs': 'UNKNOWN_NOT_A_GLOBAL_CERTIFICATE'}
    adjudication = []
    for r in r2:
        s = surface.get(r['key'], {})
        p = path.get(r['key'], {})
        reasons = []
        if s.get('proper_crossing_pairs', 0):
            reasons.append('EXACT_EXPORTED_SHELL_INTERIOR_CROSSING')
        if p.get('whole_path_status') == 'EXACT_COLLISION_COUNTERWITNESS':
            reasons.append('EXACT_EXPORTED_SHELL_PREPARATION_SWEEP_COLLISION')
        adjudication.append({'key': r['key'], 'resolution': 'PER_TOOTH', 'status': 'REJECT_EXPORTED_GEOMETRY' if reasons else 'UNKNOWN', 'reasons': reasons or [r.get('reason', 'source/shell embedding and rigorous wall/serialization enclosure not established')], 'physical_status': 'UNKNOWN_NOT_MEASURED'})
    summaries['whole_crown'] = {'requested': 18, 'certified': 0, 'geometrically_rejected': sum((r['status'] == 'REJECT_EXPORTED_GEOMETRY' for r in adjudication)), 'unknown': sum((r['status'] == 'UNKNOWN' for r in adjudication))}
    complete = len(r2) == len(r3) == len(r4) == 18
    costs = {p.stem: read(p) for p in (ROOT / 'raw').glob('*_COST.json')}
    artifacts = []
    for r in r2 + r3:
        for (k, v) in r.items():
            if k.endswith('_path') and isinstance(v, str) and Path(v).exists():
                artifacts.append({'path': v, 'bytes': Path(v).stat().st_size, 'sha256': sha(v), 'resolution': r.get('resolution', 'PER_POINT')})
    result = {'claim_type': 'capability', 'status': 'BOUNDED_DIGITAL_COMPONENTS_COMPLETE_CROWN_UNKNOWN' if complete else 'RUNNING', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'capability': 'Search an insertable preparation inside exactly locked original crown exterior; certify components and reject incomplete whole crowns.', 'answer': 'No complete crown is certified. Every frozen case receives explicit geometric results and remaining unknowns. Exact carrier union and whole-shell axial sweep are distinct from physical crown validation.', 'inherited_gate': read(ROOT / 'PREREG_R1.json')['metrics'], 'complete_digital_count': 0, 'cohort_requested': 18, 'resolution': 'PER_POINT', 'summary_resolution': 'PER_TOOTH; aggregate counts are finite frozen cohort only', 'edge_timescale': 'SIMULTANEOUS', 'external_referent': read(ROOT / 'PREREG_R1.json')['external_referent'], 'external_numerical_control': read(ROOT / 'PREREG_R1.json')['external_numerical_control'], 'summaries': summaries, 'rounds': {'R1': r1, 'R2': r2, 'R3': r3}, 'sufficiency': suff, 'controls': controls, 'dropout': {'R1': {'requested': 18, 'rejected': 18, 'fraction': 1.0, 'reasons': dict(Counter((r.get('reason', '') for r in r1)))}, 'R2': {'requested': 18, 'candidate_rejected_or_unknown': sum(('candidate_path' not in r for r in r2)), 'fraction': sum(('candidate_path' not in r for r in r2)) / 18, 'reasons': dict(Counter((r.get('reason', '') for r in r2 if 'candidate_path' not in r)))}, 'whole_crown': {'requested': 18, 'certified': 0, 'uncertified_fraction': 1.0, 'reason': 'source/whole-shell self intersection and rigorous wall/serialization enclosure missing, plus explicit sweep collision where present'}, 'external_source_screening': {'attempted': 5, 'retained': 4, 'rejected': 1, 'fraction': 0.2, 'reasons': ['S12 web PDF failed; no numeric product threshold imported']}}, 'cost': {'measured_components': costs, 'data_usage': usage(), 'fit': 0, 'legacy_acquisition': 'UNKNOWN', 'interrupted_attempt_elapsed': 'UNKNOWN; all preserved in raw logs, not counted as free', 'lab_measurement_and_instruments': 'UNKNOWN', 'technician_time_gain': 'NOT_MEASURED', '10x_gain': 'NOT_MEASURED; baseline complete_count0 has undefined ratio', 'preparation_discovery_validation_queries_fallback': 'all geometric operations included in measured component timings; historical acquisition and manual interpretation not separately timed'}, 'knowledge_debt': [{'quantity': 'true same-case prepared tooth / finish line / pulp', 'resolution': 'PER_POINT', 'replacement_measurement': 'LAB_CONTRACT.md same-object before/after and spatial seated-fit scan'}, {'quantity': 'full source/shell self-intersection and floating wall/STL enclosure', 'resolution': 'PER_SURFACE_REGION', 'replacement_operation': 'exact surface predicates plus outward-rounded metric/serialization bounds'}, {'quantity': 'product manufacture/retention/load qualification', 'resolution': 'PER_TOOTH', 'replacement_measurement': 'specified product protocol on fabricated specimen'}], 'artifacts': artifacts, 'source_license_binding': read(ROOT / 'raw/LICENSE_BINDING.json') if (ROOT / 'raw/LICENSE_BINDING.json').exists() else {'status': 'UNKNOWN'}, 'limits': read(ROOT / 'PREREG_R2.json')['limits']}
    result['rounds']['R4'] = r4
    result['case_adjudication'] = adjudication
    result['answer'] = 'Seven exported candidates are refuted by exact geometric witnesses; the other eleven frozen cases remain UNKNOWN. No complete crown is certified. Original labelled triangles are exterior truth, not empirical preparation truth.'
    dump(ROOT / 'results.json', result)
    table = ["| Case / tooth type | Exact cone | Nominal full wall, mm | The whole axial road | Complete crown |", '|---|---|---:|---|---|']
    for r in r2:
        cp = r.get('canonical_intaglio', {})
        cone = 'JA' if cp.get('exact_cone_failures') == 0 else 'UNKNOWN'
        lb = r.get('wall', {}).get('geometric_lower_mm')
        value = f'{lb:.6f}' if lb is not None else 'UNKNOWN'
        p = path.get(r['key'], {}).get('whole_path_status', 'NOT_RUN')
        table.append(f"| {r['key']} | {cone} | {value} | {p} | UNKNOWN |")
    text = "The original exterior can be held exactly while the inside is given a trialable concierge structure. Seven exported candidates are convicted by exact geometric counter-witnesses; the other eleven cases are UNKNOWNNo complete crown is certified on the18 frysta Teeth3The DS teeth. reference observations are the original annotated triangles; preparation, basal slope and material/cement scenario are virtual.\n\n"
    text += f"R1: {summaries['R1']['source_closure_unknown']}/18 source closures failed; five voxel exports broke the deposit cone and lacked working opening. R2: {summaries['R2']['exact_carrier_cone_pass']}/18 exakta konunioner, {summaries['R2']['closed_shells_nominal']}/18 nominally closed crown shells and {summaries['R2']['nominal_wall_lower_pass']}/18 floating point limits above0,5mm. R3: {summaries['R3']['exact_axis_clear']} entire surfaces are exactly free from the finite virtual preparation20mm axial road, {summaries['R3']['exact_axis_collision']} These partial outcomes are not a conjuncture for a functioning crown.\n\n"
    text += '\n'.join(table) + '\n\n![Figur](figures/demo.png)\n\n'
    text += f"R4: {summaries['R4']['exact_crossing_pairs']} rational witnesses show strict intersections between original surface and margin bridge in{summaries['R4']['exact_crossing_cases']} These five cases are separate from the two axial collision cases.49 witnesses are convicted when the point is displaced0,2mm with unchanged barycentric weights. Same49 HiGHS comparisons0 differences. Only previously positive couples are retested; the absence of counter-witnesses is not a global approval gate.\n\n"
    text += "Each  local   distance  trial  PER_POINT  before  PER_TOOTH -tabellen. The wall boundary is a geometric centroid/radial boundary over all's original facets, but its floating point enclosure is missing. Precise conunion applies to rational carriers; STL -quantization and the original/marginal full self-intersection are separate UNKNOWN . The exported two-bearer family can maintain very small volume and is not the smallest preparation . 0,5 mm is a technical scenario, no general IFU .\n\n"
    text += f"Satisfaction test: exact volume difference{suff['exact_volume_difference']}mm³, identitetsfel i samplat histogram{suff['sampled_histogram_identity_error_mm']}mm. Downstream local wall minimum separates{suff['minimum_sampled_difference_mm']:.9f}mm and the validation gate changes. The minimum additions for the validation gate are the verified minimum of the whole wall. A completely real thickness histogram can't have any other minimum; the original conflicting formulation has been corrected openly. The first volume sample missed machine exact identity and is preserved as FAILED_V1.\n\n"
    text += "Control line: same sources , carriers and triangles tested with explicit box control, direct surface control and SciPy HiGHS . No algorithm or time superior unit is claimed. 0,2 mm-bulan rejects same continuous wall gate; each control family has an injected error that is rejected.\n\n"
    text += "Physical continuation requires measured preparation , real margin, same-case powder/tissue information under such requirement, and product manufacturing/last qualification. See LAB_CONTRACT.md , README_DEMO.md , all PREREG and preserved bugs in COMMANDS.md . Status PENDING_INDEPENDENT_REVIEW.\n"
    (ROOT / 'RESULTS.md').write_text(text)
    figure(result)
    return result

def figure(result):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    rr = result['rounds']['R2']
    fig = plt.figure(figsize=(17, 5), layout='constrained')
    ax = [fig.add_subplot(1, 3, 1), fig.add_subplot(1, 3, 2), fig.add_subplot(1, 3, 3, projection='3d')]
    values = [r.get('wall', {}).get('geometric_lower_mm', np.nan) for r in rr]
    ax[0].bar(np.arange(len(rr)), values, color='#387892')
    ax[0].axhline(0.5, color='#ae4937', linestyle='--', label='0.5mm technical gate')
    ax[0].set(xlabel='Frozen case/tooth index (all18)', ylabel='Whole-wall floating lower, mm', title='Local wall model; rigorous rounding UNKNOWN')
    ax[0].legend()
    if any(('candidate_path' in r for r in rr)):
        r = next((r for r in rr if 'candidate_path' in r))
        d = load_np(r['candidate_path'])
        o = d['original']
        i = d['intaglio_vertices']
        s = d['vertices'][d['faces'][d['roles'] == 2]]
        ax[1].scatter(o[::max(1, len(o) // 1500)].mean(1)[:, 0], o[::max(1, len(o) // 1500)].mean(1)[:, 2], s=2, color='#777777', label='Locked original exterior')
        ax[1].scatter(s.mean(1)[:, 0], s.mean(1)[:, 2], s=4, color='#c07539', alpha=0.5, label='Constructed shoulder')
        ax[1].scatter(i[:, 0], i[:, 2], s=8, color='#1c6989', label='Exact two-core intaglio')
        ax[1].set(title=r['key'] + ' (unqualified geometry)', xlabel='Local x, mm', ylabel='Local z, mm')
        ax[1].axis('equal')
        ax[1].legend(fontsize=8)
    else:
        ax[1].text(0.1, 0.5, 'No surface candidate; source/closure UNKNOWN')
    crossing = next((r for r in result['rounds'].get('R4', []) if r.get('proper_crossing_pairs')), None)
    if crossing:
        w = next((q for q in read(crossing['exact_witness_file'])['pairs'] if q['status'] == 'EXACT_STRICT_INTERIOR_CROSSING'))
        A = np.array(w['added_triangle_mm'])
        B = np.array(w['source_triangle_mm'])
        p = w['point_mm']
        ax[2].add_collection3d(Poly3DCollection([A], facecolor='#c07539', alpha=0.5, edgecolor='#9d4b22'))
        ax[2].add_collection3d(Poly3DCollection([B], facecolor='#387892', alpha=0.5, edgecolor='#165675'))
        ax[2].scatter(*p, color='black', s=30, label='Exact interior crossing')
        both = np.r_[A, B]
        lo = both.min(0)
        hi = both.max(0)
        mid = (lo + hi) / 2
        radius = max(hi - lo) * 0.55
        ax[2].set_xlim(mid[0] - radius, mid[0] + radius)
        ax[2].set_ylim(mid[1] - radius, mid[1] + radius)
        ax[2].set_zlim(mid[2] - radius, mid[2] + radius)
        ax[2].set(title='Rational crossing: shoulder / source', xlabel='x, mm', ylabel='y, mm', zlabel='z, mm')
        ax[2].legend(fontsize=8)
    (ROOT / 'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT / 'figures/demo.png', dpi=170)
    plt.close(fig)
if __name__ == '__main__':
    main()
