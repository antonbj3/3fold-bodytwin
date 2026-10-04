from fc_common import *
import collections, platform, resource
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def run():
    rounds = {f'R{i}': read(ROOT / 'rounds' / f'R{i}.json') for i in range(1, 8)}
    audit = read(ROOT / 'raw/SOURCE_GEOMETRY_AUDIT.json')
    val = read(ROOT / 'raw/VALIDATION.json')
    exports = read(ROOT / 'raw/EXPORT_VALIDATION.json')
    counts = []
    for tag in ['R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
        for (participant, dd) in rounds[tag]['summary'].items():
            if participant.endswith('_outer'):
                continue
            for (family, ss) in dd.items():
                counts.append(dict(round=tag, participant=participant, family=family, resolution='PER_TOOTH', population_summary_resolution='POPULATION', **ss))
    datafiles = {str(p.relative_to(DATA)): {'bytes': p.stat().st_size, 'sha256': sha(p)} for p in DATA.rglob('*') if p.is_file() and '/replay/' not in str(p)}
    costs = {}
    for tag in ['R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
        costs[tag] = {'generation': read(DATA / (tag + '_predictions') / 'COST.json'), 'scoring_seconds': rounds[tag]['seconds']}
    costs.update(R1_diagnostic_seconds=rounds['R1']['seconds'], public_extraction_seconds=read(ROOT / 'FROZEN_PUBLIC_INPUTS.json')['seconds'], inherited_X11_fit_and_prior_discovery='UNKNOWN_NOT_ZERO', model_reasoning_and_source_reading_time='NOT_INSTRUMENTED', physical_measurements='NOT_RUN', fallback='Every failed construction retained, no positive substitute', threads=4, gpu_used=False, data_bytes=sum((x['bytes'] for x in datafiles.values())))
    result = dict(lane='PROOF_LANE-full-crown', claim_type='capability', status='FAILED_ANATOMICAL_GATE_WITH_DIAGNOSED_REPRESENTATION_AND_REFERENCE_SUPPORT', review_state='PENDING_INDEPENDENT_REVIEW', external_referent=REFERENT, primary_answer='The0/138 anatomical frontier has NOT been broken. Full-surface homolog and preparation-conditioned generators both failed the unchangedv5 gate. A full-prescan co-design branch generated actual closed crown/die pairs with conditional wall/insertion/film proofs, but it also failed the anatomy gate.', resolution_schema={'distances_and_wall': 'PER_POINT; p95 aggregated PER_TOOTH', 'contacts': 'PER_SURFACE_REGION', 'family_counts': 'POPULATION from6 case clusters, not independent method submissions', 'nominal_film_tool_material_requirements': 'PHENOMENOLOGICAL engineering design requirements; see measurement debts'}, per_type_table=counts, rounds=rounds, source_geometry_audit=audit, validation=val, exports=exports, costs=costs, data_manifest=datafiles, dropout_reporting='Per-round numerator/denominator and exact reasons in rows; methods/cases/tracks never pooled as independent patients', uncertainty={'source_surface': 'unmeasured scan uncertainty; X11 target FDI unvalidated in these datasets', 'cervical_margin': 'virtual scalar15percent quantile, not an annotated anatomical finish line', 'reconstruction': 'area-stratified1024 triangle-centroid probes; no continuous p95 enclosure;4096 refinement was registered for passes, none occurred', 'wall': 'outward-rounded whole-triangle Lipschitz lower bound against analytic cavity, subdivision allowance1e-8mm; physical scan/process uncertainties separate', 'cohort': 'retrospective six-case panel; no fresh author-blind or population estimate', 'physical': 'no candidate-specific fracture force, actual seated film, holder-awareCAM or retention validation'}, measurement_debts=[dict(quantity='actual cervical finish line', resolution='PER_POINT', replacement='independent annotated margin on complete pre/post-preparation scan'), dict(quantity='seated cement film', resolution='PER_SURFACE_REGION', replacement='registered dry/seated microCT or replica map on same printed specimen'), dict(quantity='real cutter/machine access', resolution='PER_POINT', replacement='actual tool+holder trajectory and part metrology; current proof only ideal ball sweep'), dict(quantity='loaded contacts/strength', resolution='PER_SURFACE_REGION', replacement='registered pressure/load and matched geometry/batch mechanical test')], edge_contracts=[dict(producer='native full-surface triangles with observed/virtual provenance', consumer='whole-crown scorer and source-conditioned generator', resolution='PER_POINT', time_scale='SIMULTANEOUS'), dict(producer='manufactured registered seated film (missing)', consumer='mechanical model', resolution='PER_SURFACE_REGION', time_scale='HANDOVER')], next_construction='Freeze independently annotated tooth/cervical boundaries and reference-support ownership; construct full-surface donor prior on patient-disjoint verified labels and test unchanged0.35mm plus spatial function. Real registered pre/post-preparation pair remains missing.', scientific_admission=False)
    dump(ROOT / 'results.json', result)
    (fig, axs) = plt.subplots(2, 2, figsize=(13, 9))
    ax = axs[0, 0]
    regs = ['margin', 'occlusal', 'proximal', 'buccolingual_axial']
    colors = ['#B96734', '#477BB8', '#57996C', '#8C66A4']
    for (j, track) in enumerate(['R3', 'R5']):
        rr = [r for r in rounds['R1']['regional_rows'] if r['track'] == track]
        ct = [sum((r['tail_count'] for r in rr if r['region'] == reg)) for reg in regs]
        left = 0
        for (k, c) in enumerate(ct):
            w = 100 * c / sum(ct)
            ax.barh(j, w, left=left, color=colors[k], label=regs[k] if j == 0 else None)
            left += w
    ax.set_yticks([0, 1], ['Inherited roofs', 'Inherited prep offsets'])
    ax.set_xlabel('Share of worst 5% probes, both directions (%)')
    ax.set_title('A  Where inherited errors occur')
    ax.legend(fontsize=7, loc='lower right')
    ax = axs[0, 1]
    labels = []
    vals = []
    for tag in ['R2', 'R3', 'R6', 'R7']:
        rr = [r for r in rounds[tag]['rows'] if r.get('kind') == 'shell' and r.get('reconstruction_p95_mm') is not None]
        labels.append(tag + ('\nfull prescan' if tag in ['R6', 'R7'] else '\nhidden form'))
        vals.append(min((r['reconstruction_p95_mm'] for r in rr)))
    ax.bar(labels, vals, color=['#477BB8', '#477BB8', '#57996C', '#57996C'])
    ax.axhline(0.35, c='black', linestyle='--', label='Fixed 0.35 mm gate')
    ax.set_ylabel('Best complete-shell p95 (mm)')
    ax.set_title('B  No anatomy gate passed')
    ax.legend(fontsize=8)
    ax = axs[1, 0]
    ax.axis('off')
    text0 = 'Two states, identical source p95 = 0.250000 mm\n\nState A: contact 2 mm²; penetration 0 mm²\nState B: contact 0 mm²; penetration 1 mm²\n\nIdentity error = 0 exactly\nBoth pass the 0.35 mm similarity gate\n\nRequired extension: spatial signed gap field'
    ax.text(0.02, 0.95, text0, va='top', fontsize=12)
    ax.set_title('C  Similarity does not decide function', loc='left')
    ax = axs[1, 1]
    r = next((r for r in exports['selected'] if r.get('path')))
    base = DATA / (r['round'] + '_predictions') / ('prescan_local_holes_die_shell' if r['round'] == 'R7' else 'prescan_signed_solid_die_shell') / r['key']
    m = npz(base / 'mesh.npz')
    v = m['vertices']
    ext = v[m['faces'][m['face_roles'] == 0]].mean(1)
    inner = v[m['faces'][m['face_roles'] == 1]].mean(1)
    ax.scatter(ext[:, 0], ext[:, 2], s=0.15, c='#477BB8', label='Research exterior')
    ax.scatter(inner[:, 0], inner[:, 2], s=0.15, c='#B96734', label='Analytic cavity')
    ax.set_aspect('equal')
    ax.set_xlabel('Native x (mm)')
    ax.set_ylabel('Native z (mm)')
    ax.set_title(f"D  Best molar specimen; p95 {r['p95_mm']:.3f}mm (FAIL)")
    ax.legend(markerscale=10, fontsize=8)
    fig.suptitle('Full crown: reconstruction failure, source-support diagnosis and research specimens', fontsize=14)
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/demo.png', dpi=180)
    fig.savefig(ROOT / 'figures/demo.svg')
    plt.close(fig)
    lines = ['# Complete anatomical crown: the0.35mm gate remains failed', '', 'Claim type: capability. PENDING_INDEPENDENT_REVIEW. Not a clinical recommendation.', '', 'The original missing-form frontier is not solved. The unchangedv5 scorer evaluated actual public-only full-surface generators. All18 known-source height-envelope oracles also failed, identifying a representation problem that persists even before prediction.', '', '| Track / method | Molar pass/requested | Premolar | Anterior | Best complete-shell p95 mm |', '|---|---:|---:|---:|---:|']
    for tag in ['R2', 'R3', 'R6', 'R7']:
        for (part, dd) in rounds[tag]['summary'].items():
            if part.endswith('_outer'):
                continue
            fields = [f"{dd[f]['passes']}/{dd[f]['requested']}" for f in ['molar_crown', 'premolar_crown', 'anterior_crown']]
            rr = [r for r in rounds[tag]['rows'] if r['participant'] == part and r.get('reconstruction_p95_mm') is not None]
            best = min((r['reconstruction_p95_mm'] for r in rr))
            lines.append('|' + tag + ' / ' + part + '|' + '|'.join(fields) + f'|{best:.6f}|')
    lines += ['', 'Fractions are POPULATION summaries of PER_TOOTH decisions from six case clusters. R2 has two methods on the same18 sites. R3 has3/18 shell abstentions. R6 andR7 each have8/18 failed constructions. R4 andR5 each refused18/18. R6/R7 use a full preoperative scan AND a co-designed preparation; they must not be added to the missing-form result.', '', 'Reference support: all18 original tooth-label surfaces are open; none has edges with>2 incident faces in the audited original meshes. Failures in virtual closure are separate from source nonmanifoldness. In each of10 R7 exported cases,100% of the worst5% forward probes are closer to virtual closure facets than to observed reference triangles. This diagnoses source support and completion error; it does not prove clinical validity or that0.35mm is universally inappropriate.', '', 'The original scan is a reconstruction referent, not a uniquely optimal crown. The exact-summary probe gives identicalp95=0.25mm, identity error0, but contact differs2mm² and interference1mm². The minimal extension is the spatial signed gap field. See [figure](figures/demo.png).', '', 'The strongest same-information rigid and affine homolog controls were actually run; neither passed. Anatomical reflection, preparation-conditioned deformation and topology repairs are preserved as failures. No algorithm-superiority claim.', '', 'Ten source-conditioned research pairs in each ofR6/R7 satisfy the actual final-mesh wall lower bound>=0.5mm against an analytic cavity. Ideal aligned insertion follows from downward-closed cavity halfspaces; ideal spherical-cutter accessibility has no holder/machine guarantee. Film0.08mm is a prescribed analytic/digital offset, not measured cement. No physical force, fatigue, retention or clinical margin is established.', '', 'Validation: ' + str(val['passed']) + '/' + str(val['count']) + ' checks pass. The two failures are the same case’s float32STL topology loss inR6/R7; they are retained. All fault-injection rejection checks pass. Selected indexed3MF pairs preserve float64 vertices and indices exactly and pass topology/unit faults. No anterior pair satisfies this particular swept-ball design class.', '', 'No independently annotated target FDI, measured pre/post preparation, loaded contact map or candidate-specific material calibration is available here. Signed-solid error lacks a rigorous source-geometry enclosure; the wall bound instead covers final exported triangles. Next: independent cervical/target annotation and a source-support-aware whole-surface dataset, followed by a patient-disjoint form prior and unchanged geometric/function gates.']
    (ROOT / 'RESULTS.md').write_text('\n'.join(lines) + '\n')
    dump(ROOT / 'raw/COSTS.json', costs)
    return result
if __name__ == '__main__':
    run()
