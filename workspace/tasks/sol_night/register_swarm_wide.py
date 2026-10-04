"Widths BodyTwin heat catalogue (tasks/free48/CATALOG.json) with new source families, atomic and additive.\n\nSame form as graflan's BT-CTX-Q### families 30/9: real source files copied to sources/<KEY>/,\na definition target BT-CTX-<KEY> in notes/RESEARCH_ACTIONS_20260926.json and a target context file. All new targets are\ndefine-only (numerical_dispatch_allowed=false). Existing keys, packages and files are not touched; backup first.\ntyped_relations is left empty: related native nodes are selected by graph-lan, not here.\nRuns under controller lock. Usage: python3 register_swarm_wide.py [--dry]\n"
import fcntl, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path

B = Path(''); A = B / 'tasks/free48'; R = B / 'results'
CATF = A / 'CATALOG.json'; RA = B / 'notes/RESEARCH_ACTIONS_20260926.json'
CTXD = B / 'notes/GRAPH_TARGET_CONTEXT_SWARMWIDE_20260930'
JOB = 'SOLNIGHT-SWARMWIDE-20260930'; BUDGET = 400_000
DRY = '--dry' in sys.argv

QS = ['Q009', 'Q014', 'Q017', 'Q019', 'Q020', 'Q021', 'Q022', 'Q026', 'Q031', 'Q036', 'Q043', 'Q044', 'Q049',
      'Q054', 'Q058', 'Q077', 'Q080', 'Q084', 'Q088', 'Q100', 'Q107', 'Q121', 'Q127', 'Q140', 'Q154', 'Q156', 'Q160']

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def title(q):
    for p in (R / f'BT-HX-{q}/BRIEF.md',):
        try:
            for l in p.read_text().splitlines():
                if l.strip(): return l.lstrip('# ').split('—', 1)[-1].strip()
        except OSError: pass
    return q

def families():
    out = []
    for q in QS:
        hx, dat = R / f'BT-HX-{q}', R / f'BT-DAT-{q}'
        files = [(f'{q}_model.py', hx / 'model.py'), (f'{q}_QUESTION.md', hx / 'inputs/QUESTION.md'),
                 (f'{q}_RESULTS.md', hx / 'RESULTS.md'), (f'{q}_PREREG.md', hx / 'PREREG.md'),
                 (f'{q}_results.json', hx / 'results.json'), ('DAT_RESULTS.md', dat / 'RESULTS.md'),
                 ('DAT_DATA_SOURCES.json', dat / 'DATA_SOURCES.json')]
        t = title(q)
        out.append(dict(key=q, target=f'BT-CTX-{q}', parent=f'BT-DAT-{q}', files=files, title=t,
            context=(f'Existing BodyTwin first-principles question and baseline: {t} Code and data reports are unaudited. '
                     'Recover the precise observable, units, input domain and old failures before expanding. Combine with '
                     'other catalog families through dimensionally defined ports when a shared variable, boundary or '
                     'consumer exists. The graph binding is the source-family target context '
                     f'BT-CTX-{q}; it does not assert native-graph coverage or biological validation.')))
    sol = R / 'LANE_AMBITIOUS_QUERY_GRAPH'
    out.append(dict(key='SOLBENCH', target='BT-CTX-GLUCOSE-HISTORY', parent='LANE_AMBITIOUS_QUERY_GRAPH', title='Sol coupled history benchmark',
        files=[('model.py', sol / 'shared/model.py'), ('native_glucose.py', sol / 'shared/inputs/native_glucose.py'),
               ('response_operator.py', sol / 'shared/response_operator.py'), ('context_operator.py', sol / 'shared/context_operator.py'),
               ('QUERY_GRAPH_RESULTS.md', sol / 'RESULTS.md'), ('HISTORY_INVERSE_RESULTS.md', R / 'LANE_AMBITIOUS_HISTORY_INVERSE/RESULTS.md'),
               ('COUPLED_RESPONSE_RESULTS.md', R / 'LANE_AMBITIOUS_COUPLED_RESPONSE/RESULTS.md'),
               ('OPEN_HYPOTHESES_COUPLED_RESPONSE.md', B / 'tasks/build_night/steer/LANE_AMBITIOUS_COUPLED_RESPONSE.md'),
               ('OPEN_HYPOTHESES_HISTORY_INVERSE.md', B / 'tasks/build_night/steer/LANE_AMBITIOUS_HISTORY_INVERSE.md'),
               ('OPEN_HYPOTHESES_QUERY_GRAPH.md', B / 'tasks/build_night/steer/LANE_AMBITIOUS_QUERY_GRAPH.md')],
        context=('Synthetic, uncalibrated coupled native glucose/insulin organ plus 32-512 finite-memory regions with shared '
                 'uncertainty rows, used by three Sol lanes hunting a scalable counterfactual/inverse operation that beats the '
                 'strongest equally informed sparse/ROM replay. Frozen snapshot 30/9 23:30. Useful swarm work: cheap decisive '
                 'sub-experiments on the OPEN_HYPOTHESES (mixed monotonicity sign patterns, numerical rank of history->future '
                 'exposure maps, per-query cost decomposition, amortization crossover), or importing a structure from another '
                 'catalog family. Strongest control must get identical information and all setup cost. No clinical or '
                 'empirical claims; everything PENDING_INDEPENDENT_REVIEW.')))
    DN = Path('local_path')
    out.append(dict(key='BIORESP', target='BT-CTX-BIORESP', parent='dental:BIORESP', title='Bone remodeling under implant load: a mechanistic physiological operator',
        files=[('REMODEL_crestal_REPORT.md', DN / 'results/REMODEL_crestal/REPORT.md'), ('analyze_remodel.py', DN / 'results/REMODEL_crestal/analyze_remodel.py'),
               ('MECHANO_sign_REPORT.md', DN / 'results/MECHANO_sign/REPORT.md'), ('ORTHO_rate_REPORT.md', DN / 'results/ORTHO_rate/REPORT.md'),
               ('dental_physics_material_procedure.md', DN / 'notes/northstar_inventory/02_physics_material_procedure.md')],
        context=('Cross-lane family from dental (anton-08, 30/9). Crestal/peri-implant bone remodeling. THREE generic mechanostat drivers are '
                 'REFUTED in dental and must not be repeated or re-banded: absolute Frost bands, the tooth-referenced mechanostat and the '
                 'ortho PDL-pressure driver; a biological-width model passes 2/5 contrasts. Needed: a mechanistic physiological operator '
                 '(cell populations, signaling, transport, perfusion, damage/repair from first principles) that changes a remodeling '
                 'decision and beats the existing 2/5 model on held-out contrasts. Consumer: dental implant design. Synthetic/public only.')))
    NB = Path('source_repository/'); M = NB / 'scripts/msk'
    SURG = [
     ('SURG_INCISION', 'Scalpel incision through layered skin/subcutis/fascia: cutting mechanics, wound gape and local transport',
      [('Q033_incision_model.py', R / 'BT-HX-Q033/model.py'), ('Q033_RESULTS.md', R / 'BT-HX-Q033/RESULTS.md'), ('Q033_results.json', R / 'BT-HX-Q033/results.json'),
       ('Q033_PREREG.md', R / 'BT-HX-Q033/PREREG.md'), ('skin_pulp_mechanics.py', M / 'skin_pulp_mechanics.py'), ('skin_pulp_mechanics_evidence.json', M / 'skin_pulp_mechanics_evidence.json'),
       ('bone_fracture_toughness_lefm.py', M / 'bone_fracture_toughness_lefm.py'), ('skin_barrier_tewl.py', M / 'skin_barrier_tewl.py')],
      ['SKIN-SUBCUTIS-DECOMPOSITION', 'SKIN-DERMAL-STIFFNESS-COLLAGEN-COMPOSITION', 'MSK-FASCIA-NETWORK', 'MODEL-BONE-FRACTURE-TOUGHNESS'],
      'Blade enters layered tissue. Resolve what the scalpel actually does from first principles: fracture/cutting energy (work of fracture + friction + far-field deformation; blade edge radius and angle), layer-specific toughness and prestress (epidermis, dermis collagen network, subcutis adipose, fascia), fibre orientation / tension lines setting gape, the cell damage zone at the cut front, severed microvessels and interstitial/evaporative transport. Q033 has an elliptical-gape + Darcy model with a cadaver closure-force anchor 4.2-6.0 N.'),
     ('SURG_COLLAGEN', 'Collagen bindings to fibres to tissue: crosslinks, fibril sliding, thermal denaturation and rupture under a blade',
      [('collagen_triple_helix_thermal_stability.py', M / 'collagen_triple_helix_thermal_stability.py'), ('collagen_triple_helix_thermal_stability.json', NB / 'data/msk_results/collagen_triple_helix_thermal_stability.json'),
       ('tendon_collagen_hierarchical_mechanics.py', M / 'tendon_collagen_hierarchical_mechanics.py'), ('tendon_collagen_hierarchical_mechanics.json', NB / 'data/msk_results/tendon_collagen_hierarchical_mechanics.json'),
       ('myofascial_transmission.py', M / 'myofascial_transmission.py')],
      ['MODEL-TENDON-COLLAGEN-MECHANICS', 'MODEL-COLLAGEN-TRIPLE-HELIX-THERMAL-STABILITY', 'MSK-FASCIA-NETWORK'],
      'Molecular-to-tissue chain for cutting: triple-helix stability and crosslink chemistry, fibril/fibre sliding and recruitment, fibre network orientation, and how bond rupture energy scales up to macroscopic tissue toughness and to thermal (electrosurgical) denaturation zones.'),
     ('SURG_HEMOSTASIS', 'Bleeding and hemostasis after a cut: severed vessel flow, platelet plug, coagulation and fibrinolysis',
      [('coagulation_hemostasis.py', M / 'coagulation_hemostasis.py'), ('platelet_hemostasis.py', M / 'platelet_hemostasis.py')],
      ['MODEL-COAGULATION-THROMBIN-GENERATION', 'MODEL-COAGULATION-FIBRINOLYSIS', 'MODEL-PLATELET-ACTIVATION-AGGREGATION', 'HEMATO-HEMOSTATIC-BALANCE'],
      'Couple the incision geometry (vessel density and calibre per layer, cut length/depth) to bleeding rate and time to hemostasis: Poiseuille/wall-shear at severed ends, platelet adhesion under shear, thrombin generation and fibrin, local pressure and tissue compliance.'),
     ('SURG_HEALING', 'From incision to closure: inflammation, cell migration, matrix deposition and scar mechanics',
      [('wound_healing_cascade.py', M / 'wound_healing_cascade.py'), ('wound_healing_cascade_results.json', NB / 'data/msk_smoketest/wound_healing_cascade/wound_healing_cascade_results.json'),
       ('acute_phase_inflammation.py', M / 'acute_phase_inflammation.py'), ('Q036_ischemia_model.py', R / 'BT-HX-Q036/model.py'), ('Q036_RESULTS.md', R / 'BT-HX-Q036/RESULTS.md')],
      ['WOUND-HEALING-CASCADE', 'MODEL-NLRP3-INFLAMMASOME'],
      'Healing as a function of the cut: damage zone and ischemic margin (Q036 state vector for reperfusion), inflammatory and migratory cell dynamics, collagen deposition/remodeling and the resulting scar stiffness and strength over time.')]
    for key, t, files, rel, ctx in SURG:
        out.append(dict(key=key, target='BT-CTX-' + key.replace('_', '-'), parent='surgical-chain-20260930', title=t, files=files, relations=rel,
            context=('Surgical high-resolution family (Anton 30/9: increase resolution where a scalpel opens tissue - materials, cells, tissues, bindings). '
                     + ctx + ' Combine with the other SURG_* families through dimensionally defined ports (cut geometry, damage zone, vessel map, time). '
                     'Code is unaudited; public literature only; strongest equally informed control; PENDING_INDEPENDENT_REVIEW.')))
    return out

def main():
    fams = families()
    with (A / 'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        cat = json.loads(CATF.read_text()); ra = json.loads(RA.read_text())
        before_cat = json.dumps(cat, sort_keys=True); have_pk = {p['id'] for p in ra['packets']}
        stamp = time.strftime('%Y%m%d_%H%M%S'); receipt = []
        if not DRY:
            shutil.copy2(CATF, A / f'CATALOG.json.bak_{stamp}_swarmwide'); shutil.copy2(RA, RA.with_name(RA.name + f'.bak_{stamp}_swarmwide'))
            CTXD.mkdir(exist_ok=True)
        template = next(p for p in ra['packets'] if p['id'] == 'BT-CTX-Q052')
        for f in fams:
            k = f['key']
            if k in cat: receipt.append({'key': k, 'state': 'already_present'}); continue
            ok = [(n, s) for n, s in f['files'] if s.is_file() and 0 < s.stat().st_size <= BUDGET]
            if not ok: receipt.append({'key': k, 'state': 'SKIP_no_files'}); continue
            if DRY: receipt.append({'key': k, 'state': 'would_add', 'files': [n for n, _ in ok], 'title': f['title']}); continue
            dst = A / 'sources' / k
            if dst.exists(): shutil.rmtree(dst)
            dst.mkdir(parents=True)
            for n, s in ok: shutil.copy2(s, dst / n)
            if f['target'] not in have_pk:
                ctx = CTXD / f'{k}.md'
                ctx.write_text(f"# Target context {k}: {f['title']}\n\nGraph target context for dispatch and feedback routing of the "
                               f"swarm-wide source family {k}, registered {stamp} by {JOB}. Definition-only; admits no science. "
                               "Related native nodes are not yet selected (graph lane).\n")
                p = json.loads(json.dumps(template)); p.update(id=f['target'], claim=f"Target context for {f['title']}",
                    source_file=str(ctx.relative_to(B)), source_line=1, source_sha256=sha(ctx), typed_relations=[{'type': 'related_context_only', 'target': 'bodytwin:' + r, 'target_source_status': 'OPEN'} for r in f.get('relations', [])],
                    topic=f['title'].lower()[:200])
                ra['packets'].append(p); have_pk.add(f['target'])
            cat[k] = {'target_id': f['target'], 'parent': f['parent'], 'context': f['context'],
                      'files': [{'name': n, 'source': str(s), 'sha256': sha(s)} for n, s in ok], 'source_directory': str(dst),
                      'graph_advice': {'task_kind': 'define', 'numerical_dispatch_allowed': False},
                      'registered_by': JOB, 'registered_at': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
            receipt.append({'key': k, 'state': 'added', 'files': [n for n, _ in ok]})
        if not DRY:
            for path, data in ((RA, ra), (CATF, cat)):
                tmp = path.with_name(path.name + '.tmp'); tmp.write_text(json.dumps(data, ensure_ascii=False, indent=1) + '\n'); os.replace(tmp, path)
            old = json.loads(before_cat); new = json.loads(CATF.read_text())
            assert all(new[k] == old[k] for k in old), "existing keys changed"
    bad = []
    if not DRY:
        for f in fams:
            if f['key'] in cat and cat[f['key']].get('registered_by') == JOB:
                r = subprocess.run([str(B / 'graph'), 'working', 'packet', '--id', f['target']], cwd=B, capture_output=True, text=True)
                if r.returncode != 0: bad.append(f['target'])
    out = {'job': JOB, 'dry': DRY, 'keys_after': len(cat), 'added': [r['key'] for r in receipt if r['state'] == 'added'],
           'packets_unresolved': bad, 'receipt': receipt}
    if not DRY: (B / 'tasks/build_night' / f'SWARMWIDE_RECEIPT_{stamp}.json').write_text(json.dumps(out, ensure_ascii=False, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != 'receipt'}, ensure_ascii=False))
    for r in receipt: print(' ', r.get('key'), r['state'], r.get('title', ''))

if __name__ == '__main__':
    main()
