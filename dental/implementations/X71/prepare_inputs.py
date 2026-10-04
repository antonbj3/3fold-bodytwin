"""Freeze small source evidence and an explicit, editable acquisition contract.

This is a source inventory, not laboratory measurement or numerical optimization.
"""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import csv, datetime, hashlib, json, re
HERE = Path(__file__).resolve().parent
DENT = HERE.parents[1]
PKG = DENT / 'results/DEMO48_PACKAGE'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def save(name, obj):
    p = HERE / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')

def main():
    if (HERE / 'INPUT_LOCK.json').exists():
        raise SystemExit('Already frozen; do not overwrite. Create a new version.')
    demos = json.loads((PKG / 'demos.json').read_text())['demos']
    sources = []
    loc = {}
    for d in demos:
        p = PKG / d['copy_dir']
        loc[d['id']] = p
        for name in ['results.json', 'README_DEMO.md']:
            if (p / name).exists():
                sources.append((d['id'], p / name))
    for (id_, rel) in [('X26', 'batch3/demos/X26'), ('X24', 'batch3/demos/X24'), ('X25', 'batch3/demos/X25'), ('X31', 'batch3/demos/X31'), ('X38', 'batch7/demos/X38')]:
        loc[id_] = PKG / rel
        for name in ['results.json', 'README_DEMO.md']:
            if (loc[id_] / name).exists():
                sources.append((id_, loc[id_] / name))
    reviews = []
    for n in range(10, 17):
        for p in sorted((DENT / f'results/LANE_XREVIEW_BATCH{n}').glob('REVIEW_*.json')):
            if p.name == 'REVIEW_INDEX.json':
                continue
            d = json.loads(p.read_text())
            reviews.append({'path': str(p), 'data': d})
            sources.append(('REVIEW', p))
            if d['job'].startswith('BT-'):
                for name in ['results.json', 'RESULTS.md']:
                    q = DENT / 'results' / d['job'] / name
                    if q.exists():
                        sources.append((d['job'], q))
                        loc[d['job']] = q.parent
    save('raw/REVIEW_RECORDS.json', reviews)
    extra = [PKG / 'MEASUREMENT_PLAN.md', PKG / 'demos.json', DENT / 'results/LANE_X26_DECIDABILITY/DECIDABILITY_TABLE.csv', loc['X59'] / 'FACIT.csv', loc['X59'] / 'raw/cunali2017.pdf', loc['X59'] / 'raw/cunali2017.txt', loc['X1B'] / 'LAB_PROTOCOL.md', loc['X47'] / 'LAB_PROTOCOL_RETENTION_AGING.md', DENT / 'notes/OLD_DENTAL_CELL_INDEX.jsonl', Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/bodytwin/scripts/physics_exp/cbct_multidevice_sensitometry_overdet.py'))]
    sources += [('CONTRACT', p) for p in extra if p.exists()]
    manifest = []
    unknown = []
    for (job, p) in sources:
        manifest.append({'job': job, 'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size})
        if p.suffix != '.json':
            continue

        def walk(o, ptr=''):
            if isinstance(o, dict):
                for (k, v) in o.items():
                    walk(v, ptr + '/' + str(k).replace('~', '~0').replace('/', '~1'))
            elif isinstance(o, list):
                for (i, v) in enumerate(o):
                    walk(v, ptr + '/' + str(i))
            elif isinstance(o, str) and re.search('\\bUNKNOWN\\b|NOT_MEASURED|NOT_PERFORMED|NOT_RUN', o, re.I):
                unknown.append({'job': job, 'path': str(p), 'pointer': ptr, 'text': o, 'source_sha256': sha(p)})
        walk(json.loads(p.read_text()))
    save('SOURCE_MANIFEST.json', manifest)
    save('raw/UNKNOWN_OCCURRENCES.json', unknown)
    decisions = []

    def add(id_, jobs, label, measurements=(), other=(), chain='crown', quantity='', level='PER_TOOTH', time='SIMULTANEOUS', criterion='Prospective gate required', scope='matched laboratory specimen'):
        paths = [str(loc[j] / 'README_DEMO.md') for j in jobs if j in loc]
        decisions.append(dict(id=id_, source_jobs=jobs, label=label, requires=list(measurements), other_blockers=list(other), chain=chain, quantity=quantity or label, resolution_level=level, timescale=time, criterion=criterion, scope=scope, evidence_locators=paths, physical_status='UNKNOWN', claim_type='information_link'))
    add('D01', ['X10', 'X55'], 'Signerad as-built-avvikelse i datumram', ['scan'], quantity='signed normal deviation, um', level='PER_POINT', criterion='Independent length/registration error included; no best-fit erasure')
    add('D02', ['X10', 'X1B'], 'Dry seating and margin decisions', ['scan', 'film'], quantity='vertical marginal height and solid clearance, um', level='PER_POINT', criterion='X1B distance bias ±30 um; X10 local bound is specimen specific')
    add('D03', ['X13', 'X49', 'GENCAD_V4'], 'Real cement film field on current crowns', ['scan', 'film'], quantity='assembled film h(x), um', level='PER_POINT', time='HANDOVER', criterion='Same-state calibrated film contrast and registered positions')
    add('D04', ['X13'], 'CAD -spacer to regional gap on held out batch', ['scan', 'film', 'spacer40'], chain='cement', level='PER_SURFACE_REGION', criterion='20 um error / original 90% coverage; 40 pilot crowns, not population certification')
    add('D05', ['X13', 'X26'], 'Setting dynamics and hydraulic film', ['scan', 'film', 'rheology'], other=['hydraulic_model_validation'], chain='cement', level='PER_POINT', time='HANDOVER', criterion='New pressure/height prediction must be frozen; static film alone insufficient')
    add('D06', ['X59'], 'Separate Replicator from Scaffolding Difference', ['film', 'bridge10'], chain='cement', level='PER_SURFACE_REGION', criterion='Four same-state observations; >=10 independent calibration and >=10 held specimens; 20 um gate')
    add('D07', ['X55', 'X51'], 'Separate the scannerbias from manufacturing on the same part', ['scan'], quantity='independently referenced signed shape, um', level='PER_POINT', criterion='Independent datum/reference required; repeats do not remove common bias')
    add('D08', ['X14', 'X61'], 'Physical sinter scale and regional QC update', ['scan'], other=['green_sinter_pair'], chain='process', level='PER_POINT', time='HANDOVER', criterion='Independent pre/post length gauge plus correlated scan uncertainty; missing green state cannot be recreated on ready crowns')
    add('D09', ['X53'], 'Physical trial of milling access witnesses', ['scan'], other=['matched_CAM_toolpath'], chain='process', level='PER_POINT', criterion='Exact tool/stock/pose witness registration; source digital replay is not milling')
    add('D10', ['X1B'], 'Break load and origin of tested crowns generated', ['fracture'], quantity='F_peak N; blinded failure mode and origin', time='HANDOVER', criterion='Trace calibration and origin; observation of tested specimens, no transfer to untested designs')
    add('D11', ['X1B'], 'Thickness-angle contrast Q on the entire 72 test', ['fracture72'], time='HANDOVER', criterion='Q95CI vs [0.351426,0.640340] and control [0.740818,1.349859]; 12/cell; incompatible origins block PASS')
    add('D12', ['X1B', 'X60', 'X26'], 'Absolute breaking load for new held out design', ['scan', 'film', 'fracture', 'contact'], other=['flaw_population', 'independent_held_batch', 'mesh_enclosure', 'geometry_gate_repair'], time='HANDOVER', criterion='Own generated design + new frozen absolute predictions; literature strength does not suffice')
    add('D13', ['X61', 'X51'], 'Real model/process alarm from the crown series', ['scan', 'film', 'fracture'], other=['sealed_empirical_likelihood'], chain='QC', criterion='Per-record score before update; prospective calibrated likelihood needed')
    add('D14', ['X47'], 'The Age Contrast of the Retention in the Current Town', ['retention24'], chain='retention', time='HANDOVER', criterion='12 water +12 aged; state-specific observations; new protocol, no clinical-year mapping')
    add('D15', ['X47'], 'Absolute retention for new batch', ['retention24'], other=['retention_held_batch'], chain='retention', time='HANDOVER', criterion='No Ti to composite transfer; matched grip/mode and independent held batch')
    add('D16', ['X18', 'X26'], 'Contact force map on matched crown/ceiling test', ['contact'], quantity='nonnegative registered patch fractions and total N', level='PER_SURFACE_REGION', criterion='100 N source roof; simultaneous interval errors incl. registration; new crown needs new freeze')
    add('D17', ['X18', 'X54'], 'Hold elastic response at measured contact', ['contact'], quantity='registered strain or force-displacement', level='PER_POINT', criterion='Own specimen geometry/material/support, independent DIC response; no interior peak inference')
    add('D18', ['X18', 'X26'], 'Physical local crown voltage under known load', ['contact', 'film'], other=['mesh_enclosure', 'validated_stress_operator'], level='PER_POINT', criterion='Conditional fraction budgets 0.1/0.2 are not universal force accuracy')
    add('D19', ['X54'], 'Tooth rests/PDL to regional power profile', ['pdl'], other=['individual_paired_anatomy'], chain='orthodontic', level='PER_SURFACE_REGION', criterion='Multiple physical forces/preloads and compliance, not aggregate population force')
    add('D20', ['X63'], 'Four supported contact forces after tightening', ['four_support'], quantity='four simultaneous reaction forces N and gaps um', chain='implant', level='PER_POINT', criterion='64 source support predictions; matched fixture and calibrated force sensors')
    add('D21', ['X63'], 'Reuse history to clamp force', ['four_support'], chain='implant', time='HANDOVER', criterion='Measured individual preload before/after specified torque events; event count is not time')
    add('D22', ['X56'], 'Distinguish local bone – implant - sliding from system movement', ['slip'], chain='implant', level='PER_SURFACE_REGION', criterion='Simultaneous local slip plus total compliance, same region/direction/load')
    add('D23', ['X56'], 'Biologically valid micromotion limit', ['slip'], other=['biological_motion_outcome'], chain='biology', level='PER_SURFACE_REGION', time='HANDOVER', criterion='50–150 um is not a calibrated local biological endpoint')
    add('D24', ['X65', 'X25'], 'Transfer fatigue from R=0.1 to new R', ['four_support'], other=['critical_root_stress', 'residual_stress', 'matched_fatigue_second_R'], chain='implant', level='PER_SURFACE_REGION', time='HANDOVER', criterion='Root-region stress, not shaft foil gauge; complete horizon and failure modes')
    add('D25', ['X24', 'X26'], 'CBCT -Attitude at held out position', ['bone'], chain='bone', level='PER_SURFACE_REGION', criterion='HA fit levels 0/0.1/0.3, 0.2 held; original 40 HU_ref gate is not modulus accuracy')
    add('D26', ['X24', 'X26'], 'Same-Region Mandibular Module', ['bone'], chain='bone', level='PER_SURFACE_REGION', criterion='Directional strain-range mechanics + density/fabric; holdout specimens')
    add('D27', ['X5', 'X58', 'X8'], 'Physical local channel edge and poseminimum', ['anatomy'], chain='anatomy', level='PER_POINT', criterion='Independent fiducial-matched anatomy; annotation revision alone is not physical truth')
    add('D28', ['X20', 'X26'], 'Physical bone -/sine limit and realized tool pose', ['anatomy'], chain='anatomy', level='PER_POINT', time='HANDOVER', criterion='Distinct maxillary stratum; test ±0.60 mm scenario, no clinical choice')
    add('D29', ['X15'], 'Trajectory to later bone response', ['anatomy'], other=['longitudinal_response'], chain='biology', level='PER_SURFACE_REGION', time='HANDOVER', criterion='Static anatomy cannot identify remodeling; week-one impossible')
    add('D30', ['X31', 'X12', 'X9'], 'Remaining dentin on the same horn/ray', ['anatomy'], chain='anatomy', level='PER_POINT', criterion='Same ray T,E,r; horn coverage and boundary error; no pulp-risk transfer')
    add('D31', ['X48', 'X19', 'X62'], 'Orthodontic force and bond-slip at signed pose', ['pdl'], chain='orthodontic', criterion='Physical multidirectional response and bond-slip, not two secants or nominal material')
    add('D32', ['X62'], 'Actual whole group of bracket pose exceedances', ['pdl'], chain='orthodontic', level='PER_ARCH', criterion='Full actual specimen poses relative frozen physical target centres; moments do not certify tails')
    add('D33', ['X64'], 'Optical response by current curved crown', ['optics'], chain='material', level='PER_POINT', criterion='Same spatial ray and optical operator; 236 source points not one mean')
    add('D34', ['X64'], 'Post curing to arrest mood on resin samples', ['optics'], chain='material', time='HANDOVER', criterion='16/60 min x axial/oblique; different material from 3Y crown72')
    add('D35', ['X50', 'X7'], 'Assessment noise floor for bite assessment', other=['blinded_independent_raters'], chain='diagnostic', criterion='Independent masked duplicate ratings; lab mechanical assay cannot supply this')
    add('D36', ['PATIENT360_R2'], 'Joint physical IOS – CBCT - frame', other=['paired_same_case_IOS_CBCT'], chain='anatomy', level='PER_POINT', criterion='Matched physical patient/frame; separate cases do not form a paired measurement')
    add('D37', ['X57'], 'Local periodontal etiology and stadium foundations', other=['independent_periodontal_truth'], chain='diagnostic', level='PER_SURFACE_REGION', criterion='Diagnostic rule correctness is not physical truth; no clinical acquisition here')
    add('D38', ['X57'], 'Lesion barrier and biological pulp function', other=['independent_lesion_and_pulp_truth'], chain='diagnostic', level='PER_POINT', criterion='Static label/deep text cannot supply biological vitality')
    add('D39', ['X45'], 'ISQ to local BIC and healing response', other=['local_BIC_longitudinal'], chain='biology', level='PER_SURFACE_REGION', time='HANDOVER', criterion='ISQ alone does not identify local tissue parameters')
    add('D40', ['X45'], "Coupon distribution to the crown's aged flaw fields", ['optics'], other=['matched_crown_surface_hazard'], chain='material', level='PER_POINT', time='HANDOVER', criterion='Coupon m and scale are not a validated crown hazard field')
    add('D41', ['X45'], 'Wear to volume/contact change', other=['wear_density_contact_path'], chain='material', time='HANDOVER', criterion='Mass is not volume without density; cumulative history retained')
    add('D42', ['X3', 'GENCAD_V2', 'GENCAD_V5', 'X42', 'PROOF_LANE', 'X60'], 'Digital geometry to physical generated crown', ['scan', 'film', 'fracture'], other=['design_matched_validation', 'geometry_gate_repair'], criterion='Crown72 is not every generated design/bridge; no blanket transfer')
    add('D43', ['X2', 'X21'], 'Geometric antagonist of physical contact force', ['contact'], other=['matched_antagonist_case'], criterion='Contact area/spherical antagonist is insufficient force information')
    add('D44', ['X4'], 'Mandible plate for physical fatigue', other=['manufactured_plate_fatigue'], chain='anatomy', time='HANDOVER', criterion='Crown specimen pool incompatible with plate fatigue')
    add('D45', ['X17'], 'Airway geometry to physiological flow', other=['airway_boundary_conditions'], chain='diagnostic', criterion='No new physiological measurements or clinical recommendations')
    add('D46', ['X22'], 'Synthetic caries to real lesion response', other=['physical_lesion_ground_truth'], chain='diagnostic', criterion='Own fixture does not establish lesion truth')
    for rr in reviews:
        j = rr['data']['job']
        if j.startswith('BT-'):
            decisions.append(dict(id='LEGACY:' + j, source_jobs=[j], label=rr['data'].get('correction') or rr['data']['review_scope'], requires=[], other_blockers=['reviewed_legacy_law_or_scope'], chain='legacy', quantity='review-dependent model obligation', resolution_level='PHENOMENOLOGICAL', timescale='HANDOVER', criterion='Repair rejected law/scope before a matched acquisition is valued', scope='reviewed legacy artifact', evidence_locators=[rr['path']], physical_status='UNKNOWN', claim_type='information_link', count_in_value=False))
    save('inputs/DECISIONS.json', decisions)
    measures = []

    def m(id_, name, supplies, hours, equipment, n, compatible, note, machine=(0, 0), new=0, prereq=()):
        measures.append(dict(id=id_, name=name, supplies=supplies, operator_hours_interval=hours, instrument_hours_interval=list(machine), equipment=equipment, specimens=n, new_specimens=new, crown72_compatible=compatible, notes=note, prerequisites=list(prereq), cost_status='PHENOMENOLOGICAL_PLANNING_ASSUMPTION', cost_replacement='Lab quote and timestamped pilot time log; equipment access and money UNKNOWN', preparation='Included in operator range; shared manufacture separately', fit='Included analysis; no physical fit performed', discovery='Source review recorded separately; labour UNKNOWN', validation='Held new batch excluded unless explicitly supplied', queries='Included export/scoring; software time measured', fallback='Separate budget; no automatic new batch'))
    m('M01', 'Registered as-built + independent length metrology, pilot 12', ['scan'], [6, 10], ['scanner', 'independent length and surface reference/CMM'], 12, True, 'Two/cell of original72; preserve IDs and all signed points. Independent surface reference required to separate scanner and manufacturing. Larger72 programme requires new cost quote.', machine=(2, 4))
    m('M02', 'Dry seating + cement film μCT, same pilot 12', ['film'], [6, 10], ['microCT', 'seating force/height fixture', 'contrast thickness reference'], 12, True, 'Matched dry and assembled state; CT visibility/dose and cleaning effect must be checked; no replica default correction.', machine=(6, 12), prereq=['M01'])
    m('M03', 'Cementrheogram + force – height – time', ['rheology'], [3, 6], ['rheometer', 'calibrated force-height logger'], 0, True, 'Uses same cement aliquots and pilot assembly. It does not validate hydraulic model by itself.', machine=(2, 4), prereq=['M02'])
    m('M04', 'Break load, stiffness and blind origin, pilot 12', ['fracture'], [6, 10], ['calibrated universal test machine', 'fractography microscope'], 12, True, 'Two/cell count toward original72 only if all interventions/protocol match; Q gate waits for12/cell.', machine=(2, 4), prereq=['M02'])
    m('M05', 'Full thickness-angle test72', ['fracture', 'fracture72'], [24, 40], ['test machine', 'fractography', 'independent die coupons'], 72, True, 'Replaces M04; full Q protocol12/cell. Independent die coupons are extra, number UNKNOWN.', machine=(12, 20), prereq=['M02'])
    m('M06', 'Retention water/ aging 24', ['retention24'], [12, 20], ['thermocycler', 'pull-off grip', 'test machine'], 24, 'CONFLICT', 'X47 proposes48 fracture +24 retention. Changes original72 Q trial; pull-off damages the endpoint state. Aging elapsed time UNKNOWN.', machine=(20, 60))
    m('M07', 'Registrerad patchkraft + oberoende elastiskt svar', ['contact'], [8, 16], ['regional force sensors', 'DIC', 'matched fixture'], 0, 'CONDITIONAL', 'Needs its own matched roof/crown preparation and force registration. Low-load tests on72 require non-damage proof and new freeze.', machine=(3, 8), new=4)
    m('M08', 'Contact forces in four-support fixture + preload history', ['four_support'], [10, 18], ['four force transducers', 'torque instrument', 'registered four-support bridge'], 0, False, 'Four-support bridge and screws are not individual crowns. Four supports are channels, not four independent specimens.', machine=(4, 8), new=3)
    m('M09', 'Local bone – implant - sliding + system compliance', ['slip'], [14, 24], ['DIC/local displacement', 'bone/implant fixture', 'force measurement'], 0, False, 'Local slip and total system motion measured simultaneously; no biological threshold certification.', machine=(4, 8), new=3)
    m('M10', 'HA + same-region bone mechanics', ['bone'], [18, 30], ['CBCT', 'HA phantom', 'directional compression', 'density/fabric reference'], 0, False, 'Separate bone coupons; final n UNKNOWN. Three pilot regions do not certify a bone law.', machine=(8, 16), new=3)
    m('M11', 'Fiducial-anatomi + realiserad pose', ['anatomy'], [20, 36], ['independent high-res/section anatomy', 'CBCT', 'pose reference'], 0, False, 'Three mandibular and three maxillary scouts; dentin/horns require extracted teeth. Anatomical access/licence UNKNOWN.', machine=(10, 20), new=6)
    m('M12', 'PDL /retainer forces + full signed pose', ['pdl'], [14, 24], ['multiaxis reaction fixture', 'pose metrology', 'bond-slip tracking'], 0, False, 'Matched orthodontic assembly; linear reciprocity must be tested. No patient force claim.', machine=(6, 12), new=3)
    m('M13', 'Optical beams + curing mechanics', ['optics'], [16, 28], ['spectrophotometer', 'resin cure', 'test machine'], 0, False, 'Resin and curved optical specimens incompatible with3Y crowns; same-point support required.', machine=(8, 16), new=8)
    m('M14', 'X13 spacer-/batchpilot40', ['spacer40'], [18, 32], ['manufacture', 'scanner', 'microCT'], 0, False, 'Two settings20 plus held setting10 and held batch10. Cannot reuse unchanged single-spacer72 as40 varied-spacer pieces.', machine=(20, 40), new=40, prereq=['M02'])
    m('M15', 'X59 four state , calibration 10 + held out 10', ['bridge10'], [12, 20], ['validated silicone-film CT', 'microscopy', 'same-state fiducials'], 0, 'CONDITIONAL', '20 independent copings. Replica/cleaning sequence may alter72 fracture state; use separate scout pieces until noninterference is shown.', machine=(12, 24), new=20, prereq=['M02'])
    save('inputs/MEASUREMENTS.json', measures)
    save('inputs/LAB_ASSUMPTIONS.json', dict(operator_hours_week=40, instrument_hours_week=40, ready_crowns_at_week_start='UNKNOWN', equipment_access='UNKNOWN', currency_cost='UNKNOWN', shared_manufacture_operator_hours=[12, 24], shared_manufacture_machine_hours=[24, 60], shared_sinter_elapsed_hours=[6, 18], cost_status='PHENOMENOLOGICAL', week_definition='One operator, five8h days. Conditional ready-specimen schedule; scratch manufacture separately tested.', chain_semantics='All explicitly modeled measured ports plus all other blockers; not clinical certification'))
    rows = list(csv.DictReader((loc['X59'] / 'FACIT.csv').open()))
    save('inputs/EXTERNAL_FACIT.json', [x for x in rows if x['source'] == 'Cunali2017'])
    idjobs = {j: [d['id'] for d in decisions if j in d['source_jobs']] for j in {x['job'] for x in unknown}}
    rejected = []
    linked = []
    for (i, x) in enumerate(unknown):
        x['occurrence_id'] = i
        p = x['pointer'].lower()
        if x['job'] in ['REVIEW', 'CONTRACT'] or any((k in p for k in ['/model/name', '/model/value', '/model_provenance', '/licenses', '/licence', '/license', '/tokens', '/full_cost', '/cost', '/environment', '/permissions', '/route_metadata'])):
            x['reason'] = 'metadata_cost_license_or_review_wording; not a counted physical decision'
            rejected.append(x)
        else:
            x['decision_obligations'] = idjobs.get(x['job'], [])
            x['mapping_status'] = 'SOURCE_LEVEL_OBLIGATIONS_NOT_ONE_TO_ONE_LEAF_IDENTIFICATION' if x['decision_obligations'] else 'UNMAPPED_RETAINED_NOT_COUNTED'
            linked.append(x)
    save('raw/UNKNOWN_LINKED.json', linked)
    save('raw/UNKNOWN_REJECTED.json', rejected)
    save('raw/INVENTORY_SCOPE.json', dict(review_files=len(reviews), package_index_entries=len(demos), source_occurrences=len(unknown), linked_occurrences=len(linked), excluded_occurrences=len(rejected), unresolved_leaf_mapping=sum((x['mapping_status'].startswith('UNMAPPED') for x in linked)), limitations='Every literal UNKNOWN retained. Source-level debts below are curated, not a complete one-to-one physical-leaf parser. Repeated outputs do not create independent decisions. All unclassified/legacy debts remain UNKNOWN and receive zero value.'))
    locks = [HERE / 'inputs' / n for n in ['DECISIONS.json', 'MEASUREMENTS.json', 'LAB_ASSUMPTIONS.json', 'EXTERNAL_FACIT.json']]
    save('INPUT_LOCK.json', dict(created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), inputs=[{'path': str(p.relative_to(HERE)), 'sha256': sha(p)} for p in locks], manifest_sha256=sha(HERE / 'SOURCE_MANIFEST.json'), snapshot_mode='Small scalar contracts; source arrays not copied'))
    print('Frozen inputs; see raw/INVENTORY_SCOPE.json. No laboratory measurements performed.')
if __name__ == '__main__':
    main()
