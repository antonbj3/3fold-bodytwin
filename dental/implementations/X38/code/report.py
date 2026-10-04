"""Generate final reader-facing reports from preserved raw results."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import datetime, hashlib, json, os, subprocess
R = Path(__file__).resolve().parents[1]

def read(p):
    return json.loads((R / p).read_text())

def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()

def save(p, x):
    (R / p).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def main():
    r = read('rounds/R3.json')
    s = read(r['statistics_path'])
    r1 = read('rounds/R1.json')
    env = os.environ.copy()
    env['PYTHONNOUSERSITE'] = '1'
    subprocess.run(['python3', str(R / 'code/plot.py')], check=True, env=env, cwd=R)
    maxerr = max((max(q['max_correspondence_mm'], q.get('published_trimesh_max_error_mm', 0)) for m in r['meshes'] for q in m['checks']))
    stats = {k: {a: b for (a, b) in v.items() if a != 'search_curve'} for (k, v) in s.items() if k in ['plan29', 'plan42', 'actual100N', 'binary', 'continuous', 'all_observation_screening']}
    result = dict(lane='X38-export-gate', claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', capability='Frozen prediction-bearing mesh transport and fixed laboratory experiment planning', export=dict(meshes=11, conversion_paths=55, transport_passes=55, refused_sources=r['refused_sources'], maximum_correspondence_mm=maxerr, resolution='PER_POINT', tolerance_mm=1e-05, region_probe=r['region_probe'], screening=r['screening'], physical_manufacture='UNKNOWN', vendor_CAM_interoperability='UNKNOWN', anatomical_region_validation='UNKNOWN', sinter_factor_calibration='UNKNOWN'), statistical_answers=stats, controls=dict(unittest_log='raw/tests_latest.log', published_geometry='trimesh imports and exports', exact_confidence='SciPy binomtest against beta/closed form', binary_power='all4900 count tables at69/group checked with SciPy fisher_exact', continuous_power='statsmodels finite outputs plus alternate conditional integral for all12 declared cases', wrong_injected_values_refute=True), prereg_scope_audit='raw/PREREG_SCOPE_AUDIT.json', interface_R4=dict(source_hash_join='checked before output creation', matching_report_verdict_preserved='FAIL', actual_X34_X18_report_on_other_mesh='REFUSED', tests_log='raw/tests_latest.log'), numerical_demo_code_lock='revisions/R3/CODE_LOCK.json', negative_results=[dict(round='R2', outcome='FAIL', reason='all 15 criteria failed:4 source facets below unchanged numerical degeneracy cutoff; R3 exports11 and refuses4', artifact='rounds/R2.json'), dict(round='R1', wrong_region_facets=r1['wrong_region_facets'], fraction=r1['region_disagreement_fraction'], resolution='PER_SURFACE_REGION', cause=r1['cause']), dict(round='STATS_R1', outcome='FAIL', reason='Published SciPy/statsmodels noncentral-t NaN at n 80 effect1.5; mixture fallback independently checked under unchanged tolerance', artifact='rounds/STATS_R1_FAILURE.json'), dict(outcome='UNAVAILABLE', reason='Installed libxml2 refused original3MF appendix XSD maxOccurs; no full XSD conformance claim', artifact='raw/SCHEMA_CHECK_FAILURE.json')], external_referent=dict(kind='published_code', locator=['https://trimesh.org/trimesh.exchange.stl.html', 'https://trimesh.org/trimesh.exchange.threemf.html', 'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html', 'https://www.statsmodels.org/stable/generated/statsmodels.stats.power.TTestIndPower.html'], compared_quantity='Imported physical triangles and fixed two-design test rejection probability', refutes_us=True), external_referents=[dict(kind='closed_form', locator='https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/exacbici.htm', compared_quantity='exact binomial one-sided failure-probability upper confidence limit', refutes_us=True), dict(kind='closed_form', locator='https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nct.html', compared_quantity='Gaussian/chi-square mixture for noncentral t power', refutes_us=True), dict(kind='published_dataset', locator='https://doi.org/10.5281/zenodo.22306621', compared_quantity='36 specimen cycle/event observations;27 complete horizon outcomes and9 early censors', refutes_us=True), dict(kind='our_own_fixture', locator='rounds/R1.json; raw/runs labelled_transport_probe', compared_quantity='Synthetic transport regions only, not physical/anatomical validation', refutes_us=True)], full_cost=dict(r['full_cost'], preserved_failed_and_preparation_cost=read('raw/FULL_COST_EXTRA.json')), run_directory=r['run_directory'], artifacts=[dict(path=p, sha256=sha(p)) for p in ['INPUT_LOCK.json', 'CODE_LOCK.json', 'FROZEN_PREDICTIONS.json', 'rounds/R1.json', 'rounds/R2.json', 'rounds/R3.json', r['statistics_path']]], debts=[dict(quantity='sinter factor', resolution='PHENOMENOLOGICAL', replacement_measurement='same-batch unsintered and sintered lengths in3 axes'), dict(quantity='continuous standardized effect/common SD', resolution='PHENOMENOLOGICAL', replacement_measurement='pilot same-protocol failure-force distribution from independent specimens'), dict(quantity='binary planning risks pA,pB', resolution='PHENOMENOLOGICAL', replacement_measurement='prespecified matched-horizon full events in each design')], edges=[dict(producer='X1b geometry and frozen prediction artifact', consumer='lab fabricated specimen identity', resolution='PER_POINT', time_scale='HANDOVER', status='transport checked; physical manufacture UNKNOWN'), dict(producer='lab region annotation', consumer='roundtrip region check', resolution='PER_SURFACE_REGION', time_scale='SIMULTANEOUS', status='UNKNOWN for actual anatomy'), dict(producer='X25 time/event observations', consumer='DENT-IF-IA-ASSEMBLY-FATIGUE', resolution='PER_TOOTH', time_scale='SIMULTANEOUS', status='PENDING_INDEPENDENT_REVIEW')], clinical_recommendation=False, scientific_admission=False)
    save('results.json', result)
    readme = f"""# Export and specimen counts for a dental research laboratory

Run **`./run_all.sh`**. Two tools preserve the link between design, frozen predictions and laboratory specimens. Export gate rereads geometry/metadata; specimen planner distinguishes risk limits from power. Python and packages in `RUNTIME_LOCK.json` must already exist. No GPU/network required.

| Question | Executable answer | Resolution |
|---|---|---|
| Are 11 of 15 X1b surfaces preserved through five file paths each? | 55/55 transport tests pass; maximum error {maxerr:.3g} mm, tolerance 1e-5 mm | PER_POINT |
| Are test regions preserved when triangles are reordered? | {r['region_probe']['mismatches']} errors in {r['region_probe']['facets']} triangles | PER_SURFACE_REGION; synthetic labels |
| What do 3 survivors at 100 N/5 million cycles establish? | Upper one-sided 95% fracture risk **63.16%** | POPULATION; source specimens PER_TOOTH |
| How many needed for risk<10% at 95% confidence? | **29** if all survive; **42 per condition** for 4 simultaneous conditions | POPULATION |
| 30% versus 10% risk, two-sided Fisher alpha 0.05 and 80% power? | **69 per design**, total 138; power 80.73% | POPULATION; assumed risks |
| Mean difference/common SD=0.8, two-sided t test, 80% power? | **26 per design**, total 52; power 80.75% | POPULATION; normality assumption |

![Results](figures/export_lab.png)

## Export a laboratory design

```bash
./exportgate export inputs/exports/D1/crown.stl --out my_export \\
  --stl-unit mm --material "3Y zirconia; verify product/batch" \\
  --sinter-factor 1 --geometry-state final_sintered \\
  --predictions inputs/X1B_FROZEN_PREDICTIONS.json
```

`my_export` contains binary STL, ASCII STL, 3MF, a JSON sidecar for each, frozen predictions and a receipt. Sinter factor explicitly means **manufacturing length/final sintered length**; the tool does not alter geometry. Example 1 is an uncalibrated transport marker. Supply laboratory validated factor/correct state; CAM operator should verify compensation is not applied twice.

For regions add `--regions my_regions.json`. File contains `{{"labels":["intaglio",...],"semantics":"laboratory independent annotation"}}`, one label per triangle in input order. X1b lacks anatomical labels; default `unclassified_surface`. Test labels in first figure panel synthetically partition the existing surface.

Save sidecar SHA256 from export receipt **before** CAM conversion. For a converted file:

```bash
./exportgate check converted.stl --sidecar my_export/model.stl.json \\
  --expected-sidecar-sha256 SAVED_SHA256 --stl-unit mm --converted
```

STL has no standardized unit, so explicit units/sidecar required. 3MF uses millimetres and embedded region sets; if converter removes metadata, recover from sidecar after unique matching of every oriented physical triangle. Changed scale, geometry, normal direction, regions, material or frozen files rejected. Remeshing changing triangulation rejected; nearest-surface labeling cannot guarantee regions here. Only one built mesh object supported. Exit 0=transport PASS, 2=REFUSED.

`--design-gate report.json` requires report `inputs` to contain current source file SHA256. Another crown X34 report rejected before export. Original report PASS/FAIL/UNKNOWN retained as separate evidence. Transport PASS is not manufacturing approval.

## Plan specimen count before the experiment

```bash
./labcount plan --risk .1 --confidence .95 --true-risk .1
./labcount plan --risk .1 --confidence .95 --family-tests 4
./labcount endpoints --data inputs/specimens.json --horizon-cycles 5000000 \\
  --configuration 'MUSLA 04020' --force-N 100
./labcount compare --endpoint binary --p-a .30 --p-b .10 --power .80
./labcount compare --endpoint continuous --effect .80 --power .80
```

`plan` uses a fixed sample plan and one-sided exact binomial bound. `--max-failures 1` can allow one fracture, requiring more specimens. 29 is a conditional minimum: if true risk is 10%, acceptance probability only **4.71%**. Confidence is not power. To improve acceptance probability for a good design, choose assumed true risk and permitted failures before testing; `--true-risk` reports acceptance probability.

`endpoints` retains early censoring as UNKNOWN, not failure/survival. Of 36 published specimens, 9 are censored early (**25% missing known horizon outcomes**, no rows discarded). Upper risk bound uses worst possible completion. 100 N group has 3 known survivors. Explicit selection rejects mixed loads/configurations and reports excluded count. No lifetimes extrapolated. Duplicate IDs and missing units/locators rejected.

`compare` assumes independent unpaired groups and a preselected two-sided test. Binary power sums all count tables using exact Fisher test. Continuous fracture loads use `effect=|meanA-meanB|/common SD`; normality/common SD are empirical assumptions. Table effects are planning examples, not MEASURED improvements. Exit 3 means requested power not reached within `--max-n`. All risk bounds use fixed plans; sequential stopping uses X25 separate locked-order tool.

## References and failures

Geometry checked against [published trimesh import/export](https://trimesh.org/trimesh.exchange.threemf.html), units/region sets against [3MF consortium specification](https://github.com/3MFConsortium/spec_core/blob/master/3MF%20Core%20Specification.md). Exact risk bound follows [NIST binomial inversion](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/exacbici.htm). Two-group controls: [SciPy Fisher](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html) and [statsmodels t power](https://www.statsmodels.org/stable/generated/statsmodels.stats.power.TTestIndPower.html). Fault injections cover incorrect scale, equal-volume shear, wrong labels/units/hashes, censoring, duplicate counts and false risk/power numbers.

R2 failed when 4 of 15 original surfaces each had one tiny positive-area triangle below numeric degeneracy threshold 1e-12 mm^2 for cross-product norm. Extra cutoff implemented before source test but absent as an explicit R2 PREREG field; R3 freezes it explicitly. Rejections establish tool policy, not actual unmanufacturability. R3 retains threshold: **11 exports, 4 rejections (26.67%)**, no repair. Four are M1/M2 die/preparation. R2 all-15 gate remains FAIL. R1 positional sidecar moved **{100 * r1['region_disagreement_fraction']:.2f}%** of labels despite unchanged geometry; construction rejected/preserved. Continuous power gave NaN in both SciPy/statsmodels at n80/effect1.5. After frozen followup use [normal/chi-square representation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nct.html), checked by reversed integration order. 15 tests pass. Published methods give same mathematical answers within frozen tolerances; no algorithmic improvement claimed.

**Limitations:** no exocad/hyperDENT/slicer or physical manufacturing executed; anatomical regions unvalidated. Sinter factor, batch dependence, correlation and load protocol must come from laboratory. Full self-intersection checks and general 3MF component import absent. Original 3MF appendix XSD could not compile under installed libxml2; full schema conformity UNKNOWN. Transport/specimen planning executable now; physical performance UNKNOWN.

Local data: selected X1b STS-3D-derived research meshes, **CC BY4.0**, [STS-Tooth3D dataset](https://doi.org/10.5281/zenodo.10597292) and [data description](https://doi.org/10.1038/s41597-024-04306-9), retain earlier anonymized provenance; no original volumes/patient information read. Preparations are virtual laboratory surfaces. [Zenodo22306621](https://doi.org/10.5281/zenodo.22306621), Nicolas-Silvente/Velasco-Ortega, **CC BY4.0**, is the only new measurement source used by specimen planner. X34 report retained as earlier local research artifact, not new dataset or release of underlying Bits2Bites data. 3MF specification is separate documentation, not dental data. Exact origins/SHA256 in `INPUT_LOCK.json`.

Next decisive test: laboratory performs ordinary CAM conversion of annotated crown, retains sidecar, checks reimport, manufactures and measures regions against already frozen predictions. No clinical recommendations.
"""
    (R / 'README_DEMO.md').write_text(readme)
    (R / 'RESULTS.md').write_text(f"""X38 delivers two executable CLI tools preserving mesh/prediction identity and answering how to size a fixed laboratory experiment. 11 of 15 existing X1b surfaces passed 55 file paths against published trimesh; maximum correspondence error {maxerr:.12g} mm (PER_POINT), frozen tolerance 1e-5 mm. Test regions recovered with 0 errors (PER_SURFACE_REGION; synthetic spatial annotation, no anatomical validation).

3/3 published survivors at 100 N and 5e6 cycles give 63.159685% upper one-sided 95% risk bound (POPULATION from PER_TOOTH). 29 complete zero-failure specimens give 9.814463% at predetermined load; 42/condition for 4 simultaneous decisions. All 36 published rows retained; 27 known horizon outcomes and 9 early censorings (25%). No censored times extrapolated. Two-sided comparisons with planned 80% power give 69/design for assumed 30% versus 10% risks and 26/design for assumed standardized mean difference 0.8. Conditional mathematical answers, not MEASURED design effects.

References: published trimesh, NIST exact binomial inversion, SciPy Fisher tables, statsmodels t power and published normal/chi-square definition of noncentral t. 36 actual observations with exact specimen locators from https://doi.org/10.5281/zenodo.22306621. All references/quantities in results.json. Wrong geometry/region/unit/risk/power injections rejected.

Negative R2: 4 of 15 sources rejected by unchanged numeric degeneracy threshold; rounds/R2.json. Negative R1: {r1['wrong_region_facets']} of {r1['faces']} regions moved by reordered facets despite identical volume. STATS_R1: library power NaN at n80/effect1.5; newly frozen mixture representation restores finite answer with independent integration order. Original code/log/gate retained. Original 3MF XSD failed compilation; full schema conformity UNKNOWN. Equally informed published statistical controls reproduce results; no method victory claimed.

Transport PASS is not manufacturing approval. Commercial CAM import, actual sinter factor, anatomical regions, self-intersection, physical fit and fracture performance UNKNOWN. README_DEMO.md/HANDOFF.md. Cost: {r['full_cost']['wall_seconds']:.3f} s wall time, RSS {r['full_cost']['peak_process_RSS_KiB'] / 1024:.1f} MiB; manual discovery UNKNOWN. PENDING_INDEPENDENT_REVIEW.
""")
    (R / 'HANDOFF.md').write_text(f"""X38 completed R1 -> R2 -> R3 and STATS_R1 -> STATS_R2. Run ./run_all.sh; exportgate/labcount standalone. Read README_DEMO.md, rounds results and tests/test_tools.py. Final actual raw run: {r['run_directory']}. Frozen hashes in PREREG*.sha256/FROZEN_PREDICTIONS.sha256/INPUT_LOCK/CODE_LOCK.

Delivered: 11/15 X1b surfaces through 55 transport paths; 0 region errors under reordered spatial test partition; STL+3MF+sidecar with explicit units, material, sinter convention, geometry state and frozen prediction copy/hash. Actual demo exports only unclassified_surface. Save sidecar hash outside mutable file before conversion. PASS means transport, not manufacturable/safe design. Bind source-hash-matched X34 report using --design-gate for new lab case; earlier negative thickness decision remains.

Statistics: 3 survivors ->63.16% one-sided95% upper risk; 29 complete zero-failure ->9.81%; 42 for 4 simultaneous conditions. At true risk10%, 29 specimens have 4.71% acceptance probability. Binary pA.3/pB.1 needs 69/group; continuous d.8 needs 26/group at alpha.05 and 80% power. Assumed population effects/SD are PHENOMENOLOGICAL debts. 36 local published raw specimens, 27 known/9 early censored. Duplicate/unit/locator checks executed.

Failure logs: rounds/R2.json and raw/demo_initial.log retain failed all-15 gate; 4 originals still rejected in R3. rounds/R1.json positional regions49.88% wrong; raw/tests_initial.log and revisions/STATS_R1 retain NaN. rounds/STATS_R1_FAILURE.json had incorrect provisional df2 stopping diagnosis; STATS_R1_DIAGNOSIS_CORRECTION identifies actual n80/d1.5. No outcome changed. raw/PREREG_SCOPE_AUDIT.json records 1e-12 cutoff absent as explicit R2 prereg field, first formally frozen in R3; 4 rejections numeric policy, not empirical unmanufacturability. PREREG_STATS_R2 freezes integral representation error budget. Original XSD compiler problem UNKNOWN.

Graph: DENT-IF-IA-ASSEMBLY-FATIGUE dispatched (raw/GRAPH_DISPATCH_RECEIPT.json). Bind specimen-to-horizon predicate at PER_TOOTH/SIMULTANEOUS; consumer aggregates. Export geometry-to-manufactured specimen identity PER_POINT/HANDOVER. Proposed missing coverage DENT-MFG-EXPORT-TRANSPORT-CONTRACT, never invented native dependency. GRAPH_FEEDBACK.json PENDING_INDEPENDENT_REVIEW. Source graphs unchanged. R4 interface requires X34 report source SHA matching export input; 15 tests pass after change. 55 geometry paths executed under preserved revisions/R3/CODE_LOCK.json; only optional report binding changed afterward.

Next construction: replace identical triangulation assumption with laboratory annotated surface field, region boundaries and measured geometry/registration uncertainty, allowing real CAM remeshing rebinding without moving thin margins across region boundaries. Freeze new regional error limits and actual vendor export before numerics. Minimum external reference: laboratory metadata/regions before/after vendor roundtrip, then same-specimen triple scan. Requires missing data, not another power derivation or generic metadata audit.
""")
    feedback = dict(lane=_release_expand('X38'), target_id='DENT-IF-IA-ASSEMBLY-FATIGUE', review_state='PENDING_INDEPENDENT_REVIEW', claim_type='capability', result_file='results/LANE_X38_EXPORT_GATE/RESULTS.md', sha256=sha('RESULTS.md'), measured_quantity='per-specimen horizon event/censor classification and fixed-plan exact risk/power', units='cycles; N; probability; specimen count', resolution='PER_TOOTH', uncertainty='Exact one-sided population limits conditional on independent fixed-plan specimens; empirical protocol and correlation UNKNOWN', population_regime='published36 titanium implant-abutment assemblies at prespecified5e6 cycles; separate assumed unpaired comparison scenarios', preregistered_gate={'prereg': 'results/LANE_X38_EXPORT_GATE/PREREG_STATS.json', 'sha256': sha('PREREG_STATS.json'), 'risk_tolerance': 1e-10, 'power_tolerance': 1e-08}, baseline='SciPy binomtest/fisher_exact; statsmodels power; same-information verification, no method race', outcome='CAPABILITY_DELIVERED;VENDOR_AND_PHYSICAL_VALIDATION_UNKNOWN', negative_result=True, external_referent=result['external_referents'][2], edges=result['edges'], coverage_proposal={'id': 'DENT-MFG-EXPORT-TRANSPORT-CONTRACT', 'scope': 'geometry/prediction transport into physical lab specimen; missing native coverage', 'status': 'PROPOSED_ONLY'})
    save('GRAPH_FEEDBACK.json', feedback)
    save('CURRENT_WORK_STATE.json', dict(lane='X38-export-gate', phase='DELIVERED_PENDING_INDEPENDENT_REVIEW', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), latest_gate='55/55 transport paths;15/15 independent/fault tests; fixed risk and power answers; preserved R1/STATS_R1 failures', next_operation='Independent coordinator review, then real annotated vendor CAM roundtrip and region-preserving remesh representation; physical test inputs absent', claim_type='capability'))
    print('README_DEMO, RESULTS, results.json, figure, HANDOFF and PENDING feedback written.')
if __name__ == '__main__':
    main()
