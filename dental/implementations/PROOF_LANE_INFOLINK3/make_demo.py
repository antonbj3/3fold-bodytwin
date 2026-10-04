"""Bound existing observations to dental consumers; no network or source writes."""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import csv, io, json, hashlib, time, datetime, zipfile, re, contextlib, copy, resource
from pathlib import Path
import numpy as np, pandas as pd, openpyxl
from scipy.stats import rankdata, spearmanr, mannwhitneyu
H = Path(__file__).resolve().parent
ROOT = H.parent.parent
COLD = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/New Volume/coldstore/render_match_benchmarks'))
PSP = COLD / 'ti64-lpbf-psp'
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_INFOLINK3'))
F = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/cad-to-simulation-F/scripts/physics_exp'))
B = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/bodytwin/scripts/physics_exp'))
G = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/cad-to-simulation-G/scripts/physics_exp'))
OUT = H / 'raw'
OUT.mkdir(exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)
T0 = time.perf_counter()
sources = {}
facts = {}
metrics = {}
controls = {}
legacy = []
attrition = {}
observed = {}
old = [json.loads(l) for l in (ROOT / 'notes/expansion/old_cells_20261002.jsonl').read_text().splitlines()]
prereg = json.loads((H / 'PREREG_INFOLINK3_R1.json').read_text())

def save(p, d):
    p.write_text(json.dumps(d, indent=2, ensure_ascii=False, allow_nan=False, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)))

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def source(p):
    p = Path(p)
    if str(p) not in sources:
        sources[str(p)] = {'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)}
    return str(p)

def state(stage, next_op):
    save(H / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-infolink3', phase=stage, updated_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), latest_gate=stage, next_operation=next_op))

def witness(lid, p, kind, locator, value, unit, quantity, **kw):
    level = kw.pop('resolution_level', 'PER_POINT')
    facts.setdefault(lid, []).append(dict(path=source(p), kind=kind, locator=locator, original_value=value, unit=unit, quantity=quantity, resolution_level=level, **kw))

def corr(x, y):
    return float(spearmanr(x, y).statistic)

def partial(x, y, z):
    X = np.column_stack([np.ones(len(x)), rankdata(z)])
    a = rankdata(x)
    b = rankdata(y)
    return float(np.corrcoef(a - X @ np.linalg.lstsq(X, a, rcond=None)[0], b - X @ np.linalg.lstsq(X, b, rcond=None)[0])[0, 1])

def bootstrap_partial(x, y, z, seed=3):
    rng = np.random.default_rng(seed)
    a = []
    for i in range(1000):
        ix = rng.integers(0, len(x), len(x))
        a.append(partial(np.array(x)[ix], np.array(y)[ix], np.array(z)[ix]))
    return [float(v) for v in np.percentile(a, [2.5, 97.5])]

def load_legacy(p, run=True, rebind=None):
    """Run read code in lane, redirect outputs; no package or original-tree edits."""
    p = Path(p)
    txt = p.read_text()
    source(p)
    maps = {_release_expand('@DENTAL_EXTERNAL_ROOT@/datasets/render_match_benchmarks'): str(COLD), str(COLD / 'mendeley-pla-lattice'): str(COLD / 'lattice-compression-mendeley')}
    maps.update(rebind or {})
    for (a, b) in maps.items():
        txt = txt.replace(a, b)
    if not run:
        txt = txt.split('\ndef main():')[0]
    ns = {'__name__': '__main__' if run else 'legacy_import', '__file__': str(p)}
    log = io.StringIO()
    start = time.perf_counter()
    code = 0
    (H / 'legacy/reports/probes').mkdir(parents=True, exist_ok=True)
    prev = Path.cwd()
    os.chdir(H / 'legacy')
    try:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            exec(compile(txt, str(p), 'exec'), ns)
    except SystemExit as e:
        code = e.code
    except Exception as e:
        code = 1
        log.write(repr(e))
    finally:
        os.chdir(prev)
    path = H / 'legacy' / (p.stem + '.txt')
    path.write_text(log.getvalue())
    legacy.append(dict(path=str(p), sha256=sha(p), run=run, rebindings=maps, exit_code=code, seconds=time.perf_counter() - start, log=str(path)))
    save(H / 'LEGACY_REPLAY.json', legacy)
    if code:
        print('LEGACY FAILURE', p.name, code, flush=True)
    return ns

def fatigue():
    state('fatigue_originals', 'Join specimens, stages and fracture sites')
    p = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/cad-to-simulation-F/data/external/mdpi_data11040081_ti64_pbf_fatigue/data-4162961-supplementary-Table S1.csv'))
    source(p)
    d = pd.read_csv(p)
    d['source_line'] = np.arange(len(d)) + 2
    print('FATIGUE COLUMNS', list(d.columns), flush=True)
    a = d[d.DefectMeasured == True].copy()
    rt = a.RetestFlag == True
    a['StressEff'] = np.where(rt, a.RetestSigmaMax_MPa, a.SigmaMax_MPa)
    a['CyclesEff'] = np.where(rt, a.Retest_Cycles, a.Cycles_InitialTest)
    idcol = 'physical_specimen_key'
    a[idcol] = a.BuildID.astype(str) + ':' + a.SpecimenPosition.astype(str)
    a['test_region_key'] = a[idcol] + ':' + a.TestRegion.astype(str)
    raw = a.replace({np.nan: None}).to_dict('records')
    save(OUT / 'FATIGUE_ROWS.json', raw)
    for lid in ['L01', 'L02', 'L03']:
        observed[lid] = {'original_path': str(p), 'source_rows': 'raw/FATIGUE_ROWS.json', 'source_row_key': 'source_line', 'identity_axes': ['BuildID', 'SpecimenPosition', 'TestRegion', 'InitDefectID'], 'count': len(raw)}
    ar = a[rt].iloc[0]
    for (c, u) in [('Cycles_InitialTest', 'cycle'), ('Retest_Cycles', 'cycle'), ('SigmaMax_MPa', 'MPa'), ('RetestSigmaMax_MPa', 'MPa')]:
        witness('L01', p, 'csv', {'line': int(ar.source_line), 'column': c}, ar[c], u, c)
    for (lid, col, u) in [('L02', 'RootAreaEff_um', 'um'), ('L03', 'DefectPositionClass', 'category')]:
        witness(lid, p, 'csv', {'line': int(a.iloc[0].source_line), 'column': col}, a.iloc[0][col], u, col)
    metrics['L01'] = {'total_rows': len(d), 'characterized_site_rows': len(a), 'unique_specimens': int(a[idcol].nunique()), 'unique_test_regions': int(a.test_region_key.nunique()), 'retest_site_rows': int(rt.sum()), 'retested_unique_specimens': int(a[rt][idcol].nunique()), 'identifier_column': 'BuildID + SpecimenPosition; then TestRegion; then InitDefectID', 'initial_cycles_all_retests': sorted(a[rt].Cycles_InitialTest.unique().tolist()), 'terminal_over_initial_cycle_ratio_median': float(np.median(a[rt].CyclesEff / a[rt].Cycles_InitialTest))}
    controls['L01'] = {'baseline': 'Initial cycles only', 'misassigned_fracture_site_rows': int((a.CyclesEff != a.Cycles_InitialTest).sum()), 'fault_probe': 'Replacing terminal cycle with initial cycle must break source-stage pairing.'}
    partials = {}
    rawrho = {}
    for (mode, sub) in [('all', a), *list(a.groupby('TestMode'))]:
        partials[mode] = {'n_site_rows': len(sub), 'partial_rank_correlation': partial(sub.RootAreaEff_um, np.log10(sub.CyclesEff), sub.StressEff), 'resolution_level': 'POPULATION'}
        rawrho[mode] = corr(sub.RootAreaEff_um, np.log10(sub.CyclesEff))
    unique = a.drop_duplicates(idcol, keep='first')
    metrics['L02'] = {'by_mode': partials, 'raw_correlations': rawrho, 'size_stress_correlation': corr(a.RootAreaEff_um, a.StressEff), 'first_site_per_specimen_partial': partial(unique.RootAreaEff_um, np.log10(unique.CyclesEff), unique.StressEff), 'uncertainty': 'Descriptive associations; repeated sites and retest preconditioning preclude naive independent causal inference.'}
    controls['L02'] = {'baseline': 'No initiating-defect area', 'missing_site_covariates': len(a), 'raw_vs_stress_adjusted': {'raw': rawrho['all'], 'adjusted': partials['all']['partial_rank_correlation']}}
    X = np.column_stack([np.ones(len(a)), rankdata(a.RootAreaEff_um), rankdata(a.StressEff)])
    y = rankdata(np.log10(a.CyclesEff))
    a['resid'] = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    cls = a.groupby('DefectPositionClass').resid.agg(['count', 'median'])
    emb = a[a.DefectPositionClass == 'Embedded'].resid
    surf = a[a.DefectPositionClass.isin(['Near_Surface', 'Emergent'])].resid
    metrics['L03'] = {'classes': cls.to_dict('index'), 'embedded_minus_other_median_residual_rank': float(emb.median() - surf.median()), 'one_sided_mannwhitney_p_descriptive_only': float(mannwhitneyu(emb, surf, alternative='greater').pvalue), 'caveat': 'Rank adjustment not exact matching. Class residual comparison is descriptive; no causal surface factor.'}
    controls['L03'] = {'baseline': 'Ignore position class', 'unrepresented_class_count': len(cls), 'within_same_data_class_residuals': cls['median'].to_dict()}
    attrition['fatigue'] = {'input_rows': len(d), 'characterized_rows': len(a), 'excluded_no_characterization': len(d) - len(a), 'duplicates_as_extra_sites': int(len(a) - a.test_region_key.nunique()), 'repeated_regions_and_sites_within_physical_specimens': int(len(a) - a[idcol].nunique()), 'not_removed_from_raw': True}
    for f in ['wave1094_f_ti64_pbf_fatigue_murakami_defect_size_confound_controlled.py', 'wave1096_f_ti64_pbf_fatigue_surface_vs_embedded_murakami_extension.py']:
        load_legacy(F / f)
    return a

def curves():
    state('plastic_curves', 'Replay constitutive benchmark on original stress-strain zip')
    zpath = PSP / 'stress-strain_curves.zip'
    source(zpath)
    folder = DATA / 'curves'
    folder.mkdir(exist_ok=True)
    with zipfile.ZipFile(zpath) as z:
        for nm in ['flat set 1.xlsx', 'round set 1.xlsx']:
            info = z.getinfo(nm)
            assert info.file_size < 50000000
            dest = folder / nm
            if not dest.exists():
                dest.write_bytes(z.read(nm))
            source(dest)
    p = folder / 'flat set 1.xlsx'
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    sh = wb.sheetnames[0]
    rows = list(wb[sh].iter_rows(values_only=True))
    save(OUT / 'CURVE_HEADER.json', {'zip': str(zpath), 'member': p.name, 'sheet': sh, 'header_rows': rows[:6]})
    for (r, row) in enumerate(rows, 1):
        if r >= 4 and isinstance(row[0], (float, int)) and isinstance(row[1], (float, int)) and (row[0] > 0) and (row[1] > 0):
            witness('L04', p, 'xlsx', {'sheet': sh, 'cell': f'B{r}'}, row[1], 'MPa', 'engineering_stress', zip_origin=str(zpath) + '!' + p.name)
            break
    ns = load_legacy(F / 'wave714_ti64_constitutive_hollomon_render_match_superres.py', rebind={'/tmp/external-source-1000/-home-anton/7972141c-3b6d-42ad-a34b-1c77ca7419fe/scratchpad/ti64_ss': str(folder)})
    r = json.loads((H / 'legacy/reports/probes/wave714_ti64_constitutive_v0.json').read_text())
    metrics['L04'] = r | {'resolution_level': 'POPULATION', 'legacy_hardening_gate_hollomon_rmse_below_15_MPa': bool(r['holl_rmse_median'] < 15), 'caveats': ['E fallback=110000 MPa is constitutive closure, not a measurement.', 'True-stress conversion valid only before necking under stated assumptions.', 'A linear hardening comparator is historical, not the information_link baseline.']}
    controls['L04'] = {'baseline': 'No plastic curve / elastic only', 'measured_curve_count': r['n_curves'], 'historical_hardening_benchmark': {'Hollomon_RMSE_MPa': r['holl_rmse_median'], 'linear_RMSE_MPa': r['lin_rmse_median']}}
    observed['L04'] = [{'source_file': str(folder / n), 'all_sheets': pd.ExcelFile(folder / n).sheet_names, 'columns': ['engineering strain [mm/mm]', 'engineering stress [MPa]']} for n in ['flat set 1.xlsx', 'round set 1.xlsx']]

def ebsd():
    state('ebsd_measurements', 'Read maps; distinguish c-axis components from full grains')
    paths = {name: COLD / 'lpbf-ti64-ebsd/EBSD' / rel for (name, rel) in [('Small_Gaussian', 'Small Gaussian/Gaussian_1.ctf'), ('Tailored_beam', 'Tailored_beam/Tailored_beam.ctf')]}
    entries = []
    for (name, p) in paths.items():
        source(p)
        with p.open() as f:
            lines = []
            for l in f:
                lines.append(l.rstrip('\n'))
                if l.startswith('Phase\tX\tY\t'):
                    break
        hdr = len(lines)
        arr = np.loadtxt(p, skiprows=hdr, usecols=(0, 1, 2, 5, 6))
        alpha = arr[:, 0] == 2
        a = arr[alpha]
        p1 = np.radians(a[:, 3])
        P = np.radians(a[:, 4])
        v = np.column_stack((np.sin(P) * np.sin(p1), -np.sin(P) * np.cos(p1), np.cos(P)))
        T = v.T @ v / len(v)
        ev = np.linalg.eigvalsh(T)
        line = int(np.where(alpha)[0][0]) + hdr + 1
        witness('L05', p, 'ctf', {'line': line, 'column_index': 5}, float(a[0, 3]), 'degree', 'Euler_phi1', map=name)
        ln = next((i for (i, l) in enumerate(lines, 1) if l.startswith('XStep')))
        witness('L06', p, 'text_field', {'line': ln, 'field_index': 1}, float(lines[ln - 1].split()[1]), 'um', 'map_step', map=name)
        entries.append({'map': name, 'path': str(p), 'header': lines, 'rows': len(arr), 'alpha_rows': len(a), 'excluded_non_alpha_rows': int((~alpha).sum()), 'orientation_tensor': T.tolist(), 'eigenvalues': ev.tolist(), 'resolution_level': 'POPULATION', 'tensor_population': 'alpha-phase pixels, area weighting; not independent grains', 'anisotropy_frobenius_from_I_over_3': float(np.linalg.norm(T - np.eye(3) / 3))})
    save(OUT / 'EBSD_MAPS.json', entries)
    observed['L05'] = [{'path': e['path'], 'source_columns': ['Phase', 'X', 'Y', 'Euler1', 'Euler2', 'Euler3'], 'rows': e['rows'], 'original_header': e['header']} for e in entries]
    metrics['L05'] = {'maps': entries, 'uncertainty': 'Two maps, one per condition; no between-build replication. Anisotropic texture does not alone identify anisotropic E.'}
    controls['L05'] = {'baseline': 'Isotropic T=I/3', 'baseline_lambda1': 1 / 3, 'observed_lambda1': [e['eigenvalues'][-1] for e in entries]}
    load_legacy(F / 'wave791_lpbf_ti64_ebsd_texture_qc_reproduction.py')
    load_legacy(F / 'wave696_lpbf_ti64_grainsize_two_sided_censoring.py')
    r = json.loads((H / 'legacy/reports/probes/wave696_lpbf_grainsize_censoring_v0.json').read_text())
    metrics['L06'] = {'maps': r, 'segmentation_status': 'PHENOMENOLOGICAL c-axis-connected components at 10 degrees; rotation about c ignored; these are not validated full-misorientation grains.', 'replacement_measurement': 'Full HCP symmetry-aware misorientation segmentation and a larger field of view at calibrated step.'}
    controls['L06'] = {'baseline': 'No boundary/resolution flags', 'flagged_component_fractions': {k: {q: v[q] for q in ['frac_le3px', 'frac_edge', 'frac_top5_on_edge']} for (k, v) in r.items()}}
    observed['L06'] = {'original_maps': observed['L05'], 'derived_components': 'legacy/reports/probes/wave696_lpbf_grainsize_censoring_v0.json', 'derivation_status': 'PHENOMENOLOGICAL; c-axis components are not full-misorientation grains'}
    attrition['ebsd'] = [{k: e[k] for k in ['map', 'rows', 'alpha_rows', 'excluded_non_alpha_rows']} for e in entries]

def pores_grains():
    state('pore_and_grain_joins', 'Read process, pore and feature rows with units')
    p = PSP / 'Ti-6Al-4V_PSP_feature_table.xlsx'
    source(p)
    ws = openpyxl.load_workbook(p, read_only=True, data_only=True).active
    tab = {}
    for row in ws.iter_rows(min_row=7):
        try:
            n = int(row[1].value)
        except (TypeError, ValueError):
            continue
        tab[n] = {'row': row[1].row, 'P': row[2].value, 'V': row[3].value, 'porosity': row[7].value, 'sph': row[21].value, 'grain_mean': row[27].value, 'AR': row[29].value, 'ductility': row[39].value, 'uts': row[33].value}
    raw = []
    rejected = []
    total = 0
    bad = 0
    for f in sorted((PSP / 'pores').glob('sample*.csv')):
        source(f)
        n = int(re.search('sample\\s*(\\d+)', f.name).group(1))
        rows = list(csv.reader(f.open(encoding='utf-8-sig')))
        header = rows[1]
        if not ('EqDiameter (mm' in header[2] and header[4] == 'Sphericity'):
            rejected.append({'sample': n, 'reason': 'Diameter unit missing and fifth column labelled sorted, not Sphericity', 'uninterpreted_data_rows': len(rows) - 2, 'header': header})
            continue
        parsed = []
        for (line, r) in enumerate(rows[2:], 3):
            total += 1
            try:
                dia = float(r[2])
                sph = float(r[4])
                idx = float(r[0])
                vol = float(r[1])
                area = float(r[3])
            except (ValueError, IndexError):
                bad += 1
                continue
            if not (dia > 0 and 0 < sph <= 1):
                bad += 1
                continue
            parsed.append((line, idx, dia * 1000, sph, vol, area))
        a = np.array(parsed)
        if len(a) < 50 or n not in tab:
            rejected.append({'sample': n, 'reason': '<50 valid pores or no process set'})
            continue
        maxrow = a[a[:, 2].argmax()]
        med = float(np.median(a[:, 3]))
        t = tab[n]
        raw.append(dict(sample=n, path=str(f), n_pores=len(a), header=header, max_pore_line=int(maxrow[0]), max_diameter_um=float(maxrow[2]), max_diameter_original_mm=float(maxrow[2] / 1000), sphericity_of_max=float(maxrow[3]), median_sphericity=med, q999_um=float(np.percentile(a[:, 2], 99.9)), **t))
        if n == 1:
            witness('L07', f, 'csv_index', {'line': int(maxrow[0]), 'column_index': 2}, float(maxrow[2] / 1000), 'mm', 'largest_observed_equivalent_diameter')
            witness('L07', p, 'xlsx', {'sheet': ws.title, 'cell': f"AN{t['row']}"}, t['ductility'], 'percent', 'elongation_at_failure', resolution_level='POPULATION', native_resolution='coupon mean within process set')
            witness('L08', f, 'csv_index', {'line': int(maxrow[0]), 'column_index': 4}, float(maxrow[3]), '1', 'sphericity_of_largest_observed_pore')
            witness('L09', p, 'xlsx', {'sheet': ws.title, 'cell': f"C{t['row']}"}, t['P'], 'W', 'laser_power')
            witness('L09', f, 'csv_index', {'line': int(maxrow[0]), 'column_index': 4}, float(maxrow[3]), '1', 'sphericity_of_largest_observed_pore')
    save(OUT / 'PORE_PROCESS_JOINS.json', raw)
    data = pd.DataFrame(raw)
    valid = data[['max_diameter_um', 'porosity', 'ductility', 'q999_um']].notna().all(axis=1)
    a = data[valid]
    metrics['L07'] = {'n_process_sets': len(a), 'rho_max_diameter_ductility': corr(a.max_diameter_um, a.ductility), 'rho_q999_ductility': corr(a.q999_um, a.ductility), 'rho_mean_porosity_ductility': corr(a.porosity, a.ductility), 'partial_max_given_porosity': partial(a.max_diameter_um, a.ductility, a.porosity), 'bootstrap_process_set_CI95': bootstrap_partial(a.max_diameter_um, a.ductility, a.porosity), 'resolution_level': 'POPULATION', 'scope': 'Process-set averages joined to pore-file sample IDs, not proven same tensile coupon. No fatigue response in this dataset.'}
    controls['L07'] = {'baseline': 'Mean porosity only', 'retained_3D_process_sets': len(a), 'partial_observation_added': metrics['L07']['partial_max_given_porosity']}
    metrics['L08'] = {'n_process_sets': len(data), 'n_valid_pore_rows': sum(data.n_pores), 'fraction_largest_less_spherical_than_median': float(np.mean(data.sphericity_of_max < data.median_sphericity)), 'median_largest_sphericity': float(data.sphericity_of_max.median()), 'median_sample_median_sphericity': float(data.median_sphericity.median()), 'resolution_level': 'POPULATION', 'prohibited_inference': 'No measured Kt or fatigue-initiation assignment; rho squared is not a percent independence claim.'}
    controls['L08'] = {'baseline': 'All pores spherical', 'baseline_sphericity': 1, 'median_absolute_error_at_largest_pore': float(np.median(abs(1 - data.sphericity_of_max)))}
    led = data.P * 1000 / data.V
    metrics['L09'] = {'n_process_sets': len(data), 'rho_LED_median_shape': corr(led, data.median_sphericity), 'rho_LED_largest_shape': corr(led, data.sphericity_of_max), 'median_shape_substitution_abs_error': float(np.median(abs(data.median_sphericity - data.sphericity_of_max))), 'resolution_level': 'POPULATION', 'unknown': 'No calibrated mapping from coupon pore shape to dental lattice strength; non-significance is not proof of process blindness.'}
    controls['L09'] = {'baseline': 'Largest shape equals sample median shape', 'substitution_median_absolute_error': metrics['L09']['median_shape_substitution_abs_error']}
    zpath = PSP / 'prior-beta_grain_morphology.zip'
    source(zpath)
    gr = []
    with zipfile.ZipFile(zpath) as z:
        for nm in z.namelist():
            n = int(re.search('_F(\\d+)', nm).group(1))
            txt = z.read(nm).decode()
            r = list(csv.reader(io.StringIO(txt)))
            assert r[0][2] == 'Equivalent Diameter (mm)'
            vals = np.array([float(a[2]) * 1000 for a in r[1:] if len(a) > 2 and float(a[2]) > 0])
            t = tab[n]
            gr.append(dict(sample=n, zip_member=nm, header=r[0], n_features=len(vals), feature_mean_um=float(vals.mean()), table_mean_um=t['grain_mean'], ratio=float(t['grain_mean'] / vals.mean()), AR=t['AR'], table_row=t['row']))
            if n == 1:
                witness('L10', zpath, 'zip_csv_index', {'member': nm, 'line': 2, 'column_index': 2}, float(r[1][2]), 'mm', 'feature_equivalent_diameter')
                witness('L10', p, 'xlsx', {'sheet': ws.title, 'cell': f"AB{t['row']}"}, t['grain_mean'], 'um', 'PSP_grain_equivalent_diameter_mean', resolution_level='POPULATION', native_resolution='grain population mean within process set')
    save(OUT / 'GRAIN_DEFINITION_PAIRS.json', gr)
    g = pd.DataFrame(gr)
    metrics['L10'] = {'n_process_sets': len(g), 'mean_table_to_feature_ratio': float(g.ratio.mean()), 'sd_ratio_population': float(g.ratio.std(ddof=0)), 'spearman_means': corr(g.feature_mean_um, g.table_mean_um), 'resolution_level': 'POPULATION', 'gap': 'Both raw sources report equivalent diameter; local convention evidence does not establish that one is a 3D measurement. No certified 2D-to-3D conversion.'}
    controls['L10'] = {'baseline': 'Diameters interchangeable, ratio=1', 'mean_ratio': float(g.ratio.mean()), 'relative_scale_discrepancy': float(abs(g.ratio.mean() - 1))}
    attrition['pores'] = {'files': len(data) + len(rejected), 'accepted_files': len(data), 'rejected_files': rejected, 'input_rows': total, 'invalid_or_nonphysical_rows': bad, 'ductility_incomplete_process_sets': int((~valid).sum()), 'missing_pore_files_among_42_process_sets': len(tab) - len(data)}
    observed['L07'] = [{'path': r['path'], 'sample': r['sample'], 'columns': r['header'], 'valid_rows': r['n_pores'], 'matched_PSP_row': r['row']} for r in raw]
    observed['L08'] = observed['L07']
    observed['L09'] = observed['L07']
    observed['L10'] = gr
    load_legacy(F / 'p5_f_ti64_weakest_link_3d_pore_resolves_abstain_seed_proof.py')
    load_legacy(F / 'p5_f_murakami_sqrt_area_incomplete_observable_fatigue_killer_systematically_irregular_sphericity_independent.py')
    load_legacy(F / 'p5_f_pore_sphericity_aggregate_process_controlled_but_fatigue_killer_shape_process_blind.py')
    folder = DATA / 'grain_features'
    folder.mkdir(exist_ok=True)
    with zipfile.ZipFile(zpath) as z:
        for nm in z.namelist():
            p = folder / nm
            if not p.exists():
                p.write_bytes(z.read(nm))
    load_legacy(F / 'p5_f_ti64_grain_2d3d_disagreement_stereology_reconcile_seed.py', rebind={'/tmp/external-source-1000/-home-anton/7972141c-3b6d-42ad-a34b-1c77ca7419fe/scratchpad/grains': str(folder)})

def lattice():
    state('lattice_curves', 'Test lower-infill calibration and first-peak representation')
    p = COLD / 'lattice-compression-mendeley/Compressivedata.xlsx'
    source(p)
    df = pd.read_excel(p, 'Resumen', header=1)
    df['source_row'] = np.arange(len(df)) + 3
    rows = []
    for (_, r) in df.iterrows():
        nm = str(r.get('Nombre', '')).strip()
        m = re.match('([TGH])(\\d+)', nm)
        if not m:
            continue
        try:
            ym = float(r['Young Mod'])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(ym) or ym <= 0:
            continue
        rows.append({'specimen': nm, 'topology': m[1], 'infill_proxy': int(m[2]), 'E_native_unit': ym, 'source_row': int(r.source_row)})
    save(OUT / 'LATTICE_MODULUS_ROWS.json', rows)
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    sh = wb['Resumen']
    headers = list(sh.iter_rows(min_row=2, max_row=2, values_only=True))[0]
    ci = list(headers).index('Young Mod') + 1
    first = rows[0]
    cell = openpyxl.utils.get_column_letter(ci) + str(first['source_row'])
    witness('L11', p, 'xlsx', {'sheet': 'Resumen', 'cell': cell}, first['E_native_unit'], 'UNKNOWN', 'reported_Young_modulus')
    res = []
    for topo in ['T', 'G', 'H']:
        a = [r for r in rows if r['topology'] == topo]
        d = np.array([r['infill_proxy'] for r in a])
        E = np.array([r['E_native_unit'] for r in a])
        hi = d.max()
        tr = d < hi
        n = 1 if topo == 'T' else 2
        C = np.median(E[tr] / d[tr] ** n)
        pred = C * hi ** n
        meas = float(E[~tr].mean())
        res.append({'topology': topo, 'n_lower_infill': int(tr.sum()), 'n_test': int((~tr).sum()), 'n_assumed': n, 'test_infill_proxy': int(hi), 'pred_native_unit': float(pred), 'measured_test_mean_native_unit': meas, 'relative_error': float(abs(pred - meas) / meas), 'passes_legacy_10percent_gate': bool(abs(pred - meas) / meas < 0.1), 'resolution_level': 'POPULATION'})
    metrics['L11'] = {'by_topology': res, 'caveats': ['Infill number from specimen name is a density proxy, not measured as-built relative density.', 'Historical free-fit uses highest-density outcomes when fitting n, so its error is not held-out.', 'No numerical PLA-to-Ti64 parameter transfer.']}
    controls['L11'] = {'baseline': 'Binary exponent n=1/2', 'actual_errors': [r['relative_error'] for r in res]}
    load_legacy(B / 'i1_lattice_engine_predicts_modulus.py')
    q = G / 'g396_dseed_lattice_buckling_load_nonidentifiable_from_elastic_regime_imperfection_knockdown_abstain_overdet_mechanical.py'
    ns = load_legacy(q, run=False)
    data = ns['load_all']()
    rlist = []
    excluded = []
    for (name, (e, s, art)) in data.items():
        a = ns['analyse'](e, s)
        if a:
            rlist.append(a | {'sheet': name, 'negative_stress_artifact': bool(art), 'resolution_level': 'PER_POINT', 'native_resolution': 'specimen stress-strain peak'})
        else:
            excluded.append(name)
    save(OUT / 'LATTICE_FIRST_PEAKS.json', rlist)
    sn = rlist[0]['sheet']
    w = wb[sn]
    rr = list(w.iter_rows(values_only=True))
    hi = list(rr[1])
    col = next((i for (i, c) in enumerate(hi) if 'Stress' in str(c)))
    for (k, r) in enumerate(rr[2:], 3):
        if isinstance(r[col], (float, int)) and r[col] > 0:
            witness('L12', p, 'xlsx', {'sheet': sn, 'cell': openpyxl.utils.get_column_letter(col + 1) + str(k)}, r[col], 'MPa', 'raw_compression_stress')
            break
    rr0 = np.array([r['resid'] for r in rlist])
    clean = np.array([r['resid'] for r in rlist if not r['negative_stress_artifact']])
    metrics['L12'] = {'n_analyzed_curves': len(rlist), 'fraction_negative_projection_residual': float(np.mean(rr0 < 0)), 'median_abs_projection_residual': float(np.median(abs(rr0))), 'median_abs_residual_excluding_negative_stress_artifacts': float(np.median(abs(clean))), 'legacy_A1_ge_80percent_negative': bool(np.mean(rr0 < 0) >= 0.8), 'legacy_A2_median_abs_gt_10percent': bool(np.median(abs(rr0)) > 0.1), 'resolution_level': 'POPULATION', 'gap': 'Projection uses observed peak strain, so it is a retrospective diagnostic. First peak is algorithmically extracted, not independently confirmed buckling; no nonidentifiability theorem.'}
    controls['L12'] = {'baseline': 'Elastic projection at observed first-peak strain', 'median_relative_error': metrics['L12']['median_abs_projection_residual']}
    attrition['lattice'] = {'summary_rows_retained': len(rows), 'curve_sheets_loaded': len(data), 'curves_retained_first_peak': len(rlist), 'excluded_first_peak_gate': excluded, 'negative_stress_artifacts_retained_and_flagged': [r['sheet'] for r in rlist if r['negative_stress_artifact']]}
    observed['L11'] = rows
    observed['L12'] = rlist

def read_fact(w):
    p = Path(w['path'])
    loc = w['locator']
    k = w['kind']
    if k == 'xlsx':
        return openpyxl.load_workbook(p, read_only=True, data_only=True)[loc['sheet']][loc['cell']].value
    if k in ['csv', 'csv_index']:
        with p.open(encoding='utf-8-sig', newline='') as f:
            rr = list(csv.reader(f))
        col = rr[0].index(loc['column']) if k == 'csv' else loc['column_index']
        val = rr[loc['line'] - 1][col]
    elif k == 'zip_csv_index':
        with zipfile.ZipFile(p) as z:
            rr = list(csv.reader(io.StringIO(z.read(loc['member']).decode())))
        val = rr[loc['line'] - 1][loc['column_index']]
    else:
        with p.open() as f:
            for _ in range(loc['line']):
                line = f.readline()
        val = line.split()[loc.get('column_index', loc.get('field_index'))]
    return val if isinstance(w['original_value'], str) else float(val)

def verify_value(w):
    actual = read_fact(w)
    return actual == w['original_value']
CELL_MAP = {'L01': 'wave1094_f', 'L02': 'wave1094_f', 'L03': 'wave1096_f', 'L04': 'wave714_ti64', 'L05': 'wave791_lpbf', 'L06': 'wave696_lpbf', 'L07': 'p5_f_ti64_weakest_link_3d', 'L08': 'p5_f_murakami_sqrt_area', 'L09': 'p5_f_pore_sphericity', 'L10': 'p5_f_ti64_grain_2d3d', 'L11': 'i1_lattice_engine_predicts_modulus', 'L12': 'g396_dseed'}

def finalize():
    state('source_and_fault_validation', 'Validate original witnesses and serialize twelve bounded links')
    results = []
    mut = []
    facit = []
    for c in prereg['candidates']:
        lid = c['id']
        fs = facts[lid]
        paths = [x for x in old if CELL_MAP[lid] in Path(x['evidence'][0]).name]
        assert len(paths) == 1, (lid, paths)
        assert lid in observed, ('missing observation registry', lid)
        oldnode = paths[0]
        cell = Path(oldnode['evidence'][0]).expanduser()
        source(cell)
        for w in fs:
            ok = verify_value(w)
            assert ok, ('witness failed', lid, w, read_fact(w))
            bad = copy.deepcopy(w)
            bad['original_value'] = bad['original_value'] + '__INJECTED_ERROR' if isinstance(bad['original_value'], str) else bad['original_value'] + max(1, abs(bad['original_value']) * 0.1)
            rejected = not verify_value(bad)
            assert rejected
            mut.append({'link': lid, 'kind': 'wrong_original_value', 'locator': w['locator'], 'valid_passed': ok, 'injected_rejected': rejected})
        first = fs[0]
        loc = first['path'] + '#' + json.dumps(first['locator'], ensure_ascii=False, sort_keys=True)
        relation = {'type': 'input_for', 'source': oldnode['id'], 'target': c['consumer'], 'target_kind': 'draft_node', 'role': c['key'], 'claim_type': 'information_link', 'resolution_level': c['resolution_level'], 'native_resolution': 'measurement point or specimen test record; explicit source row axes retained', 'time_scale': c['time_scale'], 'source_locator': loc, 'semantics_status': 'PENDING_INDEPENDENT_REVIEW', 'note': c['changed_operation']}
        row = {'id': 'DENT-PROOF_LANE3-INFOLINK-' + lid, 'type': 'INFORMATION-LINK', 'status': 'OPEN', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'claim_type': 'information_link', 'claim': c['capability'], 'source_node_id': oldnode['id'], 'source_id': oldnode['id'], 'consumer': c['consumer'], 'original_value': first['original_value'], 'unit': first['unit'], 'source_locator': loc, 'resolution_level': c['resolution_level'], 'time_scale': c['time_scale'], 'depends_on': [], 'evidence': [str(cell), *sorted(set((w['path'] for w in fs))), 'results/PROOF_LANE_INFOLINK3/results.json'], 'relations': [relation], 'source_contract': {'witnesses': fs, 'source_sha256': {w['path']: sources[w['path']]['sha256'] for w in fs}, 'complete_observation_registry': 'results/PROOF_LANE_INFOLINK3/OBSERVATION_REGISTRY.json#' + lid, 'original_rows_read_in_place': True, 'claim_boundary': c['changed_operation'], 'consumer_prediction_status': 'NOT_VALIDATED_ON_DENTAL_PARTS'}, 'derived_result': metrics[lid], 'current_practice_without_information': controls[lid], 'uncertainty': {'reported_source_measurement_sigma': 'UNKNOWN unless retained in primary file; no invented precision', 'scope': 'Benchmark process/specimen/map; derived summaries POPULATION, not universal material constants', 'scientific_review': 'PENDING_INDEPENDENT_REVIEW'}, 'new_delta_from_round2': c['key'], 'outcome': 'ESTABLISHED_DRAFT'}

        def valid_contract(x):
            return all((x[k] == v for (k, v) in {'consumer': c['consumer'], 'claim_type': 'information_link', 'resolution_level': c['resolution_level'], 'time_scale': c['time_scale']}.items())) and x['source_contract']['witnesses'] == fs
        assert valid_contract(row)
        for field in ['consumer', 'claim_type', 'resolution_level', 'time_scale']:
            bad = copy.deepcopy(row)
            bad[field] = 'INJECTED_ERROR'
            assert not valid_contract(bad)
            mut.append({'link': lid, 'kind': field, 'injected_rejected': True})
        if lid == 'L11':
            row['outcome'] = 'REJECTED_MISSING_ORIGINAL_UNIT'
            save(H / 'REJECTED_L11.json', row)
        else:
            results.append(row)
        facit.append({'id': row['id'], 'consumer': row['consumer'], 'original_value': row['original_value'], 'unit': row['unit'], 'source_locator': row['source_locator'], 'resolution_level': row['resolution_level'], 'time_scale': row['time_scale'], 'outcome': row['outcome']})
    save(H / 'SOURCE_MANIFEST.json', list(sources.values()))
    save(H / 'OBSERVATION_REGISTRY.json', observed)
    save(H / 'FACIT.json', facts)
    save(H / 'METRICS.json', metrics)
    save(H / 'CONTROLS.json', controls)
    save(H / 'ATTRITION.json', attrition)
    save(H / 'FAULT_INJECTION.json', mut)
    save(H / 'LINK_CONTRACTS.json', results)
    (H / 'information_links.preview.jsonl').write_text(''.join((json.dumps(r, ensure_ascii=False, allow_nan=False, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + '\n' for r in results)))
    with (H / 'FACIT.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=facit[0].keys())
        w.writeheader()
        w.writerows(facit)
    elapsed = time.perf_counter() - T0
    summary = {'lane': 'PROOF_LANE-infolink3', 'claim_type': 'information_link', 'status': 'PENDING_INDEPENDENT_REVIEW', 'established_links': len(results), 'candidate_attrition': {'assessed': 22, 'retained': len(results), 'rejected': 22 - len(results), 'fraction': (22 - len(results)) / 22, 'duplicate_exclusions': 6, 'unsupported_strong_interpretations': 4, 'missing_unit': 1, 'novel_edges_denominator': 16, 'novel_edges_rejected_fraction': 5 / 16, 'old_index_entries': 196, 'other_old_entries_not_adjudicated': True}, 'metrics': metrics, 'controls': controls, 'source_witnesses': sum((len(v) for v in facts.values())), 'fault_injections_rejected': len(mut), 'legacy_runs': legacy, 'cost': {'threads': 1, 'wall_seconds': elapsed, 'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'questions': 0, 'network_requests': 0, 'fit': 'Reported in metric and legacy records', 'discovery_and_reading_time': 'Not instrumented separately; agent wall time not included in demo timing.', 'fallbacks': 'Known unsupported interpretations rejected, source fields retained.', 'data_bytes': sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))}, 'data_artifacts': [{'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in DATA.rglob('*') if p.is_file()]}
    save(H / 'results.json', summary)
    save(H / 'VALIDATION.json', {'source_values_passed': sum((len(v) for v in facts.values())), 'fault_injections_rejected': len(mut), 'all_pass': True, 'prereg_sha256': sha(H / 'PREREG_INFOLINK3_R1.json'), 'preview_sha256': sha(H / 'information_links.preview.jsonl')})
    print('DONE', len(results), 'links', len(mut), 'injections; seconds', round(elapsed, 2), flush=True)
    state('eleven_links_validated_one_candidate_rejected', 'Build bounded reports, install new expansion file and run working build')
if __name__ == '__main__':
    for f in ['PREREG_INFOLINK3_R1.json', 'FROZEN_PREDICTIONS.json']:
        assert sha(H / f) == (H / (f + '.sha256')).read_text().strip(), f + ' drift'
    fatigue()
    curves()
    ebsd()
    pores_grains()
    lattice()
    finalize()
