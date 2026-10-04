"""Supplement the frozen R2 observations with independently inspected primary cells.

Do not interpret failure loads, fitted Weibull shapes, or local gaps as one quantity.
"""
from evidence_catalog import *
from decimal import Decimal

def converted_table(*args, canonical_value, conversion, canonical_sd=None, **kwargs):
    args = list(args)
    cell = clean(list(ET.parse(S / args[1]).getroot().find('.//table-wrap[@id="' + args[2] + '"]').findall('.//tr')[args[3]])[args[4]])
    transformed = cell.replace(',', '') if conversion == 'comma_thousands' else re.sub('(?<=\\d) (?=\\d{3}(?:\\D|$))', '', cell)
    assert Decimal(num(transformed)[0]) == Decimal(canonical_value), (args[0], cell, canonical_value)
    args[5] = num(cell)[0]
    table(*args, **kwargs)
    o = OBS[-1]
    o['value'] = str(canonical_value)
    o['reference']['numeric_conversion'] = conversion
    if canonical_sd is not None:
        o['dispersion']['value'] = str(canonical_sd)
for (row, group, v, sd) in [(3, 'control', '1400', '298'), (5, 'monolithic_1.0', '2340', '327'), (6, 'CADon_1.5', '2680', '253'), (7, 'CADon_1.0', '2560', '280')]:
    converted_table('CADON_' + group, 'PMC10742353.xml', 't1', row, 1, v[0], 'N', 'step_stress_fatigue_failure_load', '10 crowns/group; five groups; surviving crowns right censored', 'G10 supports; water; 20Hz,200N/5000cycles then +200N/10000cycles; crack inspection; limit2800N/135000cycles', 'PER_TOOTH', ['D-E-K32', 'D-E-K36', 'D-E-K49'], sd=sd, stat='reported_Kaplan_Meier_step_endpoint_mean', context={'group': group, 'censoring': 'limit2800N'}, canonical_value=v, conversion='space_thousands')
table('CADON_m', 'PMC10742353.xml', 't1', 6, 3, '22.41', '1', 'weibull_m', '10 CAD-on1.5mm crowns; seven surviving crowns at end test', 'Maximum-likelihood fit of step-stress failure load; same wet stepped protocol as CADON_CADon_1.5', 'PER_TOOTH', ['D-E-K32', 'D-E-K36'], 'weibull_m', kind='published_dataset', stat='reported_fitted_Weibull_shape', context={'group': 'CADon_1.5', 'CI95_printed': '11.74-37.97'})
for (row, (group, mean, sd, m, s0)) in enumerate([('C_before', '757', '81', '11.2', '791'), ('Z_before', '891', '115', '9.0', '940'), ('D_before', '835', '101', '9.6', '878'), ('C_after', '1077', '115', '12.4', '1120'), ('Z_after', '1126', '116', '12.8', '1170'), ('D_after', '1322', '212', '7.3', '1408')], 1):
    prefix = 3 if row == 1 else 2 if row == 4 else 1
    protocol = 'ISO6872:2008 biaxial discs; dry room temperature;1mm/min; ' + group + ' polishing/sintering order'
    ctx = {'group': group, 'test': 'biaxial', 'preparation': group.split('_')[1]}
    table('STAW_MEAN_' + group, 'PMC5456702.xml', 'materials-09-00180-t005', row, prefix + 2, mean, 'MPa', 'flexural_strength_mean', 'n40/group; six biaxial groups; other test methods excluded', protocol, 'PER_TOOTH', ['D-E-K32', 'D-E-K36'], sd=sd, sd_col=prefix, context=ctx)
    for (q, col, v, u, var) in [('m', prefix, m, '1', 'weibull_m'), ('sigma0', prefix + 2, s0, 'MPa', 'weibull_sigma0')]:
        table('STAW_' + q + '_' + group, 'PMC5456702.xml', 'materials-09-00180-t006', row, col, v, u, 'weibull_m' if q == 'm' else 'weibull_characteristic_strength', 'n40/group; same biaxial observations as Table5', protocol, 'PER_TOOTH', ['D-E-K32', 'D-E-K36'], var, kind='published_dataset', stat='reported_fitted_Weibull_parameter', context=ctx)
        OBS[-1]['dispersion'] = {'kind': 'reported_CI95', 'printed_interval': clean(list(ET.parse(S / 'PMC5456702.xml').getroot().find('.//table-wrap[@id="materials-09-00180-t006"]').findall('.//tr')[row])[col + 1])}
for (row, (mat, v, sd, m)) in enumerate([('e.max', '372.68', '24.10', '14.93'), ('Empress', '137.37', '10.97', '12.26'), ('Suprinity', '428.48', '12.39', '34.48'), ('Zenostar', '1033.90', '58.78', '16.64'), ('CopraSmile', '621.79', '30.09', '20.32')], 1):
    proto = 'ISO6872:2015 three-point flexure of CAD/CAM bars; material ' + mat
    table('CUREUS_MEAN_' + mat, 'PMC10064933.xml', 'TAB2', row, 1, v, 'MPa', 'flexural_strength_mean', 'n10/material', proto, 'PER_TOOTH', ['D-E-K32', 'D-E-K36'], sd=sd, sd_col=2, context={'material': mat})
    table('CUREUS_m_' + mat, 'PMC10064933.xml', 'TAB2', row, 3, m, '1', 'weibull_m', 'n10/material', proto, 'PER_TOOTH', ['D-E-K32', 'D-E-K36'], 'weibull_m', kind='published_dataset', stat='reported_fitted_Weibull_shape', context={'material': mat})
for (row, (mat, label, v, sd)) in enumerate([('Enamic', 1, '1545.04', '331.74'), ('Enamic', 2, '1707.09', '289.31'), ('Enamic', 3, '1204.96', '130.50'), ('Lava', 1, '1378.25', '232.76'), ('Lava', 2, '2222.74', '320.36'), ('Lava', 3, '1383.84', '208.54')], 1):
    table('ABUTMENT_' + mat + '_' + str(label), 'PMC10286427.xml', 'Tab2', row, 2, v, 'N', 'fracture_mean', '48 crowns; n8/group; lower right first molar implant crown', '10000thermocycles5/55C; titanium abutment;35Ncm;50umcement;45degree cusp loading1mm/min; ' + ('prefabricated abutment,1mm occlusal/3mm proximal' if label == 3 else f'custom abutment,{label}mm reserved crown space'), 'PER_TOOTH', ['D-E-K36', 'D-E-K49'], 'fracture_mean', sd=sd, sd_col=3, context={'material': mat, 'source_group_label_mm': label, 'abutment': 'prefabricated' if label == 3 else 'custom', 'uniform_thickness_sweep': False})
micro = [('prior', 2, 'UL', ['0.3', '1.3', '4.5', '4.5', '0.7', '0.5', '6', '1.5']), ('prior', 3, 'LL', ['3', '1.2', '5.5', '5.3', '1.4', '0.2', '6', '0.8']), ('prior', 4, 'LR', ['2.2', '4', '5.4', '4', '2', '0.8', '8', '6']), ('prior', 5, 'UR', ['0.3', '0.7', '0.2', '1.9', '0.2', '0.3', '0.7', '9']), ('100N_90deg', 10, 'UL', ['11', '9.5', '38', '40.5', '13', '6', '29', '32'])]
for (condition, row, position, values) in micro:
    offset = 2 if row in [2, 10] else 1
    for (j, v) in enumerate(values):
        system = ['MA', 'MS', 'NO', 'BE'][j // 2]
        rep = j % 2 + 1
        printed = num(clean(list(ET.parse(S / 'PMC10976688.xml').getroot().find('.//table-wrap[@id="Tab2"]').findall('.//tr')[row])[offset + j]))[0]
        assert Decimal(printed) == Decimal(v)
        table(f'MICRO_{condition}_{position}_{system}{rep}', 'PMC10976688.xml', 'Tab2', row, offset + j, printed, 'um', 'loaded_interface_microgap', 'Two implant assemblies/system; this cell is one local observation, not group mean', 'BESSY-II phase-contrast radioscopy; four trapezoid edges; manufacturer30Ncm/NO35Ncm; ' + condition, 'PER_POINT', ['D-E-K20', 'D-E-K38'], 'assembly_gap' if condition == 'prior' else None, stat='reported_local_gap', context={'system': system, 'assembly': rep, 'position': position, 'load_condition': condition})
        OBS[-1]['resolution_level'] = 'PER_POINT'
        OBS[-1]['dispersion'] = {'kind': 'source_stated_method_uncertainty', 'value': None, 'scope': '50% for gaps up to2um;2um for larger gaps; author cites Zabler. Not a verified confidence interval.'}
for (row, (group, v, sd)) in enumerate([('SLA', '63.74', '13.61'), ('Nanoblast', '62.83', '9.91')], 1):
    table('RABBIT_BIC_' + group, 'PMC8395172.xml', 'ijms-22-08507-t002', 10, row, v, 'percent', 'bic', '10 rabbits enrolled; Table2 nine retained animals; paired femoral implants', '12weeks rabbit femur histomorphometry; ' + group, 'PER_SURFACE_REGION', ['D-E-K19', 'D-E-K31'], 'bic', sd=sd, context={'week': 12, 'surface': group, 'uncertainty_conflict': 'SLA prose reportsSD10.89,Table2 reports13.61; Table2 preserved, calibration withheld'})
for (row, (F, cycles)) in enumerate([('417.70', '79423'), ('331.60', '280841'), ('265.30', '487551'), ('223.40', '1731994'), ('199.00', '783023'), ('167.50', '3351847'), ('132.60', '5000000')], 1):
    table('BARBED_F_' + str(row), 'PMC10054258.xml', 'materials-16-02228-t002', row, 0, F, 'N', 'fatigue_test_maximum_load', 'Experimental barbed implant assemblies; per-load table observation; replicate n UNKNOWN', 'ISO14801:2017 flexion-compression;9-10Hz; failure/deformation or5millioncycle cutoff', 'PER_TOOTH', ['D-E-K33', 'D-E-K49'], stat='reported_load_level', context={'row': row, 'cycles': cycles, 'runout': row == 7})
    OBS[-1]['resolution_level'] = 'PER_TOOTH'
    converted_table('BARBED_CYCLES_' + str(row), 'PMC10054258.xml', 'materials-16-02228-t002', row, 1, cycles if row != 1 else '79', 'cycles', 'cycles_to_failure_or_runout', 'Same assemblies as corresponding BARBED_F row; replicate n UNKNOWN', 'ISO14801:2017;9-10Hz; finite5millioncycle test; does not establish infinite life', 'PER_TOOTH', ['D-E-K33', 'D-E-K49'], stat='right_censored_runout' if row == 7 else 'reported_cycles_to_failure', context={'row': row, 'Fmax_N': F, 'runout': row == 7}, canonical_value=cycles, conversion='comma_thousands')
    OBS[-1]['resolution_level'] = 'PER_TOOTH'
for (row, diam, loss, v, sd) in [(2, '3.3', '1.5', '382.1', '59.2'), (3, '3.3', '3.0', '347.0', '35.7'), (4, '3.3', '4.5', '315.9', '30.9'), (6, '3.8', '1.5', '531.4', '36.2'), (7, '3.8', '3.0', '514.5', '40.8'), (8, '3.8', '4.5', '477.9', '26.3'), (10, '4.3', '1.5', '710.1', '38.2'), (11, '4.3', '3.0', '697.9', '65.2'), (12, '4.3', '4.5', '662.2', '45.9')]:
    table('DIAMETER_' + diam + '_' + loss, 'PMC10560161.xml', 'Tab1', row, 1, v, 'N', 'fracture_mean', '90Conelog assemblies;n10/diameter-embedding group', '1.2millioncycles50N30deg+10000thermocycles;then static30deg0.5mm/min; failure>20%force drop or2mmdeflection', 'PER_TOOTH', ['D-E-K33', 'D-E-K49'], 'fracture_mean', sd=sd, sd_col=2, context={'diameter_mm': diam, 'simulated_bone_loss_mm': loss, 'endpoint': 'post-aging static failure; not fatigue limit'})
for (label, anchor, v, sd) in [('ZPrime', '8.24 ± 3,21 MPa', '8.24', '3,21'), ('Clearfil', '4.60 ± 2.21 MPa', '4.60', '2.21')]:
    abstract('PRIMER_' + label, 'PMID24918652.json', anchor, v, 'MPa', 'cement_zirconia_shear_bond_strength', '132zirconia specimens,four systems; per-systemn UNKNOWN in abstract', 'Respective primer/resin system; cement cylinder3mm height/diameter; shear bond test; ' + label, 'PER_SURFACE_REGION', ['D-E-K38', 'D-E-K44', 'D-E-K49'], sd=sd, dispersion='reported_plus_minus_unspecified_in_abstract')
    OBS[-1]['dispersion']['value'] = sd.replace(',', '.')
    OBS[-1]['reference']['dispersion_numeric_conversion'] = 'decimal comma:3,21->3.21' if label == 'ZPrime' else None
for (system, v, sd) in [('Standard', '181.6', '60.0'), ('EsthetiCone', '291.3', '41.2'), ('MirusCone', '456.5', '44.0'), ('3iTitanium', '369.7', '32.9'), ('CeraOne', '643.4', '143.1'), ('GoldCylinder', '536.3', '68.6'), ('TiAdapt', '556.9', '145.6')]:
    abstract('TAN_PRELOAD_' + system, 'PMID11432656.json', v + ' +/- ' + sd + ' N', v, 'N', 'axial_preload', 'Seven common hex-top abutment systems; per-systemn UNKNOWN in abstract', 'Calibrated strain-gauged abutment load cells; manufacturer-specific recommended torque; low/high electronic driver speed; ' + system, 'PER_TOOTH', ['D-E-K20', 'D-E-K33', 'D-E-K38'], 'preload', sd=sd, dispersion='reported_plus_minus_unspecified_in_abstract')
    OBS[-1]['context'] = {'system': system, 'torque_Ncm': 'UNKNOWN; varies by manufacturer', 'speed': 'overall mean across settings'}
for (row, (label, v, sd)) in enumerate([('one_water', '-4.27', '0.94'), ('three_water', '-4.66', '2.90'), ('four_water', '-5.03', '1.08'), ('one_dry', '4.43', '3.30'), ('three_dry', '5.13', '3.27'), ('four_dry', '2.87', '2.97')], 1):
    table('CHUA_' + label, 'PMC6717145.xml', 'tbl2', row, 1, v, 'K', 'pulp_temperature_change', 'Extracted human premolars; eight groove cuts/port-condition; same handpieces reused; cuts not eight independent handpieces', '1/3/4-port different handpiece models;1.0±0.2N;30s cutting/30s rest/15s cutting;new bur/cut;37±0.6C pulp baseline;22.4C coolant', 'PER_POINT', ['D-E-K23', 'D-E-PULP-HEAT', 'D-E-K49'], sd=sd, context={'port_condition': label})
if __name__ == '__main__':
    (P / 'VERIFIED_OBSERVATIONS.jsonl').write_text(''.join((json.dumps(x, ensure_ascii=False) + '\n' for x in OBS)))
    print('primary numeric facets', len(OBS), 'sources', len({o['reference']['doi'] or o['reference']['url'] for o in OBS}))
