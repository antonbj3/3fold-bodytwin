"""A typed thermal HANDOVER and a domain-restricted measured injury query."""
import json
import time
import numpy as np
from common import ROOT, DENTAL, sha, load, dump, state, verify_freeze

def region_query(regions, protocol_id, bone_type, population, width_limit_mm=0.2):
    if population != 'ovine lumbar vertebra in vitro':
        return dict(status='UNKNOWN', reason='population outside observed source domain')
    selected = [r for r in regions if r['protocol_id'] == protocol_id and r['bone_type'] == bone_type]
    if len(selected) != 1:
        return dict(status='UNKNOWN', reason='no exact protocol/bone-region measurement')
    r = selected[0]
    return dict(status='NOMINAL_ONLY', nominal_meets_observed_benchmark=r['width_mm'] <= width_limit_mm + 1e-12, width_mm=r['width_mm'], width_limit_mm=width_limit_mm, robust_meets_benchmark='UNKNOWN', uncertainty=r['uncertainty'], locator=r['locator'], resolution=r['resolution'], clinical_safety='UNKNOWN', note='0.2 mm is an observed least-injury comparator, not a safe biological limit')

def retained_port(cem, solid):
    if cem.shape != solid.shape or not np.isfinite(cem).all() or (cem < 0).any():
        raise ValueError('Invalid dose field')
    out = np.where(solid, cem, 0.0)
    if out[~solid].any():
        raise ValueError('Removed tissue cannot be an initial healing state')
    return out

def validate_port(cem, solid):
    return cem.shape == solid.shape and bool(np.isfinite(cem).all()) and (not bool((cem < 0).any())) and (not bool(cem[~solid].any()))

def main():
    start = time.perf_counter()
    verify_freeze('PREREG_R2.json')
    verify_freeze('FROZEN_PREDICTIONS_R2.json')
    path = DENTAL / 'results/K3_drilling/thermal_fields.npz'
    source = np.load(path)
    regions = load('raw/ud_injury_regions.json')
    summaries = []
    fields = {}
    volume_errors = []
    roundtrip = []
    for key in source.files:
        if not key.endswith('_cem'):
            continue
        ident = key[:-4]
        cem = source[key].astype(float)
        solid = source[ident + '_solid']
        r = source[ident + '_r']
        z = source[ident + '_z']
        h = float(np.diff(r).mean())
        dose = retained_port(cem, solid)
        flags = (dose >= 16) & solid
        cells = 2 * np.pi * r[:, None] * h * h * np.ones((1, len(z)))
        volume = float(cells[flags].sum())
        (ii, jj) = np.nonzero(flags)
        independent = float(sum((2 * np.pi * r[int(i)] * h * h for (i, j) in zip(ii, jj))))
        volume_errors.append(abs(volume - independent))
        length = float(ident.split('_L')[1])
        Dfinal = 3.5
        wall = Dfinal / 2
        cone_length = wall / np.tan(np.radians(59.0))
        axial = (z > 0) & (z < length + 1 - cone_length)
        radial = r >= wall
        damage = flags & radial[:, None] & axial[None, :]
        widths = np.zeros(len(z))
        for j in np.flatnonzero(axial):
            hit = np.flatnonzero(damage[:, j])
            widths[j] = max(0.0, float(r[hit[-1]] + h / 2 - wall)) if len(hit) else 0.0
        fields[ident + '_retained_dose_min'] = dose.astype(np.float32)
        fields[ident + '_solid'] = solid
        fields[ident + '_r_mm'] = r
        fields[ident + '_z_mm'] = z
        summaries.append(dict(id=ident, source=str(path), source_sha256=sha(path), thermal_dose_resolution='PER_POINT', grid_spacing_mm=h, archival_precision='originalK3float32,CEMcappedat1e30min;screening16minunchanged,cappedmagnitudesnotexact', timescale='HANDOVER', fast_time='seconds during drilling and cooling', slow_time='days/weeks: consumer unresolved', dose_threshold_CEM43_min=16, dose_threshold_status='PHENOMENOLOGICAL: mapped single rabbit anchor, not calibrated transient injury law', retained_voxel_count=int(solid.sum()), removed_voxel_count=int((~solid).sum()), threshold_voxel_count=int(flags.sum()), threshold_volume_mm3=volume, lateral_width_max_mm=float(widths.max()), width_quantization_mm=h, cortical_lateral_width_max_mm=float(widths[z < 1.5].max()), depth_profile=dict(z_mm=z[axial].tolist(), width_mm=widths[axial].tolist(), resolution='PER_SURFACE_REGION'), biological_initial_condition=dict(viability='UNKNOWN', microcracks='UNKNOWN', osteogenic_competence='UNKNOWN'), predicted_healing_delay_weeks='UNKNOWN', predicted_ISQ='UNKNOWN', calibration_status='FAIL: K3 only3/16 thermal anchors withinfactor2;31% bandcoverage', safe_drilling_parameters='UNKNOWN', bone_type='ASSUMED human cortical1.5mm+BVT V0.31 model', injury_data_join='REJECTED: independent ovine histology belongs to different protocol/population'))
    outpath = ROOT / 'raw/retained_thermal_ports.npz'
    np.savez_compressed(outpath, **fields)
    reloaded = np.load(outpath)
    for (k, v) in fields.items():
        roundtrip.append(np.array_equal(v, reloaded[k]))
    fake = fields['central_L10.0_retained_dose_min'].copy()
    solid = fields['central_L10.0_solid']
    fake[np.argwhere(~solid)[0][0], np.argwhere(~solid)[0][1]] = 100
    removed_rejected = not validate_port(fake, solid)
    domain_rejected = region_query(regions, 'UD_2000_60', 'cortical', 'human mandible')['status'] == 'UNKNOWN'
    missing_rejected = region_query(regions, 'CD_2000_60', 'cortical', 'ovine lumbar vertebra in vitro')['status'] == 'UNKNOWN'
    fault_regions = [dict(r) for r in regions]
    for row in fault_regions:
        if row['protocol_id'] == 'UD_2000_60':
            row['width_mm'] *= 1000
    radius_rejected = not region_query(fault_regions, 'UD_2000_60', 'cortical', 'ovine lumbar vertebra in vitro')['nominal_meets_observed_benchmark']

    def volume_control(reported, reference):
        return abs(reported - reference) <= 1e-09
    volume_rejected = not volume_control(summaries[0]['threshold_volume_mm3'] + 1.0, summaries[0]['threshold_volume_mm3'])
    queries = []
    for ident in ['CD_500_60', 'UD_1200_30', 'UD_2000_60']:
        for bone in ['cortical', 'cancellous']:
            queries.append(dict(protocol_id=ident, bone_type=bone, **region_query(regions, ident, bone, 'ovine lumbar vertebra in vitro')))
    checks = dict(removed_tissue_rejected=removed_rejected, wrong_species_rejected=domain_rejected, missing_numeric_protocol_rejected=missing_rejected, radius_times1000_rejected=radius_rejected, altered_volume_rejected=volume_rejected, roundtrip=all(roundtrip), volume_parity=max(volume_errors) <= 1e-09)
    assert all(checks.values()), checks
    dump('THERMAL_HANDOVER_STATES.json', summaries)
    dump('REGION_PARAMETER_QUERY.json', queries)
    out = dict(round='R2', claim_type='information_link', outcome='Regional observed benchmark query and thermal HANDOVER implemented; human safety and healing prediction remain UNKNOWN', external_referent=load('PREREG_R2.json')['external_referent'], gates={'G1_roundtrip': checks['roundtrip'], 'G2_nominal_benchmark': 'UD2000/60 meets0.2mm for both tissue types; robust UNKNOWN', 'G3_domain_guard': domain_rejected and missing_rejected, 'G4_faults': all(checks.values()), 'G5_volume_parity': checks['volume_parity']}, checks=checks, max_volume_parity_error_mm3=max(volume_errors), resolution='PER_POINT dose; PER_SURFACE_REGION observed widths', timescale='HANDOVER', measured_protocol_answer=dict(protocol_id='UD_2000_60', rpm=2000, feed_mm_min=60, diameter_mm=4, ultrasound_frequency_Hz=20000, ultrasound_amplitude_um=20, cooling='NOT_REPORTED', sequence='single4mm drill', measured_mean_force_N=4.05, force_SD_N=0.3, reported_injury_width_mm=dict(cortical=0.2, cancellous=0.2), comparator='published least-injury observation0.2mm;not a safe limit', uncertainty='width SD missing'), thermal_states_count=len(summaries), attrition={'simulated_thermal_fields': {'candidates': 28, 'kept': 28, 'rejected': 0}, 'empirical_joins_to_K3': {'candidates': 6, 'kept': 0, 'rejected': 6, 'fraction': 1.0, 'reason': 'different source population/protocol and no co-located pairing'}, 'robust_observed_injury_constraints': {'candidates': 6, 'kept': 0, 'rejected': 6, 'fraction': 1.0, 'reason': 'no injury width uncertainty'}, 'nominal_least_injury_regions': {'candidates': 6, 'kept': 2, 'rejected': 4, 'fraction': 2 / 3, 'reason': 'observed width exceeds0.2mm comparator'}}, full_cost=dict(fit_s=0, validation_compute_s=time.perf_counter() - start, output_bytes=outpath.stat().st_size, preparation='localK3 fields reused without solving/copying source tree', queries=28 + 6, fallback='OOD/uncertainty/clinical healing -> UNKNOWN'), artifacts={'retained_dose_fields': {'path': str(outpath), 'sha256': sha(outpath), 'bytes': outpath.stat().st_size}})
    dump('results_R2.json', out)
    text = '# R2 handoff\n28 K3 fields now export only surviving tissue dose with coordinates, spatial resolution and failed empirical source-calibration status.\n16 CEM43 min is a screening closure, not a calibrated viability fraction. Lateral widths exclude the conical tip/deeper axis.\nIndependent volume and serialization checks pass; faults in removed-tissue dose, domain, width units and volume fail.\nUD2000rpm/60mmmin (4mm bit,20kHz,20um) has published approximate0.2mm injury in both ovine tissue types.\nThis meets the published least-injury observation nominally, not a biological safety limit; width uncertainty and cooling are missing.\nAll6 empirical joins to K3 rejected: different species/protocol. All6 robust width constraints remain UNKNOWN.\nNext construction: measured viable-cell fraction AND osteogenic competence after drilling. Local canine bone-chip study\nPMC10631188 has region/group-specific10s temperature and3day/14day outcomes; these endpoints can supply an honest HANDOVER\nwithout inventing ISQ or human healing delay. Preserve region/specimen/time identity.\n'
    (ROOT / 'HANDOFF_R2.md').write_text(text)
    (ROOT / 'HANDOFF.md').write_text(text)
    state('R2_DECIDED', out['gates'], 'Freeze R3: measured viability and competence HANDOVER in canine bone chips')
    print(out['outcome'])
if __name__ == '__main__':
    main()
