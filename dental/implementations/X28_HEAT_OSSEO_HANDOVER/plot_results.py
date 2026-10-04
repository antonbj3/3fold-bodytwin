import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from common import ROOT, load

def main():
    rows = load('raw/ud_measurements.json')
    regions = load('raw/ud_injury_regions.json')
    dog = load('raw/canine_group_measurements.json')
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    (fig, axes) = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    ax = axes[0, 0]
    for (meth, col, label) in [('CD', '#2563a5', 'Conventional'), ('UD', '#d47621', 'Ultrasound')]:
        group = [r for r in rows if r['method'] == meth]
        ax.errorbar([r['temp_peak_C'] for r in group], [r['force_N'] for r in group], xerr=[r['temp_SD_C'] for r in group], yerr=[r['force_SD_N'] for r in group], fmt='o', color=col, label=label, capsize=2)
    ax.set(xlabel='Measured peak at distant probe (°C), mean ± SD', ylabel='Mean drilling force (N), ± SD', title='Ovine ex vivo: hotter ultrasound, lower force')
    ax.legend()
    ax.text(0.02, 0.03, 'Probe PER_POINT; force process summary\nHeliyon 2024, Table 2', transform=ax.transAxes, fontsize=8)
    ax = axes[0, 1]
    x = np.arange(3)
    names = ['CD 500 / 60', 'UD 1200 / 30', 'UD 2000 / 60']
    for (bone, dx, col) in [('cortical', -0.18, '#85443e'), ('cancellous', 0.18, '#e3b284')]:
        vals = [r['width_mm'] for r in regions if r['bone_type'] == bone]
        ax.bar(x + dx, vals, width=0.35, color=col, label=bone)
    ax.axhline(0.2, color='gray', ls=':', lw=1)
    ax.set(xticks=x, xticklabels=names, ylabel='Reported injury width (mm)', title='Regional histology: three numeric protocols')
    ax.legend()
    ax.text(0.03, 0.78, 'Approximate observations, no SD\nPER_SURFACE_REGION\n0.2 mm is an observed comparator, not a safe limit', transform=ax.transAxes, fontsize=8, bbox=dict(facecolor='white', alpha=0.93, edgecolor='none'))
    ax.text(0.03, 0.03, 'rpm / feed (mm/min); method + settings differ\nHeliyon 2024, p0105/p0110', transform=ax.transAxes, fontsize=8, bbox=dict(facecolor='white', alpha=0.93, edgecolor='none'))
    ax = axes[1, 0]
    x = np.arange(4)
    ax.bar(x, [r['viability_day3_pct']['mean'] for r in dog], color='#507f84')
    ax.errorbar(x, [r['viability_day3_pct']['mean'] for r in dog], yerr=[r['viability_day3_pct']['SD'] for r in dog], fmt='none', color='#222', capsize=3)
    ax.set(xticks=x, xticklabels=[r['group'] for r in dog], ylim=(0, 110), ylabel='Assay day 3 viable cells (%), mean ± SD', title='Canine harvested chips: biological state after drilling')
    ax.text(0.03, 0.05, 'POPULATION, n=4 reported group samples\n94.13% (500 rpm) vs 56.54% (1000 rpm)\nBMC Oral Health 2023, Fig. 5e / Sec21', transform=ax.transAxes, fontsize=8, bbox=dict(facecolor='white', alpha=0.93, edgecolor='none'))
    ax = axes[1, 1]
    ax.errorbar(x, [r['cell_yield_day14_x1e4_per_g']['mean'] for r in dog], yerr=[r['cell_yield_day14_x1e4_per_g']['SD'] for r in dog], fmt='o-', color='#4255a3', capsize=3)
    ax.set(xticks=x, xticklabels=[r['group'] for r in dog], ylabel='14-day outgrowth cell yield (×10⁴ / g), mean ± SD', title='HANDOVER reaches a measured two-week assay')
    ax.text(0.03, 0.08, '500/1000 rpm yield ratio = 3.483 (group means)\nPOPULATION; culture outcome, not implant healing\nBMC Oral Health 2023, Fig. 5f / Sec21', transform=ax.transAxes, fontsize=8)
    fig.suptitle('Measured process → tissue information: distinct populations and supports', fontsize=14)
    (ROOT / 'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT / 'figures/X28_measured_handover.png', dpi=160)
    fig.savefig(ROOT / 'figures/X28_measured_handover.svg')
    plt.close(fig)
    states = load('THERMAL_HANDOVER_STATES.json')
    state = next((s for s in states if s['id'] == 'irr_off_L10.0'))
    fields = np.load(ROOT / 'raw/retained_thermal_ports.npz')
    dose = fields['irr_off_L10.0_retained_dose_min']
    solid = fields['irr_off_L10.0_solid']
    r = fields['irr_off_L10.0_r_mm']
    z = fields['irr_off_L10.0_z_mm']
    (fig, ax) = plt.subplots(figsize=(8, 4.5), layout='constrained')
    colors = np.ma.array(np.log10(np.maximum(dose, 1e-06)), mask=~solid)
    im = ax.pcolormesh(r - 1.75, z, colors.T, shading='nearest', cmap='magma', vmin=-6, vmax=3)
    if ((dose >= 16) & solid).any():
        ax.contour(r - 1.75, z, np.where(solid, dose, np.nan).T, levels=[16], colors='cyan', linewidths=1)
    ax.set(xlim=(-1.75, 2), ylim=(14, 0), xlabel='Radius minus final cylindrical wall (mm)', ylabel='Depth (mm)', title='Existing K3 dry-drilling dose: retained cells only (simulation)')
    fig.colorbar(im, ax=ax, label='log10(CEM43 min); cyan = 16 min screening contour')
    ax.text(0.02, 0.97, 'PER_POINT, grid 0.1 mm\nThermal source validation FAILED\nDose is not measured viability or healing time', transform=ax.transAxes, va='top', bbox=dict(facecolor='white', alpha=0.9, edgecolor='none'), fontsize=9)
    fig.savefig(ROOT / 'figures/X28_K3_thermal_port.png', dpi=160)
    plt.close(fig)
if __name__ == '__main__':
    main()
