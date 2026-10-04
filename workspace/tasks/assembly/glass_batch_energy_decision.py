"""Does the carbonate endotherm change the glass furnace budget enough to change a decision?

WHY THIS NODE. The project graph's own priority tool ranks MAKE-GLASS-CHAIN highest among open
actionable nodes (score 8): Swedish quartz sand to melt to form to polish to a measurable lens,
judged by the raytrace cell. Its measured v0 chain (reports/probes/glass_chain_v0.json) carries an
isokom furnace target of 1322.0259 C, a Pt-crucible margin of +446.97 C, and an energy budget of
1.3645 MJ per kg -- and that budget declares its own open branch in writing: "carbonate-decomposition
endotherm (Na2CO3/CaCO3 -> oxide + CO2) and furnace radiative losses".

That branch is not a missing measurement. It is standard thermochemistry on a composition the chain
already states, so it can be closed here. The decision is whether closing it changes anything: if
the carbonate term is a few percent, the v0 budget stands and the open branch was bookkeeping. If it
is comparable to the sensible heat, every energy, cost and furnace-sizing statement downstream of
this node moves, and the chain's first real number was wrong by a factor.

THE CONTROL is the sensible-heat-only budget already in the chain, 1.3645 MJ/kg. It is the number a
furnace designer would use from the v0 report. It is scored, not dismissed.

THE FALSIFIER, stated before the run: if the carbonate endotherm is below 10 percent of the sensible
heat, the open branch does not matter for any decision at this stage and should be closed as
negligible rather than left open.

WHAT IS NOT CLAIMED. Radiative and flue losses are the other half of the declared branch and are a
furnace-geometry question, not a thermochemical one; they are left open and named. Batch-to-glass
yield (CO2 leaves as gas, so a kg of batch is not a kg of glass) is computed explicitly rather than
assumed. PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import os

# Target composition from the chain's own input_spec, mol percent, Fluegel (2007) Table 1
# measured soda-lime-silica container glass.
COMPOSITION_MOL_PCT = {'SiO2': 74.41, 'Al2O3': 0.75, 'Na2O': 12.9, 'K2O': 0.19,
                       'MgO': 0.3, 'CaO': 11.27, 'Fe2O3': 0.16, 'TiO2': 0.01, 'SO3': 0.01}
MOLAR_MASS = {'SiO2': 60.083, 'Al2O3': 101.961, 'Na2O': 61.979, 'K2O': 94.196,
              'MgO': 40.304, 'CaO': 56.077, 'Fe2O3': 159.688, 'TiO2': 79.866, 'SO3': 80.063}
# Standard decomposition enthalpies, carbonate -> oxide + CO2, kJ per mol of oxide formed.
# Values are the textbook standard-state figures; each is tagged with what it rests on.
# Reconciled against the staged probe (reports/probes/glass_batch_decomp_v0.json) with the graph
# lane on 2026-10-04, and the four terms do NOT share an evidence status. The probe is a
# TWO-carbonate batch: its reaction list holds Na2CO3 and CaCO3 only, and K2CO3 and MgCO3 occur zero
# times in it. The gap between its 1.018232 and this file's first figure of 1.037191 decomposed to
# zero: K2O +0.012396, MgO +0.005872, Na2O enthalpy choice +0.003683, CaO enthalpy choice -0.000935,
# and the mass basis -0.002084, summing to the observed 0.018932 exactly.
#
# What remains is a BATCH question, not a thermochemistry one, and the two extra terms stand on
# different ground. The chain's own raw-material line reads "quartz sand + soda ash (Na2CO3->Na2O) +
# limestone/DOLOMITE (CaCO3->CaO) + minor Al2O3/K2O/MgO/Fe2O3/TiO2/SO3 fining agents".
#
#   MgO  -- dolomite is named as a raw material and dolomite is CaMg(CO3)2, so its magnesium enters
#           as a carbonate and its endotherm belongs. The probe names the mineral and omits its
#           magnesium, so the probe's basis is incomplete as written. READING OF THE SPEC.
#   K2O  -- potassium appears only inside the fining-agent list and no raw-material line names
#           potash. This term therefore assumes a potassium carbonate source. ASSUMPTION.
#
# So the defensible range on the spec as written is 1.0183 to 1.0242 MJ/kg, and anything above that
# rests on the potash assumption. Both are reported, labelled, and never summed silently.
DECOMPOSITION_KJ_PER_MOL = {
    'CaO': (178.3, 'CaCO3 -> CaO + CO2, standard formation enthalpies; probe uses 178.80, JANAF 178.8'),
    'Na2O': (321.0, 'Na2CO3 -> Na2O + CO2, standard formation enthalpies; probe uses 319.28'),
    'MgO': (117.9, 'MgCO3 -> MgO + CO2 via the DOLOMITE named in the raw materials; reading of the spec'),
    'K2O': (393.0, 'K2CO3 -> K2O + CO2; ASSUMES a potash source that no raw-material line declares'),
}
EVIDENCE_STATUS = {'Na2O': 'declared', 'CaO': 'declared', 'MgO': 'reading_of_the_spec',
                   'K2O': 'assumption_undeclared_source'}
PROBE_TWO_CARBONATE_MJ_PER_KG = 1.018232
PROBE_PLUS_DOLOMITE_MJ_PER_KG = 1.0242
PROBE_PROCESS_CO2_G_PER_KG = 176.582
CO2_MOLAR_MASS = 44.009
SENSIBLE_HEAT_MJ_PER_KG = 1.3644712286060199     # the chain's own v0 number
FURNACE_T_C = 1322.0258827053422                 # the chain's own isokom target


def main() -> None:
    # One mole of glass "formula units" on the stated mol percentages.
    mass_per_mole_glass = sum(COMPOSITION_MOL_PCT[k] / 100.0 * MOLAR_MASS[k]
                              for k in COMPOSITION_MOL_PCT)
    moles_glass_per_kg = 1000.0 / mass_per_mole_glass

    rows, total_kJ, co2_g = [], 0.0, 0.0
    for oxide, (dh, basis) in DECOMPOSITION_KJ_PER_MOL.items():
        frac = COMPOSITION_MOL_PCT.get(oxide, 0.0) / 100.0
        n = frac * moles_glass_per_kg          # mol of that oxide per kg of GLASS
        kJ = n * dh
        total_kJ += kJ
        co2_g += n * CO2_MOLAR_MASS
        rows.append({'oxide': oxide, 'mol_per_kg_glass': n, 'dH_kJ_per_mol': dh,
                     'MJ_per_kg_glass': kJ / 1000.0, 'basis': basis})

    carbonate_MJ = total_kJ / 1000.0
    share = carbonate_MJ / SENSIBLE_HEAT_MJ_PER_KG
    # A kg of glass needs more than a kg of batch, because the CO2 leaves.
    batch_per_glass = 1.0 + co2_g / 1000.0

    print(f'mass per mole of stated glass units : {mass_per_mole_glass:.4f} g')
    print(f'moles of glass units per kg         : {moles_glass_per_kg:.4f}')
    print()
    for r in sorted(rows, key=lambda x: -x['MJ_per_kg_glass']):
        print(f"  {r['oxide']:<5} {r['mol_per_kg_glass']:7.4f} mol/kg x "
              f"{r['dH_kJ_per_mol']:6.1f} kJ/mol = {r['MJ_per_kg_glass']:7.4f} MJ/kg")
    print()
    print(f'carbonate decomposition endotherm   : {carbonate_MJ:.4f} MJ/kg glass')
    print(f'the chain v0 sensible heat          : {SENSIBLE_HEAT_MJ_PER_KG:.4f} MJ/kg')
    print(f'carbonate as a share of sensible    : {share * 100:.1f} %')
    print(f'corrected thermochemical total      : '
          f'{carbonate_MJ + SENSIBLE_HEAT_MJ_PER_KG:.4f} MJ/kg, a factor '
          f'{(carbonate_MJ + SENSIBLE_HEAT_MJ_PER_KG) / SENSIBLE_HEAT_MJ_PER_KG:.3f} on the control')
    print(f'CO2 released                        : {co2_g:.1f} g per kg of glass, so '
          f'{batch_per_glass:.4f} kg of batch per kg of glass')
    declared = sum(r['MJ_per_kg_glass'] for r in rows
                   if EVIDENCE_STATUS[r['oxide']] == 'declared')
    with_dolomite = declared + sum(r['MJ_per_kg_glass'] for r in rows
                                   if EVIDENCE_STATUS[r['oxide']] == 'reading_of_the_spec')
    print('\nTHREE READINGS, never summed silently:')
    print(f'  two declared carbonates only   : {declared:.4f} MJ/kg  '
          f'(staged probe: {PROBE_TWO_CARBONATE_MJ_PER_KG:.6f})')
    print(f'  plus dolomite magnesium        : {with_dolomite:.4f} MJ/kg  '
          f'(graph lane: {PROBE_PLUS_DOLOMITE_MJ_PER_KG})')
    print(f'  plus assumed potash            : {carbonate_MJ:.4f} MJ/kg  '
          'rests on a source no raw-material line declares')
    print(f'  defensible range on the spec   : '
          f'{PROBE_TWO_CARBONATE_MJ_PER_KG:.4f} to {PROBE_PLUS_DOLOMITE_MJ_PER_KG:.4f} MJ/kg')
    share_low = PROBE_TWO_CARBONATE_MJ_PER_KG / SENSIBLE_HEAT_MJ_PER_KG
    share_high = PROBE_PLUS_DOLOMITE_MJ_PER_KG / SENSIBLE_HEAT_MJ_PER_KG
    print(f'  share of sensible heat, range  : {share_low * 100:.1f} to {share_high * 100:.1f} %')
    verdict = ('CHANGES THE DECISION: the open branch is not bookkeeping, and it holds across the '
               'whole defensible range, not only on the potash assumption'
               if share_low > 0.10 else
               'negligible: close the open branch rather than leaving it open')
    print(f'\nverdict: {verdict}')
    print('still open and named: furnace radiative and flue losses, which are geometry, not')
    print('thermochemistry, and the heat of mixing of the oxides into the melt.')

    out = {
        'node': 'MAKE-GLASS-CHAIN',
        'control_sensible_heat_MJ_per_kg': SENSIBLE_HEAT_MJ_PER_KG,
        'carbonate_endotherm_MJ_per_kg': carbonate_MJ,
        'share_of_control': share,
        'corrected_total_MJ_per_kg': carbonate_MJ + SENSIBLE_HEAT_MJ_PER_KG,
        'factor_on_control': (carbonate_MJ + SENSIBLE_HEAT_MJ_PER_KG) / SENSIBLE_HEAT_MJ_PER_KG,
        'process_co2_g_per_kg_glass': co2_g,
        'batch_kg_per_kg_glass': batch_per_glass,
        'furnace_isokom_target_C': FURNACE_T_C,
        'per_oxide': rows,
        'evidence_status_per_oxide': EVIDENCE_STATUS,
        'defensible_range_MJ_per_kg': [PROBE_TWO_CARBONATE_MJ_PER_KG,
                                       PROBE_PLUS_DOLOMITE_MJ_PER_KG],
        'process_co2_g_per_kg_glass_this_file': None,
        'process_co2_g_per_kg_glass_probe': PROBE_PROCESS_CO2_G_PER_KG,
        'verdict': verdict,
        'still_open': ['furnace radiative and flue losses', 'heat of mixing into the melt'],
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    d = 'results/ASSEMBLY_GLASS_BATCH_ENERGY'
    os.makedirs(d, exist_ok=True)
    with open(d + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('wrote', d + '/decision.json')


if __name__ == '__main__':
    main()
