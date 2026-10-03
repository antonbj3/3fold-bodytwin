BT-HX-Q140

Status: PASS against the frozen criterion. The primary prediction is an initial renal clearance of 518.642 mL/min, within the prespecified 357–663 mL/min interval.

Reference (source): Graham GG, Punt J, Arora M, Day RO, Doogue MP, Duong JK, et al. (2011), *Clinical Pharmacokinetics* 50:81–98, DOI `10.2165/11534750-000000000-00000`. PubMed page `https://pubmed.ncbi.nlm.nih.gov/21241070/` reports population mean metformin renal clearance 510 ± 130 mL/min. The transporter context is Gong et al. (2012), DOI `10.1097/FPC.0b013e3283559b22`, PMC page `https://pmc.ncbi.nlm.nih.gov/articles/PMC3651676/`.

Derived model output: at the initial plasma concentration (0.002491 mg/L), filtration is 100.000, secretion is 418.642, and reabsorption is 0.000 mL/min, giving 518.642 mL/min. Plasma concentrations at 0, 5, 12, and 24 h are 0.002491, 0.001113, 0.000348, and 0.0000467 mg/L. The 24 h urinary mass is 489.741 mg and urine volume is 2056.320 mL. Gross cumulative filtered, secreted, and reabsorbed amounts are 89.976, 406.016, and 5.341 mg; these are transport flows, not three independent final urine pools. Remaining plasma, lumen, and cell mass is 9.335, 0.910, and 0.015 mg; total mass is 500.000 mg.

Mechanistic interpretation (hypothesis): filtration is necessary for a renal-clearance floor; secretion is necessary here because clearance is well above fu,GFR and the filtration-only null is 100.000 mL/min and 252.566 mg urine at 24 h. Reabsorption is needed only if observed urinary mass is below the filtered load; its human metformin component remains UNKNOWN, not measured zero.

Sensitivity (±50%): secretory OCT2/MATE capacity gives 289.249–755.021 mL/min (rank 1); GFR gives 468.642–568.642 mL/min (rank 2); proximal water reabsorption does not change initial clearance and changes 24 h urine mass by 0.0384% (rank 3). The analytic filtration-only test, pH/dimensional checks, nonnegativity, and mass-conservation test all pass.

Next step: measure iohexol/inulin GFR, unbound metformin, timed urine concentration and volume, segmental pH/flow, and transporter abundance or MATE inhibition; replace effective geometry with nephron-segment geometry. `results.json` contains every numeric output and the complete parameter table.