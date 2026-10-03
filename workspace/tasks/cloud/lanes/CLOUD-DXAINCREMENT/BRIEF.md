# CLOUD-DXAINCREMENT: does DXA add predictive value beyond demographics?

## Public input
`nhanes_public.zip` contains eight unmodified CDC XPT files, 1999–2002 DEMO/BMX/MSX/DXX. Verify hashes in `MANIFEST.sha256`. Official source: `https://wwwn.cdc.gov/nchs/data/nhanes/public/<1999 or 2001>/datafiles/<filename>`, plus CDC's DXA guide `https://wwwn.cdc.gov/Nchs/Nhanes/Dxa/Dxa.aspx`. `a359.py` is an audit target and setup helper; inspect and correct it independently. No private data or code is included.

## Baseline and counterexamples
Strongest observed baseline from first local check: 1999–2000-fit mass+height+sex+age model yielded 25.33 Nm holdout RMSE in 2001–02, n=1485, under an assumed fixed 0.30 m lever arm. Leg lean × height alone was worse at 28.23 Nm; adding sex+age yielded 24.91 Nm, a small apparent gain. These are unreviewed hints. The register's reported 19.1% DXA gain used a weaker mass-only baseline. Counterexamples: sex residuals, target scaling by height, DXA imputation uncertainty, cohort exclusion and mass/lean collinearity may change the ranking.

## Work
Reproduce an age ≥50, same-person cross-cycle holdout. Compare mass+height+sex+age against the **nested** model adding bilateral leg lean (and, separately, leg lean × height), with a fixed target of measured peak force in N. Only after that convert with fixed 0.30 m for Nm illustration. Use the CDC strength and DXA codebooks to check units, flags and exclusions, or mark unresolved fields. Freeze feature sets and analysis before computing holdout results. Evaluate the 5 DXA imputations separately and pool estimates and uncertainty in a defensible way; report per-imputation variability, paired bootstrap uncertainty, and sex/age subgroups. Include a calibration plot table, explicit missingness audit and code runnable from the bundle. Distinguish prediction from causal/clinical claims.

Deliver RESULTS.md, results.json, small code and final FILE blocks ≤200 KB. Commit locally. No Git push, private repos, subagents, email, publishing or paid extra usage.
