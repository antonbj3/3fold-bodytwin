# CLOUD-A359DATA: independent empirical NHANES audit

## Input and boundary
`nhanes_public.zip` contains eight unmodified public CDC XPT files. Verify every SHA-256 in `MANIFEST.sha256` after extraction. Their official source URL is `https://wwwn.cdc.gov/nchs/data/nhanes/public/<1999 or 2001>/datafiles/<filename>`. CDC documents the DXA multiple-imputation files at `https://wwwn.cdc.gov/Nchs/Nhanes/Dxa/Dxa.aspx`. No private data or project code is included. `a359.py` is a first-pass analysis from CLOUD-A359 with only the CDC paths and DXA variable names corrected; audit it, do not trust it.

## Strongest baseline and counterexamples
The register claims held-out 2001–02 n≈1485 and leg-lean × height RMSE 28.7 Nm versus mass-only 36.5 Nm (19.1% gain). A local check on these public files and `a359.py` found n=1485 but 28.23 versus 32.49 Nm (13.1% gain) at an **assumed fixed** 0.30 m lever arm. The stronger mass+height+sex+age baseline scored 25.33 Nm. A height-dependent target gave 18.4% gain, creating potential circularity. These numbers are hints to challenge, not validation.

## Work
Independently verify the CDC codebooks, XPT variables, units, knee-strength exclusions, DXA imputation index, and provenance. Run the supplied script, reproduce or refute the local numbers, and inspect leakage, folds, imputation pooling, survey weights, uncertainty and the physical validity of the force-to-torque conversion. Implement a corrected analysis if needed. Compare on exactly the same people: mass-only, mass+height+sex+age, leg lean × height and leg lean × height+sex+age. Compute fixed-lever and height-lever sensitivity and paired bootstrap CIs. Explain any discrepancy with the register without fitting to its target numbers. Report n and missingness. Keep derivation, empirical evidence and assumptions separate.

Deliver `RESULTS.md`, `results.json`, runnable small code, and any tiny tabular summary. Commit locally. Print all result files verbatim in final FILE blocks (≤200 KB total). Do not use private repos, subagents, paid extra usage, Git push, publishing or email. Report any inaccessible codebook as UNKNOWN; never invent provenance.
