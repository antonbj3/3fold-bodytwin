### Q140 — Drugs through the kidney's different flow paths

**Question:** When are filtration, secretion and reabsorption needed separately to predict plasma time course and excreted amount?

**Resolution target:** Glomerular barrier, segmented nephrons, apical/basolateral transport and local pH/water flow.

**Physical basis:** Species-specific transport, reaction and binding with separate amounts, volumes and recipients for each flow.

**Mechanism:** Follow filtration, tubular secretion and reabsorption with named segments and transporters. Water flow, local ionisation and energy-requiring transport affect the net outflow to urine.

**Mathematical starting point:** `CL_R=f_u,p GFR+CL_sec−CL_reabs, where f_u,p is a fraction 0–1 and GFR, CL_sec, CL_reabs and CL_R are mL/min. The relation applies with steady mass flows and without unspecified metabolic conversion in the tubules.`

**Detailed reference to develop:** Vessel/tissue/cell-resolved transport with molecule-specific properties; PBPK is a comparison level with traceable regional aggregation.

**Laws/parameters to determine:** Substance-specific free/bound partition, enzyme/transporter kinetics, regional perfusion and geometry.

**Outputs:** Total renal clearance, mL/min; Urinary masses per substance, mg; Secretion and reabsorption contributions, mL/min

**Reuse:** `src/bodytwin/cells/renal/renal_filtration.py`

**Next work/data:** Operation: filtered, secreted and reabsorbed mass with held-out molecule and segment.

**Comparison:** Calibrated renal PK with filtration, secretion and reabsorption on the same free plasma and urine data.

**Test that can reject the proposal:** Reject the secretion/reabsorption layer if the contributions are not identifiable or miss held-out urinary mass.

**Cross-links to test:** Q015, Q019, Q020. Source entry points: S03.

