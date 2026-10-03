### Q154 — BBB influx, efflux, transcytosis and metabolism

**Question:** How is a defined substance divided between passive diffusion, transporter influx, efflux, transcytosis and endothelial metabolism?

**Resolution target:** Membrane patch/vesicle to capillary segment; milliseconds to hours.

**Physical basis:** Species and charge balance across two membranes and explicit binding/vesicle bookkeeping.

**Mechanism:** Keep luminal membrane, endothelial cell and abluminal membrane as separate mass pools; receptor-directed transcytosis is a testable route and not synonymous with permeability.

**Mathematical starting point:** `dN_E,s/dt=J_lum,in,s−J_lum,efflux,s−J_abl,out,s−J_met,s. N_E [mol], J [mol/s]; each metabolite has its own identity and outflow.`

**Detailed reference to develop:** 3D endothelial cell with separate membrane and vesicle events; compare held substance-specific transport and metabolite profile.

**Laws/parameters to determine:** Transporter expression and kinetics, vesicle traffic, enzyme rate, competition and intracellular binding.

**Output:** Intact substance and metabolite on the brain side [mol/s]; Efflux fraction [1]; Endothelial retention [mol]

**Reuse:** `src/bodytwin/cells/nervous/blood_brain_barrier.py`

**Next work/data:** Identify substrate and cell type before introducing active transport; couple Q153's perfusion and Q155's alternative barrier.

**Comparison:** Passive two-compartment permeability.

**Test that can falsify the proposal:** Active routes cannot be distinguished from a passive model by independent directed fluxes and metabolite measurements.

**Cross-connections to test:** Q009, Q135, Q153, Q155. Source inputs: S13.

