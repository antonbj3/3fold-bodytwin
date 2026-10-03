# Baseline and reuse before interpreting R2

| Parent | Reuse | Limit and R2 change |
|---|---|---|
| INCISION R1 `incision_r1.py` | Triangular mirrored domain, prescribed-cut energy balance, source conventions | Scalar opening becomes three coupled two-component finite-strain fields. No new Schur/ROM/certifier development. |
| INCISION R1 `nonlinear_exterior_r1.py` | Honest global nonlinear control and fallback cost precedent | New vector constitutive model uses full sparse Newton. No local-update gain is claimed. |
| BINDINGS R1 `pullout` and `MECHANICAL_RESULTS_R1.json` | Positive triangular bridge law, stress-constrained pullout work | Integrate work only up to blade arrival; distinguish 12 nominal fibril work from150 ceiling and from unknown chemical cleavage. |
| Q033 | Elliptical gap/closure and exudate source context | R2 gap is solved, not prescribed elliptic geometry; no exudate rerun or blood relabelling. |
| Q049 | Thin-layer/transfer relevance and separate closure provenance | Its cartilage moduli/permeabilities are not skin data and are not transferred. |
| Private `skin_pulp_mechanics.py` | Read-only constitutive/scope review | Confined palmar indentation is not incision; none of its fitted parameters becomes dermal fracture work. |
| Private collagen/skin/wound modules, MAP inventory and imported graph | Read-only presence/provenance check | Molecular thermal constants and tissue moduli do not identify mode-specific cutting work. Canonical code/graph untouched. |
| Pissarenko Table1, existing BINDINGS table/LOO | Exact skin tear observations already available | Recheck exact frozen table means, not a new independent experiment. No new independent skin force/gap dataset acquired. |

Phenomenological control receives distinct fitted work per observable, identical layer mechanics and every independent acquisition. Its zero training error is not heldout validation. The bridge traction law and positive inverse are conventional same-information controls and match exactly; strongest-control outcome TIE.
