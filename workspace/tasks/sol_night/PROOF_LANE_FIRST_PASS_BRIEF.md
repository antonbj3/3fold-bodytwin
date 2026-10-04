# PROOF_LANE_FIRST_PASS — can first-passage friction be identified from second-passage friction?

## The measured state, verified today
The edge `T-E18-passage_number-vs-needle_forward_friction` in `CONSTRAINT_NETS.json` stood UNKNOWN with
the justification that **no digitized series separates forward friction on first from repeated passage in
the same hole**. That premise failed: Fukushima, Saito and Naemura, *Procedia CIRP* 5:265–269 (2013),
doi `10.1016/j.procir.2013.01.052`, insert the needle **twice in the same specimen** and use the second
passage's force as the friction measurement. I looked up the record at the publisher: title, authors, volume
and pages match. Abdullah et al., *Med Phys* 39(6Pt3):3612 (2012), doi `10.1118/1.4734666`, reportedly
report ten cycles with forward force per cycle — check it yourself.

What is missing, however, is a **first-passage** number: the passage-resolved sources measure the plateau.

## The question
The force balance during insertion is usually written `F_total = F_cut + F_fric + F_stiff`. On the second
passage in the same hole `F_cut ≈ 0` is assumed, so the second passage's force is read as friction. Decide
**exactly** under which conditions the subtraction `F_cut = F_1 − F_2` identifies the first passage's
cutting force and friction separately, and under which conditions it does not.

I want the conditions as propositions, not as arguments, and for every condition a **measurable quantity**
that decides whether it holds in a given dataset. In particular:

- what is required for friction on the first passage to be the same as on the second, given that
  the hole already exists on the second passage and that the tissue has had time to relax;
- if they are unequal, which ratio or difference is identifiable from `F_1` and `F_2` alone, and
  which is not;
- which **third** measurement makes the difference identifiable, if it exists.

## Strongest control
The assumption `F_fric,1 = F_fric,2` outright, meaning that passage number does not matter. It is the
current edge's assumption. State what is lost by using it, in force or as a share of
`F_total`, with the sources' own numbers.

## Falsifier
If identifiability requires a quantity that none of the published sources measures, the answer is
that the chain is limited by the measurement protocol and not by the model — say so plainly, with the name of that
quantity, and explain why it cannot be eliminated.

## Warned source of error
A report last night called `µ = 0,4` "konstanten motorn uses". It is **not ours** — grep over
`tasks/` gives no friction coefficient in any BodyTwin model. 0,4 is in an external ocular
trocar FE study and lies 1,69–1,88 σ above the µ it itself cites (0,295 ± 0,056 static,
0,255 ± 0,086 dynamic, Urrea et al., *JMBBM* 56:98–105, 2016). Treat 0,4 as external material.

## Rules
Status `PENDING_INDEPENDENT_REVIEW`, no statement about biological validation, and no named individuals. State every published source with DOI or PMID, volume and pages.
