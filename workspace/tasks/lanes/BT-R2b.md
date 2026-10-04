# BT-R2b — reclassification of gates according to AU8 + C1b's counter test in pytest
First read tasks/lanes/_PREAMBLE.md (binding). Mandate to drive without asking. DOCUMENTS: results/AU8/AUDIT.md (R2: 79/82 cells generally class B; at least tcell, mapk, hpa have checks against closed form/convergence without
assert = class A according to BT-R2's PREREG), results/BT-R2 (classification, patches), results/BT-C1b (9 counter tests run in script but not in pytest). TASKS: (1) In clone external_mount, new branch fix/gates-exit-2 from fix/gates-exit:
loop through all 79 B cells against PREREG's definition, reclass with file:rad justification, fix new class A cases with regression test (falls before, passes after); format-patch to results/BT-R2b/patches/; commit messages
factual, English, without AI attribution/private names. (2) Copy results/BT-C1b/compose to results/BT-R2b/c1b_tests/ and add the 9 countertests as pytest cases; drive. Output only results/BT-R2b/. RESULTS.md last.
