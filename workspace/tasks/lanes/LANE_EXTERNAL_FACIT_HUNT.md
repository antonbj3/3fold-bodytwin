# LANE_EXTERNAL_FACIT_HUNT

Anton 1/10: the swarm should hunt breakthroughs against ground truth that does not come from us. The graph lane measured: 82 % of 400 swarm jobs compare only against our own setup, 8 % have an external locator. research_value now gives +16…+30 for external ground truth, and planners require the external_referent field — but planners do not know WHICH external ground truths exist per family. Results directory `results/LANE_EXTERNAL_FACIT_HUNT/`.

## Uppgift

For each BodyTwin source family in `tasks/free48/CATALOG.json` (43; start with the most used: Q052, Q168, Q115, Q012, Q146, Q005, Q013, Q090, MITOSTRESS, IMMUNITY, SURG_*, Q019, Q088, Q080, Q100):
1. Read the family’s model and question; determine which observable it predicts (quantity, unit, regime).
2. Find 1–3 **published, publicly accessible ground truths** for that very observable: independent measurements (table/figure in a primary paper), public datasets (GEO, PhysioNet, Zenodo, figshare …), someone else’s open code as control, or closed expressions. Check that the locator is actually accessible and the number is there (quote table/figure and value).
3. State for each: kind, locator (DOI/URL/arXiv/PMID), compared_quantity (exact number/statement with unit), regime/species, and whether the ground truth could show that the family’s current model is wrong (refutes_us).

## Leverans

`FACIT_INDEX.json`: family → list of external_referent objects (schema according to external_research_path) + a short justification. The coordinator puts it in the planners’ inputs so new proposals can point to actual ground truth. Families where no external ground truth exists: write that and why (e.g. pure method question). Never invent numbers or locators; an honest "missing" is the right answer when it is missing.
