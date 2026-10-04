# CX-READERPASS — rewrite the Overview for a READER, not as a logbook (a side-by-side version; the original is NOT touched)

Anton about `external_research_path`: "'compared with 5.9–8.6 ms for our previous implementation measured alongside it' — this isn't a logbook."

## Task
Write `external_research_path`. Do NOT edit the original or the html/pdf files.
- Every paragraph: what the reader learns (the idea/result) → the evidence (a number) → at most one limitation. Remove process narration: "measured alongside it", "repeated in two runs over three timed positions", "our previous implementation", "was run", "passed comparisons", "remains to be performed", run IDs, and the order of events.
- Keep EVERY number and every source reference [Sxx] that carries the claim. Round only where the precision is meaningless to the reader (e.g. 0.41–0.58 ms → about 0.5 ms). Never change the meaning, never strengthen a claim, and keep limitations that change how the result may be used.
- No trust rhetoric, no adjectives about the system, no analogies. Verb + object + number. Keep the headings as they are (the headings are Anton's).
- Example of the right level: "Implant contact is computed in about 0.5 ms per evaluation and matches an independent double-precision calculation to within 0.0004 N. [S19, S20]"
- Afterwards: write `READERPASS_DIFF.md` next to it, a table of paragraph → the number(s) in the original → the number(s) in v2 → the same meaning (yes/no). Every row must say yes.

The work is text only, with no computation. lane runner has full permissions to write the two new files in FOR_JOHN/; everything else is read-only. No publishing, no emails.

## ADDENDUM (Anton, takes precedence): NO comparisons with our own earlier versions
"You don't need to tell them that it used to take 8 ms; that's irrelevant." Remove ALL comparisons with our own earlier implementations/runs: "than before", "previous", "improved from", "N× faster" against ourselves, and the history of fixes. State only what the system does NOW. Comparisons are allowed only against external references the reader knows (e.g. the reference model, published methods, measured implant force). READERPASS_DIFF.md marks the removed comparisons as "removed (internal history)".
