# CX-OOSCHECK — ADVERSARIAL check of A1778 (frozen hybrid beats N1g by 23 % on unseen non-gait). Try to BREAK the result

Background: A1778 CX-HYBRIDOOS. Hybrid (c), learned on gait: 0.536 vs N1g 0.698 BW on 41 non-gait trials, 4/4 persons, and the time-shift control loses. Suspected weakness: N1g's k is transferred from gait and is known to be too low in squat/chair rise (A1419). The hybrid could therefore win only because N1g gets the level wrong, not because of the model's shape.

## Tasks (PREREG.md + sha256 first; you are the adversary, look for the error)
1. **The fairest N1g:** k learned LOPO on NON-GAIT (the 3 other persons' non-gait frames), per activity group and pooled. Also N1g with k fitted on the held-out person's own non-gait (an upper bound, cheating on purpose, reported only as a ceiling).
2. **Level-only control:** the hybrid's level with the shape removed = w·N1g + (1−w)·(mean of lo_q + c). If that wins as much, the shape was not what helped.
3. **The parameter-free law** (A1398) on the same frames.
4. **Code audit of results/CX-HYBRIDOOS:**
   - leakage: was anything learned on non-gait or on the held-out person's implant force?
   - the mask: does the hybrid lose frames where it would have done badly (1,638 missing, mostly JW)? Score N1g and the hybrid on EXACTLY the same frames, and report what happens if the missing frames get N1g as a fallback.
   - phase definition, BW, units.
5. **Verdict:** does (c) beat the fairest N1g (1) for ≥ 3/4 persons with ≥ 5 % person-median improvement, AND beat the level-only control (2)? Report the numbers.

Runtime ≤ 45 min, 2 threads under bigmem.lock, at most 200 MB of intermediate files. Deliver RESULTS.md starting with `# CX-OOSCHECK` and results.json. Internal data stays local. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
