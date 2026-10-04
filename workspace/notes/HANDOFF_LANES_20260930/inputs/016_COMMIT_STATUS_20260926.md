# Commit status checked 2026-09-26

Checked actual worktrees and separate git directories, rather than assuming every workspace is a standalone repository.

| Location | Observed version state |
|---|---|
| Field development worktree `3fold_staging/3fold-field-engine-lane` | Clean; latest `a72439e`, 24 September. The staging main is older. |
| Motion staging and public clone | Clean; direct-download/papers commits from 22 September. No claim that later private experiments are integrated. |
| Graph `projects/graph_workspace/pub` | Tracked files clean; ten untracked research briefs. |
| BodyTwin lane | Separate git directory `research/bodytwin_lane_git`, work tree `projects/3fold-workspaces/bodytwin`; latest `9f777f3` on 25 September. 1,063 tracked working changes at this check, including changed/generated/deleted artifacts; these require scoped reconciliation. |
| Original `projects/bodytwin` | A different, long-lived repository with substantial accumulated working changes. Not a measure of this session's new work; do not sweep into a broad commit. |
| Dental workspace | No ordinary git root found. Its other scoped delivery repositories need a separate ownership-aware inventory before any broad commit. |
| Current autonomous research directory | Was outside git. Finished mathematical/computational packages are now eligible for the explicit local archive maintained by `version_research.py`. Original packages remain preserved. |

The version archive records our finished research packages, code, commands, reports and small data with their source manifests. Large data remain in place with hashes. Ongoing experiments are distinguished from completed snapshots. This closes the missing-commit gap for our research packets; it does not resolve unrelated BodyTwin/dental working changes or constitute independent scientific review.
