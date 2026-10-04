# Research overview — last 48 hours

Updated: 2026-09-26T20:50:33.140635+00:00. Window from 2026-09-24T20:50:33.140635+00:00.

**Scope:** the named active workspaces below. This is an inventory of files, deliveries and version state, not a full scientific review or a claim that all other sessions' work has been found. Modification time does not show when an experiment started or finished.

## Status key

- **Ongoing:** no completed delivery receipt.
- **Delivered:** report/package exists; a running process may have finished.
- **Version-saved:** the report's bytes match a commit or a version-saved package copy.
- **Reviewed:** requires explicit independent review; not inferred from a report or commit.
- **Integrated:** requires a concrete implementation/commit in the receiving project; experimental gains are not enough.

## Repos and commits

| Workspace | Latest commit | Modified tracked files | Untracked entries |
|---|---|---:|---:|
| Field integration | a72439e23ed5 | 0 | 0 |
| Motion staging | f679239a69cb | 0 | 0 |
| Graph | 3928720730b3 | 0 | 10 |
| BodyTwin lane | 9f777f34b7c0 | 1063 | 9150 |
| BodyTwin original | 2c42ec52f9b3 | 2560 | 3418 |
| Dental patent lane | 6d14200907dc | 0 | 1 |
| Research archive | 68d70d62bbdc | 0 | 0 |
| Graph/bounded-flat-readout-20260922 | 1d3e9c852fb1 | 0 | 0 |
| Graph/compute-training-plan-20260922 | e1de65d14fce | 0 | 0 |
| Graph/conditional-flat-readout-20260922 | 0364295797be | 0 | 0 |
| Graph/graph-coordinator-continuum-20260922 | e3714b24dc4b | 0 | 0 |
| Graph/graph-coordinator-fast-probe-20260923 | c9dde8d64aeb | 0 | 0 |
| Graph/graph-coordinator-federation-20260922 | 073b6ebbd87f | 0 | 0 |
| Graph/graph-coordinator-oed-batching-20260923 | b6a975836f52 | 0 | 0 |
| Graph/graph-ds-publication-prep-20260924 | 1b8e211fae81 | 0 | 0 |
| Graph/graph-night-integration-20260923 | 268535857109 | 0 | 0 |
| Graph/inference-training-20260922 | 73e76dd83a60 | 0 | 0 |
| Graph/native-inference-bridge-20260922 | 0992f3976187 | 0 | 0 |
| Graph/predictive-actions-20260922 | 1851d01a1cb3 | 0 | 0 |

BodyTwin's lane repo has a separate git directory; the usual empty .git folder in the workspace is not its version history. Working changes must not be merged blindly across owners/experiments. The field integration tree is newer than staging main.

## Inventoried results

| Area | Reports found | Modified within 48 h |
|---|---:|---:|
| Field/motion private experiments | 2245 | 1130 |
| BodyTwin research | 3738 | 1588 |
| Dental research | 165 | 91 |
| Dental graph-bound report snapshots | 0 | 0 |
| Field integration reports | 0 | 0 |
| Graph reports | 0 | 0 |
| Earlier recovery research | 88 | 0 |
| PROOF_LANE_FIVE_20260926 | 10 | 10 |
| PROOF_LANE_WAVE2_20260926 | 7 | 7 |
| PROOF_LANE_WAVE3_20260926 | 8 | 8 |
| PROOF_LANE_WAVE4_20260926 | 0 | 0 |
| LOCAL_WAVE_20260926 | 19 | 19 |

A total of 2762 recently modified report files. A lane can have multiple report types; this is not the number of new or approved research results.

## Traceable packages

The full index is in [INDEX.json](INDEX.json): each report has a path, hash, nearby code, any completion receipt, process outcome and matching against a version-saved report.

### Latest mathematics and hardware tracks

| Package | Delivery marker | Version copy | Report |
|---|---|---|---|
| PROOF_LANE_FIVE_20260926/UNCERTAINTY | ongoing/report exists | yes | [UNCERTAINTY](external_research_path) |
| PROOF_LANE_FIVE_20260926/GRAPH_WILD | done | yes | [GRAPH_WILD](external_research_path) |
| PROOF_LANE_FIVE_20260926/MOTION_WILD | done | yes | [MOTION_WILD](external_research_path) |
| PROOF_LANE_FIVE_20260926/OPERATORS | ongoing/report exists | yes | [OPERATORS](external_research_path) |
| PROOF_LANE_FIVE_20260926/DYNAMICS | ongoing/report exists | yes | [DYNAMICS](external_research_path) |
| PROOF_LANE_FIVE_20260926/CAD_NTOP | done | yes | [CAD_NTOP](external_research_path) |
| PROOF_LANE_FIVE_20260926/OBSERVABILITY | done | yes | [OBSERVABILITY](external_research_path) |
| PROOF_LANE_FIVE_20260926/RSI_LOOP | done | yes | [RSI_LOOP](external_research_path) |
| PROOF_LANE_FIVE_20260926/OBJECTIVES | done | yes | [OBJECTIVES](external_research_path) |
| PROOF_LANE_FIVE_20260926/FIELD_WILD | done | yes | [FIELD_WILD](external_research_path) |
| PROOF_LANE_WAVE2_20260926/SYMMETRY | done | yes | [SYMMETRY](external_research_path) |
| PROOF_LANE_WAVE2_20260926/CONTACT_FRONTIER | done | yes | [CONTACT_FRONTIER](external_research_path) |
| PROOF_LANE_WAVE2_20260926/MATERIAL_DESIGN | done | yes | [MATERIAL_DESIGN](external_research_path) |
| PROOF_LANE_WAVE2_20260926/FOUNDATIONS | done | yes | [FOUNDATIONS](external_research_path) |
| PROOF_LANE_WAVE2_20260926/GRAPH_INTELLIGENCE | done | yes | [GRAPH_INTELLIGENCE](external_research_path) |
| PROOF_LANE_WAVE2_20260926/RECURSIVE_DISCOVERY | done | yes | [RECURSIVE_DISCOVERY](external_research_path) |
| PROOF_LANE_WAVE2_20260926/FIELD_STATE | done | yes | [FIELD_STATE](external_research_path) |
| PROOF_LANE_WAVE3_20260926/DISCOVERY_GEOMETRY | done | yes | [DISCOVERY_GEOMETRY](external_research_path) |
| PROOF_LANE_WAVE3_20260926/MIRROR_DUALITY | done | yes | [MIRROR_DUALITY](external_research_path) |
| PROOF_LANE_WAVE3_20260926/OBJECTIVE_GEOMETRY | ongoing/report exists | not yet | [OBJECTIVE_GEOMETRY](external_research_path) |
| PROOF_LANE_WAVE3_20260926/REPRESENTATION_TRANSPORT | done | yes | [REPRESENTATION_TRANSPORT](external_research_path) |
| PROOF_LANE_WAVE3_20260926/CONTACT_GEOMETRY | done | yes | [CONTACT_GEOMETRY](external_research_path) |
| PROOF_LANE_WAVE3_20260926/WHOLE_SPACE | done | yes | [WHOLE_SPACE](external_research_path) |
| PROOF_LANE_WAVE3_20260926/HIDDEN_BRANCHES | ongoing/report exists | not yet | [HIDDEN_BRANCHES](external_research_path) |
| PROOF_LANE_WAVE3_20260926/MATERIAL_REACHABILITY | done | yes | [MATERIAL_REACHABILITY](external_research_path) |

## Remaining organization work

1. Split the BodyTwin lane's working changes by completed experiment and responsible session; distinguish collector moves/generated files from code changes. No broad commit of the entire workspace has been made.
2. Dental's different scoped repos and original data need continued version mapping; this inventory links the workspace's reports.
3. Connect the largest scientific claims to independent review decisions, and then to the receiving code commit when integration actually happens.
4. Keep raw data on the data disk; large artifacts need hash/manifest and backup rather than everything being put in Git. The current local version archive is not a separate backup.
5. The new wave's deliveries get their own package copies and local commits through the final check. No experiments are automatically classified as integrated.
