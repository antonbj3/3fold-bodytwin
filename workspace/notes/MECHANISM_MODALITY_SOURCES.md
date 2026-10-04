# Mechanism, motion, video, splats and camera twin — entry points

Local paths checked 2026-09-22. This is an orientation for a targeted inventory. No models have been rerun and no general assessment of functionality or transferability has been made here.

## Mechanism and its relationship to BodyTwin

- `local_path` is a separate directory with code, data and graph.
- `source_repository/` also contains inherited code from this work.
- Three samples were byte-identical between the directories: `scripts/football_tracking/build_football_tracking_v0.py`, `scripts/video_index/football_scene_verifier.py`, `scripts/physics_exp/p20_camera_twin.py`.
- `data/MECHANISM_ANCHOR_GRAPH.json`, however, was not byte-identical between them. This does not establish where the differences lie or which graph is better updated.

Inventory unique assets, versions, reports, dependencies and negative results before transfer. Avoid importing the same evidence twice. Make a targeted comparison; the entire Mechanism graph must not automatically overwrite BodyTwin's graph.

Older `docs/MECHANISM_HARDENED_CONVENTIONS.md`, `docs/MECHANISM_STARTUP.md` and related documents describe historical workflows. Read them to understand the source projects. Do not run old waves, do not write to source graphs and do not reintroduce old model attributions in commits. The current session assignment and Anton's instructions govern the new work.

## Concrete search entry points

| Area | Entry point |
|---|---|
| Football/motion tracking | `local_path` |
| Video clips and scene checking | `local_path`, and `data/video_index_tables/` |
| Camera twin | `local_path` |
| Splat representation | `local_path`, `splat_quality_cert.py` and `selftest_splat_world.py` |
| Splat experiments | `local_path`, `l_splat_prospective.py`, `l_certified_splat.py` and related reports |
| Camera twin in later CS work | `local_path` and `graf3_kameratvilling_v1_assets/` |
| Football-track review | `local_path` and targeted experiments in `scripts/local/` |
| Older data mapping | `external_research_path` — services, licences and access need current checking before a new download |

Additional directories exist under `local_path`, `local_path` and `local_path`. Their content and relevance remain to be mapped. Directory names do not imply a finished product or an independent implementation.

## Transfer to the new BodyTwin work

Follow the entire chain from raw data to extracted quantities and onwards to a computation result. For motion: investigate identity over time, occlusion, camera motion, geometric scale and time base before results are used in biomechanics. For splats: define which task the representation should support and which separate measures can verify it. Link the camera model's assumptions to computed uncertainty.

Save the inventory in `notes/MECHANISM_TRANSFER_INVENTORY.md`: source/hash, counterpart in BodyTwin, unique change, experimental support, scope limits and recommended use. Make adaptations in your own workspace with traceability to the original and coordinate with other active sessions.
