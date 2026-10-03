# Reproduktion

Run [surgical_chain/README.md](surgical_chain/README.md). All test/run destinations must be new. The reference is `r1/nominal_complete_v4`; earlier `scenario_nominal`, `scenario_seal_memory_v2` and `scenario_final_v3` are preserved. An initial default run result ended up under an extra relative path inside the lane; `run_chain.sh` now preserves caller-CWD and adds the package to PYTHONPATH. No external destination was changed.

`prepare_r1.py`, `build_port_table_r1.py`, `ingest_response_r4.py` and `validate_constructed_chain_r1.py` are single-use freezers with exclusive outputs, not idempotent replay commands. Use the package's CLI for a new run. Graphcli in the lane uses frozen canonical receipt code with workspace root directed to `graph_review_workspace/`; the packet is retrieved read-only from the workspace's `graph working packet`.

The matplotlib figure can be reproduced with `plot_r1.py --out <new lane-local basename>`; it reads frozen raw fields without model fitting. Laboratory acquisition has not been performed. No skill/plugin, agent, queue or cloud is used.
