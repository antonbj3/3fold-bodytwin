# Reproduction of R3

Choose new output destinations. Experiments and tests verify frozen PREREG. Two CPU threads; no network traffic or external HP/abundance consumer in the forward chain.

```bash
PYTHONDONTWRITEBYTECODE=1 python results/LANE_SURGICAL_SYNTHESIS/experiment_r3.py --out results/LANE_SURGICAL_SYNTHESIS/r3/experiment_REPLAY
PYTHONDONTWRITEBYTECODE=1 python results/LANE_SURGICAL_SYNTHESIS/tests_r3.py --out results/LANE_SURGICAL_SYNTHESIS/r3/tests_REPLAY
results/LANE_SURGICAL_SYNTHESIS/run_chain.sh --inventory shared --mode scenario --config results/LANE_SURGICAL_SYNTHESIS/r3/CONFIG_ALL_R3.json --out results/LANE_SURGICAL_SYNTHESIS/r3/chain_REPLAY
results/LANE_SURGICAL_SYNTHESIS/run_chain.sh --resume results/LANE_SURGICAL_SYNTHESIS/r3/experiment_v2/all_R3/checkpoint21.npz --out results/LANE_SURGICAL_SYNTHESIS/r3/restart_REPLAY
```

The CLI now uses shared inventory as default; evidence-mode leaves native outcomes null. `surgical_chain.shared.run` is the new Python entry point. `surgical_chain.chain.run(..., inventory='legacy')` (historical API default) preserves the exact R1/R2 construction for rerunning old diagnostics. R1 package files before this round are copied to r3/legacy_package/; all old rawfiles and negative results remain. The experiment's first path error is in EXPERIMENT_v1.log and its code snapshot; experiment_v2 applies. For reproduction from another working directory, the package must be in PYTHONPATH according to run_chain.sh.
