#!/usr/bin/env bash
# Local git for the BodyTwin lane's own work in the shared workspace. The git directory is OUTSIDE the workspace
# (no .git in the shared directory, other sessions' files are unaffected). Only explicitly added paths are tracked.
exec git --git-dir=external_research_path --work-tree= "$@"
