AGENTER & MOLN

Starta: ~/bin/agent-status
Q closes the view; D shows details; R requests a new measurement. Update every 15 seconds.
The view uses curses: the same terminal area, no ongoing scroll.

Reads local processes, lane_runner active work turns in this session and
the clouds’ verified launcher PIDs, the whole systemd job’s memory/CPU,
model policy, actual shutdown timers and existing lease plan.
The coordinator is reported separately. Queued jobs do not count as active.
Other lane_runner sessions’ built-in agents do not count in the Sol column.
Upon connection error, older values are retained with age and clear marking.
The clouds’ clock difference does not affect status age measurement.

Resource slots use the existing launch_reserved.snapshot directly,
without creating locks or changing files. Model concurrency limits,
cooldown, queue eligibility and atomic final check govern actual start.
An occupied job may wait for a model without using much CPU.
Local CLI CPU/RAM includes the main process, the cloud job all subprocesses.

latest.json contains the latest read values, no command lines/API keys.
lane_runner runtime SQLite is opened read-only; schema may change upon upgrade.
The status view only reads. Automatic refill is in the clouds’ existing
shared launch_reserved guard, with a budget from completed whole jobs.
It selects 1024/1280/1500/2048/2560/3072/4096 MiB with 40% margin above measured peak,
at least three comparable successful jobs before a smaller budget. Unknown:1500.
The first15s reserves the whole limit. Then measured consumption is used
plus at least256MiB growth margin or twice observed growth over60s.
Missing or old measurements reserve the whole limit. This allows
overreservation of individual maximum limits; simultaneous peaks are not guaranteed
by the forecast. Individual and shared hard limits remain.
Shared MemoryHigh is below slice-max with the existing reserve.
CPU>=85%, or CPU>=65% with high CPU-PSI, stops refill.
Memory pressure and new OOM also stop it. CPU-PSI alone does not block available
host cores. Recovery requires60s and two fresh calm measurements.
Running jobs’ memory limits are not shrunk. Only clean inactive file cache
outside child jobs may reduce the parent’s fixed reservation.
ExecStopPost saves RESOURCE_USAGE.json and persistent resource_history.
Memory errors at 2048 MiB or larger still require review before a new attempt.
Old manual positive queue caps are replaced by resource control; cap0 still pauses.
The old dental cap12 applies only when automatic control is off.
Local free jobs start separately in the user’s systemd, with3072MiB limit,
12GiB host reserve and four reserved cores. Older jobs remain running
and get a conservative extra growth reserve. The local limit is controlled automatically.

research-resource-share.service maintains about25% swarm_worker among SIMULTANEOUSLY
active workers (including local), rounded up to whole agents. Targets are written
every20s. Provider backoff, swarm_worker quota, account, cap and expiry apply.
Cat is returned automatically by the existing router when provider backoff expires.
The previous health timer does not change adaptive memory limits or queue caps.

One-time text: agent-status --once
Details in text: agent-status --once --details
Raw data: agent-status --json

RESEARCH DIRECTION
research_value.py gives shared prioritization to the clouds and local The swarm jobs.
It reads the decision in JOB.json, sources and parent; the category name alone is
insufficient. Target distribution80% construction/calibration,15% decisive counter-testing/review,
5% report support when all these types exist. Empty groups are filled with other work.
Families, parents and projects are spread. All original queue records are preserved.
Scoring is an uncalibrated planning heuristic, not an assessment of scientific
truth or a guarantee of actual information value.
New planners get VALUE_BACKLOG.json and the same offensive instructions as
the workers. FOLLOWUPS get a research_value contract; repeated report checks
are pushed back and at most two siblings in the same parent/family/consumer are admitted
per proposal file. Raw proposals and negative results are preserved.
Existing runs finish normally. Existing graphdispatch/feedback,
models, resourcegate and leases apply. VALUE_QUEUE_PREVIEW.json shows the next
prioritization; RESEARCH_VALUE_DEPLOY.json saves the verification.
