#!/bin/bash
# The swarm keepalive without coordinator: restart drivers, stop hung OVH agents, refill the queue. Runs from a systemd timer.
[ $(date +%s) -ge 1790948885 ] && exit 0
cd "$(dirname "$0")/.." || exit 1
systemctl --user is-active -q bt-queue-ovh || systemd-run --user --collect --unit=bt-queue-ovh --slice=jobs-bodytwin.slice -p MemoryMax=300M --working-directory=$PWD /bin/bash $PWD/tasks/lanes/ovh_agents/bt_queue_ovh.sh
systemctl --user is-active -q bt-queue-driver || systemd-run --user --collect --unit=bt-queue-driver --slice=jobs-bodytwin.slice -p MemoryMax=200M --working-directory=$PWD /bin/bash $PWD/tasks/lanes/bt_queue2.sh
timeout 60 ssh -i ~/.ssh/hunt_20260923 -o UserKnownHostsFile=external_research_path -o BatchMode=yes -n ubuntu@51.77.110.4 'now=$(date +%s); for u in $(systemctl list-units --no-legend "agent-*" | awk "{print \$1}"); do d=$(systemctl show -p WorkingDirectory --value $u); [ $((now-$(stat -c %Y $d/agent.log 2>/dev/null||echo $now))) -gt 1500 ] && sudo systemctl stop $u; done'
n=$(python3 tasks/queue_supply.py)
if [ $n -lt 40 ] && true; then python3 tasks/refill_v2.py 80; fi
echo "[$(date +%T)] keepalive runnable=$n ovh=$(ls tasks/lanes/ovh_agents/running | wc -l) local=$(systemctl --user list-units --no-legend 'bt-job-*' | wc -l)" >> tasks/lanes/bunny_keepalive.log
# Free result-driven planner replaces automatic paid lane runner maintenance.
n2=$(python3 tasks/queue_supply.py)
if [ "$n2" -lt 40 ]; then python3 tasks/free_controller.py; fi
# (2) Cap: raise it while the success rate is good and OVH has memory; lower it on many failures.
ok=$(grep fetched tasks/lanes/ovh_agents/queue_ovh.log | tail -40 | grep -c results=yes); cap=$(cat tasks/lanes/ovh_agents/cap)
av=$(timeout 30 ssh -i ~/.ssh/hunt_20260923 -o UserKnownHostsFile=external_research_path -o BatchMode=yes -n ubuntu@51.77.110.4 'free -g|sed -n 2p' | awk '{print $7}')
if [ "${av:-0}" -gt 8 ] && [ $ok -ge 24 ] && [ $cap -lt 22 ]; then echo $((cap+5 > 22 ? 22 : cap+5)) > tasks/lanes/ovh_agents/cap; fi
if [ $ok -lt 12 ] && [ $cap -gt 10 ]; then echo $((cap-5)) > tasks/lanes/ovh_agents/cap; fi
