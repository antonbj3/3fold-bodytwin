#!/bin/bash
# lane runner mode from the session logs: the weekly window AND the credits.
# 2/10 21:30 (anton-5f): the old version only read "primary".used_percent and reported
# "91 % (week)", which I read as a cap and almost closed the eight lane wave on. The same file carries
# "credits":{"has_credits":true,"unlimited":false,"balance":"62500"} — Sun does NOT stop at 100 % by
# the weekly window, it continues on credits. A measure that shows half the truth is worse than nothing,
# because it reads like a check. Fifth instrument today with that characteristic.
f=$(ls -t ~/.lane_runner/sessions/*/*/*/*.jsonl 2>/dev/null | head -1)
l=$(grep -o '"rate_limits":{.*' "$f" 2>/dev/null | tail -1)
p=$(echo "$l" | grep -o '"used_percent":[0-9.]*' | head -1 | cut -d: -f2)
r=$(echo "$l" | grep -o '"resets_at":[0-9]*' | head -1 | cut -d: -f2)
bal=$(echo "$l" | grep -o '"balance":"[0-9]*"' | head -1 | grep -o '[0-9]*')
has=$(echo "$l" | grep -o '"has_credits":[a-z]*' | head -1 | cut -d: -f2)
unl=$(echo "$l" | grep -o '"unlimited":[a-z]*' | head -1 | cut -d: -f2)
spend=$(echo "$l" | grep -o '"spend_control_reached":[a-z]*' | head -1 | cut -d: -f2)
echo "lane_runner weekly window ${p:-?}% (resets $( [ -n "$r" ] && date -d @$r '+%a %d/%m %H:%M' || echo '?' ))"
echo "lane_runner credits: has_credits=${has:-?} balance=${bal:-?} unlimited=${unl:-?} spend_control_reached=${spend:-?}"
if [ "$has" = true ] && [ "${bal:-0}" -gt 0 ] 2>/dev/null; then
  echo "=> THE WEEKLY WINDOW IS NO STOP: credits available, keep running."
else
  echo "=> no credits visible in the log — check before counting the window as a stop."
fi
