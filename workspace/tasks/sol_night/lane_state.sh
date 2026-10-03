#!/bin/bash
# Lane state compared against slots.txt, NOT against list-units.
# 3/10 18:4x (anton-5f): two lanes that had finished real rounds (365,956 tokens, 19/19 numeric checks)
# sat idle and appeared in NEITHER of my checks, because a finished unit is garbage-collected and then
# falls out of `list-units` entirely. A silent exit can only be found against the slot list.
cd "$(dirname "$0")/../.." || exit 1
D=tasks/build_night; ST=$D/state
grep -v '^#' "$D/slots.txt" | awk 'NF{print $1}' | sed 's/.*>//' | sort -u > /tmp/.slots.$$
systemctl --user list-units 'bt-solnight-*' --state=active --no-legend 2>/dev/null \
  | sed 's/.*bt-solnight-//;s/\.service.*//' | sort > /tmp/.active.$$
printf 'slots %s, aktiva %s\n' "$(wc -l < /tmp/.slots.$$)" "$(wc -l < /tmp/.active.$$)"
while read -r L; do
  now=$(date +%s); back=$(cat "$ST/$L.backoff_until" 2>/dev/null || echo 0)
  printf '  %-34s go=%s rund=%s %s\n' "$L" \
    "$([ -f "$ST/$L.go" ] && echo ja || echo NEJ)" \
    "$(cat "$ST/$L.round" 2>/dev/null || echo -)" \
    "$([ "$now" -lt "$back" ] && echo "backoff $(( (back-now)/60 )) min" || echo '')"
done < <(comm -23 /tmp/.slots.$$ /tmp/.active.$$)
rm -f /tmp/.slots.$$ /tmp/.active.$$
