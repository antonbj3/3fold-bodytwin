"""Explicit Go fallback; only after BOTH free models enter cooldown.

Stops at reserve_worker_weekly_used_pct >= reserve_worker_weekly_stop_pct so a few percent
of the weekly reserve_worker allowance is preserved for manual setup work.
"""
def choose(policy,counts,now):
    if now>=policy.get('reserve_worker_reserve_until',0):return None
    if not all(now<policy.get(m+'_backoff_until',0) for m in ('swarm','free_worker')):return None
    if now<policy.get('reserve_worker_backoff_until',0):return None
    if policy.get('reserve_worker_weekly_used_pct',0)>=policy.get('reserve_worker_weekly_stop_pct',101):return None
    if sum(v for (a,m),v in counts.items() if m=='reserve_worker')>=policy.get('reserve_worker_reserve_cap',0):return None
    return policy.get('reserve_worker_reserve_account','C')
