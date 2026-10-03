"""Explicit Go fallback after free admission is unavailable.

Stops at reserve_worker_weekly_used_pct >= reserve_worker_weekly_stop_pct so a few percent
of the weekly reserve_worker allowance is preserved for manual setup work.
"""
def choose(policy,counts,now,*,free_unavailable=False):
    if now>=policy.get('reserve_worker_reserve_until',0):return None
    if not free_unavailable and not all(now<policy.get(m+'_backoff_until',0) for m in ('swarm','free_worker')):return None
    if now<policy.get('reserve_worker_backoff_until',0):return None
    if policy.get('reserve_worker_weekly_used_pct',0)>=min(95,policy.get('reserve_worker_weekly_stop_pct',95)):return None
    if sum(v for (a,m),v in counts.items() if m=='reserve_worker')>=policy.get('reserve_worker_reserve_cap',0):return None
    return policy.get('reserve_worker_reserve_account','C')
