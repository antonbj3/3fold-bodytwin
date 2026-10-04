#!/usr/bin/env python3
"""Estimate weekly reserve_worker spend from job logs and update the cloud guard.

The provider exposes usage only in its console, so this derives an estimate from
the `cost` field opencode writes into every job's agent.log over a rolling
window, per account, and writes reserve_worker_weekly_used_pct to each owned host's
FREE_MODEL_POLICY.json so fallback_reserve.choose stops reserve_worker at the threshold.

Conservatism: the cloud console counts a fixed calendar week; a rolling window
of window_days is >= that, so the estimate stops the account a little early,
which preserves at least the requested margin.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INFRA = Path('external_research_path')
STATE = ROOT / 'reserve_worker_budget_state.json'
CONF = ROOT / 'reserve_worker_budget.json'
SSH = ['ssh', '-i', 'local_config_path/hunt_20260923', '-o', 'BatchMode=yes',
       '-o', 'ConnectTimeout=15', '-o', 'UserKnownHostsFile=' + str(INFRA / 'known_hosts')]
HOSTS = ['ubuntu@51.77.110.4', 'root@212.147.226.179']
SCAN_ROOTS = [Path('external_mount')]
COST_RE = re.compile(r'"cost":([0-9.]+)')

DEFAULT_CONF = {
    'window_days': 7,
    'stop_pct': 95,
    'weekly_budget_usd': {'A': 30.0, 'B': 30.0, 'C': 30.0},
    'note': 'reserve_worker V4.1 Flash weekly allowance = 50%% of monthly limit: Go $30, Go Plus $60',
}


def load_json(path: Path, default):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def save_json(path: Path, data) -> None:
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(path)


def parse_log(path: Path) -> tuple[float, str]:
    try:
        text = path.read_text(errors='replace')
    except OSError:
        return 0.0, '?'
    cost = sum(float(x) for x in COST_RE.findall(text))
    account = '?'
    route = path.parent / 'MODEL_ROUTE.json'
    if route.exists():
        try:
            account = json.loads(route.read_text()).get('account', '?')
        except (OSError, ValueError):
            pass
    return cost, account


def local_spend(window_days: int) -> tuple[dict, dict]:
    cutoff = time.time() - window_days * 86400
    cache = load_json(STATE, {}).get('files', {})
    fresh: dict = {}
    spend: dict = {}
    seen = 0
    for root in SCAN_ROOTS:
        if not root.is_dir():
            continue
        for dirpath, _dirs, files in os.walk(root):
            if 'agent.log' not in files:
                continue
            path = Path(dirpath) / 'agent.log'
            try:
                st = path.stat()
            except OSError:
                continue
            if st.st_mtime < cutoff:
                continue
            seen += 1
            key = str(path)
            row = cache.get(key)
            if not row or row.get('mtime') != st.st_mtime or row.get('size') != st.st_size:
                cost, account = parse_log(path)
                row = {'mtime': st.st_mtime, 'size': st.st_size, 'cost': cost, 'account': account}
            fresh[key] = row
            if row['cost'] > 0:
                spend[row['account']] = spend.get(row['account'], 0.0) + row['cost']
    save_json(STATE, {'updated': time.time(), 'files': fresh})
    return spend, {'logs_in_window': seen, 'cached': len(fresh)}


REMOTE = r'''
import json,os,re
from pathlib import Path
cost_re=re.compile(r'"cost":([0-9.]+)')
spend={}; n=0
for route in Path('/opt/agents/jobs').glob('*/MODEL_ROUTE.json'):
    d=route.parent
    log=d/'agent.log'
    if not log.exists(): continue
    try: txt=log.read_text(errors='replace')
    except OSError: continue
    c=sum(float(x) for x in cost_re.findall(txt))
    if c<=0: continue
    try: acct=json.loads(route.read_text()).get('account','?')
    except (OSError,ValueError): acct='?'
    spend[acct]=spend.get(acct,0.0)+c; n+=1
print(json.dumps({'spend':spend,'jobs_with_cost':n}))
'''


def remote_spend(host: str) -> tuple[dict, int]:
    try:
        out = subprocess.run(SSH + [host, 'python3 -'], input=REMOTE, text=True,
                             capture_output=True, timeout=120).stdout
        data = json.loads(out)
        return data.get('spend', {}), data.get('jobs_with_cost', 0)
    except (subprocess.SubprocessError, ValueError):
        return {}, 0


def host_reserve_account(host: str) -> str | None:
    code = ("import json;print(json.load(open('/opt/agents/FREE_MODEL_POLICY.json'))"
            ".get('reserve_worker_reserve_account','C'))")
    try:
        out = subprocess.run(SSH + [host, 'python3 -c ' + shlex.quote(code)], text=True,
                             capture_output=True, timeout=60).stdout.strip()
        return out or None
    except subprocess.SubprocessError:
        return None


def push(host: str, account: str, used_pct: float, conf: dict) -> bool:
    payload = {'used': round(used_pct, 2), 'stop': conf['stop_pct'],
               'budget': conf['weekly_budget_usd'].get(account, 30.0)}
    code = ('import json,fcntl,time,os\n'
            'P=' + json.dumps(payload) + '\n'
            'root="/opt/agents";pol=root+"/FREE_MODEL_POLICY.json"\n'
            'lock=open(root+"/rate_policy.lock","a");fcntl.flock(lock,fcntl.LOCK_EX)\n'
            'd=json.load(open(pol))\n'
            'd["reserve_worker_weekly_used_pct"]=P["used"]\n'
            'd["reserve_worker_weekly_stop_pct"]=P["stop"]\n'
            'd["reserve_worker_weekly_budget_usd"]=P["budget"]\n'
            'd["reserve_worker_weekly_updated"]=time.time()\n'
            'tmp=pol+".tmp";json.dump(d,open(tmp,"w"),indent=2);open(tmp,"a").write("\\n");os.replace(tmp,pol)\n'
            'print("ok")\n')
    try:
        r = subprocess.run(SSH + [host, 'sudo -u ubuntu python3 -'], input=code, text=True,
                           capture_output=True, timeout=60)
        return r.returncode == 0
    except subprocess.SubprocessError:
        return False


def main() -> int:
    conf = load_json(CONF, DEFAULT_CONF)
    for key, value in DEFAULT_CONF.items():
        conf.setdefault(key, value)
    save_json(CONF, conf)
    budget = conf['weekly_budget_usd']

    spend, local_info = local_spend(conf['window_days'])
    remote_info = {}
    for host in HOSTS:
        rs, n = remote_spend(host)
        remote_info[host] = n
        for acct, cost in rs.items():
            spend[acct] = spend.get(acct, 0.0) + cost

    result = {'time': time.time(), 'window_days': conf['window_days'],
              'stop_pct': conf['stop_pct'], 'spend_usd': spend,
              'used_pct': {}, 'local': local_info, 'remote_jobs_with_cost': remote_info}
    for host in HOSTS:
        acct = host_reserve_account(host)
        if acct is None:
            continue
        used = min(100.0, spend.get(acct, 0.0) / max(1e-9, budget.get(acct, 30.0)) * 100.0)
        result['used_pct'][host] = {'account': acct, 'used_pct': round(used, 2),
                                    'budget_usd': budget.get(acct, 30.0)}
        push(host, acct, used, conf)
    save_json(ROOT / 'reserve_worker_budget_last.json', result)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
