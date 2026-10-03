#!/usr/bin/env python3
"""Shared per-host admission: reserve full cgroup limits, serialize start calls."""
import argparse
import fcntl
import json
import os
import pwd
from pathlib import Path
import subprocess
import sys

MIB = 1024**2

def snapshot(agent_mib, slice_headroom, host_headroom, host_cap):
    raw = subprocess.check_output(['systemctl', 'show', 'research.slice', '-p', 'MemoryMax', '-p', 'MemoryCurrent'], text=True)
    props = dict(line.split('=', 1) for line in raw.splitlines())
    # 1/10 (anton-5f, Field's finding): an INACTIVE research.slice gives MemoryCurrent='[not set]' and
    # MemoryMax='infinity'. int() raised ValueError -> snapshot failed -> 0 slots to ALL drivers,
    # silently. UpCloud stood idle 14:23-21:30 from exactly this. An inactive slice uses 0 bytes, so
    # '[not set]' is read as 0; 'infinity' has no cap and falls back on the host's own memory.
    def _bytes(name, fallback):
        raw_value = props.get(name, '').strip()
        if raw_value in ('[not set]', '', 'infinity'):
            return fallback
        try:
            return int(raw_value)
        except ValueError:
            return fallback
    host_total = os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE')
    limit, current = _bytes('MemoryMax', host_total), _bytes('MemoryCurrent', 0)
    base = Path('/sys/fs/cgroup/research.slice')
    reserved = child_current = total = 0
    if not base.is_dir():
        # The slice does not exist until it has been started; no children to sum, the capacity is the whole cap.
        return dict(slots=max(0, (limit - slice_headroom*MIB)//(agent_mib*MIB)), total=0,
                    slice_max=limit, slice_current=0, slice_inactive=True,
                    agent_mib=agent_mib, slice_headroom_mib=slice_headroom,
                    host_headroom_mib=host_headroom, host_cap=host_cap)
    for child in base.iterdir():
        if not child.is_dir():
            continue
        try:
            used = int((child/'memory.current').read_text())
            maximum = int((child/'memory.max').read_text())
        except FileNotFoundError:
            continue  # unit exited while taking the snapshot
        # Unbounded children fail closed (int('max') raises ValueError).
        reserved += max(used, maximum)
        child_current += used
        total += child.name.startswith('agent-')
    reserved += max(0, current-child_current)
    available = next(int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:'))
    slots = max(0, min((limit-reserved-slice_headroom*MIB)//(agent_mib*MIB),
                       (limit-current-slice_headroom*MIB)//(agent_mib*MIB),
                       (available-host_headroom*MIB)//(agent_mib*MIB), host_cap-total))
    if __import__('shutil').disk_usage('/opt/agents').free < 2*1024**3:
        slots = 0
    return dict(slots=slots, total=total, slice_max=limit, slice_current=current,
                reserved=reserved, mem_available=available, agent_mib=agent_mib,
                slice_headroom_mib=slice_headroom, host_headroom_mib=host_headroom, host_cap=host_cap)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--agent-mib', type=int, default=1500)
    p.add_argument('--slice-headroom-mib', type=int, default=2048)
    p.add_argument('--host-headroom-mib', type=int, default=4096)
    p.add_argument('--host-cap', type=int, default=17)
    p.add_argument('--json', action='store_true')
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    if a.agent_mib <= 0 or min(a.slice_headroom_mib, a.host_headroom_mib, a.host_cap) < 0:
        return 64
    # This file and its lock are always created/opened as ubuntu.
    original_uid = os.geteuid()
    if original_uid == 0:
        os.seteuid(pwd.getpwnam('ubuntu').pw_uid)
    try:
        lock_file = Path('/opt/agents/capacity_admission.lock').open('a')
    finally:
        os.seteuid(original_uid)
    with lock_file as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            s = snapshot(a.agent_mib, a.slice_headroom_mib, a.host_headroom_mib, a.host_cap)
        except (OSError, ValueError, KeyError, subprocess.SubprocessError) as e:
            print('capacity snapshot failed: '+str(e), file=sys.stderr)
            return 75
        command = a.command[1:] if a.command[:1] == ['--'] else a.command
        if not command:
            print(json.dumps(s) if a.json else s['slots'])
            return 0
        if s['slots'] < 1:
            print('capacity admission deferred: '+json.dumps(s), file=sys.stderr)
            return 75
        # Synchronous systemd-run registers the unit before releasing the lock.
        return subprocess.run(command).returncode

if __name__ == '__main__':
    raise SystemExit(main())
