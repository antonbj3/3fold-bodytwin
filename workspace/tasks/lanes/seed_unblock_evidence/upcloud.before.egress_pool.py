"""Per-(model, IP) egress pool: spread free-model jobs over several egress IPs.

The provider rate limit is per egress IP (observed by paired route probes on a
fixed account/model/task: OVH direct rate-limits while the home IP passes). Each
(model, egress) is therefore an independent bucket with its own backoff; a job is
assigned by weighted least-load so a limited IP is avoided without pausing the
whole model. Enable by creating /opt/agents/EGRESS_POOL.json.
"""
from pathlib import Path
import fcntl
import hashlib
import json
import math
import os
import shutil
import time

ROOT = Path('/opt/agents')
POOL = ROOT / 'EGRESS_POOL.json'
STATE = ROOT / 'EGRESS_POOL_STATE.json'
LOCK = ROOT / 'EGRESS_POOL.lock'


def _load(path, default):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default



def _model_key(model):
    name = model.rsplit('/', 1)[-1]
    return {'space-swarm-free': 'swarm', 'swarm_worker-2.5-preview-free': 'swarm_worker'}.get(name, name)


def _normalise_state(state):
    # Runtime records short names while apply_route supplies provider model IDs.
    # Preserve all accumulated counts and the latest backoff when merging aliases.
    for field in ('counts', 'backoff'):
        merged = {}
        for key, value in state.get(field, {}).items():
            model, separator, egress = key.partition('|')
            canonical = _model_key(model) + separator + egress
            if field == 'counts':
                merged[canonical] = merged.get(canonical, 0) + value
            else:
                merged[canonical] = max(merged.get(canonical, 0), value)
        state[field] = merged
    return state


def _healthy(entry, now):
    health_file = entry.get('health_file')
    if not health_file:
        return True
    health = _load(Path(health_file), {})
    try:
        return health.get('ready') is True and float(health.get('valid_until', 0)) > now
    except (TypeError, ValueError):
        return False


def enabled():
    pool = _load(POOL, None)
    return bool(pool and pool.get('enabled'))


def _write_state(state):
    tmp = STATE.with_suffix('.tmp')
    tmp.write_text(json.dumps(state, indent=2) + '\n')
    tmp.replace(STATE)


def _active_counts():
    """Live concurrent jobs per egress, from active slot markers + job routes."""
    counts = {}
    active = ROOT / 'active_models'
    if not active.is_dir():
        return counts
    for p in active.glob('*.json'):
        try:
            row = json.loads(p.read_text())
            os.kill(row['pid'], 0)
        except (OSError, ValueError, KeyError):
            continue
        route = ROOT / 'jobs' / str(row.get('job', '')) / 'EGRESS_ROUTE.json'
        try:
            mode = json.loads(route.read_text()).get('mode')
        except (OSError, ValueError):
            continue
        if mode:
            counts[mode] = counts.get(mode, 0) + 1
    return counts


def pick(model, directory=None):
    """Return (egress_name, proxy) for this job, or None when the pool is off
    or every egress is in backoff (caller then keeps its default route)."""
    pool = _load(POOL, None)
    if not pool or not pool.get('enabled'):
        return None
    model = _model_key(model)
    now = time.time()
    with LOCK.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        state = _load(STATE, {'window_start': 0, 'counts': {}, 'backoff': {}})
        if now - state.get('window_start', 0) >= pool.get('window_seconds', 3600):
            state = {'window_start': now, 'counts': {}, 'backoff': state.get('backoff', {})}
        state = _normalise_state(state)
        backoff = state.setdefault('backoff', {})
        counts = state.setdefault('counts', {})
        names = list(pool.get('egress', {}))
        live = _active_counts()
        eligible = []
        for name in names:
            if not _healthy(pool['egress'][name], now):
                backoff[f'{model}|{name}'] = max(backoff.get(f'{model}|{name}', 0), now + pool.get('backoff_seconds', 900))
                continue
            if backoff.get(f'{model}|{name}', 0) > now:
                continue
            cap = int(pool['egress'][name].get('cap', 0) or 0)
            if cap and live.get(name, 0) >= cap:
                continue
            eligible.append(name)
        if not eligible:
            _write_state(state)
            return None

        def load(name):
            weight = float(pool['egress'][name].get('weight', 1) or 1)
            if not math.isfinite(weight) or weight <= 0:
                weight = 1.0
            return counts.get(f'{model}|{name}', 0) / max(0.000001, weight)

        chosen = min(eligible, key=lambda n: (load(n), n))
        counts[f'{model}|{chosen}'] = counts.get(f'{model}|{chosen}', 0) + 1
        _write_state(state)
    # Already-running admission waiters may hold the older apply_route() in
    # memory. Prime their exact existing runtime cache too, without restarting
    # jobs or changing account stores. MODEL_ROUTE is written before environment().
    if directory and pool['egress'][chosen].get('model_catalog'):
        metadata = _load(Path(directory) / 'MODEL_ROUTE.json', {})
        account = metadata.get('account')
        if account in ('A', 'B', 'C', 'D'):
            tag = hashlib.sha256((account + str(Path(directory).resolve())).encode()).hexdigest()[:24]
            try:
                prepare_cache(chosen, {'XDG_CACHE_HOME': str(ROOT / 'runtime' / tag / 'cache')}, model)
            except Exception:
                record(model, chosen)
                return pick(model, directory)
    return chosen, pool['egress'][chosen].get('proxy')


def record(model, egress):
    """Mark (model, egress) as rate-limited. Return True if ALL egresses for the
    model are now in backoff (the caller should then pause the model itself)."""
    pool = _load(POOL, None)
    if not pool or not pool.get('enabled') or not egress:
        return False
    model = _model_key(model)
    now = time.time()
    with LOCK.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        state = _load(STATE, {'window_start': now, 'counts': {}, 'backoff': {}})
        state = _normalise_state(state)
        backoff = state.setdefault('backoff', {})
        backoff[f'{model}|{egress}'] = now + pool.get('backoff_seconds', 900)
        _write_state(state)
        names = list(pool.get('egress', {}))
        return all(backoff.get(f'{model}|{n}', 0) > now or not _healthy(pool['egress'][n], now) for n in names)


def prepare_cache(egress, env, model):
    """Prime a cold OpenCode runtime with public metadata for a slow new route.

    OpenCode may start from its embedded catalog before a background fetch
    completes. Never copy auth, provider settings, or another job's cache.
    """
    pool = _load(POOL, {})
    source = pool.get('egress', {}).get(egress, {}).get('model_catalog')
    if not source:
        return
    catalog = Path(source)
    models = json.loads(catalog.read_text()).get('opencode', {}).get('models', {})
    model_id = model.rsplit('/', 1)[-1]
    model_id = {'swarm_worker': 'swarm_worker-2.5-preview-free', 'swarm': 'space-swarm-free'}.get(model_id, model_id)
    if model_id not in models:
        raise ValueError('Requested free model missing from public catalog')
    cache = Path(env['XDG_CACHE_HOME']) / 'opencode' / 'models.json'
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix('.warp.tmp')
    shutil.copyfile(catalog, tmp)
    tmp.replace(cache)
