#!/usr/bin/env python3
"""Collect cloud lane artifacts from authenticated session events."""
import datetime as dt
import fcntl
import json
from pathlib import Path
import re
import sys
import time
import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT.parents[1] / 'results'
SESSION_RE = re.compile(r'session_[A-Za-z0-9]+\Z')
FILE_RE = re.compile(r'=====FILE ([A-Za-z0-9_./-]+)=====\s*\n(.*?)\n?=====END FILE=====', re.S)

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def api():
    # These are the same login files and endpoints as inspect_cloud.py. Never log responses.
    config = json.loads(Path('local_config_path').read_text())
    auth = json.loads(Path('local_config_path/.credentials.json').read_text())
    headers = {'Authorization': 'Bearer ' + auth['coordinatorAiOauth']['accessToken'],
               'anthropic-beta': 'ccr-byoc-2025-07-29',
               'x-organization-uuid': config['oauthAccount']['organizationUuid'],
               'anthropic-version': '2023-06-01', 'Accept': 'application/json'}
    def get(path):
        req = urllib.request.Request('https://api.anthropic.com' + path, headers=headers)
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=30) as response:
                    return json.load(response)
            except urllib.error.HTTPError as error:
                if error.code != 429 or attempt == 2:
                    raise
                time.sleep((5, 15)[attempt])
    return get

def event_texts(events):
    items = events.get('data', events.get('events', [])) if isinstance(events, dict) else events
    for event in items if isinstance(items, list) else []:
        message = event.get('message', {})
        if not isinstance(message, dict) or message.get('role') != 'assistant':
            continue
        content = message.get('content', [])
        if isinstance(content, list):
            yield message.get('model'), '\n'.join(part.get('text', '') for part in content
                                                   if isinstance(part, dict) and part.get('type') == 'text')

def extract(text):
    files = {}
    for name, content in FILE_RE.findall(text):
        rel = Path(name)
        if rel.is_absolute() or '..' in rel.parts or name.startswith('.'):
            continue
        if not (name in ('RESULTS.md', 'results.json') or
                (rel.suffix in ('.py', '.md', '.json', '.csv', '.txt') and len(rel.parts) <= 3)):
            continue
        files[name] = content.rstrip('\n') + '\n'
    if sum(len(x.encode()) for x in files.values()) > 200_000:
        return {}, 'artifact_limit_exceeded'
    return files, None

def main():
    args = sys.argv[1:]
    no_usage = '--no-usage' in args
    if no_usage:
        args.remove('--no-usage')
    if args and (len(args) != 2 or args[0] != '--refresh-lane'):
        raise SystemExit('usage: collect_cloud.py [--no-usage] [--refresh-lane CLOUD-<ID>]')
    refresh_lane = args[1] if args else None
    # The shared usage endpoint is reserved for BodyTwin's :00/:30 slots.
    if dt.datetime.now(dt.timezone.utc).minute not in (0, 30):
        no_usage = True
    get = api()
    usage = None
    if not no_usage:
        try:
            usage = get('/api/oauth/usage')
        except urllib.error.HTTPError as error:
            if error.code != 429:
                raise
            # The shared endpoint can rate-limit; collect events without a balance.
    credit = None
    if usage is not None:
        balances = [{k: value.get(k) for k in ('limit_dollars', 'used_dollars', 'remaining_dollars',
                                               'resets_at', 'locked_reason')}
                    for value in usage.values() if isinstance(value, dict) and value.get('limit_dollars') is not None]
        credit = next((x for x in balances if x.get('limit_dollars') == 250), None)
        if credit is None:
            raise SystemExit('promotional credit not visible; no collection/launch recommendation')
        with (ROOT / 'credit_log.jsonl').open('a') as out:
            fcntl.flock(out, fcntl.LOCK_EX)
            out.write(json.dumps({'time': now(), 'balance': credit}) + '\n')
    receipts = ROOT / 'receipts.jsonl'
    if not receipts.exists():
        print(json.dumps({'balance': credit, 'sessions': []}))
        return
    summaries = []
    for line in receipts.read_text().splitlines():
        receipt = json.loads(line)
        session = receipt.get('session_id', '')
        lane = receipt.get('lane', '')
        if not SESSION_RE.fullmatch(session) or not re.fullmatch(r'CLOUD-[A-Z0-9-]+', lane):
            continue
        dest = RESULTS / lane
        previous_file = dest / 'cloud_status.json'
        if previous_file.exists():
            previous = json.loads(previous_file.read_text())
            if lane != refresh_lane and previous.get('session_status') in ('idle', 'completed') and previous.get('files'):
                summaries.append(previous)
                continue
        try:
            status = get('/v1/sessions/' + session)
            session_status = status.get('session_status', status.get('status'))
            events = {'data': []}
            if session_status in ('idle', 'completed', 'failed'):
                # Pagination (Field lane, 24/9): the endpoint returns at most ~50 events per call with has_more/last_id;
                # without the loop the last message with the FILE blocks is lost.
                page = get('/v1/sessions/' + session + '/events?limit=1000')
                for _ in range(200):
                    items = page.get('data', page.get('events', [])) if isinstance(page, dict) else page
                    events['data'].extend(items)
                    if not (isinstance(page, dict) and page.get('has_more') and page.get('last_id')):
                        break
                    page = get('/v1/sessions/' + session + '/events?limit=1000&after_id=' + page['last_id'])
                events['has_more'] = isinstance(page, dict) and bool(page.get('has_more'))
        except urllib.error.HTTPError as error:
            summaries.append({'lane': lane, 'session_id': session, 'http_error': error.code})
            continue
        texts = list(event_texts(events))
        model = status.get('session_context', {}).get('model') or status.get('model_id')
        if not model:
            model = next((m for m, _ in texts if m), None)
        files = {}
        issue = 'events_truncated' if events.get('has_more') else None
        # Merge complete marked blocks across turns. A later follow-up may return
        # only a missing code file; the most recent block for each path wins.
        for _, body in texts:
            found, problem = extract(body)
            if problem:
                issue = problem
            files.update(found)
        if previous_file.exists():
            previous = json.loads(previous_file.read_text())
            for name in previous.get('files', []):
                if name not in files and (dest / name).is_file():
                    files[name] = (dest / name).read_text()
        dest.mkdir(parents=True, exist_ok=True)
        for name, content in files.items():
            path = dest / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)
        summary = {'time': now(), 'lane': lane, 'session_id': session,
                   'session_status': status.get('session_status', status.get('status')),
                   'model_id': model, 'files': sorted(files), 'collection_issue': issue,
                   'credit_remaining_dollars': (credit or {}).get('remaining_dollars')}
        (dest / 'cloud_status.json').write_text(json.dumps(summary, indent=2) + '\n')
        summaries.append(summary)
    print(json.dumps({'balance': credit, 'sessions': summaries}, indent=2))

if __name__ == '__main__':
    main()
