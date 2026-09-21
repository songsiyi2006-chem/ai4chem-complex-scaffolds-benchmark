"""Transactional single-host queue, fenced leases, receipts and hash-linked events.

SQLite is NOT a multi-host/NFS scheduler or a tamper-proof external notary.
One connection per worker. All changes use BEGIN IMMEDIATE transactions.
"""
import hashlib
import json
import math
import sqlite3
import time
from contextlib import contextmanager

ROLES = ("Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Omega")
SHARES = dict(zip(ROLES, (10, 25, 15, 25, 5, 10, 10)))
TERMINAL = {"SUCCEEDED", "BLOCKED", "FAILED"}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


class Ledger:
    def __init__(self, path, token_cap=100_000_000):
        if type(token_cap) is not int or token_cap < 0:
            raise ValueError("token_cap must be a nonnegative integer")
        self.db = sqlite3.connect(str(path), timeout=30, isolation_level=None)
        self.db.row_factory = sqlite3.Row
        self.db.execute("PRAGMA foreign_keys=ON")
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.executescript('''
        CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS tasks(
          id TEXT PRIMARY KEY, role TEXT NOT NULL, deps TEXT NOT NULL,
          status TEXT NOT NULL, owner TEXT, lease REAL, fence INTEGER NOT NULL DEFAULT 0,
          result TEXT);
        CREATE TABLE IF NOT EXISTS events(
          seq INTEGER PRIMARY KEY, body TEXT NOT NULL, previous TEXT NOT NULL, hash TEXT NOT NULL);
        CREATE TRIGGER IF NOT EXISTS immutable_events_update BEFORE UPDATE ON events
          BEGIN SELECT RAISE(ABORT, 'append-only events'); END;
        CREATE TRIGGER IF NOT EXISTS immutable_events_delete BEFORE DELETE ON events
          BEGIN SELECT RAISE(ABORT, 'append-only events'); END;
        CREATE TABLE IF NOT EXISTS tokens(
          id TEXT PRIMARY KEY, role TEXT NOT NULL, reserved INTEGER NOT NULL,
          actual INTEGER, receipt TEXT UNIQUE);
        CREATE TABLE IF NOT EXISTS blacklist(
          molecule_hash TEXT PRIMARY KEY, reason TEXT NOT NULL, reviewer TEXT NOT NULL);
        ''')
        try:
            with self.tx():
                row = self.db.execute("SELECT value FROM meta WHERE key='token_cap'").fetchone()
                if row and int(row[0]) != token_cap:
                    raise ValueError("frozen token cap cannot change on resume")
                self.db.execute("INSERT OR IGNORE INTO meta VALUES('token_cap',?)", (str(token_cap),))
        except BaseException:
            self.db.close()
            raise
        self.cap = token_cap

    @contextmanager
    def tx(self):
        self.db.execute("BEGIN IMMEDIATE")
        try:
            yield
            self.db.execute("COMMIT")
        except BaseException:
            self.db.execute("ROLLBACK")
            raise

    def _event(self, kind, payload):
        last = self.db.execute("SELECT hash FROM events ORDER BY seq DESC LIMIT 1").fetchone()
        prev = last[0] if last else "0" * 64
        body = canonical(dict(kind=kind, payload=payload, unix_time=time.time()))
        sha = hashlib.sha256((prev + body).encode()).hexdigest()
        self.db.execute("INSERT INTO events(body,previous,hash) VALUES(?,?,?)", (body, prev, sha))

    def add(self, task_id, role, deps=()):
        if role not in ROLES or task_id in deps or len(set(deps)) != len(deps):
            raise ValueError("invalid role or dependencies")
        deps = canonical(list(deps))
        with self.tx():
            old = self.db.execute("SELECT * FROM tasks WHERE id=?", (task_id,)).fetchone()
            if old:
                if old['role'] != role or old['deps'] != deps:
                    raise ValueError("task definition changed")
                return
            for parent in json.loads(deps):
                if not self.db.execute("SELECT 1 FROM tasks WHERE id=?", (parent,)).fetchone():
                    raise ValueError("add dependencies first; forward references forbidden")
            self.db.execute("INSERT INTO tasks(id,role,deps,status) VALUES(?,?,?,'PENDING')",
                            (task_id, role, deps))
            self._event("created", dict(task=task_id, role=role))

    def claim(self, owner, lease_seconds=600, now=None):
        if not owner or not math.isfinite(lease_seconds) or lease_seconds <= 0:
            raise ValueError("invalid worker/lease")
        now = time.time() if now is None else now
        with self.tx():
            rows = self.db.execute("SELECT * FROM tasks ORDER BY rowid").fetchall()
            states = {r['id']: r['status'] for r in rows}
            for r in rows:
                available = r['status'] == 'PENDING' or (r['status'] == 'RUNNING' and r['lease'] <= now)
                if not available or any(states[p] not in TERMINAL for p in json.loads(r['deps'])):
                    continue
                fence = r['fence'] + 1
                self.db.execute("UPDATE tasks SET status='RUNNING',owner=?,lease=?,fence=? WHERE id=?",
                                (owner, now + lease_seconds, fence, r['id']))
                self._event("claimed", dict(task=r['id'], role=r['role'], owner=owner, fence=fence))
                return dict(id=r['id'], role=r['role'], owner=owner, fence=fence)
        return None

    def finish(self, claim, status, result, now=None):
        if status not in TERMINAL:
            raise ValueError("not a terminal state")
        now = time.time() if now is None else now
        with self.tx():
            row = self.db.execute("SELECT * FROM tasks WHERE id=?", (claim['id'],)).fetchone()
            if (not row or row['status'] != 'RUNNING' or row['owner'] != claim['owner']
                    or row['fence'] != claim['fence'] or row['lease'] <= now):
                raise RuntimeError("stale or expired lease; result rejected")
            self.db.execute("UPDATE tasks SET status=?,result=?,lease=NULL WHERE id=?",
                            (status, canonical(result), claim['id']))
            self._event(status, dict(task=claim['id'], role=row['role'], result=result))

    def reserve(self, request_id, role, maximum_tokens):
        if role not in ROLES or type(maximum_tokens) is not int or maximum_tokens <= 0:
            raise ValueError("invalid reservation")
        with self.tx():
            old = self.db.execute("SELECT * FROM tokens WHERE id=?", (request_id,)).fetchone()
            if old:
                if old['role'] != role or old['reserved'] != maximum_tokens:
                    raise ValueError("idempotency key reused with changed request")
                return
            rows = self.db.execute("SELECT role,COALESCE(actual,reserved) AS n FROM tokens").fetchall()
            if (sum(r['n'] for r in rows) + maximum_tokens > self.cap or
                sum(r['n'] for r in rows if r['role'] == role) + maximum_tokens > self.cap * SHARES[role] // 100):
                raise RuntimeError("token ceiling exceeded")
            self.db.execute("INSERT INTO tokens(id,role,reserved) VALUES(?,?,?)",
                            (request_id, role, maximum_tokens))
            self._event("token_reserved", dict(id=request_id, role=role, maximum=maximum_tokens))

    def settle(self, request_id, receipt_id, input_tokens, output_tokens, cached_input_tokens=0):
        values = (input_tokens, output_tokens, cached_input_tokens)
        if not receipt_id or any(type(n) is not int or n < 0 for n in values) or cached_input_tokens > input_tokens:
            raise ValueError("invalid provider usage; cached input is already INCLUDED in input")
        actual = input_tokens + output_tokens
        receipt = canonical(dict(id=receipt_id, input=input_tokens, output=output_tokens, cached=cached_input_tokens))
        with self.tx():
            row = self.db.execute("SELECT * FROM tokens WHERE id=?", (request_id,)).fetchone()
            if not row:
                raise ValueError("unreserved call")
            if row['receipt'] is not None:
                if row['receipt'] != receipt:
                    raise ValueError("receipt conflict")
                return
            for r in self.db.execute("SELECT receipt FROM tokens WHERE receipt IS NOT NULL"):
                if json.loads(r[0])['id'] == receipt_id:
                    raise ValueError("provider receipt already charged to another call")
            # Charge overruns in full. Subsequent reservations are blocked by the cap.
            self.db.execute("UPDATE tokens SET actual=?,receipt=? WHERE id=?", (actual, receipt, request_id))
            self._event("token_settled", dict(id=request_id, actual=actual, overrun=actual > row['reserved']))

    def snapshot(self):
        tasks = [dict(r) for r in self.db.execute("SELECT * FROM tasks ORDER BY rowid")]
        rows = self.db.execute("SELECT actual,reserved FROM tokens").fetchall()
        return dict(tasks=tasks, token_cap=self.cap, actual_provider_tokens=sum(r['actual'] or 0 for r in rows),
                    unresolved_reserved_tokens=sum(r['reserved'] for r in rows if r['actual'] is None))

    def reject_molecule(self, molecule_hash, reason, reviewer):
        if (len(molecule_hash) != 64 or any(c not in "0123456789abcdef" for c in molecule_hash)
                or reason not in {"INVALID_VALENCE", "STEREOCHEMISTRY_MISMATCH", "REVIEWED_ROUTE_POLICY"}
                or not reviewer):
            raise ValueError("reviewed identity/route reason required; numerical failure is not chemical rejection")
        with self.tx():
            self.db.execute("INSERT INTO blacklist VALUES(?,?,?)", (molecule_hash, reason, reviewer))
            self._event("molecule_rejected", dict(molecule_hash=molecule_hash, reason=reason, reviewer=reviewer))

    def verify_events(self):
        prev = "0" * 64
        for r in self.db.execute("SELECT * FROM events ORDER BY seq"):
            if r['previous'] != prev or hashlib.sha256((prev + r['body']).encode()).hexdigest() != r['hash']:
                raise ValueError("broken event chain")
            prev = r['hash']
        return prev

    def events(self):
        self.verify_events()
        return [dict(r) for r in self.db.execute("SELECT * FROM events ORDER BY seq")]

    def close(self):
        self.db.close()
