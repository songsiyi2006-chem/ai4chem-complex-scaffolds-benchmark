"""Opt-in LLM boundary: no included provider, no API key discovery or network call.

Inject a reviewed provider in a future deployment. Unknown response usage retains
the reservation; reconcile with the same provider idempotency key before retrying.
"""
from dataclasses import dataclass
from typing import Protocol
import json
from swarm.memory_ledger import canonical, digest


@dataclass(frozen=True)
class Response:
    receipt_id: str
    input_tokens: int
    output_tokens: int
    cached_input_tokens: int
    structured_result: dict


class Provider(Protocol):
    def invoke(self, *, idempotency_key: str, role: str, artifact_references: list,
               max_total_tokens: int) -> Response:
        """Return actual provider usage, not an estimate from text length."""
        ...


def invoke_reserved(ledger, provider, request_id, role, artifact_references, worst_case_tokens):
    signature = digest(dict(role=role, references=artifact_references, cap=worst_case_tokens))
    signature_key = "provider_request:" + request_id
    response_key = "provider_response:" + request_id
    with ledger.tx():
        existing = ledger.db.execute("SELECT value FROM meta WHERE key=?", (signature_key,)).fetchone()
        if existing and existing[0] != signature:
            raise ValueError("idempotency key reused with changed artifact references")
        ledger.db.execute("INSERT OR IGNORE INTO meta VALUES(?,?)", (signature_key, signature))
    ledger.reserve(request_id, role, worst_case_tokens)
    old = ledger.db.execute("SELECT value FROM meta WHERE key=?", (response_key,)).fetchone()
    if old:
        response = Response(**json.loads(old[0]))
    else:
        # Provider-side idempotency is REQUIRED across a crash or concurrent calls.
        response = provider.invoke(idempotency_key=request_id, role=role,
                                   artifact_references=artifact_references,
                                   max_total_tokens=worst_case_tokens)
        serialized = canonical(response.__dict__)
        with ledger.tx():
            old = ledger.db.execute("SELECT value FROM meta WHERE key=?", (response_key,)).fetchone()
            if old and old[0] != serialized:
                raise ValueError("provider violated idempotency; keep reservation for reconciliation")
            ledger.db.execute("INSERT OR IGNORE INTO meta VALUES(?,?)", (response_key, serialized))
    ledger.settle(request_id, response.receipt_id, response.input_tokens,
                  response.output_tokens, response.cached_input_tokens)
    return response.structured_result
