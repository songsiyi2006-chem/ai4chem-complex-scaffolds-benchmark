"""Seven deterministic roles with restartable local task scheduling."""
import importlib.util
import importlib.metadata
import json
from pathlib import Path
import uuid
from common import ROOT, load_module, write_json, sha256, source_hashes
from swarm.memory_ledger import Ledger, ROLES

DEPENDENCIES = dict(Alpha=[], Beta=["Alpha"], Gamma=["Beta"], Delta=["Alpha"],
                    Epsilon=["Alpha"], Zeta=["Beta", "Delta"], Omega=list(ROLES[:-1]))


def preflight():
    names = ("rdkit", "numpy", "scipy", "matplotlib", "pyscf", "openmm", "openmmplumed", "pymbar", "ase")
    versions = {}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return dict(packages={n: importlib.util.find_spec(n) is not None for n in names},
                distribution_versions=versions,
                production_enabled=False, llm_provider_connected=False,
                requested_method="original omegaB97X-D/def2-TZVP",
                method_blocker="PySCF original omegaB97X-D dispersion interface not supported",
                reference="BRD4 BD2 / VHL / macroPROTAC-1 (6SIS); PROTAC, not a proven molecular glue")


def prepare(out, mode, cap):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    for name in ("reports", "results/quantum_states", "results/fep_trajectories", "figures", "attempts"):
        (out / name).mkdir(parents=True, exist_ok=True)
    reference_hashes = {p.relative_to(ROOT).as_posix(): sha256(p) for p in (ROOT / "results/reference").glob("*") if p.is_file()}
    write_json(out / "manifest.json", dict(schema=1, mode=mode, token_cap=cap, source_hashes=source_hashes(),
                                         reference_hashes=reference_hashes,
                                         execution="LOCAL_DETERMINISTIC_ROLES_NO_LLM"))
    write_json(out / "preflight.json", preflight())
    (out / "reports/error_healing_log.jsonl").touch()
    ledger = Ledger(out / "ledger.sqlite", cap)
    for role in ROLES:
        ledger.add(role, role, DEPENDENCIES[role])
    ledger.close()
    return out


def validate_resume(out):
    out = Path(out).resolve()
    spec = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    if spec["source_hashes"] != source_hashes():
        raise ValueError("source changed since prepare; create a new workspace, do not mix provenance")
    for name, expected in spec.get("reference_hashes", {}).items():
        if not (ROOT / name).is_file() or sha256(ROOT / name) != expected:
            raise ValueError("reference input changed since prepare")
    return out, spec


def handle(role, out, attempt, mode):
    if role == "Alpha":
        if mode == "plan":
            return "SUCCEEDED", dict(evidence="SPECIFICATION_ONLY", molecules_evaluated=0)
        result = load_module("01_macrocycle_generator").generate(attempt / "conformers")
        write_json(attempt / "conformers.json", result)
        return "SUCCEEDED", dict(evidence=result["evidence"], artifact=(attempt / "conformers.json").relative_to(out).as_posix(),
                                 sha256=sha256(attempt / "conformers.json"),
                                 converged=sum(c["converged"] for r in result["records"] for c in r["conformers"]))
    if role == "Beta":
        return "BLOCKED", dict(reason="Requested original omegaB97X-D not wired; no reviewed QM partition",
                               pyscf_available=preflight()["packages"]["pyscf"], qm_executed=False, ts_executed=False)
    if role == "Gamma":
        return "SUCCEEDED", dict(evidence="GATE_AUDIT", healing_attempts=0,
                                 reason="dependency/method blockers are not SCF oscillation; no retries")
    if role == "Delta":
        plan = load_module("04_qmmm_openmm_fep").plan()
        write_json(attempt / "lambda_plan.json", plan)
        return "BLOCKED", dict(reason="No validated forcefield/topology/partition or PLUMED/REST2 sampling adapter",
                               plan=(attempt / "lambda_plan.json").relative_to(out).as_posix(), trajectories=0)
    if role == "Epsilon":
        return "BLOCKED", load_module("05_retrosynthesis_check").audit_route(None)
    if role == "Zeta":
        return "BLOCKED", dict(reason="No matched binary/conditional binding free energies",
                               alpha=None, ubiquitination_probability=None)
    return "SUCCEEDED", dict(evidence="AUDIT_ONLY", scientific_claim_accepted=False,
                             reason="No candidate passes scientific evidence gates")


def work(out, worker=None):
    out, spec = validate_resume(out)
    worker = worker or uuid.uuid4().hex
    ledger = Ledger(out / "ledger.sqlite", spec["token_cap"])
    try:
        while True:
            claim = ledger.claim(worker, lease_seconds=600)
            if claim is None:
                break
            attempt = out / "attempts" / f"{claim['role']}_{claim['fence']}_{worker}"
            attempt.mkdir(parents=True, exist_ok=False)
            try:
                status, result = handle(claim["role"], out, attempt, spec["mode"])
            except Exception as exc:
                status, result = "FAILED", dict(error_type=type(exc).__name__, message=str(exc))
            ledger.finish(claim, status, result)
            print(f"{claim['role']}: {status}", flush=True)
        return ledger.snapshot()
    finally:
        ledger.close()


def finalize(out):
    out, spec = validate_resume(out)
    ledger = Ledger(out / "ledger.sqlite", spec["token_cap"])
    try:
        snapshot = ledger.snapshot()
        if any(t["status"] not in {"SUCCEEDED", "BLOCKED", "FAILED"} for t in snapshot["tasks"]):
            raise ValueError("workers still active/pending; finalize only after all terminal")
        events = ledger.events()
        snapshot["event_chain_tip"] = ledger.verify_events()
    finally:
        ledger.close()
    write_json(out / "results/status.json", snapshot)
    (out / "reports/events.jsonl").write_text("\n".join(json.dumps(e, ensure_ascii=False) for e in events) + "\n", encoding="utf-8")
    write_json(out / "results/pareto_top_leads.json", dict(evidence="NO_QUALIFYING_EVIDENCE", candidates=[],
                                                        reason="No QM/MM-FEP, validated TS, or candidate synthesis evidence"))
    write_json(out / "results/artifact_status.json", dict(
        free_energy_landscape_3d="NOT_GENERATED_NO_UNBIASED_OR_REWEIGHTED_SAMPLING",
        quantum_states="EMPTY_NO_QM_EXECUTED", fep_trajectories="EMPTY_NO_SAMPLING_EXECUTED",
        llm_provider_tokens=snapshot["actual_provider_tokens"]))
    from render import render
    render(out, snapshot, events)
    lines = ["# Phase31 执行记录", "", "本次是本地确定性角色管线，不是已连接 LLM 的蜂群。",
             "参考体系：BRD4 BD2–VHL / macroPROTAC-1 (6SIS)。",
             f"预算上限：{spec['token_cap']:,}；本管线实际提供方 Token：{snapshot['actual_provider_tokens']}。",
             "该计数不包含编写本项目的 Codex 对话消耗。", "", "| 角色 | 状态 | 证据 |", "|---|---|---|"]
    for task in snapshot["tasks"]:
        result = json.loads(task["result"])
        lines.append(f"| {task['role']} | {task['status']} | {result.get('reason', result.get('evidence', ''))} |")
    lines += ["", "没有执行 DFT/CASSCF、QM/MM、FEP、REST2、CI-NEB 或 IRC；没有合格候选。",
              "MMFF94s 图只比较同一分子的收敛构象相对能量，既非绝对张力能，也非结合自由能。",
              "无真实采样，不生成三维自由能曲面；无 SCF 执行，自愈日志为空。",
              "SQLite 仅适用于同机共享磁盘；不能部署到 NFS 充当多节点队列。"]
    (out / "reports/run_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    files = {p.relative_to(out).as_posix(): sha256(p) for p in out.rglob("*")
             if p.is_file() and p.suffix not in {".sqlite", ".sqlite-wal", ".sqlite-shm"}
             and p.name != "artifact_hashes.json"}
    write_json(out / "reports/artifact_hashes.json", files)
    return snapshot
