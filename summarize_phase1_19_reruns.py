"""Inventory actual fresh runs without interpreting exit zero as scientific success."""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parent

def summarize(campaign):
    rows=[]
    for phase in range(1,20):
        attempts=[]
        for status in sorted((campaign/f'phase{phase:02d}').glob('*/rerun_status.json')):
            rec=json.loads(status.read_text(encoding='utf-8'))
            if rec['status']=='process_failed' and rec.get('scientific_acceptance')=='pending output review':
                rec['scientific_acceptance']='not accepted: process failed or was deliberately stopped; inspect log/notes'
            if rec['status']=='running':
                started=time.mktime(time.strptime(rec['started'],'%Y%m%dT%H%M%S'))
                rec['elapsed_s']=max(0,time.time()-started)
            logs=status.parent/'run.log'
            rec['attempt_path']=str(status.parent.resolve())
            changed=[name for name,digest in rec.get('source_sha256',{}).items()
                if not (status.parent/name).is_file() or hashlib.sha256((status.parent/name).read_bytes()).hexdigest()!=digest]
            rec['regenerated_python_artifacts']=[name for name in changed if name.startswith('output_')]
            rec['snapshot_hash_mismatches']=[name for name in changed if not name.startswith('output_')]
            rec['log_sha256']=hashlib.sha256(logs.read_bytes()).hexdigest() if logs.exists() else None
            if rec['status']=='running' and logs.exists():
                tail=logs.read_text(encoding='utf-8',errors='replace').splitlines()[-5:]
                rec['memory_hold_reported']=bool(tail and '[memory] HOLD' in tail[-1])
            # Master result records only, not private file contents / binary data.
            rec['result_files']=[str(p.relative_to(status.parent)) for p in status.parent.glob('results*/*results*.json')]
            review=status.parent/'acceptance_review.json'
            if review.exists():
                rec['acceptance_review']=json.loads(review.read_text(encoding='utf-8'))
                if rec['acceptance_review'].get('status')=='specified_checks_passed_with_limitations':
                    rec['scientific_acceptance']='specified numerical/software checks passed; see limitations'
            if phase==19:
                audit_path=status.parent/'results_phase19/phase19_acceptance_audit.json'
                raw_path=status.parent/'results_phase19/phase19_results.json'
                if audit_path.exists() and raw_path.exists():
                    data=json.loads(audit_path.read_text(encoding='utf-8'))
                    if data.get('source_sha256')==hashlib.sha256(raw_path.read_bytes()).hexdigest():
                        rec['candidate_acceptance_review']=data
                        if data.get('candidates') and not any(c['accepted'] for c in data['candidates']):
                            rec['scientific_acceptance']='no_accepted_design_under_audited_gates'
            continuation=status.parent/'continuation_status.json'
            if continuation.exists():
                rec['continuation']=json.loads(continuation.read_text(encoding='utf-8'))
            audits=[]
            for audit_path in sorted(status.parent.glob('*_postaudit_*/audit.json')):
                audit=json.loads(audit_path.read_text(encoding='utf-8'))
                audits.append(dict(path=str(audit_path.resolve()),record=audit,
                                   sha256=hashlib.sha256(audit_path.read_bytes()).hexdigest()))
            if audits: rec['post_audits']=audits
            if phase==1:
                p=status.parent/'bench_results/benchmark_results.json'
                if p.exists():
                    data=json.loads(p.read_text(encoding='utf-8'))['results']
                    rec['registry_status']=[dict(id=r['id'],status=r['status'],
                        conformers=r.get('conformers',{}).get('status'),
                        accepted=r.get('conformers',{}).get('n_accepted')) for r in data]
            attempts.append(rec)
        rows.append(dict(phase=phase,attempts=attempts))
    return dict(generated=time.strftime('%Y-%m-%d %H:%M:%S'),campaign=str(campaign.resolve()),phases=rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--campaign',type=Path,default=ROOT.parent/'phase1-19-rerun-20260910')
    args=ap.parse_args()
    record=summarize(args.campaign)
    (ROOT/'PHASE1_19_RERUN_MANIFEST.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    lines=['# Phase 1–19 全新计算核验记录','',f"更新：{record['generated']}",'',
        '本记录来自独立运行目录，不复用历史结果。进程正常退出不等于科学结论通过；未完成、失败和替代模型必须单列。',
        '所有结果均为计算/模拟，不是实验验证。未执行硬件或自动关机；Git 发布状态另以实际提交记录为准。','',
        '墙钟耗时包含中断或休眠，不能直接当作有效计算耗时或用来比较性能。','',
        '| Phase | 尝试数 | 最新进程状态 | 范围 | 墙钟秒 | 结果验收 |',
        '|---|---:|---|---|---:|---|']
    for row in record['phases']:
        a=row['attempts']
        if not a:
            lines.append(f"| {row['phase']} | 0 | 尚未启动 | — | — | 未验收 |")
            continue
        r=a[-1]
        state={'running':'运行中','process_completed':'进程完成','process_failed':'失败或主动中止'}.get(r['status'],r['status'])
        if r.get('memory_hold_reported'): state='等待内存余量'
        scope='默认完整参数' if r['scope'].startswith('default/full') else '自定义启动参数，见命令'
        if r.get('fresh_atomic_parent'): scope='本轮原子检查点＋后续重算'
        elapsed=r.get('elapsed_s',0)
        verdict=r.get('scientific_acceptance','未验收')
        if verdict.startswith('specified numerical'): verdict='所列检查通过，见局限'
        elif verdict=='pending output review': verdict='待完成与复核'
        elif verdict=='acceptance checks FAILED': verdict='验收未通过'
        elif verdict.startswith('not accepted:'): verdict='未接受，见日志与修复记录'
        if 'continuation' in r:
            cont=r['continuation']
            state='本轮检查点续算中' if cont['status']=='running' else '分段续算已结束，需复核'
            scope='同批新算数据分段接续'
            verdict='待续算完成与复核' if cont['status']=='running' else cont.get('scientific_acceptance','未验收')
            elapsed+=cont.get('elapsed_s',0)
        if verdict=='sampling_incomplete_for_reaction_barriers': verdict='默认轨迹完成，产物采样不足，势垒未验收'
        if verdict=='no_accepted_design_under_audited_gates': verdict='完整流程结束，候选全部拒绝，无合格设计'
        if any(a['record'].get('status')=='running' for a in r.get('post_audits',[])):
            state='后处理复核中'
            verdict='审计进行中，尚未完成验收'
        lines.append(f"| {row['phase']} | {len(a)} | {state} | {scope} | {elapsed:.1f} | {verdict} |")
    lines+=['','## 证据与复现','',
        '- [公开清单](audit/phase1_19_rerun_20260910/manifest.json)：命令、源文件与选定结果的 SHA256、进程状态和局限。完整本地记录为 `PHASE1_19_RERUN_MANIFEST.json`，不包含在公开文件集合中。',
        '- `rerun_phase1_19_campaign.py`：串行启动器。阶段 6/7 使用 `--phase4-attempt` 指定本轮新算的几何，并记录输入哈希。',
        '- `PHASE*_RERUN_NOTES.md`：各阶段修复内容、测试覆盖及限制。',
        '- 不能将默认参数等同于研究级充分采样；具体轨迹长度、训练步数、收敛指标以各阶段结果为准。','']
    (ROOT/'PHASE1_19_RERUN_STATUS.md').write_text('\n'.join(lines),encoding='utf-8')
    print('\n'.join(lines[:28]))

if __name__=='__main__': main()
