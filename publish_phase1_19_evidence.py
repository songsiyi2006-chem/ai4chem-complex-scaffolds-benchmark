"""Publish selected rerun evidence, without replacing historical phase outputs.

Large binary trajectories/checkpoints are omitted, not claimed as included.
Text workstation paths are redacted. Source and published hashes are distinct.
Only explicit scientific result directories from this campaign are considered.
"""
import hashlib
import json
import math
from pathlib import Path

from rerun_phase1_19_campaign import ROOT, SCRIPTS
from summarize_phase1_19_reruns import summarize

CAMPAIGN = ROOT.parent / 'phase1-19-rerun-20260910'
OUT = ROOT / 'audit' / 'phase1_19_rerun_20260910'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def redact(value):
    # Longest prefixes first, including JSON-decoded Windows and POSIX spelling.
    for path, tag in ((CAMPAIGN, '$CAMPAIGN'), (ROOT, '$SOURCE_ROOT'),
                      (ROOT.parent, '$WORK_ROOT'), (Path('C:/Users/HUIWEI'), '$USER_PROFILE')):
        for form in (str(path).replace('\\', '/'), str(path).replace('/', '\\')):
            value = value.replace(form, tag)
    return value


def portable(value):
    if isinstance(value, dict):
        return {k: portable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [portable(v) for v in value]
    if isinstance(value, str):
        return redact(value)
    if isinstance(value, float) and not math.isfinite(value):
        # Explicit unavailable value in a public JSON, not a repaired measurement.
        return None
    return value


def main():
    record = summarize(CAMPAIGN)
    OUT.mkdir(parents=True, exist_ok=True)
    previous = json.loads((OUT/'manifest.json').read_text(encoding='utf-8')) if (OUT/'manifest.json').exists() else {}
    published, omitted = [], []

    def export(src, relative, role):
        raw = src.read_bytes()
        dst = OUT / relative
        dst.parent.mkdir(parents=True, exist_ok=True)
        transformed = False
        if src.suffix == '.json':
            data = (json.dumps(portable(json.loads(raw.decode('utf-8-sig'))),
                               indent=2, ensure_ascii=False, allow_nan=False)+'\n').encode('utf-8')
            transformed = data != raw
        elif src.suffix in ('.log', '.md'):
            data = redact(raw.decode('utf-8', errors='replace')).encode('utf-8')
            transformed = data != raw
        else:
            data = raw  # Executed Python snapshots and binary arrays stay byte-exact.
        dst.write_bytes(data)
        published.append(dict(path=str(relative).replace('\\', '/'), role=role,
                              source=str(src.resolve()), source_sha256=digest(raw),
                              published_sha256=digest(data), bytes=len(data),
                              text_normalized_or_redacted=transformed))

    for row in record['phases']:
        if not row['attempts']:
            continue
        attempt = row['attempts'][-1]
        phase = row['phase']
        folder = Path(attempt['attempt_path'])
        prefix = Path(f'phase{phase:02d}') / folder.name
        # Status is a point-in-time record, even when outputs are still running.
        export(folder/'rerun_status.json', prefix/'run_status.json', 'process snapshot')
        if attempt.get('continuation'):
            export(folder/'continuation_status.json', prefix/'continuation_status.json', 'continuation provenance')
        terminal = attempt['status'] != 'running'
        if attempt.get('continuation', {}).get('status') == 'running':
            terminal = False
        if not terminal:
            continue  # Never package incomplete mutable results as final output.
        candidates = {folder/'run.log', folder/SCRIPTS[phase], folder/'acceptance_review.json'}
        if phase == 1:
            candidates |= {folder/'molecules.py', folder/'bench_results/benchmark_results.json'}
        if phase == 2:
            candidates.add(folder/'export_openff_system.py')
        for sub in folder.glob('results*'):
            if sub.is_dir():
                candidates.update(sub.glob('*.json'))
                candidates.update(sub.glob('*.npz'))
        for sub in folder.glob('figures*'):
            if sub.is_dir():
                candidates.update(sub.glob('*.png'))
        candidates.update(folder.glob('continuation*.log'))
        for recovery in folder.glob('figure_recovery_*'):
            rec_path = recovery/'recovery_status.json'
            if rec_path.exists() and json.loads(rec_path.read_text(encoding='utf-8')).get('status') == 'completed':
                candidates.add(rec_path)
                candidates.update(recovery.glob('*.py'))
                candidates.update(recovery.glob('*.log'))
                candidates.update(recovery.glob('*.md'))
                candidates.update((recovery/'figures_phase5').glob('*.png'))
        for audit in folder.glob('*_postaudit_*'):
            if audit.is_dir():
                candidates.update(audit.glob('*.json'))
                candidates.update(audit.glob('verified_source.py'))
        for src in sorted(candidates):
            if not src.is_file():
                continue
            if src.stat().st_size > 5 * 1024 * 1024:
                omitted.append(dict(source=str(src), bytes=src.stat().st_size,
                                    sha256=digest(src.read_bytes()), reason='over 5 MiB selection limit'))
                continue
            export(src, prefix/src.relative_to(folder), 'selected evidence, not full trajectory archive')

    latest = sorted((ROOT/'results_rerun_audit').glob('*/regression_status.json'))
    if latest:
        for src in sorted(latest[-1].parent.iterdir()):
            if src.is_file() and src.suffix in ('.log', '.json'):
                export(src, Path('regressions')/latest[-1].parent.name/src.name, 'regression evidence')
    current_paths={p['path'] for p in published}
    for item in previous.get('published_files', []):
        if item['path'] not in current_paths:
            retained=dict(item, retained_from_publication=previous.get('generated'))
            published.append(retained)
    record.update(published_files=published, omitted_large_files=omitted,
                  publication_scope='Selected latest terminal attempt evidence plus running status snapshots; not a complete archive',
                  text_policy='Workstation paths redacted; JSON reserialized; nonfinite JSON numbers represented as null (unavailable), never recomputed. Source and published SHA256 recorded separately.',
                  source_scope='Only selected executed entry points are bundled. All original snapshot file hashes are retained; this is not a complete environment/dependency archive.')
    manifest = portable(record)
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    lines = ['# Phase 1–19 fresh rerun evidence', '', f"Snapshot: {record['generated']}", '',
             'This is an ongoing computational audit, not an all-phases-pass or experimental-validation claim.',
             '本目录是持续重算的阶段性证据，不能解读为 19 个阶段全部通过。历史根目录结果没有被自动覆盖。', '',
             'See [status](../../PHASE1_19_RERUN_STATUS.md) and [manifest](manifest.json).', '',
             'The manifest records source/published hashes, failed attempts, retained limitations and omitted large files.',
             'Selected executed scripts are included; dependencies, all historical source snapshots and large trajectories are not fully bundled.',
             'Text paths are redacted. JSON null may denote an originally nonfinite/unavailable value; it is not a replacement measurement.', '',
             '## Selected phase records', '']
    for row in manifest['phases']:
        if row['attempts']:
            a=row['attempts'][-1]
            folder=Path(a['attempt_path'].replace('\\','/')).name
            lines.append(f"- [Phase {row['phase']}](phase{row['phase']:02d}/{folder}/run_status.json): {a['status']}; {a.get('scientific_acceptance', 'not reviewed')}")
    (OUT/'README.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    for item in published:
        if digest((OUT/item['path']).read_bytes()) != item['published_sha256']:
            raise RuntimeError('Published artifact hash mismatch: '+item['path'])
    print(json.dumps(dict(published=len(published), bytes=sum(i['bytes'] for i in published),
                          omitted_large_files=len(omitted), manifest=str(OUT/'manifest.json')), indent=2))


if __name__ == '__main__':
    main()
