#!/usr/bin/env python3
"""Run organized sources in an isolated, compatible working layout (stdlib only)."""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'docs/layout_manifest.json'
NATIVE_PHASES = ROOT / 'docs/native_phases.json'


def phase_registry(manifest):
    """Keep the immutable 1-27 migration ledger separate from new native projects."""
    phases = dict(manifest['phases'])
    native = json.loads(NATIVE_PHASES.read_text(encoding='utf-8')) if NATIVE_PHASES.exists() else {}
    overlap = set(native) & set(phases)
    if overlap:
        raise ValueError(f'Native phase shadows legacy phase: {sorted(overlap)}')
    for key, phase in native.items():
        if not key.isdigit() or phase.get('execution') != 'native':
            raise ValueError(f'Invalid native phase: {key}')
        safe_path(ROOT, phase['entrypoint'])
    phases.update(native)
    return dict(sorted(phases.items(), key=lambda item: int(item[0])))


def digest(data):
    return hashlib.sha256(data).hexdigest()


def safe_path(root, relative):
    """Reject traversal/drive paths before any workspace write."""
    p = PurePosixPath(relative)
    if p.is_absolute() or '..' in p.parts or PureWindowsPath(relative).drive or '\\' in relative:
        raise ValueError(f'Unsafe mapped path: {relative}')
    result = (root / relative).resolve()
    if not result.is_relative_to(root.resolve()):
        raise ValueError(f'Path escapes root: {relative}')
    return result


def materialize(root, workspace, manifest, reuse=False):
    root, workspace = root.resolve(), workspace.resolve()
    if workspace == root or workspace in root.parents:
        raise ValueError('Workspace must not be the repository or one of its parents')
    if workspace.is_relative_to(root) and not workspace.is_relative_to(root / 'work'):
        raise ValueError('Use a new directory under work/ or outside this repository')
    fingerprint = digest(json.dumps(manifest, sort_keys=True).encode())
    marker = workspace / '.ai4chem-workspace.json'
    source_hashes = {}
    for record in manifest['files']:
        source_hashes[record['old']] = digest(safe_path(root, record['new']).read_bytes())
        safe_path(workspace, record['old'])
    if workspace.exists():
        if not reuse or not marker.is_file():
            raise ValueError('Workspace already exists; choose a new path (or --reuse for a prepared workspace)')
        previous = json.loads(marker.read_text(encoding='utf-8'))
        if previous['manifest_sha256'] != fingerprint or previous['source_hashes'] != source_hashes:
            raise ValueError('Organized sources changed; prepare a new workspace')
        for path, checksum in previous['copied_code_hashes'].items():
            if digest(safe_path(workspace, path).read_bytes()) != checksum:
                raise ValueError(f'Workspace source changed: {path}; use a new workspace')
        return workspace
    original_docs = json.loads(gzip.decompress((root / 'shared/original_document_bytes.json.gz').read_bytes()))
    workspace.mkdir(parents=True, exist_ok=False)
    copied_code = {}
    for record in manifest['files']:
        data = safe_path(root, record['new']).read_bytes()
        # Only undo link-only migration edits. Later user edits are copied as-is.
        if record['old'] in original_docs and digest(data) == record['organized_sha256']:
            data = base64.b64decode(original_docs[record['old']])
        target = safe_path(workspace, record['old'])
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)  # Independent copies: scientific output cannot mutate tracked data.
        if target.suffix == '.py':
            copied_code[record['old']] = digest(data)
    marker.write_text(json.dumps({'manifest_sha256': fingerprint, 'source_hashes': source_hashes,
                                 'copied_code_hashes': copied_code}, indent=2)+'\n', encoding='utf-8')
    return workspace


def main(argv=None):
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    argv = list(sys.argv[1:] if argv is None else argv)
    if '--' in argv:
        split = argv.index('--'); arguments = argv[split+1:]; argv = argv[:split]
    else:
        arguments = []
    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    phases = phase_registry(manifest)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', nargs='?', type=int, choices=[int(key) for key in phases])
    parser.add_argument('--list', action='store_true', help='List phase folders without importing scientific packages')
    parser.add_argument('--workspace', type=Path, help='New isolated working directory; never the source tree')
    parser.add_argument('--prepare-only', action='store_true', help='Copy sources/data without running a calculation')
    parser.add_argument('--reuse', action='store_true', help='Explicitly reuse a workspace with matching sources')
    parser.add_argument('--python', default=sys.executable, help='Existing scientific interpreter')
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--module', help='Existing phase25_27 module, e.g. phase25_27.test_analysis')
    selection.add_argument('--script', help='Original script name from layout_manifest.json, for tests/helpers')
    args = parser.parse_args(argv)
    if args.list:
        for n, phase in phases.items():
            print(f'{int(n):02d}  {phase["folder"]}  {phase["title"]}')
        return 0
    if args.phase is None:
        parser.error('phase is required; scientific arguments go after --')
    phase = phases[str(args.phase)]
    if phase.get('execution') == 'native':
        if args.module or args.script or args.prepare_only or args.reuse:
            parser.error('Native phases run in their project layout; use arguments after -- and optional --workspace for --out')
        if args.workspace and any(arg.split('=', 1)[0] in ('--out', '--output-dir') for arg in arguments):
            parser.error('Specify only --workspace or native --out, not both')
        entry = safe_path(ROOT, phase['entrypoint'])
        if not entry.is_file():
            parser.error(f'Native entrypoint is missing: {entry}')
        command = [args.python, str(entry)]
        if args.workspace:
            command += ['--out', str(args.workspace.resolve())]
        command += arguments
        env = os.environ.copy()
        env.setdefault('PYTHONUTF8', '1')
        return subprocess.call(command, cwd=ROOT, env=env)
    if args.workspace is None:
        parser.error('phase and --workspace are required; scientific arguments go after --')
    known = {r['old'] for r in manifest['files']}
    if args.module:
        if not args.module.startswith('phase25_27.') or args.module.replace('.', '/')+'.py' not in known:
            parser.error('module must name an existing phase25_27 module')
        command = [args.python, '-m', args.module, *arguments]
    else:
        entry = args.script or manifest['phases'][str(args.phase)]['entrypoint']
        if entry is None and not args.prepare_only:
            parser.error('Phases 25-27 have no complete production driver; use --prepare-only or an explicit --module')
        if entry is not None and (entry not in known or not entry.endswith('.py')):
            parser.error('script must be a mapped Python source')
        command = [args.python, entry, *arguments] if entry else None
    try:
        workspace = materialize(ROOT, args.workspace, manifest, args.reuse)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f'Workspace: {workspace}', flush=True)
    if args.prepare_only:
        print('Prepared only; no calculation submitted.')
        return 0
    env = os.environ.copy()
    env['PYTHONPATH'] = str(workspace) + (os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    env.setdefault('PYTHONUTF8', '1')
    return subprocess.call(command, cwd=workspace, env=env)


if __name__ == '__main__':
    raise SystemExit(main())
