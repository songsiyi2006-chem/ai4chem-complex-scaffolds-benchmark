"""Validate physical migration, byte preservation, phase coverage and Markdown paths."""
import ast
import base64
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
from run_phase import phase_registry

ROOT = Path(__file__).resolve().parents[1]


def markdown_links(text):
    text = re.sub(r'```.*?```|~~~.*?~~~', '', text, flags=re.S)
    text = re.sub(r'`+[^`]*`+', '', text)
    return re.findall(r'\]\(([^)\s]+)\)', text)


def main():
    manifest = json.loads((ROOT/'docs/layout_manifest.json').read_text(encoding='utf-8'))
    phases = phase_registry(manifest)
    originals = json.loads(gzip.decompress((ROOT/'shared/original_document_bytes.json.gz').read_bytes()))
    assert set(manifest['phases']) == {str(i) for i in range(1,28)}
    assert len({r['old'] for r in manifest['files']}) == len(manifest['files'])
    assert len({r['new'] for r in manifest['files']}) == len(manifest['files'])
    preserved = 0; python_files = 0
    for record in manifest['files']:
        path = ROOT/record['new']
        assert path.is_file(), f'Missing migrated file: {path}'
        if path.suffix != '.md' and record['old'] not in ('.gitattributes','.gitignore'):
            assert hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256'], f'Changed scientific bytes: {path}'
            preserved += 1
        if record['old'] in originals:
            assert hashlib.sha256(base64.b64decode(originals[record['old']])).hexdigest() == record['sha256']
        if path.suffix == '.py':
            ast.parse(path.read_bytes(), filename=str(path)); python_files += 1
        if record['old'] != record['new']:
            assert not (ROOT/record['old']).is_file(), f'Old duplicate remains: {record["old"]}'
    for phase in phases.values():
        folder=ROOT/phase['folder']
        for relative in ('README.md','run.py','reports','code'):
            assert (folder/relative).exists(), f'Phase entry missing: {folder/relative}'
    broken=[]; links=0; documents=0
    for base in ('README.md','docs','projects','shared'):
        path=ROOT/base
        paths=[path] if path.is_file() else path.rglob('*.md')
        for doc in paths:
            documents+=1
            for url in markdown_links(doc.read_text(encoding='utf-8')):
                if re.match(r'(?:[a-zA-Z][\w+.-]*:|#|/)',url):continue
                target=url.split('#')[0]
                if not (doc.parent/target).exists():
                    broken.append({'document':doc.relative_to(ROOT).as_posix(),'target':url})
                links+=1
    report={'mapped_files':len(manifest['files']), 'moved_files':sum(r['old']!=r['new'] for r in manifest['files']),
            'phases':len(phases),'unchanged_non_markdown_files':preserved,'python_sources_parsed':python_files,
            'markdown_documents':documents,'local_markdown_links':links,'broken_links':broken}
    print(json.dumps(report,ensure_ascii=False,indent=2))
    return bool(broken)


if __name__=='__main__':
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    raise SystemExit(main())
