"""Retrieve pinned original SI files, verify hashes and extract selected XYZ data.

python -m phase25_27.fetch_sources --download-dir /path/to/scratch --extract
Optional dependency for extraction: pypdf. Never republishes entire manuscripts.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import urllib.request
OUT=Path(__file__).resolve().parents[1]/'results_phase25_27'


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\n')

SOURCES=[
 ('ja800793t_si_001.pdf','10.1021/ja800793t','https://ndownloader.figshare.com/files/4627135','33eee138abfaa38d663af5f0aaf59e6598bdf6cf1b61504d495a8dd0011bd2c2',2928628),
 ('ic3c03822_si_001.pdf','10.1021/acs.inorgchem.3c03822','https://ndownloader.figshare.com/files/44610159','0eaa3aeebb04a9717dc089af2c4025df0e4b1390561346ae2cda7e3b1743763b',25250025),
 ('jp0c08646_si_001.pdf','10.1021/acs.jpca.0c08646','https://ndownloader.figshare.com/files/25563604','9ef02f5a916b4e8c841807cb032a5c4c849744ff6093745539b37ef65da517dd',13274896)]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--download-dir',type=Path,required=True);ap.add_argument('--extract',action='store_true');args=ap.parse_args()
    args.download_dir.mkdir(exist_ok=True,parents=True); records=[]
    for name,doi,url,sha,aid in SOURCES:
        path=args.download_dir/name
        if not path.exists():
            with urllib.request.urlopen(url,timeout=120) as r: data=r.read()
            if hashlib.sha256(data).hexdigest()!=sha: raise ValueError('Unrecognized source bytes; inspect version before use')
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=sha: raise ValueError(f'Corrupt/changed source {path}')
        records.append(dict(filename=name,doi=doi,url=url,sha256=sha,figshare_article_id=aid,
            source_license='CC BY-NC 4.0 according to Figshare metadata',retrieval_role='literature reference, not independent predictions'))
    dump(OUT/'source_downloads.json',records)
    if not args.extract:return
    from pypdf import PdfReader
    text='\n'.join(p.extract_text() for p in PdfReader(args.download_dir/'ja800793t_si_001.pdf').pages)
    text=text.split('E-imine:')[1].split('Buta')[0]
    rows=re.findall(r'^\s*(C|N|H)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$',text,re.M)
    if len(rows)!=28:raise ValueError('E-imine coordinate extraction failed; inspect original pages')
    out=OUT/'literature_structures';out.mkdir(exist_ok=True)
    records=[]
    def write(name,rows,doi,pages):
        xyz=str(len(rows))+'\nLiterature coordinates '+doi+' SI pages '+','.join(map(str,pages))+'; not recalculated\n'+'\n'.join(' '.join(row) for row in rows)+'\n'
        (out/(name+'.xyz')).write_bytes(xyz.encode('utf8'))
        records.append(dict(name=name,atoms=len(rows),source_doi=doi,SI_pages=pages,sha256=hashlib.sha256(xyz.encode()).hexdigest(),
          provenance='literature_geometry',license='CC BY-NC 4.0 per Figshare; retain author and DOI attribution',
          modification='selected Cartesian coordinates reformatted as XYZ; coordinates unchanged'))
    write('reference_E_imine',rows,'10.1021/ja800793t',[2,3])
    pdf=PdfReader(args.download_dir/'ic3c03822_si_001.pdf')
    for page,name,n in [(45,'N01_singlet',36),(46,'N01_triplet',36),(47,'N06_singlet',31),(48,'N06_triplet',31),(49,'N06_NiI_intermediate',30)]:
        rows=re.findall(r'\b(Ni|Cl|C|N|H)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',pdf.pages[page-1].extract_text())
        if len(rows)!=n:raise ValueError(f'Coordinate count mismatch {name}')
        write(name,rows,'10.1021/acs.inorgchem.3c03822',[page])
    dump(out/'manifest.json',records)


if __name__=='__main__':main()
