"""Review installed resources and saved actual MCP responses for Phases 25-27.

Run with the chem-ai4s extension interpreter. Read-only against its database.
No benchmark labels enter the prospective catalyst blind test.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import hashlib,json
import numpy as np
import duckdb
from rdkit import Chem
from ord_schema.proto import reaction_pb2
from pydantic_core import to_json
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results_phase25_27/tool_capability_audit'
WORKBENCH=Path(r'C:\Users\HUIWEI\.codex\tools\chem-ai4s')
def read(name):return json.loads((OUT/(name+'.json')).read_text(encoding='utf8'))
def write(name,data):(OUT/name).write_bytes((json.dumps(data,indent=2,ensure_ascii=False)+'\n').encode('utf8'))

def main():
    with duckdb.connect(str(WORKBENCH/'data/research.duckdb'),read_only=True) as db:
        rows=db.execute('SELECT reaction_id,reaction FROM ord_example').fetchall()
    types=Counter(); measurements=Counter(); selectivities=Counter(); failed=[]
    for key,raw in rows:
        try:
            r=reaction_pb2.Reaction();r.ParseFromString(raw)
            assert r.reaction_id==key
            for ident in r.identifiers:
                if reaction_pb2.ReactionIdentifier.ReactionIdentifierType.Name(ident.type)=='REACTION_TYPE':types[ident.value]+=1
            for outcome in r.outcomes:
                for product in outcome.products:
                    for measurement in product.measurements:
                        measurements[reaction_pb2.ProductMeasurement.ProductMeasurementType.Name(measurement.type)]+=1
                        if measurement.HasField('selectivity'):
                            field=measurement.selectivity.DESCRIPTOR.fields_by_name['type']
                            selectivities[field.enum_type.values_by_number[measurement.selectivity.type].name]+=1
        except Exception as exc:failed.append(dict(reaction_id=key,error=str(exc)))
    try:to_json({'reaction':rows[0][1]});serialization_error=None
    except Exception as exc:serialization_error=type(exc).__name__+': '+str(exc)
    ord_audit=dict(rows=len(rows),parsed=len(rows)-len(failed),failed=failed,reaction_types=dict(types),measurement_types=dict(measurements),selectivity_types=dict(selectivities),
       interface_diagnostic=serialization_error,workaround='Read local DuckDB BLOB and parse ord_schema Reaction protobuf; aggregate only.',blind_test_labels_used=0,
       relevance_note='Reaction type counts are a scope audit, not proof that a target literature dataset is absent from all of ORD.')
    write('ord_audit.json',ord_audit)
    inputs=json.loads((ROOT/'results_phase25_27/phase25/donor.json').read_text(encoding='utf8'))
    a=Chem.MolFromSmiles(inputs['mapped_smiles'])
    for atom in a.GetAtoms():atom.SetAtomMapNum(0)
    local=Chem.MolToSmiles(Chem.RemoveHs(a))
    pub=read('pubchem_structure')['data']['PropertyTable']['Properties'][0]
    remote=Chem.MolToSmiles(Chem.MolFromSmiles(pub['SMILES']))
    assert local==remote
    force=read('force_check');zero=force[0]['response'];minus=force[1]['response'];plus=force[2]['response']
    fd=-(plus['energy_eV']-minus['energy_eV'])/.0002
    analytic=zero['forces_eV_per_angstrom'][0][0]
    force_error=abs(fd-analytic);assert force_error<1e-5
    cu=read('copper_mp30')['data']['data'][0]['attributes']
    structure=Structure(cu['lattice_vectors'],cu['species_at_sites'],cu['cartesian_site_positions'],coords_are_cartesian=True)
    symmetry=SpacegroupAnalyzer(structure,symprec=.01)
    assert symmetry.get_space_group_number()==225
    model_catalog=json.loads((WORKBENCH/'models.json').read_text(encoding='utf8'))
    models=[]
    for name,item in model_catalog.items():
        p=Path(item['path']);digest=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
        assert digest==item['sha256']
        models.append(dict(name=name,bytes=p.stat().st_size,sha256=digest,license=item.get('license'),source=item.get('source'),
                           target_test_this_audit='Cu energy/force consistency' if name=='mace-mh-1.model' else 'file integrity only; no target-domain validation'))
    basis=json.loads(read('basis_Ni')['basis'])
    assert set(basis['elements'])=={'1','6','7','17','28'}
    summary=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),local_dataset_entries=len(read('status')['datasets']),
       pubchem_CID=pub['CID'],donor_connectivity_matches=True,Cu_database_id='mp-30',Cu_space_group=symmetry.get_space_group_symbol(),symmetry_tolerance_A=.01,
       copper_model_energy_eV=read('copper_model')['result']['energy_eV'],copper_model_atoms=4,Cu_force_error_eV_A=force_error,
       model_tests='One fixed-cell bulk configuration and one displaced configuration, finite-difference h=1e-4 Angstrom; not DFT accuracy validation.',
       basis_elements=sorted(basis['elements'],key=int),models=models,ord_rows_parsed=len(rows)-len(failed),production_acceptance_changed=False,
       constant_potential_engine_newly_available=False,multireference_SOC_NAC_engine_newly_available=False)
    write('summary.json',summary)
    counts='; '.join(f'{k}: {v}' for k,v in types.items())
    for lang in ['EN','ZH']:
        if lang=='EN':
            text=f'''# New resource applicability audit

Snapshot: {summary['timestamp_utc']}. These checks establish usable interfaces and model consistency, not completion of Phases25-27.

| Resource | Actual check | Use and remaining boundary |
|---|---|---|
|MACE-MH-1, omat_pbe|4-atom Cu energy {summary['copper_model_energy_eV']:.9f} eV; displaced-atom force error {force_error:.3e} eV/Angstrom|Candidate for Cu geometry preparation. Fixed-charge learned potential; no grand-canonical electron number, potential control or activation free energy. Symmetry-zero bulk forces do not prove equilibrium lattice constant.|
|PubChem|CID {pub['CID']} matches the mapped dimethyl Hantzsch donor after map removal|Structure identity lookup; no reaction-condition or ee validation. Name query returned 404; structure query succeeded.|
|Materials Project OPTIMADE|mp-30 returned; local symmetry Fm-3m, tolerance 0.01 Angstrom|Reference Cu structure. JARVIS query returned HTTP500. Query results are not a local full database.|
|Basis Set Exchange|def2-TZVP retrieved for Ni, Cl, C, N, H with references|Prepare Phase27 basis-convergence jobs; basis data are not a multireference/SOC/NAC engine.|
|ORD|{len(rows)-len(failed)}/{len(rows)} local protobuf reactions parsed|Reaction provenance and pipeline development. MCP JSON serialization of BLOBs fails; local read-only protobuf decoding works. No records used as blind-test labels.|
|MACE-POLAR and OrbMol-v2|Local checkpoint sizes and SHA256 verified|Candidates for molecular prescreening after task-specific DFT validation. No excited-state branching capability established. MCP molecular endpoint is capped at 100 atoms; use local Python for the complete 117-atom system.|
|QM9, MoleculeNet, Matbench|15 catalog entries total, including PDB and ORD samples|Useful for software/descriptor baselines. They do not supply the specified catalyst ee, electrode PMFs or Ni photodynamic labels.|

ORD reaction-type scope: {counts}. Product measurement counts: {dict(measurements)}; selectivity subtypes: {dict(selectivities)}. This is a local sample audit, not a search of the entire ORD archive.

Four model checkpoints passed local file-integrity checks. MACE weight ASL terms are retained; OrbMol licensing is recorded in its installation manifest. PySCF/GPU4PySCF/periodic DFT engines remain uninstalled according to the current provisioning catalog; UMA remains access-required/uninstalled. New model files add useful preprocessing options but do not add a callable OpenAI GPU allocation.

Next use: preserve the running Psi4 frequency queue; validate molecular surrogate forces/relative energies against suitable DFT before broad conformer prescreening; validate the Cu model on slab/adsorbate configurations before geometry preparation. Constant-potential solvent sampling and Ni multireference dynamics remain separate required calculations. Scientific acceptance is unchanged.

[Raw checks](summary.json) · [ORD scope](ord_audit.json) · [PubChem source response](pubchem_structure.json) · [Cu source response](copper_mp30.json) · [Basis with references](basis_Ni.json) · [Force check inputs/results](force_check.json)
'''
        else:
            text=f'''# 新工具与数据库适用性核查

快照：{summary['timestamp_utc']}。本轮确认接口可用性与模型力的一致性，不改变 Phase25–27 的科学验收状态。

|资源|本轮实际核查|可用范围与限制|
|---|---|---|
|MACE-MH-1，omat_pbe|4 原子 Cu 能量 {summary['copper_model_energy_eV']:.9f} eV；位移结构解析力与有限差分力误差 {force_error:.3e} eV/Å|可作为铜几何预处理候选；不是恒电位电子结构引擎。对称体相零力不能证明晶格常数已平衡。|
|PubChem|CID {pub['CID']} 与现有映射结构去除映射后的二甲酯型 Hantzsch 酯一致|结构身份核对；不验证反应条件或 ee。名称查询 404 后，结构查询成功。|
|Materials Project OPTIMADE|mp-30 返回成功；本地识别 Fm-3m，对称容差 0.01 Å|铜参考结构；JARVIS 本次返回 HTTP500。联网可查询不等于完整数据库已下载。|
|Basis Set Exchange|取得 Ni、Cl、C、N、H 的 def2-TZVP 及引用|可用于 Phase27 基组收敛输入；它本身不执行多参考、SOC 或 NAC 计算。|
|ORD|本地 {len(rows)-len(failed)}/{len(rows)} 条 protobuf 记录成功解析|用于反应来源与数据管线；MCP 直接序列化 BLOB 失败，可用本地只读解析替代。没有用这些记录充当盲测标签。|
|MACE-POLAR、OrbMol-v2|模型文件大小与 SHA256 核验通过|可尝试分子预筛，须先做目标体系 DFT 验证；未建立激发态分支能力。MCP 分子接口限制 100 原子，完整 117 原子体系应从本地 Python 调用。|
|QM9、MoleculeNet、Matbench|目录合计 15 项，含 PDB 与 ORD 示例|可用于描述符或软件基线；不提供本题指定 ee、恒电位势垒或 Ni 轨迹标签。|

ORD 本地样本的反应类型：{counts}。产物测量类型计数：{dict(measurements)}；选择性子类型：{dict(selectivities)}。这是样本范围核查，不代表全 ORD 检索。

4 个本地模型文件的哈希均一致。保留 MACE 权重的 ASL 使用条件；OrbMol 许可按安装清单记录。当前配置目录仍标记 PySCF/GPU4PySCF/周期 DFT 引擎未安装，UMA 仍需访问授权且未安装。工具增多确实改善结构准备和数据获取，但没有新增可调用的 OpenAI GPU 配额。

接入顺序：继续现有 Psi4 频率队列；先以合适 DFT 参考验证分子模型的力和相对能量，再扩大构象预筛；铜模型先做表面和吸附物独立检查，再用于几何预处理。恒电位显式溶剂采样与 Ni 多参考动力学仍需各自的实际计算，不能由上述模型检查替代。

[核查汇总](summary.json) · [ORD 范围](ord_audit.json) · [PubChem 原始响应](pubchem_structure.json) · [铜结构原始响应](copper_mp30.json) · [基组与引用](basis_Ni.json) · [力检查输入输出](force_check.json)
'''
        (OUT/f'REPORT_{lang}.md').write_bytes(text.encode('utf8'))
    print(json.dumps(dict(parsed_ORD=len(rows)-len(failed),reaction_types=dict(types),measurement_types=dict(measurements),selectivity_types=dict(selectivities),force_error_eV_A=force_error,CID=pub['CID'],spacegroup=summary['Cu_space_group']),ensure_ascii=False))

if __name__=='__main__':main()
