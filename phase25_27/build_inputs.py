"""Generate explicitly labelled preliminary structures and preregistered matrices.

Run from repository root: python -m phase25_27.build_inputs
These are starting geometries, never DFT minima or transition states.
"""
from pathlib import Path
import csv
import hashlib
import itertools
import json
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdMolTransforms

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results_phase25_27'
SEEDS = [104729, 130363, 155921, 181081, 205019, 230003, 254959, 280001, 305093, 330017]


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf8')


def table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf8') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def mapped(mol, offset=0):
    mol=Chem.Mol(mol)
    for atom in mol.GetAtoms(): atom.SetAtomMapNum(offset+atom.GetIdx()+1)
    return mol


def preliminary(mol, name, starts=12, axis=None):
    map_offset=((min(a.GetAtomMapNum() for a in mol.GetAtoms())-1)//200)*200
    mol=Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        if not atom.GetAtomMapNum(): atom.SetAtomMapNum(1000+map_offset+atom.GetIdx())
    params=AllChem.ETKDGv3(); params.randomSeed=SEEDS[0]; params.numThreads=2
    params.pruneRmsThresh=-1.0  # keep search multiplicity separately from degeneracy
    ids=list(AllChem.EmbedMultipleConfs(mol, numConfs=starts, params=params))
    if not ids: raise RuntimeError(f'No embeddings: {name}')
    if not AllChem.MMFFHasAllMoleculeParams(mol): raise RuntimeError(f'MMFF parameters absent: {name}')
    opt=AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=2, maxIters=2000)
    records=[]; path=OUT/'phase25'/'starting_structures'/f'{name}.sdf'; path.parent.mkdir(parents=True,exist_ok=True)
    with Chem.SDWriter(str(path)) as writer:
        writer.SetForceV3000(True)
        for cid,(flag,energy) in zip(ids,opt):
            one=Chem.Mol(mol); conf=Chem.Conformer(mol.GetConformer(cid)); one.RemoveAllConformers(); one.AddConformer(conf)
            one.SetProp('_Name',f'{name}_start{cid:03d}')
            one.SetProp('evidence_kind','force_field_starting_geometry')
            one.SetProp('MMFF_status',str(flag)); one.SetProp('MMFF_energy_kcal_mol',str(energy))
            row=dict(start=cid,converged=flag==0,MMFF_energy_kcal_mol=float(energy),physical_degeneracy=None,thermodynamic_weight=None)
            if axis:
                angle=rdMolTransforms.GetDihedralDeg(one.GetConformer(),*axis)
                row['axis_dihedral_deg']=angle; one.SetProp('axis_dihedral_atom_maps',json.dumps([mol.GetAtomWithIdx(i).GetAtomMapNum() for i in axis]))
                one.SetProp('absolute_axis_assignment','PENDING independent CIP/reference verification')
            writer.write(one); records.append(row)
    good=[r for r in records if r['converged']]
    if not good: raise RuntimeError(f'No converged MMFF start: {name}')
    best=min(good,key=lambda r:r['MMFF_energy_kcal_mol'])
    conf=mol.GetConformer(best['start']); xyz=conf.GetPositions()
    # Mirror is a coordinate operation, not a claim of experimentally assigned R_a/S_a.
    if axis:
        mirror=Chem.Mol(mol); mirror.RemoveAllConformers(); cf=Chem.Conformer(conf)
        for i,p in enumerate(xyz): cf.SetAtomPosition(i,(-float(p[0]),float(p[1]),float(p[2])))
        mirror.AddConformer(cf); mirror.SetProp('_Name',name+'_mirror'); mirror.SetProp('absolute_axis_assignment','PENDING')
        with Chem.SDWriter(str(path.with_name(name+'_mirror.sdf'))) as w:
            w.SetForceV3000(True); w.write(mirror)
    atoms=[a.GetSymbol() for a in mol.GetAtoms()]
    xyzpath=path.with_suffix('.xyz')
    xyzpath.write_text(str(len(atoms))+'\nMMFF starting geometry; not quantum validated\n'+'\n'.join(f'{a} {x:.10f} {y:.10f} {z:.10f}' for a,(x,y,z) in zip(atoms,xyz))+'\n',encoding='utf8',newline='\n')
    result=dict(requested_starts=starts,embedded=len(ids),converged=len(good),records=records,best_xyz=str(xyzpath.relative_to(ROOT)))
    if axis:
        result['best_axis_dihedral_deg']=best['axis_dihedral_deg']
        result['best_geometric_axis']='axis_positive' if best['axis_dihedral_deg']>0 else 'axis_negative'
        result['mirror_geometric_axis']='axis_negative' if best['axis_dihedral_deg']>0 else 'axis_positive'
    return result


def phase25():
    # Symmetric BINOL graph; O-bearing ortho carbons on BOTH sides of 1,1' axis.
    base=Chem.MolFromSmiles('O=P1(O)Oc2ccc3ccccc3c2-c2c(O1)ccc3ccccc23')
    assert rdMolDescriptors.CalcMolFormula(base)=='C20H13O4P'
    substituents=[('C01','phenyl','c1ccccc1','phenyl'),('C02','4-methylphenyl','c1ccc(C)cc1','phenyl'),
        ('C03','3,5-bis(trifluoromethyl)phenyl','c1cc(C(F)(F)F)cc(C(F)(F)F)c1','phenyl'),
        ('C04','2,4,6-triisopropylphenyl','c1c(C(C)C)cc(C(C)C)cc1C(C)C','phenyl'),
        ('C05','9-anthracenyl','c1c2ccccc2cc2ccccc12','fused_polyaryl'),
        ('C06','1-naphthyl','c1cccc2ccccc12','fused_polyaryl')]
    cats=[]; searches={}
    for cid,name,smiles,family in substituents:
        mol=Chem.Mol(base)
        for anchor in (5,17):
            fragment=Chem.MolFromSmiles(smiles); n=mol.GetNumAtoms(); rw=Chem.RWMol(Chem.CombineMols(mol,fragment)); rw.AddBond(anchor,n,Chem.BondType.SINGLE); mol=rw.GetMol(); Chem.SanitizeMol(mol)
        mol=mapped(mol)
        searches[cid]=preliminary(mol,cid,axis=(4,13,14,15))
        for axis_label in ('axis_positive','axis_negative'):
            cats.append(dict(catalyst_id=cid+'_'+axis_label,pair_id=cid,substituent_3_3prime=name,scaffold_family=family,
                formula=rdMolDescriptors.CalcMolFormula(mol),mapped_connectivity_smiles=Chem.MolToSmiles(mol),
                intended_axis=axis_label,absolute_axis_assignment='PENDING',charge=0,multiplicity=1,
                axis_atoms_maps=[14,15],priority_reference_maps=[5,16],substitution_maps=[6,18],
                origin='prospective_design; absolute stereochemistry not yet certified'))
        print('Phase25 generated',cid,searches[cid]['converged'],flush=True)
    # Four ketimines, with two phenyl and two fused-aryl family members.
    substrates=[('I01','acetophenone N-phenyl ketimine','CC(=Nc1ccccc1)c1ccccc1','phenyl_ketimine'),
        ('I02','4-methylacetophenone N-phenyl ketimine','CC(=Nc1ccccc1)c1ccc(C)cc1','phenyl_ketimine'),
        ('I03','1-acetylnaphthalene N-phenyl ketimine','CC(=Nc1ccccc1)c1cccc2ccccc12','naphthyl_ketimine'),
        ('I04','2-acetylnaphthalene N-phenyl ketimine','CC(=Nc1ccccc1)c1ccc2ccccc2c1','naphthyl_ketimine')]
    ims=[]
    for sid,name,smiles,family in substrates:
        for ez,enum in [('E',Chem.BondStereo.STEREOE),('Z',Chem.BondStereo.STEREOZ)]:
            mol=Chem.MolFromSmiles(smiles); bond=mol.GetBondBetweenAtoms(1,2)
            # Carbon priority: aryl atom 9 over methyl atom 0; nitrogen: phenyl atom 3 over lone pair.
            bond.SetStereoAtoms(9,3); bond.SetStereo(enum); Chem.SetDoubleBondNeighborDirections(mol)
            Chem.AssignStereochemistry(mol,force=True,cleanIt=False)
            assert bond.GetStereo()==enum
            mol=mapped(mol,200)
            searches[sid+'_'+ez]=preliminary(mol,sid+'_'+ez,starts=8)
            ims.append(dict(substrate_id=sid,name=name,family=family,E_Z=ez,charge=0,multiplicity=1,mapped_smiles=Chem.MolToSmiles(mol),imine_C_map=202,imine_N_map=203))
    donor=Chem.MolFromSmiles('COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1')
    donor=mapped(donor,400); searches['HE_dimethyl']=preliminary(donor,'HE_dimethyl',starts=8)
    conditions=[]
    for cat,sid,solv in itertools.product(cats,[s[0] for s in substrates],['toluene','dichloromethane']):
        sf=next(s[3] for s in substrates if s[0]==sid); held_c=cat['scaffold_family']=='fused_polyaryl'; held_s=sf=='naphthyl_ketimine'
        split='test_both' if held_c and held_s else ('test_catalyst' if held_c else ('test_substrate' if held_s else 'train'))
        conditions.append(dict(condition_id=f"{cat['catalyst_id']}_{sid}_{solv}",catalyst_id=cat['catalyst_id'],pair_id=cat['pair_id'],catalyst_family=cat['scaffold_family'],substrate_id=sid,substrate_family=sf,solvent=solv,temperature_K=298.15,imine_M=0.10,HE_equiv=1.20,catalyst_molpercent=5.0,split=split,condition_origin='proposed_not_literature_reproduction',signed_ee=None,absolute_product='selectivity_undetermined'))
    dump(OUT/'phase25'/'catalysts.json',cats); dump(OUT/'phase25'/'imines.json',ims); dump(OUT/'phase25'/'search_audit.json',searches)
    dump(OUT/'phase25'/'donor.json',dict(mapped_smiles=Chem.MolToSmiles(donor),name='dimethyl 2,6-dimethyl-1,4-dihydropyridine-3,5-dicarboxylate',stoichiometry='imine + HEH2 -> amine + aromatic HE; acid regenerated',mapping_status='component maps defined; transferred H identities pending path validation'))
    table(OUT/'phase25'/'conditions.csv',conditions)
    dump(OUT/'phase25'/'blind_protocol.json',dict(seeds=SEEDS,split_before_labels=True,labels_available=False,enantiomer_pairs_grouped=True,
         held_catalyst_families=['fused_polyaryl'],held_substrate_families=['naphthyl_ketimine'],
         compare=['ensemble_network','lowest_conformer_network','phase24_style_descriptor'],
         acquisition=['random','maxmin_diversity'],cost_unit='measured core_seconds on calibrated common hardware; GPU_seconds separate',
         budget_rule='same cumulative electronic-structure cost, including failed jobs and search; fit/scaling train only; no test acquisition',
         evaluation='configuration accuracy with abstention, signed effective ddG MAE, 95% interval coverage, interval width, learning-curve cost AUC',
         preregistration_status='protocol fixed before any labels; external independent custodian and label hash still required'))


def phase26():
    rows=[]
    for facet,U,cation,coverage in itertools.product(['100','111','211'],[-0.4,-0.7,-1.0,-1.3],['K','Cs'],[0.125,0.25]):
        rows.append(dict(interface_id=f'Cu{facet}_{U:.1f}_{cation}_{coverage}',facet=facet,U_V_SHE=U,temperature_K=298.15,pH_bulk=6.8,electrolyte_M=0.10,cation=cation,anion='HCO3',CO_coverage_ML=coverage,coverage_definition='CO per exposed Cu site, step/terrace counts separate',independent_initial_states=3,
           method_status='NOT_RUN',actual_potential_V=None,activation_free_energy_eV=None))
    table(OUT/'phase26'/'interface_matrix.csv',rows)
    dump(OUT/'phase26'/'protocol.json',dict(potential_reference='SHE; U_RHE=U_SHE+(kBT/e)*ln(10)*pH',
        charge_compensation='explicit K+ or Cs+ in equilibrated first layers; continuum electrolyte reservoir and countercharge; matched ionic strength',
        production_method='grand-canonical DFT + explicit-water constrained MD; code/pseudopotentials pending',
        proposed_methods=['PBE-D3(BJ)','RPBE-D3(BJ)'],water_thickness_A=[15,20,25],lateral_size_test='double interface area at fixed density and coverage',
        paths=['2*CO <-> *OCCO','*CO + *CHO <-> *COCHO','*CO + H2O + e- <-> *CHO + OH-','*CO + H2O + e- <-> *COH + OH-',
          '* + H2O + e- <-> *H + OH-','2*H <-> H2 + 2*','*H + H2O + e- <-> H2 + OH- + *'],
        required_extensions=['CO2 supply and CO2-to-CO network','downstream C2 network before reporting final-product FE'],
        path_fields=['N_e','neutral_N_e','total_cell_charge_e','electrode_charge_partition_e','actual_U_V_SHE','target_U_V_SHE','mu_e_eV','G_eV','Omega_eV','water_seed','CVs','forces'],
        CVs=['C-C distance','C-H/O-H proton coordination','water H-bond connectivity','cation coordination'],
        convergence=dict(numerical_eV=0.05,sampling_CI95_halfwidth_eV=0.05,potential_deviation_V=0.05,independent_interfaces=3),
        controls=['constant_charge','implicit_solvent','coverage','cation_positions'],
        forbidden_claims=['RDS from largest uphill reaction free energy','terminal FE from truncated C-C network']))


def phase27():
    series=[('N01','bpy','o-tolyl','Cl','untethered literature complex 1'),('N02','bpy','o-tolyl','Br','designed'),('N03','bpy','o-tolyl','I','designed'),
      ('N04','4,4-dimethyl-bpy','o-tolyl','Cl','designed'),('N05','4,4-bis(trifluoromethyl)-bpy','o-tolyl','Cl','designed'),
      ('N06','6-phenyl-bpy cyclometalated','tethered phenyl','Cl','literature complex 2')]
    dump(OUT/'phase27'/'series.json',[dict(member_id=a,ligand=b,aryl=c,halide=d,origin=e,charge=0,formal_Ni_oxidation=2,
        spin_multiplicities=[1,3,5],solvent='THF',coordinating_solvent_adducts='0 and 1 THF; exchange free energies required',structure_status='literature XYZ available for N01/N06; derivatives pending') for a,b,c,d,e in series])
    dump(OUT/'phase27'/'protocol.json',dict(excitation_windows_nm=[[380,400],[440,460],[510,530]],
        windows_origin='proposed excitation bands; no experimental labels assigned',temperature_K=298.15,independent_batches=3,
        literature_active_space='CAS(10e,9o), 15 singlet and 25 triplet roots, QD-NEVPT2; reference only, insufficient for all Ni-X/ligand channels without expansion',
        required_orbitals=['Ni 3d','occupied bpy pi and pi*','Ni-C sigma/sigma*','Ni-X sigma/sigma*','ligand donor orbitals by occupations/entanglement'],
        convergence=['active-space expansion','basis','state-average weights','root count','SOC','NAC and state overlap phases'],
        dynamics='spin-mixed surface hopping with QM/MM explicit THF; direct and fitted surfaces analysed separately',
        nuclear_dt_fs=[0.5,0.25],electronic_substep_fs=[0.02,0.01],
        sensitivity=['SOC off','alternative decoherence','frustrated-hop policy','larger active space','longer trajectories'],
        outcomes=['first_Ni_C_break','first_Ni_X_break','radical_pair','cage_recombination','solvent_capture','permanent_separation','right_censored'],
        branch_rule='bond excursion alone is never an absorbing product; persistence, recrossing and solvent/cage states required',
        CI=dict(confidence=0.95,target_halfwidth=0.05,unit='independent absorbed-photon trajectories',sequential_rule='fixed preregistered checkpoints or time-uniform intervals; cluster bootstrap by initial-condition batch'),
        observables=['absorption','state populations','lifetime survival','competing-risk cumulative incidence','branch quantum yields','TA=GSB+SE+ESA convolved with instrument response'],
        trajectories_run=0,production_status='NOT_RUN'))


if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    phase25(); phase26(); phase27()
    dump(OUT/'build_metadata.json',dict(rdkit=rdBase.rdkitVersion,seeds=SEEDS,scope='preliminary inputs only',production_calculations_complete=False))
