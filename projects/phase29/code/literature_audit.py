"""Recalculate public SI thermochemistry/electrical bookkeeping, no new DFT."""
import csv
import json
from pathlib import Path
from phase29 import ROOT, FARADAY_C_MOL, save_json, csv_write, panel_mol, proposed_product, carbonyl_partner, canonical
from rdkit import Chem
from rdkit.Chem import Descriptors

HARTREE_KCAL_MOL=627.509474

def run_public_audit(out):
    rows=list(csv.DictReader((ROOT/'data/literature/natcomm_tableS3_energies.csv').open(encoding='utf8')))
    energies={r['species']:float(r['G_Hartree']) for r in rows}
    thermal_error=max(abs(float(r['Gcorr_Hartree'])+float(r['E_Hartree'])-float(r['G_Hartree'])) for r in rows)
    reactants=energies['INT-1']+energies['radical_B1']
    b1=(energies['TS1']-reactants)*HARTREE_KCAL_MOL
    b2=(energies['TS2']-reactants)*HARTREE_KCAL_MOL
    second=energies['INT-2']+energies['TEMPO+']
    payload={'evidence_layer':'recalculation_of_published_DFT_table_not_new_electronic_structure',
             'source':'10.1038/s41467-026-71858-2 SI Table S3 ppS19-S20',
             'Gcorr_plus_E_minus_G_max_abs_Hartree':thermal_error,
             'TS1_from_INT1_plus_radicalB1_kcal_mol':b1,
             'TS2_from_INT1_plus_radicalB1_kcal_mol':b2,
             'TS2_minus_TS1_kcal_mol':b2-b1,
             'TS3_from_INT2_plus_TEMPOplus_kcal_mol':(energies['TS3']-second)*HARTREE_KCAL_MOL,
             'TS4_from_INT2_plus_TEMPOplus_kcal_mol':(energies['TS4']-second)*HARTREE_KCAL_MOL,
             'not_done':['rerun_DFT','frequency_or_IRC_revalidation','interfacial_free_energy','conversion_to_final_regioselectivity'],
             'SI_conflicts':['TableS1 heading C1 but footnote b C1 prime; retained as source ambiguity']}
    save_json(out/'results/public_DFT_table_audit.json',payload)
    jacs_sub=panel_mol({'id':'JACS_A1','mapped_smiles':'[n:1]1[cH:2][cH:3][c:4](-c2ccccc2)[cH:5][cH:6]1'})
    benzaldehyde=Chem.MolFromSmiles('[CH:100](=[O:101])[c:102]1[cH:103][cH:104][cH:105][cH:106][cH:107]1')
    jacs_prod=proposed_product(jacs_sub,2,benzaldehyde)
    nc_sub=panel_mol({'id':'NC_A1','mapped_smiles':'[n:1]1[c:2](-c2ccccc2)[cH:3][cH:4][cH:5][cH:6]1'})
    nc_prod=proposed_product(nc_sub,4)
    identity=[]
    for name,sub,partner,prod,site,blocked in [('JACS_A1_B1_C1',jacs_sub,benzaldehyde,jacs_prod,2,4),('NC_A1_B1_C1',nc_sub,carbonyl_partner(),nc_prod,4,2)]:
        identity.append({'id':name,'substrate':canonical(sub),'partner':canonical(partner),'product':canonical(prod),
                         'mapped_reaction':Chem.MolToSmiles(sub)+'.'+Chem.MolToSmiles(partner)+'>>'+Chem.MolToSmiles(prod),
                         'fixed_parent_site':site,'blocked_fixed_parent_site':blocked,
                         'product_identity_status':'source_structure_visually_audited_heavy_atom_map_reconstructed',
                         'stereochemistry':'unassigned; no enantioselectivity inferred','MW_g_mol':Descriptors.MolWt(prod)})
    csv_write(out/'data/literature/reference_reaction_identity.csv',identity)
    electrical=[]
    for name,i,t,n,y,zlist,mw,v,vol in [('JACS_model',.014,11,.0003,.85,[2,4],Descriptors.MolWt(jacs_prod),3.01,.006),('NC_model',.025,9,.0002,.81,[2],Descriptors.MolWt(nc_prod),None,.010),('NC_scaleup',.040,80,.003,.62,[2],Descriptors.MolWt(nc_prod),None,.030)]:
        q=i*t*3600; mass=n*y*mw
        electrical.append({'record':name,'evidence_layer':'derived_from_public_SI_not_own_experiment','limiting_reagent_mol':n,'isolated_major_yield_fraction':y,'current_A':i,'time_h':t,'charge_C':q,
                           'charge_F_per_mol_limiting_input':q/FARADAY_C_MOL/n,'charge_F_per_mol_product':q/FARADAY_C_MOL/(n*y),
                           'FE_z2_percent':2*FARADAY_C_MOL*n*y/q*100,'FE_z4_percent':4*FARADAY_C_MOL*n*y/q*100 if 4 in zlist else None,
                           'product_g_reconstructed':mass,'voltage_V_reported':v,'electricity_kWh_per_kg_if_voltage_constant':i*t*v/1000/(mass/1000) if v is not None else None,
                           'STY_g_L_h_using_solution_volume':mass/vol/t,'solution_volume_L':vol,
                           'energy_boundary':'cell electrical energy only; excludes heating, separation, workup and equipment'})
    csv_write(out/'data/literature/derived_electrical_metrics.csv',electrical)
    save_json(out/'results/public_source_summary.json',{'literature_controls_rows':sum(1 for _ in csv.DictReader((ROOT/'data/literature/condition_controls.csv').open(encoding='utf8'))),
              'DFT_species_rows':len(rows),'reference_identity_rows':len(identity),'electrical_rows':len(electrical),
              'same_substrate_partner_switch_demonstrated_by_these_two_models':False,
              'reason':'JACS C4 occupied pyridine+aldehyde; NC C2 occupied pyridine+ketone; products and input structures differ',
              'NC_source_spreadsheet':'Sheet1 B3:B582 contains coordinate labels/xyz; Sheet2 empty; no additional reaction-yield dataset',
              'unresolved_JACS_SI_conflicts':['A1 synthesis scheme4-bromo versus prose2-bromo','catalyst preparation5mg versus TON paragraph2.5mg']})
    return payload,electrical

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=ROOT)
    args=p.parse_args();run_public_audit(args.out)
