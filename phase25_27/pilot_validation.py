"""Independent raw-output audits for the low-cost pilot, not DFT acceptance.

Reject printed imaginary modes even if xTB's thermochemistry reports zero after
its -20 cm^-1 cutoff. Repeated printouts must not double the mode count.
"""
import json
from pathlib import Path
import re
import numpy as np
from rdkit import Chem
from scipy.constants import h,k,c,R,physical_constants
from .build_inputs import ROOT,OUT

HARTREE_KCAL=627.5094740631


def xyz(path):
    lines=Path(path).read_text().splitlines();n=int(lines[0]);atoms=[s.split() for s in lines[2:2+n]]
    if len(atoms)!=n:raise ValueError('truncated XYZ')
    return [s[0] for s in atoms],np.array([[float(x) for x in s[1:4]] for s in atoms])


def frequency_block(raw,natoms):
    blocks=[]
    for section in raw.split('projected vibrational frequencies')[1:]:
        values=[];started=False
        for line in section.splitlines()[1:]:
            if line.strip().startswith('eigval :'):
                started=True;values.extend(float(x) for x in line.split(':',1)[1].split())
            elif started:break
        if values:blocks.append(values)
    if not blocks or any(len(a)!=3*natoms for a in blocks):raise ValueError('Hessian mode count does not equal 3N')
    if any(not np.allclose(blocks[0],a,atol=.011,rtol=0) for a in blocks[1:]):raise ValueError('Inconsistent repeated frequency blocks')
    modes=np.array(blocks[-1])
    if np.max(np.abs(modes[:6]))>1:raise ValueError('Six translation/rotation modes not projected')
    return modes[6:]


def qrrho_entropy(freqs,cutoff,T=298.15):
    """Grimme entropy interpolation with B_av=1e-44 kg m2; no enthalpy change."""
    freqs=np.asarray(freqs,float)
    if np.any(freqs<=0):raise ValueError('Imaginary/zero vibrational modes cannot enter thermochemistry')
    hz=freqs*c*100;x=h*hz/(k*T)
    sho=R*(x/np.expm1(x)-np.log(-np.expm1(-x)))
    mu=h/(8*np.pi**2*hz);mu=mu*1e-44/(mu+1e-44)
    sfr=R*(.5+np.log(np.sqrt(8*np.pi**3*mu*k*T/h**2)))
    w=1/(1+(cutoff/freqs)**4)
    return float(np.sum(w*sho+(1-w)*sfr))


def dihedral(p):
    b0=p[0]-p[1];b1=p[2]-p[1];b2=p[3]-p[2];b1=b1/np.linalg.norm(b1)
    v=b0-np.dot(b0,b1)*b1;w=b2-np.dot(b2,b1)*b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1,v),w),np.dot(v,w))))


def last_number(raw,pattern):
    hits=re.findall(pattern,raw,re.I)
    return float(hits[-1]) if hits else None


def audit_job(directory,row):
    record=dict(name=row['name'],campaign=directory.parent.name,status=row['status'],method=row['method'],solvation=row['solvation'],
      acceptance_scope='GFN2-xTB preliminary local minimum only; not DFT or catalytic-path acceptance',accepted_local_minimum=False)
    if row['status'] in ['STARTED','HESSIAN_RUNNING']:
        record['review_status']='CALCULATION_IN_PROGRESS';return record
    if not (directory/'xtbopt.xyz').exists():return record
    elements,pos=xyz(directory/'xtbopt.xyz');before,_=xyz(directory/'input.xyz')
    if elements!=before:raise ValueError('Atom ordering or composition changed')
    record['atom_count']=len(elements)
    if not (directory/'hessian.out').exists():return record
    raw=(directory/'hessian.out').read_text(errors='replace')
    try:modes=frequency_block(raw,len(elements))
    except ValueError as exc:record['rejection']=str(exc);return record
    record.update(vibrational_mode_count=len(modes),lowest_frequency_cm1=float(modes.min()),
      frequencies_below_minus1_cm1=modes[modes < -1].tolist(),frequencies_between_minus1_and_plus1_cm1=modes[np.abs(modes)<=1].tolist())
    record['accepted_local_minimum']=bool(row.get('optimization_converged') and row.get('hessian_exit_code')==0 and np.all(modes>1))
    if not record['accepted_local_minimum']:record['rejection']='Failed convergence or unresolved nonpositive/near-zero vibrational mode'
    energy=last_number(raw,r'TOTAL ENERGY\s+([-+\d.]+)\s+Eh')
    rawg=last_number(raw,r'TOTAL FREE ENERGY\s+([-+\d.]+)\s+Eh')
    record.update(energy_including_ALPB_hartree=energy,raw_engine_G_hartree=rawg,
      engine_thermochemistry='298.15 K; default 50 cm^-1 rotor interpolation; gas translation with default gsolv ALPB reference; raw value not certified standard-state solution G',
      corrected_solution_G_hartree=None)
    thermo=OUT/'phase25/thermo_sensitivity/results.json'
    if thermo.exists():
        match=next((r for r in json.loads(thermo.read_text()) if r['name']==row['name'] and r['campaign']==directory.parent.name),None)
        if match and match['replay50_pass']:
            record['entropy_only_G_shift_cutoff50_to100_kcal_mol']=match['G_shift_50_to100_kcal_mol']
            record['entropy_sensitivity_method']='native xtb thermo; original geometry and Hessian, geometry-dependent inertia'
    if len(elements) in [28,31,34] and row['name'].startswith('I'):
        angle=dihedral(pos[[9,1,2,3]])
        record['imine_aryl_C_N_aryl_dihedral_deg']=angle
        record['geometric_E_Z']='E' if abs(angle)>90 else 'Z'
        record['C_N_distance_A']=float(np.linalg.norm(pos[1]-pos[2]))
        intended='E' if '_E' in row['name'] else 'Z'
        if record['geometric_E_Z']!=intended:
            record['accepted_local_minimum']=False;record['rejection']='E/Z identity changed'
    if row['name'].startswith('C01_I01'):
        starts=json.loads((OUT/'phase25/ternary_starts/manifest.json').read_text())
        meta=next((m for m in starts if m['name']==row['name']),None)
        if meta:
            distance=lambda a,b:float(np.linalg.norm(pos[meta[a]]-pos[meta[b]]))
            contacts={key:distance(a,b) for key,a,b in [
              ('acid_O_H','acid_O_index','acid_H_index'),('imine_N_acid_H','imine_N_index','acid_H_index'),
              ('donor_C_H','donor_C_index','donor_H_index'),('imine_C_donor_H','imine_C_index','donor_H_index'),
              ('imine_C_N','imine_C_index','imine_N_index')]}
            record['diagnostic_contacts_A']=contacts
            record['acid_proton_location_by_distance']='N_bound_candidate' if contacts['imine_N_acid_H']<1.25 and contacts['acid_O_H']>1.25 else ('O_bound_candidate' if contacts['acid_O_H']<1.25 else 'shared_or_unresolved')
            record['imine_geometry_E_Z']='E' if abs(dihedral(pos[[67,59,60,61]]))>90 else 'Z'
            record['imine_E_Z_input']=meta['imine_E_Z']
            record['imine_E_Z_preserved']=record['imine_geometry_E_Z']==meta['imine_E_Z']
            record['catalyst_axis_dihedral_deg']=dihedral(pos[[4,13,14,15]])
            _,startpos=xyz(OUT/'phase25/ternary_starts'/(row['name']+'.xyz'))
            record['axis_sign_preserved']=bool(dihedral(startpos[[4,13,14,15]])*record['catalyst_axis_dihedral_deg']>0)
            pt=Chem.GetPeriodicTable();ratios=[]
            for a,b,_ in meta['covalent_bonds']:
                if elements[a]=='H' or elements[b]=='H':continue
                radii=pt.GetRcovalent(pt.GetAtomicNumber(elements[a]))+pt.GetRcovalent(pt.GetAtomicNumber(elements[b]))
                ratios.append(float(np.linalg.norm(pos[a]-pos[b])/radii))
            record['largest_expected_heavy_bond_distance_over_covalent_radii']=max(ratios)
            record['heavy_connectivity_distance_screen_pass']=max(ratios)<1.45
            if not record['axis_sign_preserved'] or not record['heavy_connectivity_distance_screen_pass']:
                record['accepted_local_minimum']=False;record['rejection']='Catalyst axis or heavy-connectivity screen failed'
            record['product_CIP']='UNASSIGNED; no IRC/path validation'
    return record


def main():
    records=[]
    for ledger in sorted((OUT/'phase25/xtb_pilot').glob('*/results.json')):
        for row in json.loads(ledger.read_text()):
            records.append(audit_job(ledger.parent/row['name'],row))
    result=dict(records=records,accepted_preliminary_minima=sum(r['accepted_local_minimum'] for r in records),
      total_jobs=len(records),chemical_barriers_validated=0,selected_ee=None,
      note='This audit does not promote tight-binding Hessians to DFT acceptance; all imaginary modes are retained.')
    dest=OUT/'phase25/xtb_pilot/validation.json';dest.write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Audited',len(records),'xTB jobs;',result['accepted_preliminary_minima'],'preliminary local minima')


if __name__=='__main__':main()
