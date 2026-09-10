"""Generate the proposed Phase25 competing-network topology without energies."""
from .build_inputs import OUT,dump


def build():
    # Moieties conserved: catalyst, substrate-derived unit, donor-derived unit.
    species={'Cat':[1,0,0],'Cat2':[2,0,0],'H':[0,0,1],'Hox':[0,0,1],
             'IE':[0,1,0],'IZ':[0,1,0],'PR':[0,1,0],'PS':[0,1,0],'CH':[1,0,1],
             'CPR':[1,1,0],'CPS':[1,1,0],'BR':[1,1,1],'BS':[1,1,1]}
    for ez in ['E','Z']:
        for state in ['HB','IP']: species[f'CI{ez}_{state}']=[1,1,0]
        species['T'+ez]=[1,1,1]
    reactions=[]
    def add(name,r,p,kind):
        for i in range(3):
            assert sum(species[k][i]*v for k,v in r.items())==sum(species[k][i]*v for k,v in p.items())
        reactions.append(dict(reaction_id=name,reactants=r,products=p,reversible=True,kind=kind,
          Gts_kcal_mol=None,path_status='NOT_RUN',required_endpoints=[r,p],atom_mapping_status='PENDING full reactive hydrogen mapping'))
    add('free_EZ',{'IE':1},{'IZ':1},'E_Z_exchange')
    add('dimerization',{'Cat':2},{'Cat2':1},'aggregation')
    add('donor_binding',{'Cat':1,'H':1},{'CH':1},'association')
    for ez in ['E','Z']:
        add('imine_binding_'+ez,{'Cat':1,'I'+ez:1},{'CI'+ez+'_HB':1},'association')
        add('proton_transfer_'+ez,{'CI'+ez+'_HB':1},{'CI'+ez+'_IP':1},'protonation')
        add('ternary_via_imine_'+ez,{'CI'+ez+'_IP':1,'H':1},{'T'+ez:1},'association')
        add('ternary_via_donor_'+ez,{'CH':1,'I'+ez:1},{'T'+ez:1},'association')
        for rs in ['R','S']:
            add('hydride_'+ez+'_'+rs,{'T'+ez:1},{'B'+rs:1},'H_transfer')
            add('uncatalyzed_'+ez+'_'+rs,{'I'+ez:1,'H':1},{'P'+rs:1,'Hox':1},'uncatalyzed_control')
    for state in ['HB','IP']: add('bound_EZ_'+state,{'CIE_'+state:1},{'CIZ_'+state:1},'E_Z_exchange')
    add('ternary_EZ',{'TE':1},{'TZ':1},'E_Z_exchange')
    for rs in ['R','S']:
        add('product_release_'+rs,{'B'+rs:1},{'Cat':1,'P'+rs:1,'Hox':1},'release')
        add('product_inhibition_'+rs,{'Cat':1,'P'+rs:1},{'CP'+rs:1},'product_inhibition')
    payload=dict(species=[dict(id=k,conserved_moieties=v,G0_kcal_mol=None) for k,v in species.items()],reactions=reactions,
        status='proposed_topology_not_validated_complete_network',charge_rule='all aggregates net neutral; bound IP denotes iminium/phosphate pair',
        required_extensions=['conformational microstates of every aggregate','donor/oxidized-donor inhibition','solvent-separated ions if populated','higher catalyst aggregates if concentration study requires','proton-relay intermediates if located'],
        controls=['mirror all 3D coordinates','achiral diphenyl phosphate O=P(O)(Oc1ccccc1)Oc1ccccc1','Cat concentration zero'],
        note='No rates are assigned to uncomputed edges. Network cannot generate a chemical prediction until G/TS data and accepted pathways exist.')
    dump(OUT/'phase25'/'network.json',payload)
    print(len(species),'species,',len(reactions),'reversible edges; conserved moieties checked')


if __name__=='__main__': build()
