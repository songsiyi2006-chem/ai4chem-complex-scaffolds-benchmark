"""Scientific plots of actual pilot coordinates and audited E/Z differences."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .build_inputs import OUT
from .pilot_validation import xyz


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    out=OUT/'phase25/figures';out.mkdir(exist_ok=True)
    name='C01_I01_E_face-1_start0'
    meta=next(r for r in json.loads((OUT/'phase25/ternary_starts/manifest.json').read_text()) if r['name']==name)
    sy,pos=xyz(OUT/'phase25/xtb_pilot/ternary_toluene_01'/name/'xtbopt.xyz')
    groups=np.array(meta['components']);colors=['#557ea4','#df9733','#4f9b6c']
    fig=plt.figure(figsize=(12,5.4),layout='constrained');ax=fig.add_subplot(121,projection='3d')
    for a,b,order in meta['covalent_bonds']:
        ax.plot(*pos[[a,b]].T,color=colors[groups[a]],lw=1.2,alpha=.7)
    for g,label in enumerate(['Full C01 catalyst','E imine','Hantzsch ester']):
        mask=(groups==g)&(np.array(sy)!='H');ax.scatter(*pos[mask].T,color=colors[g],s=22,label=label,depthshade=False)
    ax.set_box_aspect(np.ptp(pos,axis=0));ax.view_init(25,30);ax.set_axis_off();ax.legend(loc='lower left',frameon=False)
    ax.set_title('117-atom optimized ternary complex',loc='left',fontweight='bold')
    ax=fig.add_subplot(122)
    atoms=[meta[k] for k in ['acid_O_index','acid_H_index','imine_N_index','imine_C_index','donor_H_index','donor_C_index']]
    labels=['Acid O','Acid H','Imine N','Imine C','Donor H','Donor C'];acolors=['#c84f45','#777777','#357cc0','#df9733','#777777','#4f9b6c']
    schematic=np.array([[0,0],[1.4,0],[2.8,0],[3.2,1.8],[1.8,3.2],[.2,3.2]])
    for point,label,color in zip(schematic,labels,acolors):
        ax.scatter(*point,s=75,color=color,zorder=3);ax.annotate(label,point,xytext=(0,14),textcoords='offset points',ha='center',fontsize=10)
    for a,b in [(0,1),(1,2),(2,3),(3,4),(4,5)]:
        ax.plot(*schematic[[a,b]].T,color='#555555',ls='--',lw=1)
        midpoint=(schematic[a]+schematic[b])/2
        ax.annotate(f'{np.linalg.norm(pos[atoms[a]]-pos[atoms[b]]):.2f} A',midpoint,xytext=(0,-16),textcoords='offset points',ha='center',fontsize=10,bbox={'facecolor':'white','edgecolor':'none','alpha':.8})
    ax.set(xlim=(-.6,3.9),ylim=(-.7,4));ax.set_aspect('equal');ax.set_axis_off()
    ax.set_title('Geometry-derived contact distances\n(schematic layout)',loc='left',fontweight='bold')
    fig.suptitle('GFN2-xTB / ALPB(toluene): a local minimum, not a transition state',fontsize=13,fontweight='bold')
    fig.text(.52,.02,'N-H 1.05 A; O...H 1.63 A. Hydride remains on donor.\nDistance diagnostics do not establish a catalytic mechanism.',fontsize=10)
    fig.savefig(out/'ternary_minimum.png',dpi=180);plt.close(fig)
    data=json.loads((OUT/'phase25/pilot_progress.json').read_text())['E_Z_comparisons']
    fig,ax=plt.subplots(figsize=(9,4.5),layout='constrained')
    for i,solvent in enumerate(['ALPB(toluene)','ALPB(ch2cl2)']):
        rows=[r for r in data if r['solvation']==solvent];x=np.arange(len(rows))+(i-.5)*.2
        a=np.array([r['Z_minus_E_local_minimum_G_cutoff50_kcal_mol'] for r in rows])
        b=np.array([r['Z_minus_E_local_minimum_G_cutoff100_kcal_mol'] for r in rows])
        color=['#35799b','#dc9149'][i]
        ax.plot(x,a,'o',color=color,label=solvent+'; cutoff 50')
        ax.plot(x,b,'_',color=color,ms=15,label=solvent+'; cutoff 100')
        for xx,aa,bb in zip(x,a,b):ax.plot([xx,xx],[aa,bb],color=color,lw=1)
    ax.axhline(0,color='#888888',lw=.8);ax.set_xticks(range(4),['I01','I02','I03','I04'])
    ax.set_ylabel('G(Z) - G(E), kcal/mol');ax.set_title('Audited local-minimum differences at the xTB level',loc='left',fontweight='bold')
    ax.legend(ncol=2,frameon=False,fontsize=9);ax.grid(axis='y',alpha=.2)
    fig.supxlabel('Lines show a rotor-cutoff sensitivity, not confidence intervals. No reaction ee is inferred.',fontsize=9)
    fig.savefig(out/'imine_EZ_sensitivity.png',dpi=180);plt.close(fig)
    print('Saved two scientific figures')


if __name__=='__main__':main()
