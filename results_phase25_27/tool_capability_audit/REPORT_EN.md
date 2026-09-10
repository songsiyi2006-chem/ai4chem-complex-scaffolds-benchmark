# New resource applicability audit

Snapshot: 2026-09-10T16:25:07.663924+00:00. These checks establish usable interfaces and model consistency, not completion of Phases25-27.

| Resource | Actual check | Use and remaining boundary |
|---|---|---|
|MACE-MH-1, omat_pbe|4-atom Cu energy -14.969970265 eV; displaced-atom force error 1.202e-09 eV/Angstrom|Candidate for Cu geometry preparation. Fixed-charge learned potential; no grand-canonical electron number, potential control or activation free energy. Symmetry-zero bulk forces do not prove equilibrium lattice constant.|
|PubChem|CID 87116 matches the mapped dimethyl Hantzsch donor after map removal|Structure identity lookup; no reaction-condition or ee validation. Name query returned 404; structure query succeeded.|
|Materials Project OPTIMADE|mp-30 returned; local symmetry Fm-3m, tolerance 0.01 Angstrom|Reference Cu structure. JARVIS query returned HTTP500. Query results are not a local full database.|
|Basis Set Exchange|def2-TZVP retrieved for Ni, Cl, C, N, H with references|Prepare Phase27 basis-convergence jobs; basis data are not a multireference/SOC/NAC engine.|
|ORD|750/750 local protobuf reactions parsed|Reaction provenance and pipeline development. MCP JSON serialization of BLOBs fails; local read-only protobuf decoding works. No records used as blind-test labels.|
|MACE-POLAR and OrbMol-v2|Local checkpoint sizes and SHA256 verified|Candidates for molecular prescreening after task-specific DFT validation. No excited-state branching capability established. MCP molecular endpoint is capped at 100 atoms; use local Python for the complete 117-atom system.|
|QM9, MoleculeNet, Matbench|15 catalog entries total, including PDB and ORD samples|Useful for software/descriptor baselines. They do not supply the specified catalyst ee, electrode PMFs or Ni photodynamic labels.|

ORD reaction-type scope: 1.3.1 [N-arylation with Ar-X] Bromo Buchwald-Hartwig amination: 262; 1.3.4 [N-arylation with Ar-X] Iodo Buchwald-Hartwig amination: 76; 0.0 [Unassigned] Unrecognized: 92; 1.3.2 [N-arylation with Ar-X] Chloro Buchwald-Hartwig amination: 290; 9.7.39 [Other functional group interconversion] Chloro to amino: 1; 1.3.7 [N-arylation with Ar-X] Chloro N-arylation: 9; 1.3.9 [N-arylation with Ar-X] Iodo N-arylation: 1; 1.3.6 [N-arylation with Ar-X] Bromo N-arylation: 16; 1.3.3 [N-arylation with Ar-X] Iodo Buchwald-Hartwig amination: 3. Product measurement counts: {'YIELD': 750}; selectivity subtypes: {}. This is a local sample audit, not a search of the entire ORD archive.

Four model checkpoints passed local file-integrity checks. MACE weight ASL terms are retained; OrbMol licensing is recorded in its installation manifest. PySCF/GPU4PySCF/periodic DFT engines remain uninstalled according to the current provisioning catalog; UMA remains access-required/uninstalled. New model files add useful preprocessing options but do not add a callable OpenAI GPU allocation.

Next use: preserve the running Psi4 frequency queue; validate molecular surrogate forces/relative energies against suitable DFT before broad conformer prescreening; validate the Cu model on slab/adsorbate configurations before geometry preparation. Constant-potential solvent sampling and Ni multireference dynamics remain separate required calculations. Scientific acceptance is unchanged.

[Raw checks](summary.json) · [ORD scope](ord_audit.json) · [PubChem source response](pubchem_structure.json) · [Cu source response](copper_mp30.json) · [Basis with references](basis_Ni.json) · [Force check inputs/results](force_check.json)
