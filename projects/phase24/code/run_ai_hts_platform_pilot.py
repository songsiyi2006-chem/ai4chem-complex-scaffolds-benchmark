#!/usr/bin/env python3
"""Phase 24: reproducible GXNU HTE SOFTWARE PILOT, not an experimental claim.

python run_ai_hts_platform_pilot.py --out . --pptx /path/to/source.pptx
python run_ai_hts_platform_pilot.py --out another_empty_directory --seed 25
python run_ai_hts_platform_pilot.py --self-test

Dependencies: requirements_phase24.txt. Opentrons is optional and isolated.
No network, robot connection, git operation or external data download occurs.
RDKit EHT is a semiempirical model, geometry/charge descriptors are proxies.
All response data in the default run originate from an explicitly synthetic
oracle. The learner receives integrated chromatograms, never oracle labels.
The generated OT-2 protocol performs WATER commissioning by default. Organic
chemistry is locked until a separate, reviewed hardware/chemical qualification.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import importlib.metadata
import io
import json
import logging
import math
import sqlite3
import sys
import time
import unittest
import warnings
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, nnls
from scipy.special import expit
from scipy.stats import norm
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator, rdEHTTools
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning

LOG = logging.getLogger('phase24')
PPT_SHA = '3345f6610d1c31e20d376ce214ff2c0b3eeb8eface7a7b7a4e0f7363261f490b'
SOURCE_NAME = '广西师范大学-人工智能药物分子合成平台论证会-260904-1.pptx'
PPT_FACTS = {
    '3-4': '无水无氧自动化、多反应并行、结构化实验记录与AI可学习数据资产。',
    '7-8': '固体称量、液体工作站、手套箱、堆栈、机械臂、LC、封膜与离心浓缩。',
    '9': '一期独立运行，模块化设计，预留软硬件扩展接口。',
    '11': 'JSON任务参数，HTTP请求，MQTT状态返回，设备协议由执行端适配。',
    '12': '性能指标导向，不指定品牌型号，公开采购。',
    '14': '列示85万元数据模块、215万元合成模块、350万元合计、300万元优惠价及19.8094万元装修；分项包含关系需另行核对。',
    '15': '公共共享平台、科研支撑、人才培养与标准化数据资源。',
}
REFS = {
    'EHT': 'https://www.rdkit.org/docs/source/rdkit.Chem.rdEHTTools.html',
    'charges': 'https://www.rdkit.org/docs/source/rdkit.Chem.rdPartialCharges.html',
    'OT2': 'https://docs.opentrons.com/python-api/pipettes/characteristics/',
    'temperature': 'https://docs.opentrons.com/python-api/modules/temperature-module/',
    'liquids': 'https://docs.opentrons.com/python-api/building-block-commands/liquids/',
    'oxalamide': 'https://sioc.cas.cn/sourcedb/cn/lw/202306/t20230621_6784633.html',
}
METALS = ['CuI', 'Cu(OAc)2', '[Rh(COD)Cl]2', '[Ir(COD)Cl]2',
          '[Ru(p-cymene)Cl2]2', 'Pd(OAc)2']
BASES = ['DBU', 'DIPEA', 'triethylamine', 'pyridine', 'DABCO',
         'K2CO3', 'KOtBu', 'NaOtBu']
SOLVENTS = ['ethanol', 'acetonitrile', 'toluene', 'DMF']
ARYLS = ['c1ccccc1', 'c1ccccc1OC', 'c1ccc(OC)cc1',
         'c1ccc(C(F)(F)F)cc1', 'c1ccc2ccccc2c1', 'c1cccs1']
FAMILIES = ['chiral_oxalamide', 'amino_alcohol', 'phosphino_amine', 'pyridyl_oxazoline']
FEATURES = ['MW', 'TPSA', 'logP', 'L_proxy_A', 'B1_proxy_A', 'B5_proxy_A',
            'donor_charge_mean', 'donor_charge_min', 'donor_charge_max',
            'charge_dipole_proxy_D', 'EHT_gap_eV', 'CIP_R_minus_S']
RT = np.array([1.8, 3.60, 4.08, 5.18, 6.20])  # calibration, minutes
RF = np.array([1.0, 1.15, 1.15, .8, 1.1])     # UV area per nmol
SIGMA = np.array([.060, .050, .052, .070, .060])
PEAKS = ['substrate', 'R', 'S', 'byproduct', 'internal_standard']
TIME = np.linspace(0, 7.2, 721)
INITIAL_NMOL, IS_NMOL = 100., 25.
VOLUMES = {'substrate_uL': 16., 'catalyst_ligand_premix_uL': 10.,
           'base_uL': 10., 'makeup_solvent_uL': 44.}


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n',
                    encoding='utf-8', newline='\n')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, records):
    if not records:
        raise ValueError('Empty table')
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader()
        w.writerows(records)


def ppt_context(path=None):
    context = dict(source=SOURCE_NAME, sha256=PPT_SHA, slides=15,
                   page_number_note='References 14/15 are printed labels on current file slides 13/14.',
                   verified_page_notes=PPT_FACTS, raw_deck_published=False,
                   extraction='embedded, previously verified page notes')
    if path:
        if sha(path) != PPT_SHA:
            raise ValueError('PPTX differs from the reviewed version. Review it before reusing page claims.')
        with zipfile.ZipFile(path) as z:
            names = sorted((n for n in z.namelist() if n.startswith('ppt/slides/slide')
                            and n.endswith('.xml')),
                           key=lambda n: int(Path(n).stem[5:]))
            context['slide_text_counts'] = [len(ET.fromstring(z.read(n)).findall(
                './/{http://schemas.openxmlformats.org/drawingml/2006/main}t')) for n in names]
            if len(names) != 15:
                raise ValueError('Unexpected slide count')
        context['extraction'] = 'source hash and all 15 slide XML text layers checked'
    return context


def ligand_smiles(family, aryl, chirality):
    ch = '@' if chirality == 0 else '@@'
    return [f'O=C(N[C{ch}H](C){aryl})C(=O)NCC',
            f'OC[C{ch}H](N){aryl}',
            f'N[C{ch}H]({aryl})CP(c1ccccc1)c1ccccc1',
            f'[C{ch}H]9({aryl})COC(=N9)c1ccccn1'][family]


def steric_envelope(mol):
    """Donor-centroid axial molecular envelope, NOT canonical Sterimol.

    L = full axial vdW extent; B1/B5 = min/max radial support over 180 angles.
    The donor centroid is a geometric convention, not an optimized metal site.
    """
    xyz = mol.GetConformer().GetPositions()
    donor = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ('N', 'P', 'O')]
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    origin = xyz[donor].mean(0)
    axis = xyz[heavy].mean(0)-origin
    if np.linalg.norm(axis) < 1e-5:
        axis = np.linalg.svd(xyz[heavy]-origin, full_matrices=False)[2][0]
    axis /= np.linalg.norm(axis)
    ref = np.eye(3)[np.argmin(np.abs(axis))]
    e1 = np.cross(axis, ref); e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    p = xyz-origin
    radii = np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    along = p@axis
    L = (along+radii).max()-(along-radii).min()
    theta = np.linspace(0, 2*np.pi, 180, endpoint=False)
    dirs = np.cos(theta)[:, None]*e1+np.sin(theta)[:, None]*e2
    supports = (p@dirs.T+radii[:, None]).max(0)
    return float(L), float(supports.min()), float(supports.max())


def build_ligands(seed):
    ligands, fingerprints = [], []
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=512,
                                                          includeChirality=True)
    for f in range(4):
        for a, aryl in enumerate(ARYLS):
            for stereo in range(2):
                mol = Chem.MolFromSmiles(ligand_smiles(f, aryl, stereo))
                if mol is None:
                    raise RuntimeError('Invalid generated ligand')
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                if not centers or any(c == '?' for _, c in centers):
                    raise RuntimeError('Missing specified chirality')
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                fingerprints.append(generator.GetFingerprint(mol))
                hmol = Chem.AddHs(mol)
                params = AllChem.ETKDGv3()
                params.randomSeed = int(seed+len(ligands)*37)
                ids = list(AllChem.EmbedMultipleConfs(hmol, numConfs=3, params=params))
                if not ids or not AllChem.MMFFHasAllMoleculeParams(hmol):
                    raise RuntimeError(f'No conformer/MMFF parameters: {smiles}')
                opt = AllChem.MMFFOptimizeMoleculeConfs(hmol, numThreads=1, maxIters=1500)
                converged = [(i, energy) for i, (status, energy) in zip(ids, opt) if status == 0]
                if not converged:
                    raise RuntimeError(f'MMFF not converged: {smiles}')
                best, energy = min(converged, key=lambda v: v[1])
                one = Chem.Mol(hmol); conf = Chem.Conformer(hmol.GetConformer(best))
                one.RemoveAllConformers(); one.AddConformer(conf, assignId=True)
                AllChem.ComputeGasteigerCharges(one, throwOnParamFailure=True)
                charges = np.array([atom.GetDoubleProp('_GasteigerCharge') for atom in one.GetAtoms()])
                if not np.isfinite(charges).all() or abs(charges.sum()) > 1e-5:
                    raise RuntimeError('Invalid neutral-ligand charges')
                donor_ids = [atom.GetIdx() for atom in one.GetAtoms() if atom.GetSymbol() in ('N','P','O')]
                q = charges[donor_ids]
                xyz = one.GetConformer().GetPositions()
                dipole = np.linalg.norm((charges[:, None]*(xyz-xyz.mean(0))).sum(0))*4.8032047
                ok, eht = rdEHTTools.RunMol(one)
                if not ok or eht.numElectrons % 2:
                    raise RuntimeError('Closed-shell extended Huckel calculation failed')
                orbitals = np.asarray(eht.GetOrbitalEnergies())
                homo = eht.numElectrons//2-1
                gap = float(orbitals[homo+1]-orbitals[homo])
                L, B1, B5 = steric_envelope(one)
                vals = [Descriptors.MolWt(mol), Descriptors.TPSA(mol), Descriptors.MolLogP(mol),
                        L, B1, B5, q.mean(), q.min(), q.max(), dipole, gap,
                        sum(1 if cip == 'R' else -1 for _, cip in centers)]
                if not np.isfinite(vals).all() or gap < 0:
                    raise RuntimeError('Invalid descriptor')
                ligands.append(dict(ligand_id=f'L{len(ligands)+1:02}', family=FAMILIES[f],
                    aryl_variant=a, stereo_variant=stereo, smiles=smiles,
                    conformer_method='ETKDGv3; lowest converged of 3 MMFF conformers',
                    MMFF_energy_kcal=float(energy), descriptors=dict(zip(FEATURES,map(float,vals))),
                    heteroatom_indices=donor_ids, heteroatom_Gasteiger_charges=q.tolist(),
                    atoms=[at.GetSymbol() for at in one.GetAtoms()], coordinates_A=xyz.tolist()))
        LOG.info('Ligands: %d/48', len(ligands))
    if len(set(x['smiles'] for x in ligands)) != 48:
        raise RuntimeError('Duplicate stereochemical structures')
    return ligands, fingerprints


def build_space(ligands):
    rows, features = [], []
    ld = StandardScaler().fit_transform([[x['descriptors'][k] for k in FEATURES] for x in ligands])
    for l in range(48):
        for m in range(6):
            for b in range(8):
                for s in range(4):
                    # Software demonstration mask only. All real chemistry remains unapproved.
                    envelope = b < 5  # avoid treating solid base suspensions as verified liquid stocks
                    rows.append(dict(condition_id=len(rows), ligand_id=ligands[l]['ligand_id'],
                        ligand_index=l, metal_index=m, base_index=b, solvent_index=s,
                        metal=METALS[m], base=BASES[b], solvent=SOLVENTS[s],
                        demo_liquid_envelope=envelope, wet_lab_approved=False,
                        requested_reaction_temperature_C=60., provenance='virtual_design'))
                    features.append(np.r_[ld[l]/np.sqrt(ld.shape[1]),
                        np.eye(6)[m], np.eye(8)[b], np.eye(4)[s],
                        np.eye(4)[FAMILIES.index(ligands[l]['family'])]])
    return rows, np.array(features), np.array([r['demo_liquid_envelope'] for r in rows])


def maxmin_seeds(X, rows, fingerprints, eligible, n, seed):
    ids = np.flatnonzero(eligible)
    chosen = [int(np.random.default_rng(seed).choice(ids))]
    li = np.array([r['ligand_index'] for r in rows])
    mind = np.full(len(X), np.inf)
    for _ in range(n-1):
        last = chosen[-1]
        d = np.linalg.norm(X-X[last], axis=1)/np.sqrt(X.shape[1])
        tan = 1-np.array(DataStructs.BulkTanimotoSimilarity(fingerprints[li[last]], fingerprints))
        mind = np.minimum(mind, .5*d+.5*tan[li])
        mind[~eligible] = -np.inf; mind[chosen] = -np.inf
        chosen.append(int(np.argmax(mind)))
    return np.array(chosen)


def synthetic_oracle(rows):
    """Fixed, disclosed toy landscape, NEVER passed to acquisition/training.

    Values are engineered software test data, not chemical predictions. This
    function has no literature fit. It does not search seeds until success.
    """
    out = []
    for r in rows:
        l, m, b, s = (r[k] for k in ('ligand_index','metal_index','base_index','solvent_index'))
        family, aryl, stereo = l//12, (l%12)//2, l%2
        quality = ([.0,.4,1.3,.7][family] + [0,.3,.7,.15,.45,-.4][aryl]
                   + [.0,.15,1.15,.7,.5,-.4][m] + [.6,.2,.45,-.15,.1,-.4,-.7,-.5][b]
                   + [.35,.15,-.2,-.1][s] + .4*math.sin((family+1)*(m+1)))
        conversion = 55+44*expit(quality-.5)
        selectivity = .88+.105*expit(quality-1)
        product_yield = conversion*selectivity
        ee_mag = 55+44.6*expit(quality-1.05+.4*(family == 2)-.35*(s == 2))
        ee_signed = ee_mag*(1 if stereo == 0 else -1)
        out.append([conversion, product_yield, ee_signed])
    return np.array(out)


def peak_basis(t, shift=0., scale=1., width=1., eta=.18):
    centers = .45 + (RT-.45)*scale + shift  # column dead time separately accounted
    sigma = SIGMA*width
    z = t[:, None]-centers
    gaussian = np.exp(-.5*(z/sigma)**2)/(sigma*np.sqrt(2*np.pi))
    gamma = sigma*np.sqrt(2*np.log(2))
    lorentz = gamma/(np.pi*(z*z+gamma*gamma))
    return (1-eta)*gaussian+eta*lorentz


def simulate_trace(truth, seed):
    rng = np.random.default_rng(seed)
    conv, y, ee = truth
    # Injection and recovery variation canceled by independently known internal standard.
    amounts = np.array([100-conv, y*(1+ee/100)/2, y*(1-ee/100)/2, conv-y, IS_NMOL])
    shift, scale = rng.normal(0,.024), rng.normal(1,.003)
    width, eta = rng.uniform(.90,1.16), rng.uniform(.10,.27)
    signal = peak_basis(TIME,shift,scale,width,eta)@(amounts*RF*rng.uniform(.90,1.10))
    signal += rng.uniform(.1,.8)+rng.uniform(-.015,.015)*TIME
    signal += rng.normal(0,.025,size=len(TIME))
    return signal


def integrate_trace(t, signal):
    """Calibrated pseudo-Voigt variable projection, no oracle/true areas as input.

    Estimate baseline from independent peak-free intervals. Jointly fit retention
    shift, time dilation, width and Gaussian/Lorentzian mixing. Nonnegative areas
    are solved by NNLS at every nonlinear step. Missing/bad traces fail QC.
    """
    t, signal = np.asarray(t,float), np.asarray(signal,float)
    if t.ndim != 1 or t.shape != signal.shape or len(t)<100 or not np.isfinite(signal).all():
        raise ValueError('Invalid chromatogram')
    if not np.all(np.diff(t)>0) or t[0]>.1 or t[-1]<7.:
        raise ValueError('Invalid retention-time coverage/order')
    baseline_mask = (t<.7) | (t>7.)
    # Baseline and peak tails solved together (unconstrained baseline, nonnegative areas).
    B = np.column_stack([np.ones(len(t)),t-t.mean()])
    Q = np.linalg.qr(B, mode='reduced')[0]
    detrended = signal-Q@(Q.T@signal)
    def solve(p):
        A = peak_basis(t,*p)
        projected = A-Q@(Q.T@A)
        areas, _ = nnls(projected, detrended)
        baseline = np.linalg.lstsq(B,signal-A@areas,rcond=None)[0]
        fitted = A@areas+B@baseline
        return areas, fitted
    fit = least_squares(lambda p: solve(p)[1]-signal, [0,1,1,.18],
                        bounds=([-.15,.98,.65,.02],[.15,1.02,1.55,.45]),
                        max_nfev=70, xtol=1e-8, ftol=1e-8, gtol=1e-7)
    areas, fitted = solve(fit.x)
    rms = float(np.sqrt(np.mean((fitted-signal)**2)))
    if areas[-1] < .1:
        raise ValueError('Missing internal standard')
    amounts = areas/RF*(IS_NMOL/(areas[-1]/RF[-1]))
    product = amounts[1]+amounts[2]
    ee = float(100*(amounts[1]-amounts[2])/product) if product>.1 else None
    conv, y = 100*(1-amounts[0]/INITIAL_NMOL), 100*product/INITIAL_NMOL
    # Conservative Gaussian baseline-width estimate; tailing still needs empirical validation.
    rs = float((RT[2]-RT[1])*fit.x[1]/(2*(SIGMA[1]+SIGMA[2])*fit.x[2]))
    mass_error = float(abs(amounts[:4].sum()-INITIAL_NMOL))
    passed = bool(fit.success and rms<.15 and rs>=1.5 and mass_error<=5
                  and ee is not None and 0<=y<=100.5 and -.5<=conv<=100.5)
    return dict(conversion_pct=float(conv), yield_pct=float(y), ee_signed_pct=ee,
                ee_abs_pct=abs(ee) if ee is not None else None,
                areas=areas.tolist(), nmol=amounts.tolist(), fit_parameters=fit.x.tolist(),
                residual_rms=rms, resolution_R_S=rs, mass_balance_error_nmol=mass_error,
                minor_enantiomer_above_demo_LOQ=bool(min(amounts[1:3])>=.05),
                qc_pass=passed), fitted


def init_db(path, rows):
    db = sqlite3.connect(path)
    db.execute('PRAGMA foreign_keys=ON')
    db.execute('CREATE TABLE IF NOT EXISTS conditions (condition_id INTEGER PRIMARY KEY, spec TEXT NOT NULL)')
    db.execute('''CREATE TABLE IF NOT EXISTS measurements (
        measurement_id TEXT PRIMARY KEY, condition_id INTEGER NOT NULL REFERENCES conditions,
        round INTEGER NOT NULL, well TEXT NOT NULL, provenance TEXT NOT NULL,
        raw_sha256 TEXT NOT NULL, yield_pct REAL, ee_abs_pct REAL, conversion_pct REAL,
        qc_pass INTEGER NOT NULL, record TEXT NOT NULL, UNIQUE(round, well, provenance))''')
    if db.execute('SELECT COUNT(*) FROM measurements').fetchone()[0]:
        db.close()
        raise FileExistsError('Existing measurements: choose a new --out directory; no data is overwritten.')
    db.executemany('INSERT OR IGNORE INTO conditions VALUES (?,?)',
                   [(r['condition_id'],json.dumps(r)) for r in rows])
    db.commit()
    return db


def ingest(db, condition_id, round_id, well, record, raw_hash, provenance='simulated'):
    if provenance not in ('simulated','experimental') or not valid_well(well):
        raise ValueError('Invalid provenance or well')
    # No duplicate replacement. Failed QC stays in the database but not in training.
    key = f'{provenance}:R{round_id}:{well}'
    with db:
        db.execute('INSERT INTO measurements VALUES (?,?,?,?,?,?,?,?,?,?,?)',
            (key,int(condition_id),int(round_id),well,provenance,raw_hash,
             record['yield_pct'],record['ee_abs_pct'],record['conversion_pct'],
             int(record['qc_pass']),json.dumps(record,allow_nan=False)))


def training_rows(db):
    return db.execute('SELECT condition_id,yield_pct,ee_abs_pct FROM measurements '
                      "WHERE qc_pass=1 AND provenance='simulated' ORDER BY rowid").fetchall()


def utility(y, ee):
    return np.minimum(np.asarray(y)/90., np.asarray(ee)/95.)


def gp_fit(X, target, seed):
    kernel = ConstantKernel(1.,(.05,5.))*Matern(length_scale=2., length_scale_bounds=(.3,10), nu=2.5)
    kernel += WhiteKernel(.002,(1e-5,.08))
    gp = GaussianProcessRegressor(kernel=kernel, alpha=1e-6, normalize_y=True,
                                  random_state=seed, n_restarts_optimizer=0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', ConvergenceWarning)
        gp.fit(X,target)
    return gp


def select_batch(X, observed, eligible, n, seed):
    """Greedy Monte Carlo JOINT qEI on bottleneck-scaled objectives.

    A scalar GP models u=min(Y/90, |ee|/95). Correlated joint posterior draws
    score improvement of the batch maximum. Candidate shortlisting (EI plus
    random coverage) and greedy selection are approximations to global qEI.
    Auxiliary GPs supply yield/ee forecasts. No hidden outcomes are accessed.
    """
    data = np.asarray(observed,float)
    train = data[:,0].astype(int); y, ee = data[:,1],data[:,2]
    u = utility(y,ee)
    gp = gp_fit(X[train],u,seed)
    pool = np.flatnonzero(eligible & ~np.isin(np.arange(len(X)),train))
    if len(pool)<n:
        raise ValueError('Insufficient unmeasured eligible conditions')
    mean, sd = gp.predict(X[pool],return_std=True)
    imp = mean-u.max(); z = imp/np.maximum(sd,1e-12)
    ei = imp*norm.cdf(z)+sd*norm.pdf(z)
    rng = np.random.default_rng(seed)
    best = pool[np.argsort(-ei,kind='stable')[:min(320,len(pool))]]
    extra = rng.choice(pool,size=min(160,len(pool)),replace=False)
    shortlist = np.unique(np.r_[best,extra])
    mu, covariance = gp.predict(X[shortlist],return_cov=True)
    covariance = (covariance+covariance.T)/2
    chol = np.linalg.cholesky(covariance+1e-8*np.eye(len(shortlist)))
    draws = mu[:,None]+chol@rng.standard_normal((len(shortlist),256))
    incumbent = np.full(256,u.max())
    chosen, gains = [], []
    for _ in range(n):
        gain = np.maximum(draws-incumbent,0).mean(1)
        gain[chosen] = -np.inf
        if np.max(gain)<=1e-12:
            # Marginal qEI saturated: deterministic posterior variance exploration.
            gain = np.diag(covariance).copy(); gain[chosen] = -np.inf
        j = int(np.argmax(gain)); chosen.append(j)
        gains.append(float(np.maximum(draws[j]-incumbent,0).mean()))
        incumbent = np.maximum(incumbent,draws[j])
    ids = shortlist[chosen]
    ygp, egp = gp_fit(X[train],y,seed), gp_fit(X[train],ee,seed)
    ym, ys = ygp.predict(X[ids],return_std=True)
    em, es = egp.predict(X[ids],return_std=True)
    predictions = {int(i):dict(predicted_yield_pct=float(a), predicted_ee_abs_pct=float(b),
                             yield_sd=float(c),ee_sd=float(d))
                   for i,a,b,c,d in zip(ids,ym,em,ys,es)}
    return ids, predictions, dict(kernel=str(gp.kernel_), shortlist=len(shortlist),
                                  mc_samples=256, marginal_qei=gains,
                                  training_ids=train.tolist())


def valid_well(name):
    return isinstance(name,str) and len(name)>=2 and name[0] in 'ABCDEFGH' and name[1:].isdigit() and 1<=int(name[1:])<=12 and name==f'{name[0]}{int(name[1:])}'


def well_names(n=96):
    if not 0<=n<=96:
        raise ValueError('Plate size out of bounds')
    return [f'{row}{col}' for col in range(1,13) for row in 'ABCDEFGH'][:n]


def make_plate(ids, rows, predictions, round_id):
    if len(set(map(int,ids))) != len(ids) or len(ids)>96:
        raise ValueError('Duplicate conditions or plate overflow')
    plate = []
    for well,i in zip(well_names(len(ids)),ids):
        r = rows[int(i)]
        plate.append(dict(round=round_id, well=well, condition_id=int(i),
            ligand_id=r['ligand_id'],metal=r['metal'],base=r['base'],solvent=r['solvent'],
            source_catalyst_plate=well,source_base_plate=well,source_solvent_plate=well,
            **VOLUMES,total_uL=sum(VOLUMES.values()),
            substrate_stock_mM=6.25, catalyst_metal_equivalent_mM=.5,ligand_stock_mM=1.,
            base_stock_mM=20., simulated_substrate_nmol=100.,
            source_premix_load_uL=50.,source_base_load_uL=50.,source_solvent_load_uL=100.,
            predicted_yield_pct=predictions.get(int(i),{}).get('predicted_yield_pct'),
            predicted_ee_abs_pct=predictions.get(int(i),{}).get('predicted_ee_abs_pct'),
            chemistry_approved=False,provenance='simulation_plan'))
    return plate


PROTOCOL_TEMPLATE = r'''"""Phase 24 OT-2 API 2.15 commissioning protocol.
Default WATER ONLY. Source plates must be pre-arrayed by a separately validated
single-channel/workstation process using the matching plate map. No automatic
solid dosing, inert atmosphere, sealing or chromatograph control is claimed.
"""
from opentrons import protocol_api
metadata = {'protocolName':'GXNU Phase24 96-well WATER commissioning',
            'author':'GXNU pilot prototype', 'apiLevel':'2.15'}
CHEMISTRY_APPROVED = False
RUN_CHEMISTRY = False
PLATE_MAP = __PLATE_MAP__

def run(ctx: protocol_api.ProtocolContext):
    if RUN_CHEMISTRY and not CHEMISTRY_APPROVED:
        raise RuntimeError('Chemical/material/atmosphere/seal risk assessment not approved')
    ctx.pause('Verify WATER-only source plates and map, tips, calibrated labware offsets. '
              'No organic stocks or pressurized hydrogen in commissioning.')
    module = ctx.load_module('temperature module gen2', '1')
    adapter = module.load_adapter('opentrons_96_well_aluminum_block')
    plate = adapter.load_labware('nest_96_wellplate_100ul_pcr_full_skirt')
    reservoir = ctx.load_labware('nest_12_reservoir_15ml', '2')
    catalyst = ctx.load_labware('nest_96_wellplate_2ml_deep', '3')
    base = ctx.load_labware('nest_96_wellplate_2ml_deep', '8')
    solvent = ctx.load_labware('nest_96_wellplate_2ml_deep', '9')
    tips20 = [ctx.load_labware('opentrons_96_tiprack_20ul', slot) for slot in ('4','5','6')]
    tips300 = [ctx.load_labware('opentrons_96_tiprack_300ul', '7')]
    p20 = ctx.load_instrument('p20_multi_gen2','left',tip_racks=tips20)
    p300 = ctx.load_instrument('p300_multi_gen2','right',tip_racks=tips300)
    module.set_temperature(25)  # dispensing at ambient; never heat an open volatile plate

    def move(pip, source, destination, volume, gap):
        assert volume+gap <= pip.max_volume
        pip.pick_up_tip()
        pip.flow_rate.aspirate = 3 if pip is p20 else 20
        pip.flow_rate.dispense = 5 if pip is p20 else 30
        pip.flow_rate.blow_out = 8 if pip is p20 else 40
        # Prewet twice into the same source before touching any destination.
        for _ in range(2):
            pip.aspirate(volume,source.bottom(1))
            pip.dispense(volume,source.bottom(1))
        pip.aspirate(volume,source.bottom(1))
        ctx.delay(seconds=1)
        pip.air_gap(gap)
        # API 2.15 tracks air as volume; explicitly discharge gap above liquid first.
        pip.dispense(gap,destination.top(-1))
        pip.dispense(volume,destination.bottom(1))
        pip.blow_out(destination.top(-1))
        pip.touch_tip(destination, radius=.7, v_offset=-1, speed=20)
        pip.drop_tip()

    assert len(PLATE_MAP)==96
    for col in range(1,13):
        # A-row references address all EIGHT channels. Never loop over 96 wells.
        dest = plate[f'A{col}']
        move(p300,solvent[f'A{col}'],dest,44,5)
        move(p20,reservoir['A1'],dest,16,2)
        move(p20,base[f'A{col}'],dest,10,2)
        move(p20,catalyst[f'A{col}'],dest,10,2)
    if RUN_CHEMISTRY:
        ctx.pause('Stop for independently validated compatible sealing, inert handling and '
                  'temperature/pressure review. Confirm before approved 60 C hold.')
        module.set_temperature(60)
        ctx.delay(minutes=30)  # illustrative hold only, not a validated reaction duration
    module.deactivate()
    ctx.comment('Dispensed 80 uL/well; 288 P20 tips and 96 P300 tips. '
                'Chromatography requires a separate validated sampling method.')
'''


def generate_protocol(path, plate):
    if len(plate)!=96 or any(p['total_uL']!=80 for p in plate):
        raise ValueError('Protocol requires exactly one 80-uL 96-well plate')
    if [p['well'] for p in plate] != well_names():
        raise ValueError('Protocol requires column-major A1..H12 layout')
    for p in plate:
        if any(p[k] != v for k,v in VOLUMES.items()):
            raise ValueError('Protocol template volumes differ from plate map')
        if any(p[k] != p['well'] for k in ['source_catalyst_plate','source_base_plate','source_solvent_plate']):
            raise ValueError('Multi-channel source plates must match destination columns')
    text = PROTOCOL_TEMPLATE.replace('__PLATE_MAP__',repr(plate))
    ast.parse(text)
    Path(path).write_text(text,encoding='utf-8',newline='\n')


def pareto_mask(y, ee):
    a = np.column_stack([y,ee])
    return np.array([not np.any(np.all(a>=p,axis=1)&np.any(a>p,axis=1)) for p in a])


def run_measurements(db, plate, truths, root, seed):
    signals, fitted, records, errors = [], [], [], []
    rawdir = root/'results_phase24'/'chromatograms'; rawdir.mkdir(exist_ok=True)
    for item in plate:
        i, r, well = item['condition_id'],item['round'],item['well']
        signal = simulate_trace(truths[i], seed+10007*i+103*r)
        raw = rawdir/f'R{r}_{well}.csv'
        with raw.open('w',encoding='utf-8',newline='') as f:
            f.write('retention_time_min,uv_mAU\n')
            np.savetxt(f,np.column_stack([TIME,signal]),delimiter=',',fmt='%.10g')
        # Read exactly the saved detector output back. No oracle input to integration.
        loaded = np.loadtxt(raw,delimiter=',',skiprows=1)
        record, fit = integrate_trace(loaded[:,0],loaded[:,1])
        record.update(condition_id=i,round=r,well=well,provenance='simulated',
                      raw_path=str(raw.relative_to(root)).replace('\\','/'))
        ingest(db,i,r,well,record,sha(raw))
        records.append(record); signals.append(signal); fitted.append(fit)
        errors.append([record['conversion_pct']-truths[i,0],record['yield_pct']-truths[i,1],
                       record['ee_signed_pct']-truths[i,2]])
    return records, np.array(signals), np.array(fitted), np.array(errors)


def figures(root, X, eligible, seeds, rounds, records, predictions, traces, fits, random_history):
    figdir = root/'figures_hts_pilot'; figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False,'axes.spines.right':False})
    LOG.info('t-SNE: fitting all %d conditions (visualization only)', len(X))
    tsne = TSNE(n_components=2,perplexity=40,init='pca',learning_rate='auto',
                random_state=24,max_iter=750)
    projection = tsne.fit_transform(X)
    np.savez_compressed(root/'results_phase24'/'design_tsne.npz',condition_id=np.arange(len(X)),xy=projection)
    fig,ax = plt.subplots(figsize=(9,7))
    ax.scatter(*projection[~eligible].T,s=3,c='#cccccc',label='Outside demo liquid envelope')
    ax.scatter(*projection[eligible].T,s=3,c='#abc5cf',label='Demo envelope (NOT safety approval)')
    colors = ['#215e8c','#df8a24','#a5373a']
    ax.scatter(*projection[seeds].T,s=50,marker='*',c='black',label='24 MaxMin seeds',zorder=5)
    centroids = [projection[seeds].mean(0)]
    for r,ids in enumerate(rounds):
        ax.scatter(*projection[ids].T,s=15,c=colors[r],label=f'Batch {r+1}: 96')
        centroids.append(projection[ids].mean(0))
    centroids = np.array(centroids)
    ax.plot(*centroids.T,'k--',lw=1,label='Batch centroid path (embedding only)')
    ax.set(title='9,216 virtual conditions: chemistry-aware sampling',xlabel='t-SNE 1',ylabel='t-SNE 2')
    ax.legend(fontsize=8,loc='best'); fig.tight_layout()
    fig.savefig(figdir/'fig1_chemical_design_space_tsne.png',dpi=300); plt.close(fig)

    final = [r for r in records if r['round']==3]
    fig,axs=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for ax,key,title,pred in zip(axs.ravel(),['yield_pct','ee_abs_pct']*2,
                                 ['GP yield forecast','GP |ee| forecast','HPLC-derived yield','HPLC-derived |ee|'],
                                 [True,True,False,False]):
        v=np.full((8,12),np.nan)
        for r in final:
            value = predictions[r['condition_id']]['predicted_'+('yield_pct' if key=='yield_pct' else 'ee_abs_pct')] if pred else r[key]
            v[ord(r['well'][0])-65,int(r['well'][1:])-1]=value
        im=ax.imshow(v,vmin=0,vmax=100,cmap='cividis',aspect='auto')
        ax.set(xticks=range(12),xticklabels=range(1,13),yticks=range(8),yticklabels=list('ABCDEFGH'),title=title)
        for (i,j),value in np.ndenumerate(v):
            ax.text(j,i,f'{value:.0f}',ha='center',va='center',fontsize=6,color='white' if value<50 else 'black')
        fig.colorbar(im,ax=ax,label='%')
    fig.suptitle('Round 3 | simulated measurements, no physical experiment')
    fig.savefig(figdir/'fig2_96well_plate_heatmap_yield_ee.png',dpi=300); plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    counts,best=[],[]
    for r in range(4):
        seen=[m for m in records if m['round']<=r and m['qc_pass']]
        counts.append(len(seen));best.append(float(utility([m['yield_pct'] for m in seen],[m['ee_abs_pct'] for m in seen]).max()))
    axs[0].plot(counts,best,'o-',color='#215e8c',label='Joint qEI, integrated HPLC')
    baseline=np.array(random_history)
    axs[0].plot(counts,baseline.mean(0),'s--',color='gray',label='Random, 5 seeds (oracle benchmark)')
    axs[0].fill_between(counts,baseline.min(0),baseline.max(0),color='gray',alpha=.15,label='Random min/max')
    axs[0].axhline(1,color='black',ls=':',label='Simultaneous 90% yield / 95% |ee|')
    axs[0].set(xlabel='Conditions evaluated (24 + 3 x 96)',ylabel='Best min(yield/90, |ee|/95)',title='Observed convergence; success is not imposed')
    axs[0].legend(fontsize=7)
    for r,c in zip(range(4),['black']+colors):
        a=[m for m in records if m['round']==r and m['qc_pass']]
        axs[1].scatter([m['yield_pct'] for m in a],[m['ee_abs_pct'] for m in a],s=16,c=c,label=f'R{r}')
    good=[m for m in records if m['qc_pass']]
    y=np.array([m['yield_pct'] for m in good]); e=np.array([m['ee_abs_pct'] for m in good])
    front=pareto_mask(y,e); order=np.argsort(y[front])
    axs[1].plot(y[front][order],e[front][order],'k-',lw=1)
    axs[1].axvline(90,color='gray',ls=':');axs[1].axhline(95,color='gray',ls=':')
    axs[1].set(xlabel='Product yield (%)',ylabel='|ee| (%)',title='Measured synthetic Pareto front');axs[1].legend()
    fig.savefig(figdir/'fig3_active_learning_convergence_pareto.png',dpi=300);plt.close(fig)

    fig,axs=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for ax,j in zip(axs.ravel(),[0,31,63,95]):
        r=final[j]
        ax.plot(TIME,traces[j],c='#215e8c',lw=.8,label='Synthetic detector')
        ax.plot(TIME,fits[j],'--',c='#ab442e',lw=.8,label='Fitted total')
        basis=peak_basis(TIME,*r['fit_parameters'])
        for k,pk in enumerate(PEAKS):
            ax.fill_between(TIME,0,basis[:,k]*r['areas'][k],alpha=.18,label=pk)
        ax.set(title=f"{r['well']} | yield {r['yield_pct']:.1f}%, |ee| {r['ee_abs_pct']:.1f}%, Rs {r['resolution_R_S']:.2f}",
               xlabel='Retention time (min)',ylabel='UV response (mAU)')
    axs[0,0].legend(fontsize=6,ncol=3)
    fig.suptitle('Simulated HPLC with calibrated internal standard and peak integration')
    fig.savefig(figdir/'fig4_automated_hplc_96well_trace.png',dpi=300);plt.close(fig)
    return float(tsne.kl_divergence_)


def briefs_and_report(root, metrics):
    m=metrics; best=m['best_joint_condition']
    zh=f'''# 广西师范大学 AI 辅助高通量不对称催化平台试点建议书

## 第一页：建设定位与已交付能力

本建议依据《{SOURCE_NAME}》当前15页材料。文中14/15按页内原标号，对应当前文件第13/14页。项目定位为可安装运行的软件与数据模板，供一期平台联调使用，尚未通过化学实验及设备验收。PPT第3–4页提出的无水无氧操作、多反应并行、全程数字化，在本试点中分别落实为人工审核边界、96孔任务编排、可追溯数据库。第12页要求性能指标导向与品牌开放，OT-2仅为可替换适配器，不构成采购指定。

本次构建48个手性结构、6种金属前体、8种碱和4种溶剂，共9,216个虚拟条件。48个结构由24种连接结构及其对映体组成，不代表48种已合成或已采购配体。描述符包括TPSA、供体原子Gasteiger电荷、几何空间包络、点电荷偶极矩近似和扩展Hückel能隙。模型不能替代DFT、配位构象计算或配体实验表征。

24孔MaxMin初始批次后，GP联合后验Monte Carlo qEI依次推荐三块96孔板。模拟色谱经拟合积分后进入SQLite，模型不读取隐藏真值。总计{m['measurements']}条模拟记录，QC通过{m['qc_pass_count']}条。三轮后同一条件的产物收率为{best['yield_pct']:.2f}%，|ee|为{best['ee_abs_pct']:.2f}%；联合达标条件数为{m['joint_hits']}。是否达标是计算输出，不构成真实催化效果或方法优于随机搜索的保证。

建议一期验收先覆盖文件格式、孔位映射、单位、日志、幂等性与数据导出，再进入有水示踪和真实反应。软件原型不替代总控、手套箱、机械臂、密封及LC厂商系统。

<div style="page-break-after: always;"></div>

## 第二页：接口、实施与人员职责

依PPT第11页，前端提交JSON任务，后端通过HTTP下发，经设备适配器执行后以MQTT回传状态。本试点提供本地任务与事件契约，不宣称已部署HTTP服务、MQTT代理或厂商控制接口。任务包含唯一ID、孔位、试剂批次接口、体积、审批状态、原始色谱校验值及QC结果。

OT-2采用P20/P300八通道、温控模块、12槽储液槽及三块预排源板。八个通道同步动作，不把96个条件错误实现为96次独立通道操作。预排配体/金属混合液、碱与溶剂分别占一块源板，由已验证的单通道工作站或人工准备。每孔80 µL，整板288支P20吸头和96支P300吸头，默认仅运行水替代液、25 ℃加样。60 ℃有机反应必须另行确认材料相容性、密封、惰性气氛和热/压风险。此OT-2布置不执行氢气加压氢化。

建议0–2周完成数据字典、接口仿真和单元测试，3–4周开展水/染料称量验证与孔位追踪，5–8周由导师选择单一低风险、已有方法的基准反应并验证手性色谱，后续再推进盲测、重复与跨批次评估。上述为建议时间表，非已承诺交付周期。

本科生可在培训和监督下承担库存标识、协议单元测试、原始数据整理、标曲/空白数据检查与软件回归。仪器首次通电、气路、密封、危化品和反应审批由有资质人员负责，不能用软件“通过”代替实验安全审批。

PPT第14页给出350万元合计、300万元优惠价和19.8094万元装修，分项与包含关系需预算责任人书面确认。本试点不补算差额、不新增采购报价、不将OT-2成本写入既有预算。下一决策门：确定一个具体底物/产物、手性色谱方法、源液溶解性和供应商接口后，批准小规模实验验证。
'''
    en=f'''# GXNU AI-assisted high-throughput catalysis pilot: executive brief

## Page 1: Scope and demonstrated software capability

The current 15-slide GXNU platform deck dated 4 September 2026 calls for inert handling, parallel synthesis and structured research data (slides 3–4). References to 14/15 use printed labels, corresponding to file slides 13/14. This pilot supplies an executable software/data template for commissioning. It is not a commissioned synthesis platform. Slide 12 specifies performance-based, brand-neutral procurement. The OT-2 implementation is one replaceable adapter, not a purchase recommendation.

The virtual library contains 48 stereochemical structures (24 connectivities and their enantiomers), six metal precursors, eight bases and four solvents: 9,216 conditions. RDKit supplies TPSA, Gasteiger charges, geometric steric-envelope proxies, a point-charge dipole proxy and an extended-Huckel frontier gap. These are not DFT or measured ligand properties; the molecules are not certified stock items.

Twenty-four diverse seeds precede three 96-well recommendations. A GP uses correlated posterior samples for greedy joint qEI on min(yield/90, |ee|/95). Chromatographic fits, not oracle labels, enter the learning database. This run produced {m['measurements']} synthetic records, {m['qc_pass_count']} passing QC, and {m['joint_hits']} simultaneous target hits. The best joint condition yielded {best['yield_pct']:.2f}% with {best['ee_abs_pct']:.2f}% |ee| in the simulation. Success is reported as observed and is not guaranteed in real chemistry.

Deliverables include a standalone script, SQLite database, raw detector traces, plate/source maps, an API 2.15 protocol, four 300-dpi figures and an auditable metrics manifest. Existing supplier controllers and safety interlocks remain responsible for real equipment.

<div style="page-break-after: always;"></div>

## Page 2: Integration and decision gates

Slide 11 describes JSON task specifications, HTTP requests and MQTT state returns. Local schemas and example events demonstrate this contract; no live server, broker or vendor link has been commissioned. Preserve task IDs, well positions, provenance, units and raw-file hashes across adapters. Reject duplicates and quarantine QC failures before model updates.

Eight-channel pipettes actuate all eight tips together. Three pre-arrayed source plates enable different chemistry in every destination well while a 12-channel reservoir supplies common substrate. Preparation of those source plates is a separate, reviewed workflow. Dispensing totals 80 uL/well, 288 P20 tips and 96 P300 tips. Default execution uses water at 25 C. A 60 C chemistry branch requires independent approval of materials, atmosphere, sealing and thermal/pressure risks. No pressurized hydrogen handling is implemented.

Proposed gates: weeks 0–2 for schemas and software tests; weeks 3–4 for gravimetric water/dye commissioning; weeks 5–8 for a supervisor-selected benchmark and chiral HPLC validation. These are planning assumptions. Students can assist with labels, database checks, simulations and supervised calibration records. Qualified staff retain authority over utilities, gas systems, chemical handling and experiment approval.

Slide 14 lists a 3.50-million-CNY total, 3.00-million-CNY discounted price and 198,094 CNY of renovation. The inclusion relationships require written reconciliation. No price gap or incremental OT-2 budget is inferred here. The next approval requires a defined substrate/product, calibrated chiral method, stock-solubility evidence and documented vendor interfaces.
'''
    report=f'''# Phase 24 广西师范大学人工智能药物分子合成平台试点综合报告

## 1. 结论与证据等级

本项目已经实现并运行软件闭环：分子结构构建、描述符、主动学习选样、孔板任务编排、协议生成、模拟HPLC、积分入库和再训练。当前证据等级为可复现计算与协议仿真，物理设备及化学实验尚未验收。所有输出的实验结果字段带有`provenance=simulated`；不能在论文或简历中表述为实测产率、不对称催化新方法或自主实验室建成。

本报告与中英文两页式决策摘要分别面向技术核查和平台决策。Markdown中的显式分页用于排版建议，实际页数取决于渲染器。

## 2. PPT需求逐项对应

来源为用户提供的《{SOURCE_NAME}》，SHA-256：`{PPT_SHA}`。文件在处理过程中由16页更新为15页，删除可视化展示页后重新核对当前全部15页文字层，核心建设要求和预算文字不变。以下第14/15页按页内原标号引用，对应当前文件第13/14页。对设备图片不作型号识别或空间尺寸推断。原始PPT未随本报告公开上传。

| PPT页码 | 原始要求 | 本次实现 | 留给真实平台的工作 |
|---|---|---|---|
| 3–4 | 无水无氧、多条件并行、数字记录 | 96孔映射、显式审批字段、数据库 | 手套箱联锁、气氛与溶剂验证 |
| 7–8 | 固液投料、机械臂、LC、封膜与浓缩 | 液体任务与源板配置、色谱数据结构 | 固体称量、机械臂与密封驱动 |
| 9 | 模块化与扩展接口 | 单文件可运行、设备协议独立导出 | 总控及厂商SDK集成 |
| 11 | JSON / HTTP / MQTT | 任务与状态事件样例、幂等数据主键 | 网络服务、鉴权、消息重试和故障恢复 |
| 12 | 性能采购、不指定品牌 | OT-2作为参考适配器 | 正式技术参数与公开采购 |
| 14 | 经费与建设预算 | 引用已列合计并提示核对 | 确认折扣、装修及分项包含关系 |
| 15 | 公共共享、科研与人才培养 | 开放数据格式、分层学生任务 | 平台运营、预约与权限制度 |

PPT第7页的用电/占地总值仅为原方案估算，不作为本代码的电气或建筑设计依据。第14页的85、215、350、300及19.8094万元不能直接任意相加；软件原型不声称补全预算中的缺项。

## 3. 化学任务与48配体的边界

使用手性草酰胺、氨基醇、膦氮、吡啶噁唑啉四类结构，每类6种芳基变化及2种对映体，共48个立体化学唯一SMILES、24种连接结构。RDKit核验价态、明确的手性中心和3D构象。结构合法不等于可合成、可采购或能催化指定反应。六种金属前体及八种碱是虚拟变量目录，不表示任意金属/配体/碱组合相容。

选题受草酰胺配体促进铜催化偶联的研究启发。上海有机所的2017年论文记录列有伍海波与马大为等作者，支持该研究背景，但不证明本48个手性衍生结构或模拟ee与其论文一致。[原始机构论文记录]({REFS['oxalamide']})。把一般Ullmann偶联直接写成高ee反应并不成立。真实先导实验必须先定义底物、产物、手性来源与对照，再选择已有先例的偶联或转移氢化模型；本试点没有执行气体氢化。

描述符方法：ETKDGv3生成3个构象，经MMFF优化后取最低能已收敛构象；这不是全构象搜索。L为供体质心轴的整体范德华包络长度，B1/B5为径向支持函数极值，均明确为Sterimol-like近似，不是规范金属配位轴Sterimol。TPSA为拓扑量。Gasteiger电荷使用N/P/O原子，偶极矩用全部显式氢的点电荷位置加权求和并换算Debye，属于经验电荷近似。[RDKit电荷文档]({REFS['charges']})。轨道能隙来自RDKit的扩展Hückel接口，未计算DFT、溶剂化金属配合物轨道或实验光谱。[EHT文档]({REFS['EHT']})。CIP编码单独保留，避免把无手性描述符当成对映选择性信息。

## 4. 设计空间与主动学习算法

条件数48×6×8×4=9,216。供体/空间/电子特征标准化，金属、碱、溶剂和结构家族采用独热编码。初始24个条件以混合Tanimoto和欧氏距离贪心MaxMin选取；t-SNE只用于事后可视化，不参与选样。演示液体掩码排除三个未验证固体碱槽位后有{m['demo_envelope_count']}个候选。该掩码仅是数据/液体工作流演示范围，所有条件的`wet_lab_approved`仍为false，不是化学安全认证。

模型训练只使用QC通过的积分结果。两个目标由`u=min(Y/90,|ee|/95)`形成瓶颈效用，u≥1才表示同一条件同时达标。采用Matern GP和白噪声项，先按单点EI及随机覆盖形成候选短名单，再从联合后验协方差生成256次相关抽样，逐步最大化批次最大值的期望改进。批量qEI为贪心蒙特卡罗近似，不是独立EI排序，也不是严格多目标EHVI；效用压缩可能牺牲Pareto面覆盖。边际qEI耗尽时按后验方差探索，代码保留每步边际增益。另两个GP提供收率/ee预测。

固定运行预算为24个初始点+3×96=312个条件。数据产生、模型预测、事后隐藏真值误差分别保存，算法看不到尚未测量的响应。对照为同一初始种子集合和相同样本预算的5条随机轨迹，其终点使用无测量噪声oracle，仅作玩具景观诊断，不能据此证明统计显著优越性。需要后续跨种子、多噪声水平及真实历史数据回放评测。

## 5. 液体自动化的可执行范围

按列A1–A12引用八通道，同时服务A–H行。三块96孔预排源板分别承载金属/配体预混液、碱液和补足溶剂；12槽储液槽的A1为共同底物。这解决12槽容不下48配体及全部离散试剂的问题，但源板准备仍依赖单通道工作站/人工与条码核验，程序没有虚构这部分已自动化。[八通道约束]({REFS['OT2']})。

每孔模拟体积为底物16、预混10、碱10、补溶剂44 µL，共80 µL。示例底物6.25 mM给出100 nmol，金属浓度按金属原子当量0.5 mM给出5 mol%，二聚前体应按两个金属中心换算；配体1 mM给出10 mol%，碱20 mM给出2当量。这是体积/单位测试夹具，真实 stock 浓度、载体溶剂、溶解性和反应条件均需单独验证，不是实验处方。源板每孔预加载50/50/100 µL，保留40/40/56 µL余量；储液槽底物至少准备2.6 mL，实际死体积需称量确认。

甲板1温控、2储液槽、3/8/9源板、4/5/6为P20吸头、7为P300吸头。P20最大18 µL（含气隙），P300最大49 µL。预润洗在接触目标前回原源孔，吸液后等待并加气隙，气隙先在目标孔上方排出，随后低速排液并触壁；每次转移更换吸头，整板384支。这些参数是待验证起点，不能保证挥发性有机溶剂准确度。[液体控制API]({REFS['liquids']})。

API 2.15使用独立温控适配器加载。默认水替代、25 ℃，无化学物料及60 ℃加热。化学分支需要显式审批，封膜/惰性/热压控制在外部完成后才允许60 ℃；温控模块并不提供密封、气体输送、搅拌或反应容器耐压认证。[温控API]({REFS['temperature']})。未将机器人连接至真实设备。厂商模拟器结果另见`results_phase24/opentrons_validation.json`（若未生成则仅有静态验证）。

## 6. HPLC、质量控制与数据资产

模拟检测包括底物、R、S、副产物和内标。0.45 min柱死时间独立于保留因子时标，加入保留时间偏移、流速导致的时间缩放、Gaussian/Lorentzian混合展宽、进样变化、基线与噪声。通过非线性时间/峰宽拟合及非负最小二乘求面积，同时拟合线性基线。不是用模拟真值作为积分答案。

已知相对响应因子与25 nmol内标换算定量。转化率=1−残余底物/初始底物，产物收率=(R+S)/初始底物，ee=(R−S)/(R+S)，三个量分别记录。副产物不并入目标产物。QC检查内标、积分非负、保留区间、拟合残差、R/S分离度与物料衡算。微量对映体LOQ为示例阈值，需真实方法学确定。预设峰库仅覆盖本示例；未知杂质、共洗脱及峰归属需要标准品、手性柱方法与必要的LC-MS核验。UV峰出现顺序不能直接证明绝对构型。

SQLite采用condition_id外键和round/well/provenance唯一约束，重复写入报错，失败QC可入库但不进入训练。原始CSV逐条计算SHA-256，保留原始时间序列和积分参数。模拟与实验来源必须隔离；默认模型查询仅训练模拟来源。当前本地事件不等于联网消息服务，实际HTTP/MQTT仍需鉴权、重试、去重、取消及故障恢复设计。

## 7. 本次可复核结果

| 指标 | 计算值 |
|---|---:|
| 结构/条件数 | 48 / 9,216 |
| 初始+三轮测量 | {m['measurements']} |
| QC通过 | {m['qc_pass_count']} |
| 同时满足Y≥90%、|ee|≥95%的模拟条件 | {m['joint_hits']} |
| 最佳联合条件ID | {best['condition_id']} |
| 该条件产物收率 / |ee| | {best['yield_pct']:.4f}% / {best['ee_abs_pct']:.4f}% |
| 转化率积分RMSE | {m['hplc_rmse_pct'][0]:.5f} 百分点 |
| 产物收率积分RMSE | {m['hplc_rmse_pct'][1]:.5f} 百分点 |
| 有符号ee积分RMSE | {m['hplc_rmse_pct'][2]:.5f} 百分点 |
| 最小R/S分离度 | {m['minimum_resolution']:.3f} |

误差低主要因为合成色谱与积分器共享已知峰模型，不表示真实LC精度。本次初始24个条件中已经有{m['initial_joint_hits']}个达标，因此不能把首次达标归功于主动学习；图3只展示后续最优效用的改善。三轮是否成功仅对固定合成景观与种子成立，未设置保证成功的种子搜索或事后数值替换。应同时查看未达标条件、随机基线、原始色谱及训练ID，不能只摘取最佳结果。

## 8. 实施、学生职责和验收建议

软件验收首先固定SMILES与试剂批次字典，检查单位、孔位、Tip/源液预算、重复任务和异常色谱处理。之后使用水/染料进行称量校准和96孔条码盲测。化学启动前由导师确定一个已有先例的模型反应，并开展空白、无金属、无配体、阳性对照及至少三次独立重复。正式96孔板建议保留8–16孔给对照/重复，将探索条件数相应降低；本次96个独立候选是软件压力测试，不是推荐的最终实验板设计。

本科生适合参与数据清洗、API模拟、校准记录、错误注入测试和结果可视化，工作须经培训和监督。气路、电气、压力、密封及危险化学品授权由专业人员负责。基础准确度建议验收目标（非已达成指标）：水转移偏差≤5%、CV≤3%、孔位追踪100%正确、已知对映体混合物ee偏差≤2个百分点；最终阈值应按仪器量程和研究目的批准。

明确下一阶段缺项：真实substrate/product及方法学、可采购配体及溶解度、材料/气氛/安全认证、源板自动配制、厂商HTTP/MQTT集成、用户权限与备份，以及真实三轮盲测。本原型不宣称AGI、化学发现、无人值守连续运行或临床药物成果。

## 9. 复现与文件索引

```sh
python -m pip install -r requirements_phase24.txt
python run_ai_hts_platform_pilot.py --out .
python run_ai_hts_platform_pilot.py --self-test
# 重复实验请使用新的输出目录，以保护原始数据库
python run_ai_hts_platform_pilot.py --out repeat_seed25 --seed 25
```

主结果：`results_phase24/metrics.json`、`master_measurements.sqlite`、`measurements.json`、`plate_round1/2/3.csv`、`ligands.json`、`conditions.csv`和`chromatograms/`。隐藏真值仅保存在`simulation_truth_not_training.npz`供误差评估。协议为`output_hts_96well_screening.py`（第三轮），各轮协议在结果目录中。图像位于`figures_hts_pilot/`。源文件哈希、依赖版本和产物哈希列于`manifest.json`。
'''
    for name,body in [('GXNU_AI_HTS_PLATFORM_PROPOSAL_ZH.md',zh),
                      ('GXNU_AI_HTS_PLATFORM_PROPOSAL_EN.md',en),
                      ('PHASE24_GXNU_AI_HTS_REPORT_ZH.md',report)]:
        (root/name).write_text(body,encoding='utf-8',newline='\n')


class PilotTests(unittest.TestCase):
    def test_library_smiles(self):
        structures=set()
        for f in range(4):
            for a in ARYLS:
                for s in range(2):
                    mol=Chem.MolFromSmiles(ligand_smiles(f,a,s))
                    self.assertIsNotNone(mol)
                    centers=Chem.FindMolChiralCenters(mol,includeUnassigned=True)
                    self.assertTrue(centers);self.assertNotIn('?',dict(centers).values())
                    structures.add(Chem.MolToSmiles(mol))
        self.assertEqual(len(structures),48)

    def test_plate_and_multichannel_geometry(self):
        names=well_names()
        self.assertEqual(len(set(names)),96)
        self.assertEqual(names[:8],[f'{r}1' for r in 'ABCDEFGH'])
        self.assertEqual(names[-1],'H12')
        self.assertFalse(valid_well('A13')); self.assertFalse(valid_well('I1'))
        self.assertEqual(sum(VOLUMES.values()),80)
        self.assertLessEqual(16+2,20);self.assertLessEqual(44+5,300)
        self.assertEqual(12*3*8,288)

    def test_hplc_recovery(self):
        for truth in ([97.,93.,96.],[65.,42.,-80.],[80.,70.,0.]):
            signal=simulate_trace(truth,91)
            r,_=integrate_trace(TIME,signal)
            self.assertTrue(r['qc_pass'])
            np.testing.assert_allclose([r['conversion_pct'],r['yield_pct'],r['ee_signed_pct']],truth,atol=.15)
            self.assertLessEqual(r['yield_pct'],r['conversion_pct']+.1)

    def test_bad_chromatograms(self):
        with self.assertRaises(ValueError): integrate_trace(TIME,np.zeros(len(TIME)))
        with self.assertRaises(ValueError): integrate_trace(TIME[::-1],np.ones(len(TIME)))
        signal=simulate_trace([95,90,95],91);signal[0]=np.nan
        with self.assertRaises(ValueError): integrate_trace(TIME,signal)

    def test_pareto_and_joint_target(self):
        np.testing.assert_equal(pareto_mask([90,80,95],[96,80,92]),[True,False,True])
        self.assertLess(utility([99],[90])[0],1)
        self.assertGreaterEqual(utility([90],[95])[0],1)

    def test_database_duplicate_and_provenance(self):
        db=init_db(':memory:',[dict(condition_id=0)])
        r=dict(yield_pct=91.,ee_abs_pct=96.,conversion_pct=95.,qc_pass=True)
        ingest(db,0,1,'A1',r,'abc')
        with self.assertRaises(sqlite3.IntegrityError): ingest(db,0,1,'A1',r,'abc')
        ingest(db,0,1,'A1',r,'def','experimental')
        self.assertEqual(len(training_rows(db)),1)
        with self.assertRaises(sqlite3.IntegrityError): ingest(db,99,1,'A2',r,'ghi')
        db.close()

    def test_qei_no_repeat(self):
        X=np.random.default_rng(4).normal(size=(50,4))
        obs=[(i,50+i,70+i) for i in range(10)]
        ids,_,audit=select_batch(X,obs,np.ones(50,bool),8,24)
        self.assertEqual(len(set(ids)),8)
        self.assertFalse(set(ids)&set(range(10)))
        self.assertEqual(audit['training_ids'],list(range(10)))

    def test_protocol_compiles_and_water_default(self):
        text=PROTOCOL_TEMPLATE.replace('__PLATE_MAP__','[]')
        tree=ast.parse(text)
        constants={n.targets[0].id:ast.literal_eval(n.value) for n in tree.body
                   if isinstance(n,ast.Assign) and isinstance(n.value,ast.Constant)}
        self.assertFalse(constants['RUN_CHEMISTRY'])
        self.assertFalse(constants['CHEMISTRY_APPROVED'])

    def test_metrics_are_json_serializable(self):
        r,_=integrate_trace(TIME,simulate_trace([97,93,96],24))
        hits=sum(v['yield_pct']>=90 and v['ee_abs_pct']>=95 for v in [r])
        json.dumps(dict(joint_hits=hits,record=r),allow_nan=False)

    def test_protocol_rejects_misaligned_source(self):
        plate=[dict(well=w,total_uL=80,source_catalyst_plate=w,
                    source_base_plate=w,source_solvent_plate=w,**VOLUMES) for w in well_names()]
        plate[0]['source_base_plate']='B1'
        with self.assertRaises(ValueError): generate_protocol('unused.py',plate)


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out',type=Path,default=Path('.'))
    ap.add_argument('--seed',type=int,default=24)
    ap.add_argument('--pptx',type=Path)
    ap.add_argument('--self-test',action='store_true')
    args=ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO,format='[phase24 %(asctime)s] %(message)s',datefmt='%H:%M:%S')
    if args.self_test:
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(PilotTests)
        result=unittest.TextTestRunner(verbosity=2).run(suite)
        return 0 if result.wasSuccessful() else 1
    if not 0<=args.seed<2**30:
        ap.error('--seed must lie in [0, 2**30)')
    start=time.perf_counter(); root=args.out.resolve(); root.mkdir(parents=True,exist_ok=True)
    res=root/'results_phase24';res.mkdir(exist_ok=True)
    if (res/'master_measurements.sqlite').exists():
        raise FileExistsError('This output directory already contains a database. Use a new --out.')
    dump(res/'pptx_context.json',ppt_context(args.pptx))
    LOG.info('SIMULATION ONLY: no robot connection or real chemistry execution')
    ligands,fps=build_ligands(args.seed); dump(res/'ligands.json',ligands)
    rows,X,envelope=build_space(ligands);write_csv(res/'conditions.csv',rows)
    assert len(rows)==9216
    db=init_db(res/'master_measurements.sqlite',rows)
    truths=synthetic_oracle(rows)
    np.savez_compressed(res/'simulation_truth_not_training.npz',condition_id=np.arange(len(rows)),
                        conversion_yield_ee=truths)
    seeds=maxmin_seeds(X,rows,fps,envelope,24,args.seed)
    plate=make_plate(seeds,rows,{},0);write_csv(res/'plate_round0.csv',plate)
    all_records,_,_,errors=run_measurements(db,plate,truths,root,args.seed)
    error_blocks=[errors];round_ids=[];events=[]
    for r in range(1,4):
        obs=training_rows(db)
        if len(obs)<10: raise RuntimeError('Insufficient QC-passing initial measurements')
        ids,predictions,audit=select_batch(X,obs,envelope,96,args.seed+r)
        plate=make_plate(ids,rows,predictions,r);write_csv(res/f'plate_round{r}.csv',plate)
        dump(res/f'task_round{r}.json',dict(schema_version='1.0',task_id=f'phase24-seed{args.seed}-round{r}',
            state='simulation_only',chemistry_approved=False,well_volume_unit='uL',
            stock_concentration_unit='mM',inventory_lot_validation='required before physical work',
            substrate_identity='unassigned synthetic benchmark; supervisor must define',
            source_preparation='external validated pre-arraying required',wells=plate))
        dump(res/f'acquisition_round{r}.json',dict(**audit,predictions=predictions,selected_ids=ids.tolist()))
        generate_protocol(res/f'ot2_round{r}.py',plate)
        records,traces,fits,errors=run_measurements(db,plate,truths,root,args.seed)
        all_records.extend(records);error_blocks.append(errors);round_ids.append(ids)
        events.append(dict(task_id=f'phase24-seed{args.seed}-round{r}',
                           event='analysis.completed',transport='local JSON; HTTP/MQTT adapter contract only',
                           provenance='simulated',n_conditions=96,n_qc_pass=sum(x['qc_pass'] for x in records)))
        good=[x for x in all_records if x['qc_pass']]
        hits=sum(x['yield_pct']>=90 and x['ee_abs_pct']>=95 for x in good)
        LOG.info('Round %d complete: total %d; QC %d; joint target hits %d',r,len(all_records),len(good),hits)
    generate_protocol(root/'output_hts_96well_screening.py',plate)
    dump(res/'measurements.json',all_records);dump(res/'workflow_events.json',events)
    good=[x for x in all_records if x['qc_pass']]
    best=max(good,key=lambda x:utility(x['yield_pct'],x['ee_abs_pct']))
    baseline=[]
    for s in range(5):
        remaining=np.setdiff1d(np.flatnonzero(envelope),seeds)
        random_ids=np.random.default_rng(args.seed+100+s).choice(remaining,288,replace=False)
        seen=np.r_[seeds,random_ids]
        baseline.append([float(utility(truths[seen[:n],1],np.abs(truths[seen[:n],2])).max()) for n in [24,120,216,312]])
    metrics=dict(phase=24,seed=args.seed,provenance='simulated',ligands=48,conditions=9216,
                 demo_envelope_count=int(envelope.sum()),measurements=len(all_records),
                 qc_pass_count=len(good),joint_hits=int(sum(x['yield_pct']>=90 and x['ee_abs_pct']>=95 for x in good)),
                 initial_joint_hits=int(sum(x['round']==0 and x['yield_pct']>=90 and x['ee_abs_pct']>=95 for x in good)),
                 best_joint_condition={k:best[k] for k in ['condition_id','round','well','yield_pct','ee_abs_pct','conversion_pct']},
                 hplc_rmse_pct=np.sqrt(np.mean(np.vstack(error_blocks)**2,axis=0)).tolist(),
                 minimum_resolution=min(x['resolution_R_S'] for x in good),
                 random_baseline_utility=baseline,initial_condition_ids=seeds.tolist(),
                 physical_hardware_validated=False,chemistry_experiment_performed=False,
                 protocol_tips=dict(P20=288,P300=96),per_well_uL=80.)
    metrics['tsne_KL']=figures(root,X,envelope,seeds,round_ids,all_records,predictions,traces,fits,baseline)
    metrics['wall_s']=time.perf_counter()-start
    dump(res/'metrics.json',metrics);briefs_and_report(root,metrics)
    db.close()
    outputs=[p for folder in [res,root/'figures_hts_pilot'] for p in folder.rglob('*') if p.is_file()]
    outputs += [root/x for x in ['output_hts_96well_screening.py','GXNU_AI_HTS_PLATFORM_PROPOSAL_ZH.md',
                                'GXNU_AI_HTS_PLATFORM_PROPOSAL_EN.md','PHASE24_GXNU_AI_HTS_REPORT_ZH.md']]
    versions={x:importlib.metadata.version(x) for x in ['numpy','scipy','matplotlib','scikit-learn']}
    versions['rdkit']=rdBase.rdkitVersion
    dump(res/'manifest.json',dict(seed=args.seed,python=sys.version,versions=versions,
        script_sha256=sha(__file__),pptx_sha256=PPT_SHA,
        outputs={str(p.relative_to(root)).replace('\\','/'):sha(p) for p in outputs}))
    LOG.info('COMPLETE: %s',json.dumps(metrics,ensure_ascii=False))
    return 0


if __name__=='__main__':
    sys.exit(main())
