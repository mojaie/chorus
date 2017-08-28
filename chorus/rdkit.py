#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit import DataStructs

from chorus.model.graphmol import Compound
from chorus.model.atom import Atom, atom_number
from chorus.model.bond import Bond
from chorus import molutil

# TODO: stereo, aromatic etc


def to_rdmol(mol):
    """Convert molecule to RDMol"""
    rwmol = Chem.RWMol(Chem.MolFromSmiles(''))
    key_to_idx = {}
    bond_type = {1: Chem.BondType.SINGLE,
                 2: Chem.BondType.DOUBLE,
                 3: Chem.BondType.TRIPLE}
    conf = Chem.Conformer(rwmol.GetNumAtoms())
    for k, a in mol.atoms_iter():
        i = rwmol.AddAtom(Chem.Atom(atom_number(a.symbol)))
        key_to_idx[k] = i
        conf.SetAtomPosition(i, a.coords)
    rwmol.AddConformer(conf)
    for u, v, b in mol.bonds_iter():
        ui = key_to_idx[u]
        vi = key_to_idx[v]
        rwmol.AddBond(ui, vi, bond_type[b.order])
    Chem.GetSSSR(rwmol)  # Ring recognition is required for fingerprint
    rwmol.UpdatePropertyCache(strict=False)
    return rwmol.GetMol()


def from_rdmol(rdmol, assign_descriptor=True):
    """Convert RDMol to molecule"""
    mol = Compound()
    conf = rdmol.GetConformer()
    Chem.Kekulize(rdmol)
    for atom in rdmol.GetAtoms():
        key = atom.GetIdx()
        a = Atom(atom.GetSymbol())
        a.coords = conf.GetAtomPosition(key)
        mol.add_atom(key, a)
    for bond in rdmol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()
        b = Bond()
        b.order = int(bond.GetBondTypeAsDouble())
        mol.add_bond(u, v, b)
    if assign_descriptor:
        molutil.assign_descriptors(mol)
    return mol


def morgan_sim(mol1, mol2, radius=2, digit=3):
    """Calculate morgan fingerprint similarity by using RDKit
    radius=2 roughly equivalent to ECFP4
    """
    rdmol1 = to_rdmol(mol1)
    rdmol2 = to_rdmol(mol2)
    fp1 = AllChem.GetMorganFingerprint(rdmol1, radius)
    fp2 = AllChem.GetMorganFingerprint(rdmol2, radius)
    return round(DataStructs.DiceSimilarity(fp1, fp2), digit)


def morgan_dist(mol1, mol2, radius=2, digit=3):
    return round(1 - morgan_sim(mol1, mol2, radius, digit), digit)


def fmcs(mol1, mol2, timeout=2, digit=3):
    rdmol1 = to_rdmol(mol1)
    rdmol2 = to_rdmol(mol2)
    edges1 = rdmol1.GetNumBonds()
    edges2 = rdmol2.GetNumBonds()
    mcs = rdFMCS.FindMCS((rdmol1, rdmol2), timeout=timeout)
    # Jaccard-Tanimoto coefficient
    try:
        sim = round(mcs.numBonds / (edges1 + edges2 - mcs.numBonds), digit)
    except ZeroDivisionError:
        sim = 0
    return {
        "mcs_edges": mcs.numBonds,
        "similarity": sim,
        "canceled": mcs.canceled
    }
