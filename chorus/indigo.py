#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

from indigo import Indigo

from chorus.model.graphmol import Compound
from chorus.model.atom import Atom
from chorus.model.bond import Bond
from chorus import molutil

idg = Indigo()
# TODO: stereo, aromatic etc


def to_real_mol(mol):
    """Convert molecule to indigo real molecule"""
    real = idg.createMolecule()
    return indigo_mol(mol, real)


def to_query_mol(mol):
    """Convert molecule to indigo real molecule"""
    q = idg.createQueryMolecule()
    return indigo_mol(mol, q)


def indigo_mol(mol, idgmol):
    key_to_obj = {}
    for k, a in mol.atoms_iter():
        aobj = idgmol.addAtom(a.symbol)
        if a.coords is not None:
            crds = list(a.coords) + [0] * (3 - len(a.coords))
            aobj.setXYZ(*crds)
        key_to_obj[k] = aobj
    for u, v, b in mol.bonds_iter():
        uobj = key_to_obj[u]
        vobj = key_to_obj[v]
        uobj.addBond(vobj, b.order)
    return idgmol


def from_indigo(idgmol, assign_descriptor=True):
    """Convert RDMol to molecule"""
    mol = Compound()
    for atom in idgmol.iterateAtoms():
        key = atom.index()
        a = Atom(atom.symbol())
        a.coords = list(atom.xyz())
        mol.add_atom(key, a)
    for bond in idgmol.iterateBonds():
        u = bond.source()
        v = bond.destination()
        b = Bond()
        b.order = int(bond.bondOrder())
        mol.add_bond(u.index(), v.index(), b)
    if assign_descriptor:
        molutil.assign_descriptors(mol)
    return mol


def fingerprint_similarity(mol1, mol2):
    """Calculate Indigo fingerprint similarity
    """
    idmol1 = to_real_mol(mol1)
    idmol2 = to_real_mol(mol2)
    fp1 = idmol1.fingerprint("sim")
    fp2 = idmol2.fingerprint("sim")
    return round(idg.similarity(fp1, fp2, "tanimoto"), 2)


def fingerprint_dist(mol1, mol2):
    return round(1 - fingerprint_similarity(mol1, mol2), 2)
