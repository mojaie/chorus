#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

from collections import Counter
import itertools
import pickle


from chorus.model.graphmol import Compound
from chorus.model.atom import atom_number
from chorus import descriptor, remover, topology, wclogp


def clone(mol):
    """ Returns deepcopy of the object.
    Pickle and unpickle is faster than deepcopy.
    TODO: this query function should be renamed to "cloned"
    """
    return pickle.loads(pickle.dumps(mol, protocol=4))


def assign_descriptors(mol):
    topology.recognize(mol)
    topology.minify_ring(mol)
    descriptor.assign_valence(mol)
    descriptor.assign_rotatable(mol)
    descriptor.assign_aromatic(mol)


def make_Hs_implicit(original_mol, keep_stereo=True):
    """Return molecule that explicit hydrogens removed
    TODO: this query function should be renamed to "explicitHs_removed"
    """
    mol = clone(original_mol)
    mol.descriptors.clear()  # Reset descriptor
    to_remove = set()
    for i, nbrs in mol.neighbors_iter():
        if mol.atom(i).symbol == "H":  # do not remove H2
            continue
        for nbr, bond in nbrs.items():
            if mol.atom(nbr).symbol == "H" and \
                    not (bond.order == 1 and bond.type and keep_stereo):
                to_remove.add(nbr)
    for r in to_remove:
        mol.remove_atom(r)
    assign_descriptors(mol)
    return mol


def null_molecule(assign_descs=True):
    mol = Compound()
    if assign_descs:
        assign_descriptors(mol)
    return mol


def largest_graph(mol):
    """Return a molecule which has largest graph in the compound
    Passing single molecule object will results as same as molutil.clone
    """
    mol.require("Valence")
    mol.require("Topology")
    m = clone(mol)  # Avoid modification of original object
    if m.isolated:
        for k in itertools.chain.from_iterable(m.isolated):
            m.remove_atom(k)
    return m


def mw(mol, ndigits=2):
    """Return standard molecular weight
    :param ndigits: number of digits
    """
    mol.require("Valence")
    return round(sum(a.mw() for _, a in mol.atoms_iter()), ndigits)


def mw_wo_sw(mol, ndigits=2):
    """Molecular weight without salt and water
    :param ndigits: number of digits
    """
    cp = clone(mol)  # Avoid modification of original object
    remover.remove_water(cp)
    remover.remove_salt(cp)
    return round(sum(a.mw() for _, a in cp.atoms_iter()), ndigits)


def charge(mol):
    """Total charge"""
    return sum(a.charge for _, a in mol.atoms_iter())


def non_hydrogen_count(mol):
    """Non-hydrogen atom count """
    return sum(1 for _, a in mol.atoms_iter() if a.symbol != "H")


def H_donor_count(mol):
    """Hydrogen bond donor count """
    mol.require("Valence")
    return sum(1 for _, a in mol.atoms_iter() if a.H_donor)


def H_acceptor_count(mol):
    """Hydrogen bond acceptor count """
    mol.require("Valence")
    return sum(1 for _, a in mol.atoms_iter() if a.H_acceptor)


def rotatable_count(mol):
    """Rotatable bond count """
    mol.require("Rotatable")
    return sum(1 for _, _, b in mol.bonds_iter() if b.rotatable)


def rule_of_five_violation(mol):
    """Lipinski's rule of five violation count """
    v = 0
    if mw(mol) > 500:
        v += 1
    if H_donor_count(mol) > 5:
        v += 1
    if H_acceptor_count(mol) > 10:
        v += 1
    try:
        if wclogp.wclogp(mol) > 5:
            v += 1
    except TypeError:  # N/A
        v += 1
    return v


def composition(mol):
    """Molecular composition in dict format
    (ex. Glucose {'C': 6, 'H': 12, 'O': 6}).
    """
    mol.require("Valence")
    c = Counter()
    for _, a in mol.atoms_iter():
        c += a.composition()
    return c


def available_key(mol):
    """Eldest key +1 """
    return max(mol.key_set() | {1})


def mols_iter(mol):
    keyset = mol.key_set()
    for iso in mol.isolated:
        keyset -= set(iso)
        yield iso
    yield list(keyset)


def formula(mol):
    """Chemical formula.
    Atoms should be arranged in order of C, H and other atoms.
    Molecules should be arranged in order of length of formula text.
    """
    mol.require("Valence")
    mol.require("Topology")
    total_cntr = Counter()
    for m in sorted(mols_iter(mol), key=len, reverse=True):
        cntr = Counter()
        for i in m:
            cntr += mol.atom(i).composition()
        text = []
        Cs = cntr.pop("C", 0)
        if Cs:
            text.append("C")
            if Cs > 1:
                text.append(str(Cs))
        Hs = cntr.pop("H", 0)
        if Hs:
            text.append("H")
            if Hs > 1:
                text.append(str(Hs))
        heteros = sorted(cntr.items(), key=lambda x: atom_number(x[0]))
        for k, v in heteros:
            text.append(k)
            if v > 1:
                text.append(str(v))
        total_cntr["".join(text)] += 1
    total = sorted(total_cntr.items(), key=lambda x: len(x[0]), reverse=True)
    total_text = []
    for k, v in total:
        if v > 1:
            total_text.append(str(v) + k)
        else:
            total_text.append(k)
    return ".".join(total_text)
