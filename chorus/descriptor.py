#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

from chorus.util import iterator


def assign_valence(mol):
    """Assign pi electron and hydrogens"""
    for u, v, bond in mol.bonds_iter():
        if bond.order == 2:
            mol.atom(u).pi = 1
            mol.atom(v).pi = 1
            if mol.atom(u).symbol == "O" and not mol.atom(u).charge:
                mol.atom(v).carbonyl_C = 1
            if mol.atom(v).symbol == "O" and not mol.atom(v).charge:
                mol.atom(u).carbonyl_C = 1
        elif bond.order == 3:
            mol.atom(u).pi = mol.atom(v).pi = 2
    max_nbr = {"C": 4, "Si": 4, "N": 3, "P": 3, "As": 3,
               "O": 2, "S": 2, "Se": 2, "F": 1, "Cl": 1, "Br": 1, "I": 1}
    for i, nbrs in mol.neighbors_iter():
        atom = mol.atom(i)
        if len(nbrs) == 2 and all(bond.order == 2 for bond in nbrs.values()):
            atom.pi = 2  # sp (allene, ketene)
        if atom.symbol in max_nbr:
            h_cnt = max_nbr[atom.symbol] - len(nbrs) - atom.pi + atom.charge
            if h_cnt > 0:
                mol.atom(i).add_hydrogen(h_cnt)
    mol.descriptors.add("Valence")


def assign_rotatable(mol):
    """Assign rotatable bonds (single and not terminal and not in ring)"""
    mol.require("Valence")
    mol.require("Topology")
    bs = set([(u, v) for u, v, _ in mol.bonds_iter()])
    for r in mol.rings:
        for a, b in iterator.consecutive(r + [r[0]], 2):
            bs -= {(a, b)}
            bs -= {(b, a)}
    for u, v in bs:
        unc = mol.neighbor_count(u)
        vnc = mol.neighbor_count(v)
        bc = mol.bond(u, v).order
        if bc == 1 and (unc > 1 and vnc > 1):
            mol.bond(u, v).rotatable = True
    mol.descriptors.add("Rotatable")


def assign_aromatic(mol):
    """Assign aromatic ring
    sp2 atom:
    pi=1 -> +1
    N, O, S, C- -> +2
    >C=O, B, C+ -> 0
    sp3 atom -> not aromatic
    sum of the score satisfies 4n+2 -> aromatic
    """
    mol.require("Valence")
    mol.require("MinifiedRing")
    for ring in mol.rings:
        pi_cnt = 0
        for r in ring:
            if mol.atom(r).pi == 0:
                if mol.atom(r).symbol == "C":
                    if mol.atom(r).charge == 1:
                        pass
                    elif mol.atom(r).charge == -1:
                        pi_cnt += 2
                    else:
                        break
                elif mol.atom(r).charge == 0:
                    if mol.atom(r).symbol in ("N", "O", "S"):
                        pi_cnt += 2
                    elif mol.atom(r).symbol == "B":
                        pass
                    else:
                        break
                else:
                    break
            elif mol.atom(r).pi == 1:
                if mol.atom(r).carbonyl_C:
                    pass
                else:
                    pi_cnt += 1
            else:
                break
        else:
            if pi_cnt % 4 == 2:
                for u, v in iterator.consecutive(ring + [ring[0]], 2):
                    mol.atom(u).aromatic = True
                    mol.bond(u, v).aromatic = True
    mol.descriptors.add("Aromatic")


def assign_charge(mol, force_recalc=False):
    """Assign charges in physiological condition"""
    # TODO: not implemented yet
    mol.require("Aromatic")
    for i, nbrs in mol.neighbors_iter():
        atom = mol.atom(i)
        nbrcnt = len(nbrs)
        if atom.symbol == "N":
            if not atom.pi:
                # non-conjugated amines are anion
                mol.atom(i).charge_phys = 1
            elif nbrcnt == 1 and atom.pi == 2:
                # amidine, guanidine are conjugated cation
                ni = list(nbrs.keys())[0]
                conj = False
                sp2n = None
                for nni, nnb in mol.neighbors(ni).items():
                    if mol.atom(nni).symbol == "N" and nnb.order == 2 \
                            and not mol.atom(nni).aromatic:
                        mol.atom(nni).charge_conj = 1
                        conj = True
                    elif mol.atom(nni).symbol == "N" and nni != i:
                        sp2n = nni
                if conj:
                    mol.atom(i).charge_phys = 1
                    if sp2n is not None:
                        mol.atom(sp2n).charge_conj = 1
        elif atom.symbol == "O" and nbrcnt == 1 and atom.pi == 2:
            # oxoacid are conjugated anion
            ni = list(nbrs.keys())[0]
            conj = False
            if mol.atom(ni).symbol == "N":
                mol.atom(i).n_oxide = True
                mol.atom(ni).n_oxide = True
            for nni, nnb in mol.neighbors(ni).items():
                if mol.atom(nni).symbol in ("O", "S") \
                        and nnb.order == 2 and not mol.atom(ni).n_oxide:
                    mol.atom(nni).charge_conj = -1
                    conj = True
            if conj:
                mol.atom(i).charge_phys = -1
        elif atom.symbol == "S" and nbrcnt == 1:
            # thiophenols are anion
            ni = list(nbrs.keys())[0]
            if mol.atom(ni).aromatic:
                mol.atom(i).charge_phys = -1
    mol.charge_assigned = True
    mol.descriptors.add("Phys_charge")
