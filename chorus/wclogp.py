#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import os
import yaml

""" Import Wildman-Crippen logP parameters """
with open(os.path.join(os.path.dirname(__file__), "wclogp.yaml")) as file:
    DATA = yaml.load(file.read())


def wclogp(mol):
    return wcparams(mol, "logP")


def wcmr(mol):
    return wcparams(mol, "MR")


def wcparams(mol, type_):
    """
    Calculate Wildman-Crippen logP
    Wildman S., Crippen G., Prediction of Physicochemical Parameters by Atomic
    Contribution, J. Chem. Inf. Model. 39 (1999) 868-873
    """
    try:
        assign_wctype(mol)
    except:
        return "N/A"
    scores = []
    for i, atom in mol.atoms_iter():
        scores.append(float(DATA[type_][atom.wctype]))
        if atom.H_count:
            if atom.symbol == "C":
                scores.append(DATA[type_]["H1"] * atom.H_count)
            elif atom.symbol == "N":
                scores.append(DATA[type_]["H3"] * atom.H_count)
            elif atom.wctype == "O2a":
                scores.append(DATA[type_]["H4"] * atom.H_count)
            else:
                scores.append(DATA[type_]["H2"] * atom.H_count)
    return round(sum(scores), 2)


def assign_wctype(mol):
    """ Assign atom type of Wildman-Crippen logP (Wildman and Crippen 1999) """
    mol.require("Aromatic")
    p_block = ["Al", "B", "Si", "Ga", "Ge", "As", "Se", "Sn", "Te", "Pb",
               "Ne", "Ar", "Kr", "Xe", "Rn"]
    d_block = ["Fe", "Co", "Cu", "Zn", "Tc", "Cd", "Pt", "Au", "Hg", "Gd"]
    for i, atom in mol.atoms_iter():
        nbrs = mol.neighbors(i)
        nbrcnt = len(nbrs)
        nbratms = [mol.atom(i) for i, _ in nbrs.items()]
        isnbrarm = any([mol.atom(i).aromatic for i, _ in nbrs.items()])

        if atom.symbol == "C":
            if atom.aromatic:  # Aromatic carbon
                ex = [(mol.atom(i), b) for i, b in nbrs.items()
                      if not b.aromatic]
                if ex:
                    sbst, sbstb = ex[0]
                    cnosx = {"C": "C21", "N": "C22", "O": "C23", "S": "C24",
                             "F": "C14", "Cl": "C15", "Br": "C16", "I": "C17"}
                    if sbst.symbol not in cnosx:
                        atom.wctype = "C13"
                    elif sbst.aromatic:
                        atom.wctype = "C20"
                    elif sbst.symbol in ("C", "N", "O") and sbstb.order == 2:
                        atom.wctype = "C25"
                    else:
                        atom.wctype = cnosx[sbst.symbol]
                else:
                    if nbrcnt == 2:  # No substituent
                        atom.wctype = "C18"
                    elif nbrcnt == 3:  # Bridgehead
                        atom.wctype = "C19"
            elif atom.pi == 0:  # Aliphatic carbon
                hcnopsx = ["H", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
                if not all([a.symbol in hcnopsx for a in nbratms]):
                    atom.wctype = "C27"
                elif isnbrarm:  # Adjacent to aromatic
                    if atom.H_count == 3:
                        if nbratms[0].symbol == "C":
                            atom.wctype = "C8"
                        else:
                            atom.wctype = "C9"
                    else:
                        numh = {0: "C12", 1: "C11", 2: "C10"}
                        atom.wctype = numh[atom.H_count]
                elif all([a.symbol == "C" for a in nbratms]):
                    if atom.H_count >= 2:
                        atom.wctype = "C1"
                    else:
                        atom.wctype = "C2"
                else:
                    if atom.H_count >= 2:
                        atom.wctype = "C3"
                    else:
                        atom.wctype = "C4"
            elif atom.pi == 1:
                dbi = [i for i, b in nbrs.items() if b.order == 2][0]
                dba = mol.atom(dbi)
                if dba.symbol == "C":
                    if isnbrarm:  # C=C adjacent to aromatic
                        atom.wctype = "C26"
                    else:  # C=C
                        atom.wctype = "C6"
                else:  # C=(!C)
                    atom.wctype = "C5"
                    # Overwrite carbonyl oxygens
                    if dba.symbol == "O":
                        if all(a.symbol != "C" for a in nbratms):
                            mol.atom(dbi).wctype = "O11"
                        elif any(a.aromatic for a in nbratms):
                            mol.atom(dbi).wctype = "O10"
                        else:
                            mol.atom(dbi).wctype = "O9"
            elif atom.pi == 2:
                od = max([b.order for _, b in nbrs.items()])
                if od == 3:  # Alkyne
                    atom.wctype = "C7"
                elif od == 2:  # Allene
                    atom.wctype = "C6"
            if atom.wctype is None:
                atom.wctype = "CS"

        elif atom.symbol == "N":
            if atom.charge == 0:
                if atom.aromatic:
                    atom.wctype = "N11"
                elif atom.pi == 0:  # amine
                    amines = {
                        0: {0: "N7", 1: "N2", 2: "N1"},
                        1: {0: "N8", 1: "N4", 2: "N3"}
                    }
                    atom.wctype = amines[isnbrarm][atom.H_count]
                elif atom.pi == 1:  # imine
                    if atom.H_count:
                        atom.wctype = "N5"
                    else:
                        atom.wctype = "N6"
                elif atom.pi == 2:  # nitrile
                    atom.wctype = "N9"
            elif atom.charge > 0:
                if atom.aromatic:
                    atom.wctype = "N12"
                elif atom.pi == 0:  # ammonium
                    if atom.H_count:
                        atom.wctype = "N10"
                    else:
                        atom.wctype = "N13"
                elif atom.pi == 1:  # iminium
                    atom.wctype = "N13"
                elif atom.pi == 2:
                    if any(a.charge < 1 for a in nbratms):  # 1,3-dipole
                        if all(a.symbol == "N" for a in nbratms):  # azide
                            atom.wctype = "N14"
                        else:  # diazo
                            atom.wctype = "N13"
                    else:  # nitrilium
                        atom.wctype = "N14"
            else:  # nitride
                atom.wctype = "N14"
            if atom.wctype is None:
                atom.wctype = "NS"

        elif atom.symbol == "O":
            if atom.wctype:  # skip if already assigned by carbonyl C
                continue
            if atom.aromatic:
                atom.wctype = "O1"
            elif atom.pi == 0 and atom.charge == 0:
                if atom.H_count == 1:
                    if nbratms[0].symbol in ("O", "S"):  # ROOH, RSOH
                        atom.wctype = "O2a"
                    elif nbratms[0].carbonyl_C:  # Carboxylic acid
                        atom.wctype = "O2a"
                    else:  # Alcohol
                        atom.wctype = "O2"
                elif atom.H_count == 2:  # Water
                    atom.wctype = "O2"
                elif isnbrarm:  # Ether adjacent to aromatic
                    atom.wctype = "O4"
                else:  # Ether
                    atom.wctype = "O3"
            elif nbrcnt == 1:
                nbr = nbratms[0]
                if nbr.aromatic:  # O=Aromatic
                    atom.wctype = "O8"
                elif nbr.symbol in ("N", "O"):  # N-oxide, O2, NO
                    atom.wctype = "O5"
                elif nbr.symbol == "S":  # S-oxide
                    atom.wctype = "O6"
                elif nbr.carbonyl_C and atom.charge < 0:  # Carboxylate anion
                    atom.wctype = "O12"
                else:  # ex. Phosphine oxide
                    atom.wctype = "O7"
            if atom.wctype is None:
                atom.wctype = "OS"

        elif atom.symbol in ("F", "Cl", "Br", "I"):
            if atom.charge:
                atom.wctype = "Hal"
            else:
                atom.wctype = atom.symbol
        elif atom.symbol == "P":
            atom.wctype = "P"
        elif atom.symbol == "S":
            if atom.aromatic:
                atom.wctype = "S3"
            elif atom.charge:
                atom.wctype = "S2"
            else:
                atom.wctype = "S1"
        elif atom.symbol in p_block:
            atom.wctype = "Me1"
        elif atom.symbol in d_block:
            atom.wctype = "Me2"
        elif atom.symbol == "H":  # Explicit hydrogen
            nbr = nbratms[0]
            if nbr.symbol == "C":
                atom.wctype = "H1"
            elif nbr.symbol == "N":
                atom.wctype = "H3"
            else:
                atom.wctype = "H2"
        else:
            raise ValueError("Undefined atom symbol: {}".format(atom.symbol))
    mol.descriptors.add("Wildman-Crippen")
