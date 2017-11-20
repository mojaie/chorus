#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

def assign_type(mol, force_recalc=False):
    """ PATTY [Bush et al. J. Inf. Comput. Sci 33 (1993) 756-762]
    TODO: not yet implemented

    1:cation 2:anion 3:donor 4:acceptor
    5:polar 6:hydrophobe 7:others
    """
    if "PATTY" in mol.descriptors and not force_recalc:
        return
    mol.require("Phys_charge")
    for i, atom in mol.atoms_iter():
        # default is 7 (others)
        nbrcnt = mol.neighbor_count(i)
        if atom.charge > 0 or atom.charge_phys > 0 or \
                atom.charge_conj > 0 and not atom.n_oxide:
            atom.type = 1  # cation
        elif atom.charge < 0 or atom.charge_phys < 0 or \
                atom.charge_conj < 0 and not atom.n_oxide:
            atom.type = 2  # anion
        elif atom.symbol == "N":
            if nbrcnt in (1, 2):
                if atom.pi == 2:
                    atom.type = 3  # donor
                elif atom.pi == 1:
                    atom.type = 4  # acceptor
        elif atom.symbol == "O":
            if nbrcnt == 1 and not atom.pi:
                atom.type = 5  # polar
            else:
                atom.type = 4  # acceptor
        elif atom.symbol in ("C", "Si", "S", "Se", "P", "As"):
            ewg = False
            for n, bond in mol.neighbors(i).items():
                natom = mol.atom(n)
                if natom.symbol in ("N", "O", "S") and atom.pi \
                        and not (natom.pi == 2 and mol.neighbor_count(n) == 3):
                    # the sp2 adjacent to neg (but not conj tert amine) is 7
                    ewg = True
                    break
            if not ewg:
                atom.type = 6  # hydrophobes
        elif atom.symbol in ("F", "Cl", "Br", "I") and nbrcnt == 1:
            atom.type = 6  # typical halogens are hydrophobic
    mol.descriptors.add("PATTY")


# TODO:
# Conjugated amine is not h_acceptor ?
# if atom.symbol == "N" and not atom.pi:
#    if any(mol.atom(nbr).pi for nbr in nbrs.keys()):
#        mol.atom(i).h_acceptor = False
# N of triphenylamine is sp2 and thus is not H acceptor
# Diphenylamine is very weak base (pKb=0.79)
# Phenylamine is as strong base as pyridine (pKb=4-5)
# Conjugated O can be sp2? (ex. O in diphenylether is H acceptor?)


# TODO: -O-, -S-, =O, =S (H-acceptor)
# TODO: -OH, -SH (H-donor, H-acceptor)
# TODO: non-conjugated -NH2, -NH- (H-donor, H-acceptor)
# TODO: non-conjugated -N< (H-acceptor)
# TODO: conjugated -NH2, -NH- (H-donor)
# TODO: conjugated -N< (None)
# TODO: primary imine(H-donor, H-acceptor)
# TODO: secondary imine, aromatic N (H-acceptor)
# TODO: halogen(wildcard: skip evaluation)

# TODO: nitro, diazo, nitrile, isocyan, azide, nitrile oxide(center: N)


"""

- assign_pi
double bond -> atom.pi = 1
triple bond -> atom.pi = 2
=C= -> atom.pi = 2

- assign_valence
N, O with hydrogen -> atom.h_donor = True
F, O, N (except sp3 N adj to multiple bond) -> atom.h_acceptor = True

- assign_aromatic (not yet implemented)
double_bond -> +2, pi=0 and lone_pair -> +2
pi=0 -> abort(not aromatic)
4n+2: aromatic=True
TODO: dipole=-1 -> +2

- electron withdrawing group (dipole=-1) (not yet implemented)
>C=O, >C=S, >S=O, >P=O, -NO2, -C#N
ammonium, sulfonate(-M)
p-nitrile, p-pyridyl(-M)
trifluoro, trichrolo(-E?)
electron donating group (dipole=+1)
lone_pair
-I?
amine


- assign_charge (not yet implemented)
assign charge in physiological condition
evaluate atom.charge_phys, atom.h_donor, atom.h_acceptor, atom.pi
non-conjugated amine (N and lone_pair) -> cation
amidine, guanidine (=NH adj to dipole=1) -> cation
oxoacids C(=O)O, P(=O)(O)O (OH adj to dipole=-1) -> anion
thiophenol (S adj to aromatic) -> anion
1,3-dipole (h_acceptor, none, h_acceptor)
X=N=X formula(wrong) -> ok
X=[N+]-[X-] formula -> neutralize charge_phys

# nitro -[N+]([O-])=O  -[N+](=O)[O-]
# azide -[N-]-[N+]#N  -N=[N+]=[N-]
# diazo　>[C-]-[N+]#N  >C=[N+]=[N-]
# nitrile oxide -[C-]=[N+]=O  -C#[N+][O-]
N, S, P, As-oxide
∋N=O(wrong) -> ok
∋[N+][O-] -> neutralize charge_phys
∋P=O -> ok
∋[P+][O-] -> neutralize charge_phys

# N-oxide(quart amine but one nbr is O-)
# P,As-oxide(pi=1 and nbr=4 but pi nbr is O-)

special case
6-membered ring interaction
# 2-(pyridin-2-yl)acetaldehyde(pyrophthalone)
# 1,3-diketone (acetylaceton)
# TODO: sulfoneamide, phosphoneamide (SorP=N)
# TODO: imidazole,  diaminopyrimidine (basic hetero ring)
# TODO: tetrazole, thiazole, muscimol (acidic hetero rings)
"""
