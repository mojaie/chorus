#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#


def remove_salt(mol):
    to_be_removed = []
    for n, atom in mol.atoms_iter():
        if not mol.neighbor_count(n) and atom.charge:
            to_be_removed.append(n)
    for n in to_be_removed:
        mol.remove_atom(n)


def remove_water(mol):
    to_be_removed = []
    for n, atom in mol.atoms_iter():
        if not mol.neighbor_count(n) and atom.symbol == "O":
            to_be_removed.append(n)
    for n in to_be_removed:
        mol.remove_atom(n)


def remove_coordinated_metal(mol):
    # TODO: should be renamed to "remove_heavy_metals"
    heavy_metals = ("Fe", "Co")  # TODO: list transition metals
    to_be_removed = []
    for n, atom in mol.atoms_iter():
        if mol.neighbor_count(n) > 4 and atom.symbol in heavy_metals:
            to_be_removed.append(n)
    for n in to_be_removed:
        mol.remove_atom(n)


def remove_sugar(mol):
    pass
