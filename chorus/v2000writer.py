#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import chorus


def atom_block(mol):
    idx_table = {}
    lines = []
    for i, (k, atom) in enumerate(mol.atoms_iter()):
        atom_line = []
        for ax in atom.coords:
            int_, frac = str(ax).split(".")
            atom_line.append("{:>5}".format(int_))
            atom_line.append("{:0<4}".format(frac))
        atom_line.append("{:<3}".format(atom.symbol))
        lines.append(
            "{}.{}{}.{}{}.{} {} 0  0  0  0  0  0  0  0  0  0  0  0".format(
                *atom_line)
        )
        idx_table[k] = i + 1
    return lines, idx_table


def bond_block(mol, idx_table):
    stereo_conv = {
        1: {0: 0, 1: 1, 3: 4, 2: 6},
        2: {0: 0, 1: 1, 3: 3, 2: 6}
    }
    lines = []
    for u, v, b in mol.bonds_iter():
        bond_line = []
        bond_line.append("{:>3}".format(idx_table[u]))
        bond_line.append("{:>3}".format(idx_table[v]))
        bond_line.append("{:>3}".format(b.order))
        try:
            bond_line.append("{:>3}".format(stereo_conv[b.order][b.type]))
        except KeyError:
            bond_line.append("{:>3}".format(0))
        lines.append("{}{}{}{}  0  0  0".format(*bond_line))
    return lines


def prop_block(mol, idx_table):
    lines = []
    chg_line = []
    rad_line = []
    iso_line = []
    for k, atom in mol.atoms_iter():
        idx = "{:>3}".format(idx_table[k])
        if atom.charge:
            chg = "{:>3}".format(atom.charge)
            chg_line.append(" {} {}".format(idx, chg))
        if atom.multi != 1:
            rad = "{:>3}".format(atom.multi)
            rad_line.append(" {} {}".format(idx, rad))
        if atom.mass is not None:
            iso = "{:>3}".format(atom.mass)
            iso_line.append(" {} {}".format(idx, iso))
    if len(chg_line):
        chg_count = "{:>3}".format(len(chg_line))
        lines.append("M  CHG{}{}".format(chg_count, "".join(chg_line)))
    if len(rad_line):
        rad_count = "{:>3}".format(len(rad_line))
        lines.append("M  RAD{}{}".format(rad_count, "".join(rad_line)))
    if len(iso_line):
        iso_count = "{:>3}".format(len(iso_line))
        lines.append("M  ISO{}{}".format(iso_count, "".join(iso_line)))
    return lines


def data_block(mol):
    lines = []
    for k, v in mol.data.items():
        lines.append("> <{}>".format(k))
        lines.append(str(v))
        lines.append("")
    lines.append("$$$$")
    return lines


def mol_block(mol, sdfile=True):
    lines = [
        "",
        "Chorus version {}".format(chorus.VERSION),
        ""
    ]
    chiral_flag = 0
    count_line = [
        "{:>3}".format(mol.atom_count()),
        "{:>3}".format(mol.bond_count()),
        "{:>3}".format(chiral_flag)
    ]
    lines.append("{}{}  0  0{}  0  0  0  0  0999 V2000".format(*count_line))
    atoms, idx_table = atom_block(mol)
    lines.extend(atoms)
    bonds = bond_block(mol, idx_table)
    if bonds:
        lines.extend(bonds)
    props = prop_block(mol, idx_table)
    if props:
        lines.extend(props)
    lines.append("M  END")
    if sdfile:
        lines.extend(data_block(mol))
    lines.append("")
    return "\n".join(lines)


def mols_to_text(mols):
    """Converts molecules to the SDFile format text

    Args:
        mols: list of molecule objects

    Returns:
        SDFile text
    """
    return "".join(mol_block(mol) for mol in mols)


def mols_to_file(mols, path):
    """Save molecules to the SDFile format file

    Args:
        mols: list of molecule objects
        path: file path to save
    """
    with open(path, 'w') as f:
        f.write(mols_to_text(mols))
