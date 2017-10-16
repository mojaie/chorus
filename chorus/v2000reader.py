#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import traceback
import re

from chorus.model.atom import Atom
from chorus.model.bond import Bond
from chorus.model.graphmol import Compound
from chorus import molutil
import chorus.util.text as tx


def inspect(lines):
    """Inspect SDFile list of string

    Returns:
        tuple: (data label list, number of records)
    """
    labels = set()
    count = 0
    exp = re.compile(r">.*?<([\w ]+)>")  # Space should be accepted
    valid = False
    for line in lines:
        if line.startswith("M  END\n"):
            valid = True
        elif line.startswith("$$$$"):
            count += 1
            valid = False
        else:
            result = exp.match(line)
            if result:
                labels.add(result.group(1))
    if valid:
        count += 1
    return list(labels), count


def inspect_text(text):
    # Lazy line splitter. More efficient memory usage than str.split.
    exp = re.compile(r"[^\n]*\n|.")
    lines = (x.group(0) for x in re.finditer(exp, text))
    return inspect(lines)


def inspect_file(path):
    """Inspect SDFile structure

    Returns:
        tuple: (data label list, number of records)
    """
    with open(path, 'rb') as f:
        labels, count = inspect(tx.decode(line) for line in f)
    return labels, count


def optional_data(lines):
    """Parse SDFile data part into dict"""
    data = {}
    exp = re.compile(r">.*?<([\w ]+)>")  # Space should be accepted
    for i, line in enumerate(lines):
        result = exp.match(line)
        if result:
            data[result.group(1)] = lines[i + 1]
    return data


def atoms(lines):
    """Parse atom block into atom objects

    Returns:
        dict: networkx nodes
    """
    # Convert sdf style charge to actual charge
    conv_charge_table = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}
    results = {}
    for i, line in enumerate(lines):
        symbol = line[31:34].rstrip()
        try:
            atom = Atom(symbol)
        except KeyError:
            raise ValueError(symbol)
        xpos = float(line[0:10])
        ypos = float(line[10:20])
        zpos = float(line[20:30])
        atom.coords = (xpos, ypos, zpos)
        atom.mass_diff = int(line[34:37])
        old_sdf_charge = int(line[37:40])
        atom.charge = conv_charge_table[old_sdf_charge]
        if old_sdf_charge == 4:
            atom.radical = 1
        # atom.stereo_flag = int(line[40:43])  # Not used
        # valence = int(line[46:49])
        # if valence:
        #     atom.valence = valence
        results[i + 1] = {"atom": atom}
    return results


def bonds(lines, atoms):
    """Parse bond block into bond objects

    Returns:
        dict: networkx adjacency dict
    """
    # Convert sdf style stereobond (see chem.model.bond.Bond)
    conv_stereo_table = {0: 0, 1: 1, 3: 3, 4: 3, 6: 2}
    results = {a: {} for a in atoms}
    for line in lines:
        bond = Bond()
        first = int(line[0:3])
        second = int(line[3:6])
        if first > second:
            bond.is_lower_first = 0
        order = int(line[6:9])
        if order < 4:
            bond.order = order
        bond.type = conv_stereo_table[int(line[9:12])]
        results[first][second] = {"bond": bond}
        results[second][first] = {"bond": bond}
    return results


def properties(lines):
    """Parse properties block

    Returns:
        dict: {property_type: (atom_index, value)}
    """
    results = {}
    for i, line in enumerate(lines):
        type_ = line[3:6]
        if type_ not in ["CHG", "RAD", "ISO"]:
            continue  # Other properties are not supported yet
        count = int(line[6:9])
        results[type_] = []
        for j in range(count):
            idx = int(line[10 + j * 8: 13 + j * 8])
            val = int(line[14 + j * 8: 17 + j * 8])
            results[type_].append((idx, val))
    return results


def add_properties(props, mol):
    """apply properties to the molecule object

    Returns:
        None (alter molecule object directly)
    """
    if not props:
        return
    # The properties supersedes all charge and radical values in the atom block
    for _, atom in mol.atoms_iter():
        atom.charge = 0
        atom.multi = 1
        atom.mass = None
    for prop in props.get("CHG", []):
        mol.atom(prop[0]).charge = prop[1]
    for prop in props.get("RAD", []):
        mol.atom(prop[0]).multi = prop[1]
    for prop in props.get("ISO", []):
        mol.atom(prop[0]).mass = prop[1]


def molecule(lines):
    """Parse molfile part into molecule object

    Args:
        lines (list): lines of molfile part

    Raises:
        ValueError: Symbol not defined in periodictable.yaml
                    (Polymer expression not supported yet)
    """
    count_line = lines[3]
    num_atoms = int(count_line[0:3])
    num_bonds = int(count_line[3:6])
    # chiral_flag = int(count_line[12:15])  # Not used
    # num_prop = int(count_line[30:33])  # "No longer supported"
    compound = Compound()
    compound.graph.node = atoms(lines[4: num_atoms+4])
    compound.graph.adj = bonds(lines[num_atoms+4: num_atoms+num_bonds+4],
                               compound.graph.node.keys())
    compound.graph.edge = compound.graph.adj
    props = properties(lines[num_atoms+num_bonds+4:])
    add_properties(props, compound)
    return compound


def mol_supplier(lines, no_halt, assign_descriptors):
    """Yields molecules generated from CTAB text

    Args:
        lines (iterable): CTAB text lines
        no_halt (boolean):
            True: shows warning messages for invalid format and go on.
            False: throws an exception for it and stop parsing.
        assign_descriptors (boolean):
            if True, default descriptors are automatically assigned.
    """
    def sdf_block(lns):
        mol = []
        opt = []
        is_mol = True
        for line in lns:
            if line.startswith("$$$$"):
                yield mol[:], opt[:]
                is_mol = True
                mol.clear()
                opt.clear()
            elif line.startswith("M  END"):
                is_mol = False
            elif is_mol:
                mol.append(line.rstrip())
            else:
                opt.append(line.rstrip())
        if mol:
            yield mol, opt

    for i, (mol, opt) in enumerate(sdf_block(lines)):
        try:
            c = molecule(mol)
            if assign_descriptors:
                molutil.assign_descriptors(c)
        except ValueError as err:
            if no_halt:
                print("Unsupported symbol: {} (#{} in v2000reader)".format(
                      err, i + 1))
                c = molutil.null_molecule(assign_descriptors)
            else:
                raise ValueError("Unsupported symbol: {}".format(err))
        except:
            if no_halt:
                print("Unexpected error (#{} in v2000reader)".format(i + 1))
                c = molutil.null_molecule(assign_descriptors)
                c.data = optional_data(opt)
                yield c
                continue
            else:
                print(traceback.format_exc())
        c.data = optional_data(opt)
        yield c


def mols_from_text(text, no_halt=True, assign_descriptors=True):
    """Returns molecules generated from sdfile text

    Throws:
        StopIteration: if the text does not have molecule
        ValueError: if Unsupported symbol is found
    """
    if isinstance(text, bytes):
        t = tx.decode(text)
    else:
        t = text
    # Lazy line splitter. More efficient memory usage than str.split.
    exp = re.compile(r"[^\n]*\n|.")
    sp = (x.group(0) for x in re.finditer(exp, t))
    for c in mol_supplier(sp, no_halt, assign_descriptors):
        yield c


def mol_from_text(text, assign_descriptors=True):
    """Parse CTAB text and return first one as a Compound object.

    Throws:
        StopIteration: if the text does not have molecule
        ValueError: if Unsupported symbol is found
    """
    cg = mols_from_text(text, False, assign_descriptors)
    return next(cg)


def mols_from_file(path, no_halt=True, assign_descriptors=True):
    """Compound supplier from CTAB text file (.mol, .sdf)"""
    with open(path, 'rb') as f:
        fd = (tx.decode(line) for line in f)
        for c in mol_supplier(fd, no_halt, assign_descriptors):
            yield c


def mol_from_file(path, assign_descriptors=True):
    """Parse CTAB file and return first one as a Compound object."""
    cs = mols_from_file(path, False, assign_descriptors)
    return next(cs)
