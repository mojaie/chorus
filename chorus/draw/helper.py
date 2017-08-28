#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import itertools
import math
import statistics

from chorus.util import geometry
from chorus.util import iterator
from chorus import molutil


def display_terminal_carbon(mol):
    """Set visible=True to the terminal carbon atoms.
    """
    for i, a in mol.atoms_iter():
        if mol.neighbor_count(i) == 1:
            a.visible = True


def equalize_terminal_double_bond(mol):
    """Show equalized double bond if it is connected to terminal atom.
    """
    for i, a in mol.atoms_iter():
        if mol.neighbor_count(i) == 1:
            nb = list(mol.neighbors(i).values())[0]
            if nb.order == 2:
                nb.type = 2


def spine_to_terminal_wedge(mol):
    """Arrange stereo wedge direction from spine to terminal atom
    """
    for i, a in mol.atoms_iter():
        if mol.neighbor_count(i) == 1:
            ni, nb = list(mol.neighbors(i).items())[0]
            if nb.order == 1 and nb.type in (1, 2) \
                    and ni > i != nb.is_lower_first:
                nb.is_lower_first = not nb.is_lower_first
                nb.type = {1: 2, 2: 1}[nb.type]


def format_ring_double_bond(mol):
    """Set double bonds around the ring.
    """
    mol.require("Topology")
    mol.require("ScaleAndCenter")
    for r in sorted(mol.rings, key=len, reverse=True):
        vertices = [mol.atom(n).coords for n in r]
        try:
            if geometry.is_clockwise(vertices):
                cpath = iterator.consecutive(itertools.cycle(r), 2)
            else:
                cpath = iterator.consecutive(itertools.cycle(reversed(r)), 2)
        except ValueError:
            continue
        for _ in r:
            u, v = next(cpath)
            b = mol.bond(u, v)
            if b.order == 2:
                b.type = int((u > v) == b.is_lower_first)


def scale_and_center(mol):
    """Center and Scale molecule 2D coordinates.
    This method changes mol coordinates directly to center but not scale.
    This method returns width, height and MLB(median length of bond)
    and scaling will be done by drawer method with these values.

    Returns:
        width: float
        height: float
        mlb: median length of bond
    """
    cnt = mol.atom_count()
    if cnt < 2:
        mol.size2d = (0, 0, 1)
        mol.descriptors.add("ScaleAndCenter")
        return
    xs = []
    ys = []
    for _, atom in mol.atoms_iter():
        xs.append(atom.coords[0])
        ys.append(atom.coords[1])
    xmin, xmax = (min(xs), max(xs))
    ymin, ymax = (min(ys), max(ys))
    width = xmax - xmin
    height = ymax - ymin
    x_offset = width / 2 + xmin
    y_offset = height / 2 + ymin
    dists = []
    for u, v, _ in mol.bonds_iter():
        dists.append(geometry.distance(mol.atom(u).coords, mol.atom(v).coords))
    if not len(dists):  # No connection
        mlb = math.sqrt(max([width, height]) / cnt)  # empirical
    elif not max(dists):  # All connected atoms are overlapped
        mol.size2d = (0, 0, 1)
        mol.descriptors.add("ScaleAndCenter")
        return
    else:
        mlb = statistics.median(dists)
    # Centering
    for _, atom in mol.atoms_iter():
        atom.coords = (atom.coords[0] - x_offset, atom.coords[1] - y_offset)
    mol.size2d = (width, height, mlb)
    mol.descriptors.add("ScaleAndCenter")


def ready_to_draw(mol):
    """Shortcut function to prepare molecule to draw.
    Overwrite this function for customized appearance.
    It is recommended to clone the molecule before draw
    because all the methods above are destructive.
    """
    copied = molutil.clone(mol)
    # display_terminal_carbon(mol)
    equalize_terminal_double_bond(copied)
    # spine_to_terminal_wedge(copied)
    scale_and_center(copied)
    format_ring_double_bond(copied)
    return copied
