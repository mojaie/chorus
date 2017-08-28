#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

""" SMILES parser module (smilessupplier.py)"""

# TODO: make it simple and eager function

import copy
import re
import logging

from chorus.model.atom import Atom
from chorus.model.bond import Bond
from chorus.model.graphmol import Compound
from chorus import molutil


logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


def element():
    atom = []
    bond = Bond()

    def flush():
        """ Flash a connection if exist
        Returns:
          (Atom(), Bond())
          if no atom, returns (None, Bond())
            (ring flag with multiple bond like =1)
          if neither atom and bond, returns None
        Raises:
          ValueError: if the atom symbol is wrong or not supported
        """
        nonlocal bond
        a = None
        logtext = []
        if atom:
            symbol = ''.join(atom)
            a = Atom(symbol.capitalize())
            if symbol in ('c', 'n', 'o', 'p', 's'):
                a.pi = 1
                a.aromatic = True
            logtext.append('Symbol:{}'.format(a.symbol))
        else:
            logtext.append('No element')
        if bond:
            logtext.append('Bond:{}'.format(bond.order))
        else:
            logtext.append('No bond')
        logger.debug(', '.join(logtext))
        result = (a, bond)  # copy
        atom.clear()  # init
        bond = Bond()  # init
        return result

    def exec_(token):
        """ Process a token and return an atom and a bond if available.
        Args:
          token: token
        Returns:
          flush()
        """
        nonlocal bond
        result = (None, None)
        if not token:
            return flush()
        if atom and token.upper() in (
                'B', 'C', 'N', 'O', 'P', 'I', 'S', 'F',
                '=', '#', '/', '\\', '.'):
            result = flush()
        if token == '.':
            bond = None
        elif token == '=':
            bond.order = 2
        elif token == '#':
            bond.order = 3
        elif token == '/':
            bond.smiles_cis_trans = 1
        elif token == '\\':
            bond.smiles_cis_trans = -1
        else:
            atom.append(token)
        logger.debug('{} is queued'.format(token))
        return result

    return exec_


def group():
    container = []
    exp = re.compile("([0-9]*)([a-zA-Z][a-z]?)H?[0-9]*([\+\-]*[0-9]*)")

    def charge_sign(text):
        """ Process charge signs.
        + -> 1, - -> -1, +3 -> 3, +++ -> 3
        Args:
          text: charge signs (+, - or number)
        Returns:
          charge count (neg -2, -1, 0, 1, 2 pos)
        """
        result = 0
        for c in text:
            cases = {'-': -1, '+': 1}
            if c in cases:
                result += cases[c]
            elif c.isdigit():
                result *= int(c)
        return result

    def flush():
        """ Flash an atom group
        explicit hydrogen on the chiral center is added if exist.
        Returns:
          (Atom, hydrogen(as Compound branch))
        Raises:
          ValueError: if the atom symbol is wrong or not supported
        """
        contents = ''.join(container)
        sym = ''
        stereo = 0
        atom = None
        hydro = None
        if '@@' in contents:
            sym = '@@'
            stereo = -1
        elif '@' in contents:
            sym = '@'
            stereo = 1
        if sym:
            # chiral
            c, h = contents.split(sym)
            atom = Atom(c.capitalize())
            if c in ('c', 'n', 'o', 'p', 's'):
                atom.pi = 1
            atom.stereo = stereo
            if h:
                mol = molecule()
                mol('H')
                hydro = mol(None)
        else:
            # charged, isotope or inorganic atoms
            m = exp.match(contents)
            if m:
                logger.debug('Group parsed: {}'.format(m.groups()))
                isotope = m.group(1)
                atom = Atom(m.group(2).capitalize())
                if m.group(2) in ('c', 'n', 'o', 'p', 's'):
                    atom.pi = 1
                atom.charge = charge_sign(m.group(3))
        container.clear()  # init
        return atom, hydro

    def exec_(token):
        """ Process a token and return a group if available.
        Args:
          token: token
        Returns:
          flush()
        """
        if not token:
            return flush()
        container.append(token)

    return exec_


def molecule():
    idx = 0
    depth = 0
    is_group = False
    cp = Compound()
    cp.add_atom(idx, Atom('Dummy'))
    elem = element()
    grp = group()
    branches = []
    ring_labels = {}

    def connect(e):
        """ connect a bond and branches to the molecule
        Args:
          e: (Atom(), Bond())
        """
        nonlocal idx
        for b, r in branches:
            ai = cp.available_idx()
            bond = b.bond(0, 1)
            b.graph.remove_node(0)
            cp.add_molecule(b, bond, idx, 1)
            for k, v in r.items():
                i, bond = v
                ring(k, bond, i + ai - 1)
        logger.debug('{} branches are merged'.format(len(branches)))
        branches.clear()
        if e[0]:
            ai = cp.available_idx()
            cp.add_atom(ai, e[0])
            if e[1]:
                cp.add_bond(idx, ai, e[1])
                logger.debug('Connected:{} {}'.format(e[0].symbol, e[1].order))
            else:
                logger.debug('Added:{}'.format(e[0].symbol))
            idx = ai
        logger.debug(cp)

    def ring(label, bond=None, i=None):
        """ add ring label to the ring_labels dict.
        if the label matched to existing label, flush and close ring
        Args:
          label: number (string type)
          bond: Bond() to close the ring
          i: index (for branch merging)
        """
        if not i:
            i = idx
        if label in ring_labels:
            u = ring_labels[label][0]
            if not bond:
                bond = ring_labels[label][1]
                if not bond:
                    bond = Bond()
            cp.add_bond(u, i, bond)
            del ring_labels[label]
            logger.debug('ring closed')
            logger.debug(ring_labels)
        else:
            ring_labels[label] = (i, bond)
            logger.debug('ring opened')
            logger.debug(ring_labels)

    def flush():
        """ Flash a molecule
        Returns:
          (Compound, ring flag dict)
        """
        nonlocal idx
        logger.debug('Flushed ({} atoms)'.format(cp.atom_count() - 1))
        logger.debug('Ring labels: {}'.format(ring_labels))
        logger.debug(cp)
        result = (copy.deepcopy(cp), copy.deepcopy(ring_labels))  # copy
        idx = 0  # init
        cp.clear()  # init
        cp.add_atom(idx, Atom('Dummy'))  # init
        ring_labels.clear()  # init
        return result

    def exec_(token):
        """ Process a token and return a compound if available.
        Args:
          token: token
        Returns:
          flush()
        """
        nonlocal depth
        nonlocal is_group
        if not token:
            connect(elem(None))
            return flush()
        logger.debug('{} is supplied'.format(token))
        if token == ':':
            return
        elif token == ')':
            depth -= 1
            if depth > 0:
                branches[-1](token)
                logger.debug('Depth {}: {} is sent to branch'.format(depth, token))
            elif not depth:
                branches[-1] = branches[-1](None)
                logger.debug('branch closed'.format(depth))
            else:
                raise ValueError('Syntax Error: unexpected symbol \")\"')
        elif token == '(':
            if depth:
                branches[-1](token)
                logger.debug('Depth {}: {} is sent to branch'.format(depth, token))
            else:
                connect(elem(None))
                branches.append(molecule())
                logger.debug('Depth {}: branch started'.format(depth))
            depth += 1
        else:
            if depth:
                branches[-1](token)
                logger.debug('Depth {}: {} is sent to branch'.format(depth, token))
            elif token == ']':
                g, h = grp(None)
                e = elem(None)
                if e[0]:
                    connect(e)
                    e = (None, Bond())
                connect((g, e[1]))
                if h:
                    branches.append(h)
                is_group = False
            elif is_group:
                    grp(token)
            elif token == '[':
                    is_group = True
            elif token.isdigit():
                e = elem(None)
                if e[0] is None:
                    ring(token, e[1])
                else:
                    connect(e)
                    ring(token)
            else:
                connect(elem(token))

    return exec_


def smiles_to_compound(smiles, assign_descriptors=True):
    """Convert SMILES text to compound object

    Raises:
        ValueError: SMILES with unsupported format
    """
    it = iter(smiles)
    mol = molecule()
    try:
        for token in it:
                mol(token)
        result, _ = mol(None)
    except KeyError as err:
        raise ValueError("Unsupported Symbol: {}".format(err))
    result.graph.remove_node(0)
    logger.debug(result)
    if assign_descriptors:
        molutil.assign_descriptors(result)
    return result
