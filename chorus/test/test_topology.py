#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest
from collections import deque, Counter

from chorus.util import debug
from chorus.smilessupplier import smiles_to_compound
from chorus import v2000reader as reader
from chorus.test.ctabprovider import CTABS
from chorus import topology


class TestTopology(unittest.TestCase):
    def test_resolve_inclusion(self):
        a = [1, 2, 3, 4, 5]
        b = [5, 6, 7, 8, 9, 1, 2, 3, 4]
        self.assertEqual(topology.resolve_inclusion(a, b)[1],
                         [5, 6, 7, 8, 9, 1])
        b = [1, 6, 7, 8, 9]
        self.assertFalse(topology.resolve_inclusion(a, b))
        b = [3, 4, 5, 6, 7, 8, 9]
        self.assertFalse(topology.resolve_inclusion(a, b))
        b = [2, 3, 4, 5, 6, 7, 8]
        self.assertEqual(topology.resolve_inclusion(a, b)[1],
                         [5, 6, 7, 8, 2, 1])
        self.assertEqual(topology.resolve_inclusion(b, a)[1],
                         [1, 2, 3, 4, 5])
        # reversed
        a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        b = [9, 8, 12, 11, 10]
        self.assertEqual(topology.resolve_inclusion(a, b)[0],
                         [11, 1, 2, 3, 4, 5, 6, 7, 8, 12])

    def equivalent_ring(self, ring, ref):
        dqr = deque(ring)
        rev = deque(reversed(ring))
        for i in range(len(ring)):
            dqr.rotate(1)
            rev.rotate(1)
            if list(dqr) == list(ref):
                return True
            if list(rev) == list(ref):
                return True
        return False

    def test_recognize(self):
        # Phenylalanin
        mol = reader.mol_from_text(CTABS["Phe"])
        self.assertTrue(
            self.equivalent_ring(mol.rings[0], [8, 11, 12, 10, 7, 6]))
        self.assertEqual(mol.scaffolds, [[0]])
        self.assertEqual(mol.isolated, [])
        # Premarin
        mol = reader.mol_from_text(CTABS["Premarin"])
        for a, b in zip(sorted(mol.rings), sorted([
            [10, 9, 8, 7, 4, 5],
            [10, 14, 13, 12, 11, 9],
            [24, 25, 26, 12, 11],
            [2, 3, 4, 5, 6, 1]
        ])):
            self.assertTrue(self.equivalent_ring(a, b))
        self.assertEqual(mol.scaffolds, [[0, 1, 2, 3]])
        self.assertEqual(mol.isolated, [[28]])
        # Pyrene
        mol = smiles_to_compound("C12=CC=C3C=CC=C4C=CC(C2=C34)=CC=C1")
        self.assertEqual([len(r) for r in mol.rings], [6, 6, 6, 6])
        self.assertEqual(mol.scaffolds, [[0, 1, 2, 3]])
        self.assertEqual(mol.isolated, [])
        # KCl
        mol = reader.mol_from_text(CTABS["KCl"])
        self.assertEqual(mol.rings, [])
        self.assertEqual(mol.scaffolds, [])
        self.assertEqual(mol.isolated, [[2]])
        # Goserelin
        mol = reader.mol_from_text(CTABS["Goserelin"])
        self.assertEqual(sorted([len(r) for r in mol.rings]),
                         [5, 5, 5, 5, 6, 6])
        self.assertEqual(
            Counter([len(s) for s in mol.scaffolds]), {1: 4, 2: 1})
        # Tetrahedrane (K4 graph)
        mol = smiles_to_compound("C12C3C1C23")
        self.assertEqual([len(r) for r in mol.rings], [3, 3, 3])

    # @debug.profile
    def test_minify_ring(self):
        # Cyanocobalamin
        mol = reader.mol_from_text(CTABS["Cyanocobalamin"])
        self.assertEqual(
            Counter([len(r) for r in mol.rings]), {5: 7, 6: 4, 19: 1})
        # Rifabutin
        # TODO: this can be [5, 5, 6, 6, 6, 25] or [5, 5, 6, 6, 6, 24]
        # mol = reader.mol_from_text(CTABS["Rifabutin"])
        # self.assertEqual(sorted([len(r) for r in mol.rings]),
        #                 [5, 5, 6, 6, 6, 24])
        # Cubane
        mol = smiles_to_compound("C12C3C4C1C5C4C3C25")
        self.assertEqual(
            Counter([len(r) for r in mol.rings]), {4: 5})
        # Pinene
        mol = smiles_to_compound("CC1(C2CCC(=C)C1C2)C")
        self.assertEqual(
            Counter([len(r) for r in mol.rings]), {4: 1, 6: 1})
        # Fullerene
        # TODO: It does not work with such a comprecated case.
        """
        mol = reader.mol_from_text(CTABS["Fullerene"])
        topology.recognize(mol)
        topology.minify_ring(mol)
        c = Counter([len(r) for r in mol.rings])
        self.assertEqual(dict(c), {5: 12, 6: 19})
        """
