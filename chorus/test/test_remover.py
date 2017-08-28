#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

from chorus import smilessupplier as ss
from chorus import remover
from chorus import v2000reader as reader
from chorus.test.ctabprovider import CTABS


class TestRemover(unittest.TestCase):
    def test_remove_water(self):
        mol = ss.smiles_to_compound("CCO.O.O")
        remover.remove_water(mol)
        self.assertEqual(len(mol), 3)

    def test_remove_salt(self):
        mol = ss.smiles_to_compound("CC[O-].[Na+]")
        remover.remove_salt(mol)
        self.assertEqual(len(mol), 3)

    def test_remove_coordinated_metal(self):
        mol = reader.mol_from_text(CTABS["Cyanocobalamin"])
        self.assertEqual(len(mol), 95)
        remover.remove_coordinated_metal(mol)
        self.assertEqual(len(mol), 94)
