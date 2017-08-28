#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import json
import unittest

from chorus.smilessupplier import smiles_to_compound
from chorus.substructure import equal, substructure, filter_


class TestSubstructure(unittest.TestCase):
    def test_filter(self):
        query = smiles_to_compound("c1ccccc1C2CCC2C")
        mol = smiles_to_compound("c1ccccc1C2CCC2")
        self.assertFalse(filter_(mol, query))
        mol = smiles_to_compound("c1ccccc1C2CC=C2C")
        self.assertFalse(filter_(mol, query))
        mol = smiles_to_compound("c1ccccc1C2CC2CC")
        self.assertFalse(filter_(mol, query))
        mol = smiles_to_compound("c1cccc2c1CC2CCC")
        self.assertFalse(filter_(mol, query))
        mol = smiles_to_compound("c1ccccc1C2CC(C)C2")
        self.assertTrue(filter_(mol, query))

    def test_equal(self):
        mol = smiles_to_compound("c1ccccc1CCC")
        mj = json.dumps(mol.jsonized())[:1000]
        query = smiles_to_compound("c1ccccc1CCC")
        self.assertTrue(equal(mol, query))
        query = smiles_to_compound("c1ccccc1CCC.O")
        qj = json.dumps(query.jsonized())[:1000]
        self.assertTrue(equal(mol, query))
        self.assertFalse(equal(mol, query, largest_only=False))
        # Non destructive
        self.assertEqual(mj, json.dumps(mol.jsonized())[:1000])
        self.assertEqual(qj, json.dumps(query.jsonized())[:1000])

    def test_substructure(self):
        mol = smiles_to_compound("c1ccccc1")
        query = smiles_to_compound("c1ccccc1CCC")
        qj = json.dumps(query.jsonized())[:1000]
        self.assertTrue(substructure(mol, query))
        mol = smiles_to_compound("c1ccccc1.CCO.O.O")
        mj = json.dumps(mol.jsonized())[:1000]
        self.assertTrue(substructure(mol, query))
        self.assertFalse(substructure(mol, query, largest_only=False))
        # Non destructive
        self.assertEqual(mj, json.dumps(mol.jsonized())[:1000])
        self.assertEqual(qj, json.dumps(query.jsonized())[:1000])
