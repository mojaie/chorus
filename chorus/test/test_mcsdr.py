#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

import networkx as nx

from chorus.demo import MOL
from chorus import v2000reader as reader
from chorus.smilessupplier import smiles_to_compound
from chorus import mcsdr
from chorus.util import debug


class TestMCS(unittest.TestCase):
    def test_reachables(self):
        g = nx.Graph([(1, 2), (1, 3), (1, 4), (2, 5), (5, 6),
                     (6, 7), (3, 8), (8, 9), (4, 9)])
        self.assertEqual(mcsdr.reachables(g, 1, 100),
                         {1: 0, 2: 1, 3: 1, 4: 1, 5: 2,
                         6: 3, 7: 4, 8: 2, 9: 2})
        self.assertEqual(mcsdr.reachables(g, 1, 2),
                         {1: 0, 2: 1, 3: 1, 4: 1, 5: 2, 8: 2, 9: 2})

    def test_mcsdr1(self):
        # TODO: pi mismatch is not acceptable
        mol1 = reader.mol_from_text(MOL["Phe"])
        mol2 = reader.mol_from_text(MOL["Arg"])
        arr1 = mcsdr.DescriptorArray(mol1)
        arr2 = mcsdr.DescriptorArray(mol2)
        self.assertEqual(mcsdr.from_array(arr1, arr2).edge_count(), 5)
        # Delta-y exchange will not occur due to distance descriptor
        mol1 = smiles_to_compound("C1OC1CCC(=O)O")
        mol2 = smiles_to_compound("CC(O)CCC(=O)O")
        arr1 = mcsdr.DescriptorArray(mol1)
        arr2 = mcsdr.DescriptorArray(mol2)
        self.assertEqual(mcsdr.from_array(arr1, arr2).edge_count(), 7)

    def test_mcsdr2(self):
        # Disconnected
        mol1 = smiles_to_compound("C1CCCC1CCCC(=O)O")
        mol2 = reader.mol_from_text(MOL["CaAcO2"])
        arr1 = mcsdr.DescriptorArray(mol1)
        arr2 = mcsdr.DescriptorArray(mol2)
        self.assertEqual(mcsdr.from_array(arr1, arr2).edge_count(), 3)
        # No line graph
        mol1 = smiles_to_compound("CO")
        mol2 = smiles_to_compound("CC")
        arr1 = mcsdr.DescriptorArray(mol1)
        arr2 = mcsdr.DescriptorArray(mol2)
        self.assertEqual(mcsdr.from_array(arr1, arr2).edge_count(), 0)
        # TODO: minimum MCS edge size is 2
        mol1 = smiles_to_compound("CCO")
        mol2 = smiles_to_compound("CCC")
        arr1 = mcsdr.DescriptorArray(mol1)
        arr2 = mcsdr.DescriptorArray(mol2)
        self.assertEqual(mcsdr.from_array(arr1, arr2).edge_count(), 0)
        # This works well
        mol1 = smiles_to_compound("CCCO")
        mol2 = smiles_to_compound("CCCC")
        arr1 = mcsdr.DescriptorArray(mol1)
        arr2 = mcsdr.DescriptorArray(mol2)
        self.assertEqual(mcsdr.from_array(arr1, arr2).edge_count(), 2)

    def test_timeout(self):
        mol = reader.mol_from_text(MOL["Buckminsterfullerene"])
        tout = mcsdr.DescriptorArray(mol, timeout=0.1)
        self.assertFalse(tout.valid)
        arr = mcsdr.DescriptorArray(mol, diameter=5, timeout=0.1)
        sim = mcsdr.from_array(arr, arr, timeout=0.1)
        self.assertGreater(sim.local_sim(), 0)
