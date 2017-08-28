#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from chorus import rdkit
except ImportError:
    raise unittest.SkipTest("RDKit is not available.")

from chorus import v2000reader as reader
from chorus.demo import MOL
from chorus import molutil


class TestRDKit(unittest.TestCase):
    def test_to_rdmol(self):
        mol = reader.mol_from_text(MOL["Phe"])
        rdmol = rdkit.to_rdmol(mol)
        self.assertEqual(rdmol.GetNumBonds(), mol.bond_count())
        self.assertAlmostEqual(Descriptors.MolWt(rdmol), molutil.mw(mol), 2)

    def test_from_rdmol(self):
        rdmol = Chem.MolFromMolBlock(MOL["Phe"])
        mol = rdkit.from_rdmol(rdmol)
        self.assertEqual(rdmol.GetNumBonds(), mol.bond_count())
        self.assertAlmostEqual(Descriptors.MolWt(rdmol), molutil.mw(mol), 2)

    def test_morgan(self):
        mol1 = reader.mol_from_text(MOL["Goserelin"])
        mol2 = reader.mol_from_text(MOL["Formestane"])
        self.assertEqual(rdkit.morgan_sim(mol1, mol2, 2), 0.107)

    def test_fmcs(self):
        mol1 = reader.mol_from_text(MOL["Phe"])
        mol2 = reader.mol_from_text(MOL["Arg"])
        self.assertEqual(rdkit.fmcs(mol1, mol2)["mcs_edges"], 7)
        mol1 = reader.mol_from_text(MOL["Lys"])
        mol2 = reader.mol_from_text(MOL["Formestane"])
        self.assertEqual(rdkit.fmcs(mol1, mol2)["similarity"], 0.194)
        # null molecule
        mol1 = reader.mol_from_text(MOL["null"])
        mol2 = reader.mol_from_text(MOL["Phe"])
        mcs = rdkit.fmcs(mol1, mol2)
        self.assertEqual(mcs["mcs_edges"], 0)
        self.assertEqual(mcs["similarity"], 0)

    @unittest.skip('This takes 1 sec.')
    def test_fmcs_timeout(self):
        mol1 = reader.mol_from_text(MOL["Docetaxel"])
        mol2 = reader.mol_from_text(MOL["Paclitaxel"])
        res = rdkit.fmcs(mol1, mol2, timeout=1)
        self.assertTrue(res["canceled"])
