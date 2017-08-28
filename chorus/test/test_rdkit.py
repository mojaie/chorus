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
from chorus.test.ctabprovider import CTABS
from chorus import molutil


class TestRDKit(unittest.TestCase):
    def test_to_rdmol(self):
        mol = reader.mol_from_text(CTABS["Phe"])
        rdmol = rdkit.to_rdmol(mol)
        self.assertEqual(rdmol.GetNumBonds(), mol.bond_count())
        self.assertAlmostEqual(Descriptors.MolWt(rdmol), molutil.mw(mol), 2)

    def test_from_rdmol(self):
        rdmol = Chem.MolFromMolBlock(CTABS["Phe"])
        mol = rdkit.from_rdmol(rdmol)
        self.assertEqual(rdmol.GetNumBonds(), mol.bond_count())
        self.assertAlmostEqual(Descriptors.MolWt(rdmol), molutil.mw(mol), 2)

    def test_morgan(self):
        mol1 = reader.mol_from_text(CTABS["Goserelin"])
        mol2 = reader.mol_from_text(CTABS["Formestane"])
        self.assertEqual(rdkit.morgan_sim(mol1, mol2, 2), 0.107)

    def test_fmcs(self):
        mol1 = reader.mol_from_text(CTABS["Phe"])
        mol2 = reader.mol_from_text(CTABS["Arg"])
        self.assertEqual(rdkit.fmcs(mol1, mol2)["mcs_edges"], 7)
        mol1 = reader.mol_from_text(CTABS["Lys"])
        mol2 = reader.mol_from_text(CTABS["Formestane"])
        self.assertEqual(rdkit.fmcs(mol1, mol2)["similarity"], 0.194)
        # null molecule
        mol1 = reader.mol_from_text(CTABS["null"])
        mol2 = reader.mol_from_text(CTABS["Phe"])
        mcs = rdkit.fmcs(mol1, mol2)
        self.assertEqual(mcs["mcs_edges"], 0)
        self.assertEqual(mcs["similarity"], 0)

    @unittest.skip('This takes 1 sec.')
    def test_fmcs_timeout(self):
        mol1 = reader.mol_from_text(CTABS["Docetaxel"])
        mol2 = reader.mol_from_text(CTABS["Paclitaxel"])
        res = rdkit.fmcs(mol1, mol2, timeout=1)
        self.assertTrue(res["canceled"])
