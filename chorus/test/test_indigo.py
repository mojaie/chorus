#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

try:
    from chorus import indigo
    from indigo import Indigo
    idg = Indigo()
except ImportError:
    raise unittest.SkipTest("Indigo is not available.")

from chorus import v2000reader as reader
from chorus.test.ctabprovider import CTABS
from chorus import molutil


class TestIndigo(unittest.TestCase):
    def test_to_real_mol(self):
        mol = reader.mol_from_text(CTABS["Phe"])
        idmol = indigo.to_real_mol(mol)
        # self.assertEqual(idmol.GetNumBonds(), mol.bond_count())
        # self.assertAlmostEqual(Descriptors.MolWt(rdmol), molutil.mw(mol), 2)

    def test_from_indigo(self):
        idmol = idg.loadMolecule(CTABS["Phe"])
        mol = indigo.from_indigo(idmol)
        # self.assertEqual(rdmol.GetNumBonds(), mol.bond_count())
        # self.assertAlmostEqual(Descriptors.MolWt(rdmol), molutil.mw(mol), 2)

    def test_fingerprint(self):
        mol1 = reader.mol_from_text(CTABS["Goserelin"])
        mol2 = reader.mol_from_text(CTABS["Goserelin"])
        self.assertEqual(indigo.fingerprint_similarity(mol1, mol2), 1)
        mol2 = reader.mol_from_text(CTABS["Formestane"])
        self.assertEqual(indigo.fingerprint_similarity(mol1, mol2), 0.23)
