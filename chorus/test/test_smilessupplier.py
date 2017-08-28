#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest
# import logging
# logging.basicConfig(level=logging.DEBUG)

from chorus.smilessupplier import smiles_to_compound

# logger = logging.getLogger('cheddar.chem.smilessupplier')


class TestSmilesSupplier(unittest.TestCase):
    """
    def debug(self, d=False):
        if d:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.WARNING)

    def setUp(self):
        self.debug()
    """

    def test_null(self):
        mol = smiles_to_compound("")
        self.assertTrue(len(mol.descriptors))  # descriptors assigned
        self.assertEqual(mol.atom_count(), 0)
        self.assertEqual(mol.bond_count(), 0)

    def test_methane(self):
        compound = smiles_to_compound("C")
        self.assertEqual(compound.atom(1).symbol, 'C')
        self.assertEqual(compound.bond_count(), 0)

    def test_ethanol(self):
        compound = smiles_to_compound("CCO")
        self.assertEqual(compound.atom_count(), 3)
        compound = smiles_to_compound("[CH3][CH2][OH]")
        self.assertEqual(compound.atom_count(), 3)

    def test_nitrogen(self):
        compound = smiles_to_compound("N#N")
        self.assertEqual(compound.atom(2).symbol, 'N')
        self.assertEqual(compound.bond(1, 2).order, 3)

    def test_vinylchloride(self):
        compound = smiles_to_compound("C=CCl")
        self.assertEqual(compound.bond(1, 2).order, 2)
        self.assertEqual(compound.atom(3).symbol, 'Cl')

    def test_t_butyl_alchol(self):
        compound = smiles_to_compound("CC(C)(C)O")
        self.assertEqual(compound.atom(5).symbol, 'O')
        self.assertCountEqual(compound.neighbors(2), [1, 3, 4, 5])

    def test_trifluoroacetic_acid(self):
        compound = smiles_to_compound("FC(F)(F)C(=O)O")
        self.assertEqual(compound.bond(5, 6).order, 2)
        self.assertCountEqual(compound.neighbors(2), [1, 3, 4, 5])

    def test_edta(self):
        compound = smiles_to_compound("C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O")
        self.assertCountEqual(compound.neighbors(5), [4, 6, 7])
        self.assertCountEqual(compound.neighbors(12), [1, 13, 17])
        self.assertCountEqual(compound.neighbors(14), [13, 15, 16])

    def test_cyclohexene(self):
        compound = smiles_to_compound("C1C=CCCC1")
        self.assertCountEqual(compound.neighbors(1), [2, 6])
        self.assertEqual(compound.bond(2, 3).order, 2)
        compound = smiles_to_compound("C=1CCCCC1")
        self.assertCountEqual(compound.neighbors(1), [2, 6])
        self.assertEqual(compound.bond(1, 6).order, 2)
        compound = smiles_to_compound("C1CCCCC=1")
        self.assertCountEqual(compound.neighbors(1), [2, 6])
        self.assertEqual(compound.bond(1, 6).order, 2)

    def test_bicyclohexyl(self):
        compound = smiles_to_compound("C1CCCCC1C1CCCCC1")
        self.assertEqual(compound.bond(1, 6).order, 1)
        self.assertEqual(compound.bond(7, 12).order, 1)
        compound = smiles_to_compound("C1CCCCC1C2CCCCC2")
        self.assertEqual(compound.bond(1, 6).order, 1)
        self.assertEqual(compound.bond(7, 12).order, 1)

    def test_naphthalene(self):
        compound = smiles_to_compound("C=1C=CC=C2C1C=CC=C2")
        self.assertEqual(compound.bond(1, 6).order, 2)
        self.assertEqual(compound.bond(5, 10).order, 1)
        compound = smiles_to_compound("C1=CC=C2C(=C1)C=CC=C2")
        self.assertEqual(compound.bond(1, 6).order, 1)
        self.assertEqual(compound.bond(4, 10).order, 1)

    def test_nacl(self):
        compound = smiles_to_compound("[Na+].[Cl-]")
        self.assertEqual(compound.atom(1).symbol, 'Na')
        self.assertEqual(compound.atom(1).charge, 1)
        self.assertEqual(compound.atom(2).symbol, 'Cl')
        self.assertEqual(compound.atom(2).charge, -1)

    def test_iron2_cation(self):
        compound = smiles_to_compound("[Fe++]")
        self.assertEqual(compound.atom(1).symbol, 'Fe')
        self.assertEqual(compound.atom(1).charge, 2)
        compound = smiles_to_compound("[Fe+2]")
        self.assertEqual(compound.atom(1).symbol, 'Fe')
        self.assertEqual(compound.atom(1).charge, 2)

    def test_copper_sulfate(self):
        compound = smiles_to_compound("[Cu+2].[O-]S(=O)(=O)[O-]")
        self.assertEqual(compound.atom(1).charge, 2)
        self.assertEqual(compound.atom(6).charge, -1)
        self.assertCountEqual(compound.neighbors(3), [2, 4, 5, 6])

    def test_alanine(self):
        compound = smiles_to_compound("N[C@H](C)C(=O)O")
        self.assertCountEqual(compound.neighbors(2), [1, 3, 4, 5])
        self.assertEqual(compound.atom(2).stereo, 1)
        compound = smiles_to_compound("N[C@@H](C)C(=O)O")
        self.assertCountEqual(compound.neighbors(2), [1, 3, 4, 5])
        self.assertEqual(compound.atom(2).stereo, -1)

    def test_benzene(self):
        compound = smiles_to_compound("c1ccccc1")
        self.assertEqual(compound.atom(1).pi, 1)
        self.assertEqual(compound.bond(1, 6).order, 1)

    def test_2_pyridone(self):
        compound = smiles_to_compound("O=c1[nH]cccc1")
        self.assertEqual(compound.atom(3).pi, 1)
        # Implicit H without chiral specification should be ignored
        self.assertCountEqual(compound.neighbors(3), [2, 4])

    def test_sodium_phenoxide(self):
        compound = smiles_to_compound("c1cc([O-].[Na+])ccc1")
        self.assertFalse(compound.neighbors(5))
        self.assertCountEqual(compound.neighbors(3), [2, 4, 6])

    def test_thiamin(self):
        compound = smiles_to_compound("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2")
        self.assertEqual(compound.atom(7).charge, 1)
        self.assertEqual(compound.bond(4, 9).order, 1)
        self.assertEqual(compound.bond(11, 17).order, 1)

    def test_unsupported_symbol(self):
        with self.assertRaises(ValueError):
            smiles_to_compound("Hoge")
