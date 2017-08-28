#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

from chorus.test.ctabprovider import CTABS
from chorus import v2000reader as reader
from chorus.smilessupplier import smiles_to_compound
from chorus import wclogp


class TestWClogP(unittest.TestCase):
    def test_wcpaper(self):
        """ Wildman and Crippen 1999, Table.2 """
        mol = smiles_to_compound("C=1C=CC=C(OC)C=1O")
        wclogp.assign_wctype(mol)
        self.assertEqual([a.wctype for _, a in mol.atoms_iter()],
                         ["C18", "C18", "C18", "C18", "C23",
                          "O4", "C3", "C23", "O2"])
        self.assertEqual(wclogp.wclogp(mol), 1.40)
        self.assertEqual(wclogp.wcmr(mol), 34.66)
        mol = smiles_to_compound("C1=CC=CC=C1C2=CC=CC=N2")
        wclogp.assign_wctype(mol)
        self.assertEqual([a.wctype for _, a in mol.atoms_iter()],
                         ["C18", "C18", "C18", "C18", "C18", "C20",
                          "C20", "C18", "C18", "C18", "C18", "N11"])
        self.assertEqual(wclogp.wclogp(mol), 2.75)
        self.assertEqual(wclogp.wcmr(mol), 49.67)
        # maybe typo in the original paper: calcd = 49.67, expt = 50.39

    def test_wctype(self):
        mol = reader.mol_from_text(CTABS["wctype_C_alip"])
        wclogp.assign_wctype(mol)
        self.assertEqual(
            [a.wctype for _, a in mol.atoms_iter()],
            ["C2", "C1", "C7", "C4", "N2", "C3", "N1", "C7", "C5", "O9",
             "C27", "Me1", "C27", "C27", "C27", "H1", "C6", "C26", "C21", "C21",
             "N12", "C18", "C21", "C21", "C11", "O3", "C12", "C1", "C1", "C1",
             "C9", "C10", "C1"])
        mol = reader.mol_from_text(CTABS["wctype_C_arom"])
        wclogp.assign_wctype(mol)
        self.assertEqual(
            [a.wctype for _, a in mol.atoms_iter()],
            ["C21", "C23", "C22", "C24", "C19", "C19", "C14", "C20", "C25", "O1",
             "O8", "O2", "N3", "S1", "C20", "C15", "C16", "C17", "C13", "C18",
             "Cl", "Br", "I", "C8", "P", "F"])
        mol = reader.mol_from_text(CTABS["wctype_N"])
        wclogp.assign_wctype(mol)
        self.assertEqual(
            [a.wctype for _, a in mol.atoms_iter()],
            ["N11", "C21", "C22", "C18", "C22", "C22", "N4", "C4", "N8", "C4",
             "N6", "N7", "N5", "C3", "N13", "C3", "C3", "C3", "C3", "N6",
             "N10", "C7", "N9", "N14", "N14"])
        mol = reader.mol_from_text(CTABS["wctype_OS"])
        wclogp.assign_wctype(mol)
        self.assertEqual(
            [a.wctype for _, a in mol.atoms_iter()],
            ["S2", "C4", "N13", "P", "C4", "O6", "O5", "O7", "C5", "O9",
             "O3", "C5", "O4", "C23", "S3", "C23", "C24", "C21", "O2a", "O11",
             "C5", "O12", "O10", "O4", "S1", "O2a", "O2a"])
