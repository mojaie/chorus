#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

from chorus.smilessupplier import smiles_to_compound
from chorus import v2000reader as reader
from chorus.test.ctabprovider import CTABS
from chorus import molutil
from chorus.draw.svg import SVG


class TestMolUtil(unittest.TestCase):
    def test_clone(self):
        # this takes 0.003 ms
        m = reader.mol_from_text(CTABS["Indinavir"])
        cp = molutil.clone(m)
        self.assertEqual(len(m), len(cp))
        # deepcopy takes 0.013 ms
        # m = reader.mol_from_text(CTABS["Indinavir"])
        # cp = copy.deepcopy(m)
        # self.assertEqual(len(m), len(cp))

    def test_nullmol(self):
        mol = molutil.null_molecule()
        self.assertTrue(len(SVG(mol).contents()))  # draw.svg should works
        obj = mol.jsonized()
        self.assertEqual(len(obj["atoms"]), 0)
        self.assertEqual(len(obj["connections"]), 0)
        self.assertTrue("Topology" in obj["descriptors"])
        self.assertTrue("Valence" in obj["descriptors"])
        self.assertTrue("Rotatable" in obj["descriptors"])

    def test_mw(self):
        mol = reader.mol_from_text(CTABS["Phe"])
        self.assertAlmostEqual(molutil.mw(mol), 165.19, 2)

    def test_composition(self):
        mol = reader.mol_from_text(CTABS["Phe"])
        self.assertEqual(molutil.composition(mol),
                         {'H': 11, 'C': 9, 'O': 2, 'N': 1})

    def test_formula(self):
        mol = reader.mol_from_text(CTABS["Phe"])
        self.assertEqual(molutil.formula(mol), "C9H11NO2")
        mol = reader.mol_from_text(CTABS["KCl"])
        self.assertEqual(molutil.formula(mol), "Cl.K")  # longer text first
        mol = smiles_to_compound("CCO.O.O")
        self.assertEqual(molutil.formula(mol), "C2H6O.2H2O")
        mol = smiles_to_compound("CCCSCC(Cl)C(O)O")
        self.assertEqual(molutil.formula(mol), "C6H13O2SCl")

    def test_make_Hs_implicit(self):
        mol = smiles_to_compound("C(H)(H)(H)C(H)(H)C(=O)H", False)
        newmol = molutil.make_Hs_implicit(mol)
        self.assertEqual(newmol.atom_count(), 4)
        self.assertEqual(molutil.mw(newmol), 58.08)

    def test_non_hydrogen_count(self):
        mol = smiles_to_compound("C(H)(H)(H)C(H)(H)C(=O)H", False)
        self.assertEqual(molutil.non_hydrogen_count(mol), 4)
