#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

from chorus.smilessupplier import smiles_to_compound
from chorus import v2000reader as reader
from chorus.test.ctabprovider import CTABS
from chorus import descriptor
from chorus import molutil


class Testdescriptor(unittest.TestCase):
    def test_assign_pai(self):
        mol = smiles_to_compound("CC=CC#C")
        self.assertEqual([a.pi for _, a in mol.atoms_iter()],
                         [0, 1, 1, 2, 2])
        mol = smiles_to_compound("C1=CC=CN1")
        self.assertEqual([a.pi for _, a in mol.atoms_iter()],
                         [1, 1, 1, 1, 0])
        mol = smiles_to_compound("CC(=O)N=N=N")
        self.assertEqual([a.pi for _, a in mol.atoms_iter()],
                         [0, 1, 1, 1, 2, 1])

    def test_assign_hydrogen(self):
        mol = smiles_to_compound("CCC(=O)N")
        self.assertTrue(mol.atom(4).H_acceptor)
        self.assertTrue(mol.atom(5).H_donor)
        self.assertTrue(mol.atom(5).H_acceptor)
        mol = smiles_to_compound("CCN(CO)CF")
        self.assertFalse(mol.atom(3).H_donor)
        self.assertTrue(mol.atom(5).H_donor)
        self.assertTrue(mol.atom(5).H_acceptor)
        self.assertTrue(mol.atom(7).H_acceptor)

    def test_assign_rotatable(self):
        mol = reader.mol_from_text(CTABS["Phe"])
        descriptor.assign_rotatable(mol)
        self.assertEqual(molutil.rotatable_count(mol), 3)
        mol = reader.mol_from_text(CTABS["KCl"])
        self.assertEqual(molutil.rotatable_count(mol), 0)
        mol = reader.mol_from_text(CTABS["Dipyridamole"])
        self.assertEqual(molutil.rotatable_count(mol), 12)
        mol = reader.mol_from_text(CTABS["Paclitaxel"])
        self.assertEqual(molutil.rotatable_count(mol), 15)

    def test_assign_aromatic(self):
        # Furan is aromatic
        mol = smiles_to_compound("O1C=CC=C1")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 1, 1, 1, 1])
        # Quinone is not aromatic
        mol = smiles_to_compound("C1(=O)C=CC(=O)C=C1")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [0, 0, 0, 0, 0, 0, 0, 0])
        # Tropone is aromatic
        mol = smiles_to_compound("C1(=O)C=CC=CC=C1")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 0, 1, 1, 1, 1, 1, 1])
        # Azepine is not aromatic
        mol = smiles_to_compound("N1C=CC=CC=C1")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [0, 0, 0, 0, 0, 0, 0])
        # 1,2-Dihydro-1,2-azaborine is aromatic
        mol = smiles_to_compound("C1=CC=CBN1")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 1, 1, 1, 1, 1])
        # Fulvene is not aromatic
        mol = smiles_to_compound("C1=CC=CC1=C")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [0, 0, 0, 0, 0, 0])
        # Cyclopropenyl cation is aromatic
        mol = smiles_to_compound("C1=C[C+]1")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 1, 1])
        # 1,3,5-Hexatriene is not aromatic
        mol = smiles_to_compound("C=CC=CC=C")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [0, 0, 0, 0, 0, 0])
        # Borazine is aromatic
        mol = smiles_to_compound("B1NBNBN1")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 1, 1, 1, 1, 1])
        # 1,2-Methylenedioxybenzene is not aromatic
        mol = smiles_to_compound("C1=CC=CC2=C1OCO2")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 1, 1, 1, 1, 1, 0, 0, 0])
        # Naphthalene is aromatic
        mol = smiles_to_compound("C=1C=CC=C2C1C=CC=C2")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        # Coumarin is aromatic
        mol = smiles_to_compound("C1=CC(=O)OC2=C1C=CC=C2")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1])
        # Pyrene is aromatic
        mol = smiles_to_compound("C12=CC=C3C=CC=C4C=CC(C2=C34)=CC=C1")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        # Caffeine is aromatic
        mol = smiles_to_compound("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0])
        # Pyromellitimide is not aromatic
        mol = smiles_to_compound("C1(=O)NC(=O)C=2C1=CC=3C(=O)NC(=O)C=3C=2")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1])
        # TODO: Azulene is not aromatic
        mol = smiles_to_compound("C=1C=CC=2C1C=CC=CC2")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        # TODO: o-Quinodimethane is aromatic
        mol = smiles_to_compound("C1=CC=CC(=C)C1(=C)")
        self.assertEqual([a.aromatic for _, a in mol.atoms_iter()],
                         [1, 1, 1, 1, 1, 0, 1, 0])

        # Discussion:
        # [10]annulene is not aromatic due to steric strain.
        # [14]annulene, [18]annulene are aromatic.
        # Cyclononatetraenyl anion may have aromaticity
        # Cyclononatetraenyl cation may have Moebius aromaticity
        # Theoretical aromaticity of Heteronines: Azo > Thio > oxo
        # Pyromellitimide ring system has 10 pi electrons but not aromatic.
        # (Its benzene moiety is of course aromatic)
    """
    def test_assign_charge(self):
        # Oxoacids
        mol = smiles_to_compound("CC(=O)O")
        descriptor.assign_charge(mol)
        self.assertEqual([a.charge_phys for _, a in mol.atoms_iter()],
                         [0, 0, 0, -1])
        self.assertEqual([a.charge_conj for _, a in mol.atoms_iter()],
                         [0, 0, -1, 0])
        mol = smiles_to_compound("CP(=O)(=O)O")
        descriptor.assign_charge(mol)
        self.assertEqual([a.charge_phys for _, a in mol.atoms_iter()],
                         [0, 0, 0, 0, -1])
        self.assertEqual([a.charge_conj for _, a in mol.atoms_iter()],
                         [0, 0, -1, -1, 0])
        # Amidine
        mol = smiles_to_compound("CC(=N)N")
        descriptor.assign_charge(mol)
        self.assertEqual([a.charge_phys for _, a in mol.atoms_iter()],
                         [0, 0, 0, 1])
        self.assertEqual([a.charge_conj for _, a in mol.atoms_iter()],
                         [0, 0, 1, 0])
        # Guanidine
        mol = smiles_to_compound("CNC(=N)N")
        descriptor.assign_charge(mol)
        self.assertEqual([a.charge_phys for _, a in mol.atoms_iter()],
                         [0, 0, 0, 0, 1])
        self.assertEqual([a.charge_conj for _, a in mol.atoms_iter()],
                         [0, 1, 0, 1, 0])
        # thiophenol
        mol = smiles_to_compound("C1(S)=CC=CC=C1")
        descriptor.assign_charge(mol)
        self.assertEqual([a.charge_phys for _, a in mol.atoms_iter()],
                         [0, -1, 0, 0, 0, 0, 0])

    def test_assign_type(self):
        # Arginine
        mol = smiles_to_compound("C(N)(C(=O)O)CCCNC(=N)N")
        descriptor.assign_type(mol)
        self.assertEqual([a.type for _, a in mol.atoms_iter()],
                         [6, 1, 7, 2, 2, 6, 6, 6, 1, 7, 1, 1])
        # Carbosone
        mol = smiles_to_compound("NC(=O)NC1=CC=C(C=C1)[As](=O)(O)O")
        descriptor.assign_type(mol)
        self.assertEqual([a.type for _, a in mol.atoms_iter()],
                         [3, 7, 4, 3, 7, 6, 6, 6, 6, 6, 7, 2, 2, 2])
        # Chloramphenicol
        mol = smiles_to_compound("ON(=O)C1=CC=C(C=C1)C(O)C(CO)NC(=O)C(Cl)Cl")
        descriptor.assign_type(mol)
        self.assertEqual([a.type for _, a in mol.atoms_iter()],
                         [4, 7, 4, 7, 6, 6, 6, 6, 6, 6, 5,
                          6, 6, 5, 3, 7, 4, 6, 6, 6])
        # TODO: Citrinin
        # mol = smiles_to_compound("OC(=O)C1=C(O)C(=COC(C)C2(C))C2=C(C)C1(=O)")
        # descriptor.assign_type(mol)
        # self.assertEqual([a.type for _, a in mol.atoms_iter()],
        #                 [2, 7, 2, 6, 7, 5, 6, 6, 4, 6, 6,
        #                  6, 6, 6, 6, 6, 7, 5])
        # Cyclophosphamide
        mol = smiles_to_compound("C1OCCNP1(=O)N(CCCl)CCCl")
        descriptor.assign_type(mol)
        self.assertEqual([a.type for _, a in mol.atoms_iter()],
                         [6, 4, 6, 6, 3, 7, 4, 7, 6, 6,
                          6, 6, 6, 6])
        # Cytidine
        mol = smiles_to_compound("OCC(O1)C(O)C(O)C1N2C=CC(N)=NC2(=O)")
        descriptor.assign_type(mol)
        self.assertEqual([a.type for _, a in mol.atoms_iter()],
                         [5, 6, 6, 4, 6, 5, 6, 5, 6, 7, 6, 6, 7,
                          3, 4, 7, 4])
        # Diazepam
        mol = smiles_to_compound("C1=CC(Cl)=CC2=C1N(C)C(=O)CN=C2C3C=CC=CC=3")
        descriptor.assign_type(mol)
        self.assertEqual([a.type for _, a in mol.atoms_iter()],
                         [6, 6, 6, 6, 6, 6, 6, 7, 6, 7, 4, 6, 4, 7,
                          6, 6, 6, 6, 6, 6])
        # Diazinon
        mol = smiles_to_compound("CC(C)C1N=C(C)C=C(N=1)OP(=S)(OCC)OCC")
        descriptor.assign_type(mol)
        self.assertEqual([a.type for _, a in mol.atoms_iter()],
                         [6, 6, 6, 6, 4, 6, 6, 6, 7, 4, 4, 7, 6,
                          4, 6, 6, 4, 6, 6])
        """
