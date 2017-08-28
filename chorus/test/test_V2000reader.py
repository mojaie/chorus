#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import os
import unittest

from chorus.test.ctabprovider import CTABS
from chorus.util import debug
import chorus.v2000reader as reader


class TestV2000Reader(unittest.TestCase):
    # @debug.profile
    @unittest.skip("")
    def test_sdf(self):
        from chorus.topology import recognize
        f = os.path.join(os.path.dirname(__file__),
                         "DrugBank_All.sdf")
        for i, mol in enumerate(reader.file_to_mols(f)):
            recognize(mol)

    @unittest.skip("It takes long time.")
    def test_option_labels(self):
        f = os.path.join(os.path.dirname(__file__),
                         "DrugBank_FDA_Approved.sdf")
        labels, count = reader.inspect_file(f)
        self.assertEqual(len(labels), 37)
        self.assertEqual(count, 1543)

    def test_parse_options(self):
        # Space should be accepted
        parsed = reader.optional_data((">  <hoge fuga>", "piyo"))
        self.assertEqual(parsed, {"hoge fuga": "piyo"})

    def test_null(self):
        mol = reader.mol_from_text(CTABS["null"])
        self.assertTrue(len(mol.descriptors))  # descriptors assigned
        self.assertEqual(mol.atom_count(), 0)
        self.assertEqual(mol.bond_count(), 0)

    def test_phe(self):
        compound = reader.mol_from_text(CTABS["Phe"])
        self.assertEqual(compound.atom_count(), 12)
        self.assertEqual(compound.bond_count(), 12)
        self.assertEqual(compound.data['GENERIC_NAME'], 'L-Phenylalanine')

    def test_phe_omit(self):
        """ omitted CTAB """
        compound = reader.mol_from_text(CTABS["Phe_omit"])
        self.assertEqual(compound.atom_count(), 12)
        self.assertEqual(compound.bond_count(), 12)
        self.assertEqual(compound.data['GENERIC_NAME'], 'L-Phenylalanine')

    def test_kcl(self):
        """ No bonds """
        compound = reader.mol_from_text(CTABS["KCl"])
        self.assertEqual(compound.atom_count(), 2)
        self.assertEqual(compound.bond_count(), 0)

    def test_colesevelam(self):
        """ Polymer Expression not supported yet """
        with self.assertRaises(ValueError):
            next(reader.mols_from_text(CTABS["Colesevelam"], False))
        with self.assertRaises(ValueError):
            reader.mol_from_text(CTABS["Colesevelam"])

    @debug.mute  # Unsupported symbol: A (#1 in v2000reader)
    def test_no_halt(self):
        mol = next(reader.mols_from_text(CTABS["Colesevelam"]))
        self.assertTrue(len(mol.descriptors))  # descriptors assigned
        self.assertEqual(mol.atom_count(), 0)
        self.assertEqual(mol.data['GENERIC_NAME'], 'Colesevelam')
