#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

from chorus.draw import helper, calc2dcoords
from chorus import smilessupplier
from chorus import v2000reader as reader
from chorus.test.ctabprovider import CTABS


class TestHelper(unittest.TestCase):
    def test_format_double_bond(self):
        m = reader.mol_from_text(CTABS["Phe"])
        helper.scale_and_center(m)
        self.assertEqual(m.bond(8, 11).type, 0)
        self.assertEqual(m.bond(12, 10).type, 0)
        self.assertEqual(m.bond(7, 6).type, 0)
        self.assertEqual(m.bond(2, 9).type, 0)
        helper.format_ring_double_bond(m)
        self.assertEqual(m.bond(8, 11).type, 1)
        self.assertEqual(m.bond(12, 10).type, 0)
        self.assertEqual(m.bond(7, 6).type, 0)
        helper.equalize_terminal_double_bond(m)
        self.assertEqual(m.bond(2, 9).type, 2)

    def test_format_stereo(self):
        m = reader.mol_from_text(CTABS["Phe"])
        m.bond(3, 5).is_lower_first = 1
        m.bond(3, 5).type = 2
        helper.spine_to_terminal_wedge(m)
        self.assertEqual(m.bond(3, 5).is_lower_first, 0)
        self.assertEqual(m.bond(3, 5).type, 1)

    def test_scale_and_center(self):
        m2 = reader.mol_from_text(CTABS["Goserelin"])
        helper.scale_and_center(m2)
        self.assertAlmostEqual(m2.size2d[0], 15.21, 2)
        self.assertAlmostEqual(m2.size2d[1], 16.06, 2)
        self.assertAlmostEqual(m2.size2d[2], 0.83, 2)
        m3 = reader.mol_from_text(CTABS["KCl"])
        helper.scale_and_center(m3)
        self.assertAlmostEqual(m3.size2d[0], 1.0, 2)
        self.assertAlmostEqual(m3.size2d[1], 0, 2)
        self.assertAlmostEqual(m3.size2d[2], 0.71, 2)
        m = smilessupplier.smiles_to_compound("[K+].[Cl-]")  # TODO: overlap
        calc2dcoords.calc2dcoords(m)
        helper.scale_and_center(m)
        self.assertAlmostEqual(m.size2d[0], 0, 2)
        self.assertAlmostEqual(m.size2d[1], 0, 2)
        self.assertAlmostEqual(m.size2d[2], 0, 2)
