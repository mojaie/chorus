#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

from chorus.util import debug
from chorus.model.graphmol import Compound
from chorus.smilessupplier import smiles_to_compound


class TestUtil(unittest.TestCase):
    def setUp(self):
        pass

    @unittest.skip("")
    def test_total_size(self):
        # TODO: memory size depends on the implementation
        # 64 bit
        self.assertEqual(debug.total_size(0), 24)
        self.assertEqual(debug.total_size([]), 64)
        self.assertEqual(debug.total_size([1, 2, 3]), 172)
        # Empty dict size: Python3.5 -> 288, Python3.6 -> 240
        self.assertEqual(debug.total_size({}), 240)
        # Python3.5 -> 471, Python3.6 -> 423
        self.assertEqual(debug.total_size({"hoge": "fuga", "piyo": 0}), 423)
        self.assertEqual(debug.total_size(((1, 2), (3, 4))), 304)
        # Python3.5 -> 2662, Python3.6 -> 2470
        self.assertEqual(debug.total_size(Compound()), 2470)
        c = smiles_to_compound("C1CCCCC1CNC=O")
        # Python3.5 -> 16530, Python3.6 -> 15330
        self.assertEqual(debug.total_size(c), 15330)
