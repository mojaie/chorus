#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

from chorus.model.atom import Atom


class TestAtom(unittest.TestCase):
    def test_atom(self):
        atom = Atom('C')
        atom.add_hydrogen(4)
        self.assertEqual(atom.composition(), {'C': 1, 'H': 4})
        self.assertAlmostEqual(atom.mw(), 16.043)
        self.assertEqual(atom.formula_html(), "CH<sub>4</sub>")
        self.assertEqual(atom.formula_html(True), "H<sub>4</sub>C")

    def test_add_hydrogen(self):
        atom = Atom('N')
        self.assertEqual(atom.H_donor, False)
        atom.add_hydrogen(3)
        self.assertEqual(atom.H_donor, True)
        atom.charge = 1
        atom.add_hydrogen(4)
        self.assertEqual(atom.composition(), {'N': 1, 'H': 4})
        self.assertAlmostEqual(atom.mw(), 18.039)
        self.assertEqual(atom.formula_html(), "NH<sub>4</sub><sup>+</sup>")
        self.assertEqual(atom.formula_html(1), "<sup>+</sup>H<sub>4</sub>N")

    def test_charge(self):
        atom = Atom('Cu')
        atom.charge = 2
        self.assertEqual(atom.formula_html(), "Cu<sup>2+</sup>")
        atom = Atom('S')
        atom.charge = -2
        self.assertEqual(atom.formula_html(), "S<sup>2–</sup>")
        # en dash, not hyphen-minus
