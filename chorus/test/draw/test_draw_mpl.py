#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

from chorus.draw.matplotlib import Matplotlib
from chorus import v2000reader as reader
from chorus.test.ctabprovider import CTABS


class TestDraw(unittest.TestCase):
    def test_draw_mol(self):
        compound = reader.mol_from_text(CTABS["demo"])
        # compound = reader.mol_from_text(CTABS["Phe"])  # Small
        # compound = reader.mol_from_text(CTABS["Goserelin"])  # Large
        # compound = reader.mol_from_text(CTABS["CyclosporinA"])  # Atypical bond len
        # compound = reader.mol_from_text(CTABS["Carbidopa"])  # Isolated components
        # compound = reader.mol_from_text(CTABS["Gadodiamide"])  # Charged
        # compound = reader.mol_from_text(CTABS["Premarin"])  # Stereo
        # compound = reader.mol_from_text(CTABS["Nitroprusside"])  # Transition metal
        # compound = reader.mol_from_text(CTABS["Fondaparinux"])  # Multi-line props
        # compound = reader.mol_from_text(CTABS["KCl"])  # No bond, width = 0
        # compound = reader.mol_from_text(CTABS["Cyanocobalamin"])
        mpl = Matplotlib(compound)
        # mpl.save("docs/_static/demo.png")

    @unittest.skip("")
    def test_lys_arg(self):
        compound = reader.mol_from_text(CTABS["Lys"])
        print("lys")
        mpl = Matplotlib(compound)
        mpl.save("_test_Lys.png")
        compound = reader.mol_from_text(CTABS["Arg"])
        print("arg")
        mpl = Matplotlib(compound)
        mpl.save("_test_Arg.png")

    def test_get_size(self):
        compound = reader.mol_from_text(CTABS["demo"])
        mpl = Matplotlib(compound)
        self.assertEqual(mpl.get_size(), (733, 733), 2)
