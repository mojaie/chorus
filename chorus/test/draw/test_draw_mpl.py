#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

from chorus.draw.matplotlib import Matplotlib
from chorus import v2000reader as reader
from chorus.demo import MOL


class TestDraw(unittest.TestCase):
    def test_draw_mol(self):
        compound = reader.mol_from_text(MOL["demo"])
        # compound = reader.mol_from_text(MOL["Phe"])  # Small
        # compound = reader.mol_from_text(MOL["Goserelin"])  # Large
        # compound = reader.mol_from_text(MOL["CyclosporinA"])  # Atypical bond len
        # compound = reader.mol_from_text(MOL["Carbidopa"])  # Isolated components
        # compound = reader.mol_from_text(MOL["Gadodiamide"])  # Charged
        # compound = reader.mol_from_text(MOL["Premarin"])  # Stereo
        # compound = reader.mol_from_text(MOL["Nitroprusside"])  # Transition metal
        # compound = reader.mol_from_text(MOL["Fondaparinux"])  # Multi-line props
        # compound = reader.mol_from_text(MOL["KCl"])  # No bond, width = 0
        # compound = reader.mol_from_text(MOL["Cyanocobalamin"])
        mpl = Matplotlib(compound)
        # mpl.save("docs/_static/demo.png")

    @unittest.skip("")
    def test_lys_arg(self):
        compound = reader.mol_from_text(MOL["Lys"])
        print("lys")
        mpl = Matplotlib(compound)
        mpl.save("_test_Lys.png")
        compound = reader.mol_from_text(MOL["Arg"])
        print("arg")
        mpl = Matplotlib(compound)
        mpl.save("_test_Arg.png")

    def test_get_size(self):
        compound = reader.mol_from_text(MOL["demo"])
        mpl = Matplotlib(compound)
        self.assertEqual(mpl.get_size(), (733, 733), 2)
