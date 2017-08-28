#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import copy
import math
import unittest

from chorus.test.ctabprovider import CTABS
from chorus import v2000reader as reader
from chorus.draw.svg import SVG
from chorus.draw import calc2dcoords
from chorus import smilessupplier


class TestCalc2DCoords(unittest.TestCase):
    def assert_coords(self, res, expected):
        for k, v in expected.items():
            for i, n in enumerate(res[k][:3]):
                self.assertAlmostEqual(n, v[i], 3)

    @unittest.skip("")
    def draw_chain(self, res):
        mol = smilessupplier.smiles_to_compound("C" * len(res))
        for k, a in mol.atoms_iter():
            mol.atom(k).coords = res[k - 1][:2]
        svg = SVG(mol)
        svg.save("_test_result.svg")

    @unittest.skip("")
    def draw_test(self, mol):
        svg = SVG(mol)
        svg.save("_test_result.svg")

    @unittest.skip("")
    def test_scaffold(self):
        scaffold = {0: [0, 0, 0, 1]}
        calc2dcoords.attach_spiro(scaffold, [0, 1, 2, 3, 4, 5], 0)
        expected = {
            0: [0, 0, 3.142], 1: [0.5, 0.866, 2.094],
            2: [1.5, 0.866, 1.047], 3: [2.0, -0.0, 0],
            4: [1.5, -0.866, 5.236], 5: [0.5, -0.866, 4.189]}
        self.assert_coords(scaffold, expected)
        scaffold2 = copy.deepcopy(scaffold)
        calc2dcoords.attach_fused(scaffold2, [0, 1, 6, 7], (0, 1))
        expected = {
            0: [0, 0, 4.189], 1: [0.5, 0.866, 1.047],
            2: [1.5, 0.866, 1.047], 3: [2.0, -0.0, 0],
            4: [1.5, -0.866, 5.236], 5: [0.5, -0.866, 4.189],
            6: [-0.366, 1.366, 1.833], 7: [-0.866, 0.5, 3.403]
            }
        self.assert_coords(scaffold2, expected)
        # self.draw_chain(scaffold2)
        scaffold3 = copy.deepcopy(scaffold)
        calc2dcoords.attach_fused(scaffold3, [5, 4, 6, 7, 8, 9, 10], (4, 5))
        # self.draw_chain(scaffold3)

    def test_bridge_angle(self):
        step, angle = calc2dcoords.bridge_angle(1, 2)
        self.assertAlmostEqual(math.degrees(step), -90, 3)
        self.assertAlmostEqual(math.degrees(angle), 90, 3)
        step, angle = calc2dcoords.bridge_angle(1 + math.sqrt(3), 2)
        self.assertAlmostEqual(math.degrees(step), -30, 3)
        self.assertAlmostEqual(math.degrees(angle), 30, 3)
        step, angle = calc2dcoords.bridge_angle(1, 3)
        self.assertAlmostEqual(math.degrees(step), -72, 3)
        self.assertAlmostEqual(math.degrees(angle), 108, 3)

    # @unittest.skip("")
    def test_coords(self):
        # mol = reader.mol_from_text(CTABS["Phe"])
        # TODO: without scaffold (Linoleic acid)
        mol = smilessupplier.smiles_to_compound("CCCCCC=CCC=CCCCCCCCC(=O)O")
        # mol = reader.mol_from_text(CTABS["KCl"])  # TODO: overlap
        # mol = reader.mol_from_text(CTABS["Carbidopa"])  # TODO: quart overlap
        # mol = reader.mol_from_text(CTABS["Formestane"])  # TODO: overlap
        # mol = reader.mol_from_text(CTABS["Ceftazidime"])  # TODO: overlap
        # mol = reader.mol_from_text(CTABS["Daunorubicin"])  # TODO: quart overlap
        # mol = reader.mol_from_text(CTABS["Paclitaxel"])  # TODO: overlap
        # mol = reader.mol_from_text(CTABS["Spinosad"])  # TODO: overlap
        calc2dcoords.calc2dcoords(mol)
        self.draw_test(mol)
