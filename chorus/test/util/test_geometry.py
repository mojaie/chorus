#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import math
import unittest

import chorus.util.geometry as gm


class TestGeometry(unittest.TestCase):
    def test_length(self):
        o = (2, 2)
        p1 = (5, 6)
        p2 = (2, 2)
        self.assertEqual(gm.distance(p1, o), 5)
        self.assertEqual(gm.distance(p2, o), 0)

    def test_rotate(self):
        o = (2, 2)
        p1 = (3, 2)
        self.assertEqual(gm.rotate(p1, math.pi / 2, o), (2, 3))

    def test_cross_product(self):
        v1 = (3, 0)
        v2 = (0, 4)
        self.assertEqual(gm.cross_product(v1, v2), 12)

    def test_dot_product(self):
        v1 = (1, 2)
        v2 = (3, 4)
        self.assertEqual(gm.dot_product(v1, v2), 11)

    def test_interior_angle(self):
        o = (0, 0)
        base = (3, 0)
        v1 = (5, 0)
        v2 = (-20, 0)
        v3 = (3, math.sqrt(3) * -3)
        v4 = (-3, math.sqrt(3) * 3)
        ev = (0, 0)
        self.assertAlmostEqual(gm.interior_angle(base, v1, o), 0)
        self.assertAlmostEqual(gm.interior_angle(base, v2, o), math.pi)
        self.assertAlmostEqual(gm.interior_angle(base, v3, o), math.pi / 3)
        self.assertAlmostEqual(gm.interior_angle(base, v4, o), math.pi / 3 * 2)
        with self.assertRaises(ValueError):
            gm.interior_angle(base, ev, o)

    def test_is_clockwise(self):
        vertices1 = [(0, 1), (1, 0), (0, -1), (-1, 0)]
        vertices2 = [(0, 1), (-1, 0), (0, -1), (1, 0)]
        vertices3 = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
        vertices4 = [(-1, 1), (1, 1), (1, 2), (2, 2), (
            2, -2), (-2, -2), (-2, 2), (-1, 2)]
        self.assertEqual(gm.is_clockwise(vertices1), True)
        self.assertEqual(gm.is_clockwise(vertices2), False)
        with self.assertRaises(ValueError):
            gm.is_clockwise(vertices3)
        self.assertEqual(gm.is_clockwise(vertices4), True)

    def test_move_seg(self):
        p1 = (1, 1)
        p2 = (2, 2)
        seg = gm.m_seg(p1, p2, math.pi / -2, math.sqrt(2) / 2)
        self.assertEqual(seg, ((1.5, 0.5), (2.5, 1.5)))

    def test_trim_seg(self):
        p1 = (1, 1)
        p2 = (5, 5)
        seg_a1 = gm.t_seg(p1, p2, 1 / 4, 1)
        seg_a2 = gm.t_seg(p1, p2, 1 / 4, 2)
        seg_c = gm.t_seg(p1, p2, 1 / 4)
        self.assertEqual(seg_a1, ((1, 1), (4, 4)))
        self.assertEqual(seg_a2, ((2, 2), (5, 5)))
        self.assertEqual(seg_c, ((1.5, 1.5), (4.5, 4.5)))
