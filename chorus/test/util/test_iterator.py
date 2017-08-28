#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

import itertools
import chorus.util.iterator as it


class TestIterator(unittest.TestCase):

    def test_consecutive(self):
        self.assertEqual("".join(list(it.consecutive("ABCDEF", 3))[3]), "DEF")
        self.assertEqual(len(list(it.consecutive("", 3))), 0)
        self.assertEqual(len(list(it.consecutive("ABCDEF", 0))), 0)
        self.assertEqual("".join(list(it.consecutive("ABCDEF", 6))[0]),
                         "ABCDEF")
        self.assertEqual(len(list(it.consecutive("ABCDEF", 7))), 0)
        # Loop
        looped = it.consecutive(itertools.cycle("ABCDEF"), 3)
        self.assertEqual(
            "".join(list(itertools.islice(looped, 4, 5))[0]), "EFA")

    def test_repeat_each(self):
        self.assertEqual("".join(list(it.repeat_each("ABCDEF", 3))),
                         "AAABBBCCCDDDEEEFFF")
        self.assertEqual(len(list(it.repeat_each("", 3))), 0)
        self.assertEqual(len(list(it.repeat_each("ABCDEF", 0))), 0)

    def test_chunk(self):
        self.assertEqual(len(list(it.chunk("ABCDEFGHI", 3))), 3)
        self.assertEqual(len(list(it.chunk("ABCDEFGHIJ", 3))), 4)
        self.assertEqual(len(list(it.chunk("", 3))), 0)
        self.assertEqual(len(list(it.chunk("ABCDEFGHIJ", 0))), 0)
        self.assertEqual(len(list(it.chunk("ABCDEFGHIJ", 10))), 1)
        self.assertEqual(len(list(it.chunk("ABCDEFGHIJ", 11))), 1)

    def test_chunk_truncate(self):
        self.assertEqual(len(list(it.chunk_truncate("ABCDEFGHI", 3))), 3)
        self.assertEqual(len(list(it.chunk_truncate("ABCDEFGHIJ", 3))), 3)
        self.assertEqual(len(list(it.chunk_truncate("", 3))), 0)
        self.assertEqual(len(list(it.chunk_truncate("ABCDEFGHIJ", 0))), 0)
        self.assertEqual(len(list(it.chunk_truncate("ABCDEFGHIJ", 10))), 1)
        self.assertEqual(len(list(it.chunk_truncate("ABCDEFGHIJ", 11))), 0)

    def test_chunk_fill(self):
        self.assertEqual(len(list(it.chunk_fill("ABCDEF", 3, 'x'))), 2)
        self.assertEqual(
            "".join(list(it.chunk_fill("ABCDEFG", 3, 'x'))[2]), "Gxx")
        self.assertEqual(len(list(it.chunk_fill("", 3, 'x'))), 0)
        self.assertEqual(len(list(it.chunk_fill("ABCDEF", 0, 'x'))), 0)
        self.assertEqual(len(list(it.chunk_fill("ABCDEF", 6, 'x'))), 1)
        self.assertEqual(
            "".join(list(it.chunk_fill("ABCDEF", 7, 'x'))[0]), "ABCDEFx")
