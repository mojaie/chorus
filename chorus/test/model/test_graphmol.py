#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import json
import pickle
import unittest

from chorus import v2000reader as reader
from chorus.test.ctabprovider import CTABS
from chorus.model.graphmol import Compound
from chorus import molutil, substructure


class TestGraphMol(unittest.TestCase):
    def test_pickle(self):
        # Compound is picklable
        m = reader.mol_from_text(CTABS["Phe"])
        dmp = pickle.dumps(m)
        m2 = pickle.loads(dmp)
        self.assertEqual(len(m), len(m2))

    def test_hide_carbon(self):
        m = reader.mol_from_text(CTABS["Phe"])
        self.assertTrue(m.atom(3).visible)
        self.assertFalse(m.atom(6).visible)

    def test_json(self):
        m = reader.mol_from_text(CTABS["Cyanocobalamin"])
        d = m.jsonized()
        j = json.dumps(d)
        # Compressed file size
        """
        from chorus.util import debug
        print("Original: {} Bytes".format(debug.total_size(m)))
        print("Dict: {} Bytes".format(debug.total_size(d)))
        print("JSON: {} Bytes".format(debug.total_size(j)))
        print("Pickled: {} Bytes".format(
            debug.total_size(pickle.dumps(m))))
        print("Pickled dict: {} Bytes".format(
            debug.total_size(pickle.dumps(d))))
        """
        # Results
        """
        Original: 144414 Bytes
        Pickled: 30033 Bytes
        Dict: 195291 Bytes
        JSON: 32373 Bytes
        Pickled dict: 22943 Bytes
        """
        m2 = Compound(json.loads(j))
        # Atom key should be integer, not string
        self.assertIsInstance(next(iter(m2.graph.node.keys())), int)
        self.assertIsInstance(next(iter(m2.graph.adj.keys())), int)
        self.assertEqual(molutil.mw(m), molutil.mw(m2))
        self.assertTrue(substructure.equal(m, m2))
