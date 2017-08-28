#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import unittest

import networkx as nx

from chorus.test.ctabprovider import CTABS
from chorus import v2000reader as reader
from chorus.smilessupplier import smiles_to_compound
from chorus import mcsdr
from chorus.draw.svg import SVG
from chorus.util import debug


class TestMCS(unittest.TestCase):
    def test_reachables(self):
        g = nx.Graph([(1, 2), (1, 3), (1, 4), (2, 5), (5, 6),
                     (6, 7), (3, 8), (8, 9), (4, 9)])
        self.assertEqual(mcsdr._reachables(g, 1, 100, 100),
                         {1: 0, 2: 1, 3: 1, 4: 1, 5: 2,
                         6: 3, 7: 4, 8: 2, 9: 2})
        self.assertEqual(mcsdr._reachables(g, 1, 2, 100),
                         {1: 0, 2: 1, 3: 1, 4: 1, 5: 2, 8: 2, 9: 2})

    def test_mcsdr1(self):
        # TODO: pi mismatch is not acceptable
        mol1 = reader.mol_from_text(CTABS["Phe"])
        mol2 = reader.mol_from_text(CTABS["Arg"])
        arr1 = mcsdr.comparison_array(mol1)
        arr2 = mcsdr.comparison_array(mol2)
        self.assertEqual(mcsdr.local_sim(arr1, arr2)["mcsdr_edges"], 5)
        # Delta-y exchange will not occur due to distance descriptor
        mol1 = smiles_to_compound("C1OC1CCC(=O)O")
        mol2 = smiles_to_compound("CC(O)CCC(=O)O")
        arr1 = mcsdr.comparison_array(mol1)
        arr2 = mcsdr.comparison_array(mol2)
        self.assertEqual(mcsdr.local_sim(arr1, arr2)["mcsdr_edges"], 7)

    def test_mcsdr2(self):
        # Disconnected
        mol1 = smiles_to_compound("C1CCCC1CCCC(=O)O")
        mol2 = reader.mol_from_text(CTABS["CaAcO2"])
        arr1 = mcsdr.comparison_array(mol1)
        arr2 = mcsdr.comparison_array(mol2)
        self.assertEqual(mcsdr.local_sim(arr1, arr2)["mcsdr_edges"], 3)
        # No line graph
        mol1 = smiles_to_compound("CO")
        mol2 = smiles_to_compound("CC")
        arr1 = mcsdr.comparison_array(mol1)
        arr2 = mcsdr.comparison_array(mol2)
        self.assertEqual(mcsdr.local_sim(arr1, arr2)["mcsdr_edges"], 0)

    @unittest.skip("")
    @debug.profile
    def test_mcsperformance(self):
        mol1 = reader.mol_from_text(CTABS["Fondaparinux"])
        mol2 = reader.mol_from_text(CTABS["Goserelin"])
        arr1 = mcsdr.comparison_array(mol1)
        arr2 = mcsdr.comparison_array(mol2)
        print(mcsdr.find_mcs(arr1, arr2))

    @unittest.skip("")
    def test_clique_dist(self):
        sq = "./datasource/DrugBank_FDA_Approved.sqlite3"
        from cheddar.data.sqliteconnection import CON
        CON.connect(sq, "DrugBank_FDA_Approved")
        mcol = CON.columns.index("Mol_Block")
        d = {}
        # small
        d["s1"] = CON.find_by("DRUGBANK_ID", "DB00120")
        d["s2"] = CON.find_by("DRUGBANK_ID", "DB00968")
        # mid
        d["m1"] = CON.find_by("DRUGBANK_ID", "DB00279")
        d["m2"] = CON.find_by("DRUGBANK_ID", "DB00451")
        # large
        d["l1"] = CON.find_by("DRUGBANK_ID", "DB00881")
        d["l2"] = CON.find_by("DRUGBANK_ID", "DB00691")
        # cyclic polypeptide
        d["c1"] = CON.find_by("DRUGBANK_ID", "DB00093")
        d["c2"] = CON.find_by("DRUGBANK_ID", "DB00035")
        # porphyrin
        d["p1"] = CON.find_by("DRUGBANK_ID", "DB00115")
        d["p2"] = CON.find_by("DRUGBANK_ID", "DB00200")
        # long chain sugar
        d["inu"] = CON.find_by("DRUGBANK_ID", "DB00638")

        d["s1"] = mcsdr.comparison_array(reader.mol_from_text(d["s1"][mcol]))
        d["s2"] = mcsdr.comparison_array(reader.mol_from_text(d["s2"][mcol]))
        print(len(d["s1"][0]), d["s1"][1], len(d["s2"][0]), d["s2"][1])
        print(mcsdr.mcs_score(d["s1"], d["s2"]))
        # print(mcs_dist(d["s1"], d["s2"]))
        d["m1"] = mcsdr.comparison_array(reader.mol_from_text(d["m1"][mcol]))
        d["m2"] = mcsdr.comparison_array(reader.mol_from_text(d["m2"][mcol]))
        print(len(d["m1"][0]), d["m1"][1], len(d["m2"][0]), d["m2"][1])
        print(mcsdr.mcs_score(d["m1"], d["m2"]))
        # print(mcs_dist(d["m1"], d["m2"]))
        d["l1"] = mcsdr.comparison_array(reader.mol_from_text(d["l1"][mcol]))
        d["l2"] = mcsdr.comparison_array(reader.mol_from_text(d["l2"][mcol]))
        print(len(d["l1"][0]), d["l1"][1], len(d["l2"][0]), d["l2"][1])
        print(mcsdr.mcs_score(d["l1"], d["l2"]))
        # print(mcs_dist(d["l1"], d["l2"]))
        d["c1"] = mcsdr.comparison_array(reader.mol_from_text(d["c1"][mcol]))
        d["c2"] = mcsdr.comparison_array(reader.mol_from_text(d["c2"][mcol]))
        print(len(d["c1"][0]), d["c1"][1], len(d["c2"][0]), d["c2"][1])
        print(mcsdr.mcs_score(d["c1"], d["c2"]))
        # print(mcs_dist(d["c1"], d["c2"]))
        dp1 = mcsdr.comparison_array(reader.mol_from_text(d["p1"][mcol]))
        dp2 = mcsdr.comparison_array(reader.mol_from_text(d["p2"][mcol]))
        print(len(dp1[0]), dp1[1], len(dp2[0]), dp2[1])
        print(mcsdr.mcs_score(dp1, dp2))
        """
        m = reader.mol_from_text(d["c1"][mcol])
        svg = SVG(m)
        svg.save("hoge.svg")
        dp1 = comparison_array(reader.mol_from_text(d["c1"][mcol]))
        dp2 = comparison_array(reader.mol_from_text(d["c2"][mcol]))
        print(len(dp1[0]), dp1[1], len(dp2[0]), dp2[1])
        # print(mcs_score(d["p1"], d["p2"]))
        s = SVG(mcs_structure(reader.mol_from_text(d["c1"][mcol]),
                              reader.mol_from_text(d["c2"][mcol])))
        s.save("fuga.svg")
        """
        # print(mcs_dist(d["p1"], d["p2"]))
        # print(len(d["inu"][0]), d["inu"][1])
