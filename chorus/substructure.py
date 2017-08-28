#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import operator
from collections import Counter

from networkx.algorithms.isomorphism.vf2userfunc import GraphMatcher

from chorus import molutil


def atom_match(n1, n2):
    a1 = n1['atom']
    a2 = n2['atom']
    if a1.symbol == a2.symbol and a1.pi == a2.pi:
        return True


def filter_(mol, query, f=operator.ne):
    mol.require("Topology")
    query.require("Topology")
    if f(molutil.composition(mol), molutil.composition(query)):
        return False
    if f(Counter(a.pi for _, a in mol.atoms_iter()),
         Counter(a.pi for _, a in query.atoms_iter())):
        return False
    if f(Counter(len(r) for r in mol.rings),
         Counter(len(r) for r in query.rings)):
        return False
    if f(Counter(len(s) for s in mol.scaffolds),
         Counter(len(s) for s in query.scaffolds)):
        return False
    return True


def equal(mol, query, largest_only=True, ignore_hydrogen=True):
    """ if mol is exactly same structure as the query, return True
    Args:
      mol: Compound
      query: Compound
    """
    m = molutil.clone(mol)
    q = molutil.clone(query)
    if largest_only:
        m = molutil.largest_graph(m)
        q = molutil.largest_graph(q)
    if ignore_hydrogen:
        m = molutil.make_Hs_implicit(m)
        q = molutil.make_Hs_implicit(q)
    if molutil.mw(m) == molutil.mw(q):
        gm = GraphMatcher(q.graph, m.graph, node_match=atom_match)
        return gm.is_isomorphic()
    return False


def substructure(mol, query, largest_only=True, ignore_hydrogen=True):
    """ if mol is a substructure of the query, return True
    Args:
      mol: Compound
      query: Compound
      largest_only: compare only largest graph molecule
    """
    def subset_filter(cnt1, cnt2):
        diff = cnt2
        diff.subtract(cnt1)
        if any(v < 0 for v in diff.values()):
            return True

    if not (len(mol) and len(query)):
        return False  # two blank molecules are not isomorphic
    m = molutil.clone(mol)
    q = molutil.clone(query)
    if largest_only:
        m = molutil.largest_graph(m)
        q = molutil.largest_graph(q)
    if ignore_hydrogen:
        m = molutil.make_Hs_implicit(m)
        q = molutil.make_Hs_implicit(q)
    if filter_(m, q, f=subset_filter):
        gm = GraphMatcher(q.graph, m.graph, node_match=atom_match)
        return gm.subgraph_is_isomorphic()
    return False


def superstructure(mol, query, largest_only=True, ignore_hydrogen=True):
    return substructure(query, mol, largest_only, ignore_hydrogen)
