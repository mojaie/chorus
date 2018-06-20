#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

""" Maximal Common Substructure with Diameter Restriction (MCS-DR)

This calculetes maximum common edge subgraph (MCES) between two molecule
by using the algorithm based on maximum clique detection in the modular
product graph of the two molecule.
Diameter restriction is used for fast and stable local MCS estimation
and for common substructure mapping

"""

import time

import networkx as nx
import numpy as np

from chorus.model.graphmol import Compound
from chorus import molutil
from chorus import remover

try:
    from chorus.cython import mcsdr
    find_cliques = mcsdr.find_cliques
    comparison_graph = mcsdr.comparison_graph
    CYTHON_AVAILABLE = True
except ImportError:
    from networkx.algorithms.clique import find_cliques
    CYTHON_AVAILABLE = False
    try:
        import numexpr as ne
        NUMEXPR_AVAILABLE = True
    except ImportError:
        NUMEXPR_AVAILABLE = False


def reachables(G, root, max_dist):
    visited = {root: 0}
    nbrs = G[root]
    for d in range(max_dist):
        if not nbrs:
            # all nodes are reachable
            return visited
        new_nbrs = {}
        for n in nbrs:
            if n not in visited:
                visited[n] = d + 1
                new_nbrs.update(G[n])
        nbrs = new_nbrs
    # exceed dist limit
    return visited


class DescriptorArray(object):
    """ MCS-DR descriptor array

    Args:
        molecule(chorus.model.graphmol.Compound):  molecule object
        diameter(int): diameter parameter of MCS-DR
        ignore_hydrogen(bool): ignore all hydrogens (implicit and explicit)
        timeout(float): timeout of maximum clique detection in seconds

    Returns:
        dict of the result. array: comparison_array, max_size: max fragment
        size of Glaph-based local similarity, int_to_node: dict of line graph
        node index to original bond in tuple of the atom indices, elapsed_time:
        elapsed time, timeout: underwent timeout or not
    """
    def __init__(self, mol, diameter=8, ignore_hydrogen=True, timeout=5):
        mol.require("Valence")
        self.diam = diameter
        self.ignoreh = ignore_hydrogen
        self.timeout = timeout
        # Results
        self.array = []
        self.max_size = 0
        self.int_to_node = {}
        self.elapsed_time = 0
        self.valid = False
        if len(mol) < 3:
            return
        start_time = time.perf_counter()
        self.mol = mol
        self.preprocess()
        # Generate line graph and reindexing
        lg = nx.line_graph(self.mol.graph)
        node_to_int = {}
        for i, ln in enumerate(lg.nodes()):
            node_to_int[ln] = i
            lg.nodes[ln]["type"] = self.node_desc(ln)
        self.int_to_node = {v: k for k, v in node_to_int.items()}
        self.graph = nx.relabel_nodes(lg, node_to_int)
        # Edges
        edges = []
        for u, v, attr in self.edge_gen():
            edges.append((u, v))
            self.array.append((u, v, self.edge_desc(attr)))
        # Max fragment size determination
        fcres = find_cliques(self.graph.nodes(), edges, timeout=timeout)
        self.max_size = len(fcres["max_clique"])
        self.elapsed_time = round(time.perf_counter() - start_time, 7)
        self.valid = not fcres["timeout"]

    def preprocess(self):
        if self.ignoreh:
            self.mol = molutil.make_Hs_implicit(self.mol)  # clone
        # Ignore salt, water
        remover.remove_salt(self.mol)
        remover.remove_water(self.mol)
        # multivalent coordinated metals notably affect the performance
        remover.remove_coordinated_metal(self.mol)

    def node_desc(self, atoms):
        """default 9 bits descriptor
        7 bits of atomic number (0-127) and 2 bits of pi electrons (0-3)
        """
        a1 = self.mol.atom(atoms[0])
        a2 = self.mol.atom(atoms[1])
        a1t = a1.number << 2 | a1.pi
        a2t = a2.number << 2 | a2.pi
        pair = sorted((a1t, a2t))
        return pair[0] << 9 | pair[1]

    def edge_gen(self):
        for u in self.graph.nodes:
            r = reachables(self.graph, u, self.diam)
            del r[u]
            for v, d in r.items():
                attr = {
                    "dist": d,
                    "umol": self.graph.nodes[u]["type"],
                    "vmol": self.graph.nodes[v]["type"]
                }
                yield (u, v, attr)

    def edge_desc(self, attr):
        """default 42 bits descriptor
        6 bits of distance descriptor (0-31)
        18 bits of bond (9 bits of atom x2) x2
        """
        return (attr["dist"] << 18 | attr["umol"]) << 18 | attr["vmol"]


def comparison_graph_py(arr1, arr2):
    """ DEPRECATED: Generate comparison graph
    Comparison graph is a modular product of molecule edges
    """
    # timeout is not implemented
    u1, v1, c1 = zip(*arr1)
    u2, v2, c2 = zip(*arr2)
    c1 = np.array(c1, dtype=int)
    c2 = np.array(c2, dtype=int)
    product = nx.Graph()
    c1 = c1[:, np.newaxis]  # transpose
    if NUMEXPR_AVAILABLE:
        m = ne.evaluate("c2 == c1")
    else:
        m = c2 == c1
    edges = []
    for x, y in zip(*np.nonzero(m)):
        edges.append({"u1": u1[x], "v1": v1[x], "u2": u2[y], "v2": v2[y]})
    # Graph.add_edges is expensive. Add adjacency dict manually.
    node = {}
    for e in edges:
        node[(e["u1"], e["u2"])] = {}
        node[(e["v1"], e["v2"])] = {}
    adj = node.copy()
    for e in edges:
        adj[(e["u1"], e["u2"])][(e["v1"], e["v2"])] = {}
        adj[(e["v1"], e["v2"])][(e["u1"], e["u2"])] = {}
    product = nx.Graph()
    product.node = node
    product.adj = adj
    return product


if not CYTHON_AVAILABLE:
    comparison_graph = comparison_graph_py


class McsdrGls(object):
    def __init__(self, arr1, arr2, timeout=10):
        # DescriptorArray
        self.max1 = arr1.max_size
        self.max2 = arr2.max_size
        self.map1 = arr1.int_to_node
        self.map2 = arr2.int_to_node
        self.arr1_time = arr1.elapsed_time
        self.arr2_time = arr2.elapsed_time
        # Results
        self.max_clique = []
        self.mod_product_time = None
        self.max_clique_time = None
        self.valid = False
        self.elapsed_time = 0
        if not (arr1.array and arr2.array) or timeout <= 0:
            return
        start_time = time.perf_counter()
        cgout = timeout / 2  # TODO: empirical
        cgres = comparison_graph(arr1.array, arr2.array, timeout=cgout)
        self.mod_product_time = round(cgres["elapsed_time"], 7)
        rest = timeout - cgres["elapsed_time"]
        if rest <= 0:
            self.elapsed_time = round(time.perf_counter() - start_time, 7)
            return
        fcres = find_cliques(
            cgres["decoder"].keys(), cgres["edges"], timeout=rest)
        self.max_clique = fcres["max_clique"]
        self.max_clique_time = round(fcres["elapsed_time"], 7)
        if arr1.valid and arr2.valid \
                and not cgres["timeout"] and not fcres["timeout"]:
            self.valid = True
        self.elapsed_time = round(time.perf_counter() - start_time, 7)

    def edge_count(self):
        return len(self.max_clique)

    def local_sim(self, digit=3):
        ecnt = self.edge_count()
        try:
            sim = ecnt / (self.max1 + self.max2 - ecnt)
        except ZeroDivisionError:
            # if both arr1 and arr2 have no edges, no common structures.
            sim = 0
        return round(sim, digit)

    def common_struct(self, mol1):
        new_mol = Compound()
        edges = [self.map1[n1] for n1, _ in self.max_clique]
        atoms = set()
        for u, v in edges:
            atoms |= {u, v}
            bond = mol1.bond(u, v)
            new_mol.add_bond(u, v, bond)
        for a in atoms:
            new_mol.add_atom(a, mol1.atom(a))
        return new_mol


def from_array(arr1, arr2, timeout=10, gls_cutoff=None, edge_cutoff=None):
    sm, bg = sorted((arr1.max_size, arr2.max_size))
    if not sm:
        return McsdrGls(arr1, arr2, timeout=0)
    if gls_cutoff is not None and gls_cutoff > sm / bg:
        return McsdrGls(arr1, arr2, timeout=0)
    if edge_cutoff is not None and edge_cutoff > sm:
        return McsdrGls(arr1, arr2, timeout=0)
    return McsdrGls(arr1, arr2, timeout)


def from_mol(mol1, mol2, diameter=8, ignore_hydrogen=True,
             timeout=10, arr_timeout=2):
    arr1 = DescriptorArray(mol1, diameter=diameter,
                           ignore_hydrogen=ignore_hydrogen,
                           timeout=arr_timeout)
    arr2 = DescriptorArray(mol2, diameter=diameter,
                           ignore_hydrogen=ignore_hydrogen,
                           timeout=arr_timeout)
    rest = timeout - (arr1.elapsed_time + arr2.elapsed_time)
    return McsdrGls(arr1, arr2, rest)
