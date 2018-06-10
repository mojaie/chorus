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


def preprocess(mol, ignore_hydrogen):
    if ignore_hydrogen:
        m = molutil.make_Hs_implicit(mol)  # clone
    # Ignore salt, water
    remover.remove_salt(m)
    remover.remove_water(m)
    # multivalent coordinated metals notably affect the performance
    remover.remove_coordinated_metal(m)
    return m


def node_desc(atom1, atom2):
    a1t = atom1.number << 2 | atom1.pi
    a2t = atom2.number << 2 | atom2.pi
    pair = sorted((a1t, a2t))
    return pair[0] << 9 | pair[1]


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


def edge_gen(G, diam):
    for u in G.nodes:
        r = reachables(G, u, diam)
        del r[u]
        for v, d in r.items():
            attr = {
                "dist": d,
                "umol": G.nodes[u]["type"],
                "vmol": G.nodes[v]["type"]
            }
            yield (u, v, attr)


def edge_desc(attr):
    return (attr["dist"] << 18 | attr["umol"]) << 18 | attr["vmol"]


def comparison_array(mol, diameter=8, ignore_hydrogen=True, timeout=5):
    """ Generate comparison array
    Comparison array consists of node pairs in the graph and a collater.
    42 bit collater
        6 bit of distance (0-31)
        18 bit of bond attribute x2
            9 bit of atom attribute x 2
                7 bit of atom number (0-127)
                2 bit of atom pi(0-3)

    Reference:
    [Sheridan, R.P. and Miller, M.D., J. Chem. Inf. Comput. Sci. 38 (1998) 915]

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
    Throws:
        ValueError: if len(mol) < 3
    """
    mol.require("Valence")
    res = {
        "array": [],
        "max_size": 0,
        "int_to_node": {},
        "elapsed_time": 0,
        "valid": False
    }
    if len(mol) < 3:
        return res
    start_time = time.perf_counter()
    mol = preprocess(mol, ignore_hydrogen)
    # Generate line graph and reindexing
    lg = nx.line_graph(mol.graph)
    node_to_int = {}
    for i, ln in enumerate(lg.nodes()):
        node_to_int[ln] = i
        lg.nodes[ln]["type"] = node_desc(mol.atom(ln[0]), mol.atom(ln[1]))
    int_to_node = {v: k for k, v in node_to_int.items()}
    g = nx.relabel_nodes(lg, node_to_int)
    # Edges
    edges = []
    for u, v, attr in edge_gen(g, diameter):
        edges.append((u, v))
        res["array"].append((u, v, edge_desc(attr)))
    # Max fragment size determination
    fcres = find_cliques(g.nodes(), edges, timeout=timeout)
    res["int_to_node"] = int_to_node
    res["max_size"] = len(fcres["max_clique"])
    res["elapsed_time"] = time.perf_counter() - start_time
    res["valid"] = not fcres["timeout"]
    return res


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
        self.max1 = arr1["max_size"]
        self.max2 = arr2["max_size"]
        self.map1 = arr1["int_to_node"]
        self.map2 = arr2["int_to_node"]
        self.max_clique = []
        self.perf = {
            "arr1_time": round(arr1["elapsed_time"], 5),
            "arr2_time": round(arr2["elapsed_time"], 5),
            "mod_product_time": None,
            "max_clique_time": None,
            "valid": False,
        }
        if not (arr1["array"] and arr2["array"]) or not timeout:
            return
        cgout = timeout / 2  # TODO: empirical
        cgres = comparison_graph(arr1["array"], arr2["array"], timeout=cgout)
        self.perf["mod_product_time"] = round(cgres["elapsed_time"], 5)
        rest = timeout - cgres["elapsed_time"]
        fcres = find_cliques(
            cgres["decoder"].keys(), cgres["edges"], timeout=rest)
        self.max_clique = fcres["max_clique"]
        self.perf["max_clique_time"] = round(fcres["elapsed_time"], 5)
        if arr1["valid"] and arr2["valid"] \
                and not cgres["timeout"] and not fcres["timeout"]:
            self.perf["valid"] = True

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
    sm, bg = sorted((arr1["max_size"], arr2["max_size"]))
    if not sm:
        return McsdrGls(arr1, arr2, timeout=0)
    if gls_cutoff is not None and gls_cutoff > sm / bg:
        return McsdrGls(arr1, arr2, timeout=0)
    if edge_cutoff is not None and edge_cutoff > sm:
        return McsdrGls(arr1, arr2, timeout=0)
    return McsdrGls(arr1, arr2, timeout)


def from_mol(mol1, mol2, diameter=8, ignore_hydrogen=True,
             timeout=10, arr_timeout=2):
    arr1 = comparison_array(mol1, diameter=diameter,
                            ignore_hydrogen=ignore_hydrogen,
                            timeout=arr_timeout)
    arr2 = comparison_array(mol2, diameter=diameter,
                            ignore_hydrogen=ignore_hydrogen,
                            timeout=arr_timeout)
    return McsdrGls(arr1, arr2, timeout)
