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


def comparison_array(molecule, diameter=8, size=30,
                     ignore_hydrogen=True, timeout=1):
    """ Generate comparison array
    Comparison array consists of node pairs in the graph and a collater.
    42 bit collater
        6 bit of distance (0-31)
        18 bit of bond attribute x2
            9 bit of atom attribute x 2
                7 bit of atom number (0-127)
                2 bit of atom pi(0-3)


    [Sheridan, R.P. and Miller, M.D., J. Chem. Inf. Comput. Sci. 38 (1998) 915]

    Using possible_path_length instead of shortest_path_length
    and find intersection of the distance set each other,
    it is possible to detect roundabout matching path.
    Therefore, exact graph isomorphism can be determined.
    However, set intersection is too costful to do
    in spite of a trivial improvement.

    Args:
        mol: Compound
        cutoff: more distant connection than the value is no longer used
                for comparison matrix graph due to performance reason.

    Returns:
        arr(list): comparison array. list of tuple (node1, node2, collater)
        max_mcs(int): maximum size of possible mcs
        int_to_node(dict): int -> node index pair reconverter dict
    Throws:
        ValueError: if len(mol) < 3
    """
    molecule.require("Valence")
    mol = molutil.clone(molecule)
    if ignore_hydrogen:
        mol = molutil.make_Hs_implicit(mol)
    # Ignore salt, water
    remover.remove_salt(mol)
    remover.remove_water(mol)
    # multivalent coordinated metals notably affect the performance
    remover.remove_coordinated_metal(mol)
    g = nx.line_graph(mol.graph)
    node_to_int = {}
    for i, e in enumerate(g.nodes()):
        node_to_int[e] = i
        a1 = mol.atom(e[0])
        a2 = mol.atom(e[1])
        a1t = a1.number << 2 | a1.pi
        a2t = a2.number << 2 | a2.pi
        pair = sorted((a1t, a2t))
        g.node[e]["type"] = pair[0] << 9 | pair[1]
    # convert node index pair to integer expression
    g = nx.relabel_nodes(g, node_to_int)
    # interger -> index pair reconverter
    int_to_node = {v: k for k, v in node_to_int.items()}
    arr = []
    matrix = nx.Graph()
    for ui, ua in g.nodes(data=True):
        r = _reachables(g, ui, diameter, size)
        for vi, d in r.items():
            if not d:
                continue
            matrix.add_edge(ui, vi)
            code = (d << 18 | ua["type"]) << 18 | g.node[vi]["type"]
            arr.append((ui, vi, code))
    cliques, elapsed = find_cliques(matrix, timeout=timeout)
    if elapsed > timeout:
        raise RuntimeError("Max fragment determination has timed out")
    max_size = len(max(cliques, key=len, default=[]))
    return arr, max_size, int_to_node, round(elapsed, 5)


def _reachables(G, root, max_dist, max_size):
    # TODO: is the BFS tree deterministic ?
    # TODO: complete exploration of the started level even if max_size exceeded
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
            if len(visited) == max_size:
                # exceed size limit
                return visited
        nbrs = new_nbrs
    # exceed dist limit
    return visited


def comparison_graph_py(arr1, arr2):
    """ Generate comparison graph
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
    def __init__(self, arr1, arr2, timeout):
        self.max1 = arr1[1]
        self.max2 = arr2[1]
        self.map1 = arr1[2]
        self.map2 = arr2[2]
        self.max_clique = []
        self.cg_elapsed = 0
        self.cl_elapsed = 0
        if (arr1[0] and arr2[0]):
            cg, elapsed = comparison_graph(arr1[0], arr2[0], timeout=timeout)
            rest = timeout - elapsed
            cliques, elapsed2 = find_cliques(cg, timeout=rest)
            self.max_clique = max(cliques, key=len, default=[])
            self.cg_elapsed = elapsed
            self.cl_elapsed = elapsed2

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

    def exec_time(self):
        return {
            "comparison_graph": round(self.cg_elapsed, 5),
            "find_cliques": round(self.cl_elapsed, 5),
            "total": round(self.cg_elapsed + self.cl_elapsed, 5)
        }


def from_array(arr1, arr2, timeout=10):
    return McsdrGls(arr1, arr2, timeout)


def from_mol(mol1, mol2, timeout=10):
    arr1 = comparison_array(mol1)
    arr2 = comparison_array(mol2)
    return from_array(arr1, arr2, timeout)
