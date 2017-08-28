# cython: boundscheck=False, wraparound=False, langage_level=3

#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import networkx as nx
from libcpp.vector cimport vector


cdef struct Edge:
    int u1
    int v1
    int u2
    int v2


def comparison_graph(arr1, arr2):
    u1a, v1a, c1a = zip(*arr1)
    u2a, v2a, c2a = zip(*arr2)
    cdef vector[int] u1 = u1a
    cdef vector[int] v1 = v1a
    cdef vector[long long] c1 = c1a
    cdef vector[int] u2 = u2a
    cdef vector[int] v2 = v2a
    cdef vector[long long] c2 = c2a
    cdef int i, j
    cdef int i_max = c1.size()
    cdef int j_max = c2.size()
    cdef Edge e
    cdef vector[Edge] edges
    for i in range(i_max):
        for j in range(j_max):
            if c1[i] == c2[j]:
                e.u1 = u1[i]
                e.v1 = v1[i]
                e.u2 = u2[j]
                e.v2 = v2[j]
                edges.push_back(e)
    # Graph.add_edges is expensive. Add adjacency dict manually.
    node = {}
    for edge in edges:
        node[(edge.u1, edge.u2)] = {}
        node[(edge.v1, edge.v2)] = {}
    adj = node.copy()
    for edge in edges:
        adj[(edge.u1, edge.u2)][(edge.v1, edge.v2)] = {}
        adj[(edge.v1, edge.v2)][(edge.u1, edge.u2)] = {}
    product = nx.Graph()
    product.node = node
    product.adj = adj
    return product


def find_cliques(graph, root=None):
    # tuple index may expensive. convert to integer
    cdef int i
    decode = {i: node for i, node in enumerate(graph)}
    encode = {node: i for i, node in decode.items()}
    adj = {}
    for node, adjs in graph.adjacency_iter():
        eadjs = set(encode[a] for a in adjs)
        adj[encode[node]] = eadjs
    # initialize
    result = []
    R = []
    P = set(adj.keys())
    X = set()
    stack = []
    if root:
        R.append(root)
        P = adj[root] - set(R)
    Pv = P - pivot(P, len(P) - 1, adj)
    if root:
        stack.append((P, X, Pv))
    # DFS loop
    cdef int n, Pcnt, pvdeg
    while True:
        if Pv:
            n = Pv.pop()
        else:
            if stack:
                P, X, Pv = stack.pop()
            else:
                break
            R.pop()
            continue
        P.remove(n)
        X.add(n)
        nbrs = adj[n]
        new_P = P & nbrs
        new_X = X & nbrs
        Pcnt = len(new_P)
        if Pcnt in (0, 1):
            if not new_X:
                if Pcnt:
                    result.append([decode[r] for r in R + [n] + list(new_P)])
                else:
                    result.append([decode[r] for r in R + [n]])
            continue
        pvnbrs = pivot(new_X, Pcnt, adj, new_P)
        pvdeg = len(pvnbrs)
        if pvdeg == Pcnt:
            continue
        pvPnbrs = pivot(new_P, Pcnt - 1, adj)
        if len(pvPnbrs) > pvdeg:
            pvnbrs = pvPnbrs
        stack.append((P, X, Pv))
        R.append(n)
        P = new_P
        X = new_X
        Pv = P - pvnbrs
    return result


cdef pivot(vs, int goal, adj, filter_=None):
    cdef int max_num = -1
    cdef int d
    cdef int v
    for v in vs:
        if filter is not None:
            a = adj[v]
        else:
            a = adj[v] & filter_
        d = len(a)
        if d > max_num:
            res = a
            max_num = d
            if d == goal:  # short cut
                break
    if max_num == -1:
        return set()
    return res
