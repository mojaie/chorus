# cython: boundscheck=False, wraparound=False, langage_level=3, profile=True

#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import cython
import networkx as nx

from libcpp.vector cimport vector
from libc.time cimport clock_t, clock, CLOCKS_PER_SEC


cdef struct Edge:
    int u1
    int v1
    int u2
    int v2


@cython.cdivision(True)
def comparison_graph(arr1, arr2, double timeout):
    cdef clock_t t0 = clock()
    cdef clock_t expire = t0 + <clock_t>(timeout * CLOCKS_PER_SEC)
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
    product = nx.Graph()
    cdef int interrupted = 0
    for i in range(i_max):
        for j in range(j_max):
            if c1[i] == c2[j]:
                e.u1 = u1[i]
                e.v1 = v1[i]
                e.u2 = u2[j]
                e.v2 = v2[j]
                edges.push_back(e)
                if clock() >= expire:
                    interrupted = 1
                    break
    if not interrupted:
        # Graph.add_edges is expensive. Add adjacency dict manually.
        node = {}
        for edge in edges:
            node[(edge.u1, edge.u2)] = {}
            node[(edge.v1, edge.v2)] = {}
        adj = node.copy()
        for edge in edges:
            adj[(edge.u1, edge.u2)][(edge.v1, edge.v2)] = {}
            adj[(edge.v1, edge.v2)][(edge.u1, edge.u2)] = {}
        product._node = node
        product._adj = adj
    cdef double elapsed = <double>(clock() - t0) / CLOCKS_PER_SEC
    return product, elapsed

"""
def find_cliques_old(graph, timeout=None):
    cdef float t0
    if timeout is not None:
        t0 = time()
    # tuple index may expensive. convert to integer
    cdef int i
    decode = {i: node for i, node in enumerate(graph)}
    encode = {node: i for i, node in decode.items()}
    adj = {}
    for node, adjs in graph.adjacency():
        eadjs = set(encode[a] for a in adjs)
        adj[encode[node]] = eadjs
    # initialize
    result = []
    R = []
    P = set(adj.keys())
    X = set()
    stack = []
    Pv = P - pivot(P, len(P) - 1, adj)
    # DFS loop
    cdef int n, Pcnt, pvdeg
    cdef float now
    while True:
        if timeout is not None:
            now = time.time()
            if now - t0 > timeout:
                return result, False
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
    return result, True


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
"""

@cython.cdivision(True)
def find_cliques(G, double timeout):
    # Based on networkX v2.1
    # TODO: slower???
    result = []
    if len(G) == 0:
        return result, 0
    cdef clock_t t0 = clock()
    cdef clock_t expire = t0 + <clock_t>(timeout * CLOCKS_PER_SEC)
    cdef int i, u, q, l
    # tuple index may expensive. convert to integer
    decode = {i: node for i, node in enumerate(G)}
    encode = {node: i for i, node in decode.items()}
    adj = {}
    for node, adjs in G.adjacency():
        eadjs = set()
        for a in adjs:
            eadjs.add(encode[a])
        adj[encode[node]] = eadjs
    Q = [None]
    subg = set(decode)
    cand = set(decode)
    u = max_adj(subg, cand, adj)
    ext_u = cand - adj[u]
    stack = []
    cdef double elapsed
    while True:
        if clock() >= expire:
            break
        if not ext_u:
            Q.pop()
            if not stack:
                break
            subg, cand, ext_u = stack.pop()
            continue
        q = ext_u.pop()
        cand.remove(q)
        Q[len(Q) - 1] = q
        adj_q = adj[q]
        subg_q = subg & adj_q
        if not subg_q:
            result.append([decode[r] for r in Q])
            continue
        cand_q = cand & adj_q
        if cand_q:
            stack.append((subg, cand, ext_u))
            Q.append(None)
            subg = subg_q
            cand = cand_q
            u = max_adj(subg, cand, adj)
            ext_u = cand - adj[u]
    elapsed = <double>(clock() - t0) / CLOCKS_PER_SEC
    return result, elapsed


cdef max_adj(subg, cand, adj):
    cdef int maxadj = -1
    cdef int s, numadj, res
    for s in subg:
        numadj = len(cand & adj[s])
        if numadj > maxadj:
            maxadj = numadj
            res = s
    return res
