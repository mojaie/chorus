# cython: boundscheck=False, wraparound=False, langage_level=3, profile=True

#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import cython

from libcpp.vector cimport vector
from libc.time cimport clock_t, clock, CLOCKS_PER_SEC


@cython.cdivision(True)
def comparison_graph(arr1, arr2, double timeout):
    cdef clock_t t0 = clock()
    cdef clock_t expire = t0 + <clock_t>(timeout * CLOCKS_PER_SEC)
    result = {
        "edges": [],
        "decoder": {},
        "elapsed_time": 0,
        "timeout": False
    }
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
    cdef tuple node1, node2
    encoder = {}
    for i in range(i_max):
        for j in range(j_max):
            if c1[i] == c2[j]:
                node1 = (u1[i], u2[j])
                if node1 not in encoder:
                    encoder[node1] = len(encoder)
                idx1 = encoder[node1]
                node2 = (v1[i], v2[j])
                if node2 not in encoder:
                    encoder[node2] = len(encoder)
                idx2 = encoder[node2]
                result["edges"].append((idx1, idx2))
                if clock() >= expire:
                    result["timeout"] = True
                    break
        else:
            continue
        break
    result["decoder"] = {v: k for k, v in encoder.items()}
    result["elapsed_time"] = <double>(clock() - t0) / CLOCKS_PER_SEC
    return result


@cython.cdivision(True)
def find_cliques(nodes, edges, double timeout):
    # Based on networkX v2.1
    result = {
        "max_clique": [],
        "elapsed_time": 0,
        "timeout": False
    }
    if len(nodes) == 0:
        return result
    cdef clock_t t0 = clock()
    cdef clock_t expire = t0 + <clock_t>(timeout * CLOCKS_PER_SEC)
    cdef int i, u, q, l
    adj = {n: set() for n in nodes}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    Q = [None]
    subg = set(nodes)
    cand = set(nodes)
    u = max_adj(subg, cand, adj)
    ext_u = cand - adj[u]
    stack = []
    cdef double elapsed
    cdef int itcount = 0
    cdef int max_q = 0
    while True:
        if not ext_u:
            Q.pop()
            if not stack:
                break
            if clock() >= expire:
                result["timeout"] = True
                break
            subg, cand, ext_u = stack.pop()
            continue
        q = ext_u.pop()
        cand.remove(q)
        Q[len(Q) - 1] = q
        adj_q = adj[q]
        subg_q = subg & adj_q
        if not subg_q:
            if len(Q) > max_q:
                max_q = len(Q)
                result["max_clique"] = Q[:]
            continue
        cand_q = cand & adj_q
        if cand_q:
            stack.append((subg, cand, ext_u))
            Q.append(None)
            subg = subg_q
            cand = cand_q
            u = max_adj(subg, cand, adj)
            ext_u = cand - adj[u]
    result["elapsed_time"] = <double>(clock() - t0) / CLOCKS_PER_SEC
    return result


cdef max_adj(subg, cand, adj):
    cdef int maxadj = -1
    cdef int s, numadj, res
    for s in subg:
        numadj = len(cand & adj[s])
        if numadj > maxadj:
            maxadj = numadj
            res = s
    return res
