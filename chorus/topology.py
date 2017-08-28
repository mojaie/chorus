#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

from collections import deque

import networkx as nx


def recognize_nx(mol):
    """ NetworkX implementation (for comparison)"""
    mol.rings = nx.cycle_basis(mol.graph)
    mol.isolated = sorted(nx.connected_components(mol.graph),
                          key=len, reverse=True)[1:]
    mol.scaffolds = [x for x in nx.biconnected_components(mol.graph)
                     if len(x) > 2]


def recognize(mol):
    """ Detect cycle basis, biconnected and isolated components (DFS).
    This will add following attribute to the molecule instance object.

    mol.ring: Cycle basis
    mol.scaffold: biconnected components
    mol.isolated: isolated components other than the largest one

    To find minimum set of rings, additionally execute topology.minify_ring.

    Reference:
        networkx cycle_basis function
    """
    g = set(i for i, _ in mol.atoms_iter())
    bccs = {}  # BiConnected Components
    isoc = []  # ISOlated Components
    while g:
        start = g.pop()
        stack = [start]
        pred = {start: start}
        used = {start: set()}
        root = {start: start}
        while stack:
            tail = stack.pop()
            for nbr in mol.neighbors(tail):
                if nbr not in used:  # New node
                    pred[nbr] = tail
                    stack.append(nbr)
                    used[nbr] = {tail}
                    root[nbr] = nbr
                elif nbr in stack:  # Cycle found
                    pn = used[nbr]
                    cyc = [nbr, tail]
                    p = pred[tail]
                    end = pred[nbr]
                    root[nbr] = root[tail] = root[end]
                    while p not in pn:  # Backtrack
                        cyc.append(p)
                        root[p] = root[end]
                        if p in bccs:  # Append scaffold to new cycle
                            if root[end] not in bccs:
                                bccs[root[end]] = []
                            bccs[root[end]].extend(bccs[p])
                            del bccs[p]
                        p = pred[p]
                    cyc.append(p)
                    if root[end] not in bccs:  # Append new cycle to scaffold
                        bccs[root[end]] = []
                    bccs[root[end]].append(cyc)
                    used[nbr].add(tail)
        isoc.append(list(pred.keys()))
        # print(pred)
        g -= set(pred)
    mol.rings = []
    mol.scaffolds = []
    for cycles in bccs.values():
        rcnt = len(mol.rings)
        mol.rings.extend(cycles)
        mol.scaffolds.append(list(range(rcnt, rcnt + len(cycles))))
    mol.isolated = list(sorted(isoc, key=len, reverse=True))[1:]
    mol.descriptors.add("Topology")


def minify_ring(mol, verbose=False):
    """ Minify ring set (similar to SSSR)
    Limitation: this can not correctly recognize minimum rings
    in the case of non-outerplanar graph.
    Note: concept of SSSR is controversial. Roughly reduce the size of
    cycle basis can help some scaffold-based analysis
    """
    mol.require("Topology")
    for cyc_idx in mol.scaffolds:
        rings = deque(sorted([mol.rings[c] for c in cyc_idx], key=len))
        minified = []
        cnt = 0
        while rings:
            # TODO: can be infinite loop ?
            cnt += 1
            if cnt > 100:
                print("Minimization failed")
                print(mol.options)
                mol.descriptors.add("MinifiedRing")
                return
            r = rings.popleft()
            init_r = r
            if verbose:
                print(len(r), "Ring:{}".format(r))
            for m in minified:
                if verbose:
                    print(len(m), "Minified:{}".format(m))
                resolved = resolve_inclusion(r, m)
                if resolved:
                    if verbose:
                        print(len(resolved[0]), len(resolved[1]),
                              "Resolved:{}".format(resolved))
                    r = resolved[0]
            if verbose:
                print(len(r), "New ring:{}\n".format(r))
            if len(r) == len(init_r):  # no longer be able to minified
                minified.append(r)
            else:  # further minification required
                rings.append(r)
        for c in cyc_idx:
            mol.rings[c] = minified.pop()
    mol.descriptors.add("MinifiedRing")


def resolve_inclusion(a, b):
    lt, bg = (a, b)
    if len(lt) > len(bg):
        lt, bg = (bg, lt)
        rev = True
    else:
        rev = False
    isec = set(lt) & set(bg)
    isize = len(isec)
    if isize == 3 and len(lt) == 4:
        pass  # TODO: cubane special case
    elif isize != len(lt) and isize <= len(lt) / 2 + 1:
        return False  # Already minimum
    lq = deque(lt)
    bq = deque(bg)
    while bq[0] not in isec or bq[-1] in isec:
        bq.rotate(1)
    while lq[0] != bq[0]:
        lq.rotate(1)
    if bq[1] != lq[1]:  # Reverse
        lq = deque(reversed(lq))
        while lq[0] != bq[0]:
            lq.rotate(1)
    if list(bq)[:isize] != list(lq)[:isize]:
        return False  # Aboid irregular minification
    bx = list(bq)[isize - 1:]
    lx = list(lq)[isize:]
    new_bg = bx + [bq[0]] + list(reversed(lx))
    if rev:
        lt, new_bg = (new_bg, lt)
    return lt, new_bg
