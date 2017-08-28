#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

from collections import deque
import math

from chorus import topology
import chorus.util.geometry as gm

# BOND_LENGTH = 1


def calc2dcoords(mol):
    """ Calculate optimal 2D coordinates of chemical structure
    """
    topology.recognize(mol)
    g = set(i for i, _ in mol.atoms_iter())
    # 1: get nodes in scaffolds
    scaffolds = []
    belongs = {}
    for i, rkeys in enumerate(sorted(mol.scaffolds, key=len)):
        scf = []
        for rkey in rkeys:
            ring = mol.rings[rkey]
            for r in ring:
                belongs[r] = i
            scf.append(ring)
            g -= set(ring)
        scaffolds.append(scf)
    # 2: traverse nodes and scaffolds
    # the node and scaffold graph should be a tree (no cycles)
    f = True
    # print(scaffolds)
    coords = {}
    while g:
        if f and scaffolds:  # largest scaffold is first
            stack = [scaffolds[-1][0][0]]
            f = False
        else:
            stack = [g.pop()]
        pred = {}
        branch = {}
        while stack:
            # print("stack: {}".format(stack))
            tail = stack.pop()
            # print("tail: {}".format(tail))
            if tail in belongs:  # scaffolds
                scf = scaffold_coords(scaffolds[belongs[tail]])
                # print(scf.keys())
                # rotate and translate
                if not coords:
                    coords = scf
                else:
                    u = coords[pred[tail]]
                    v = scf[tail]
                    op = [u[0] + math.cos(u[2]), u[1] + math.sin(u[2])]
                    translate(scf, gm.vector(v[:2], op))
                    rotate(scf, op, gm.rad(u[2] + math.pi - v[2]))
                    coords.update(scf)
                # stack nbrs of scaffold
                for k in scf.keys():
                    pred[k] = None
                    for nbr in mol.neighbors(k):
                        if nbr not in scf.keys():
                            stack.append(nbr)
                            pred[nbr] = k
            else:  # append linker
                if tail not in pred:  # isolated
                    coords[tail] = [0, 0, 0, 1]
                    continue
                p = pred[tail]
                x, y, ang, d = coords[p]
                # TODO: ring configuration
                coords[tail] = [x + math.cos(ang), y + math.sin(ang),
                                ang + d * math.pi / 3, d * -1]
                if p not in branch:
                    coords[p][2] = gm.rad(coords[p][2] + math.pi * 2 / 3 * d)
                    branch[p] = 1
                elif branch[p] == 1:
                    coords[p][2] = gm.rad(coords[p][2] + math.pi * d)
                    branch[p] += 1
                elif branch[p] == 2:
                    coords[p][2] = gm.rad(coords[p][2] + math.pi * 2 / 3 * d)
                    branch[p] += 1
                for nbr in mol.neighbors(tail):
                    if nbr not in pred:
                        stack.append(nbr)
                        pred[nbr] = tail
        g -= set(pred)
    resolve_overlap(coords)
    for i, a in mol.atoms_iter():
        mol.atom(i).coords = coords[i][:2]


def translate(vertices, operator):
    for k, v in vertices.items():
        new_pos = list(gm.translate(v[:2], operator))
        vertices[k] = new_pos + v[2:]


def rotate(vertices, center, angle):
    for k, v in vertices.items():
        new_pos = list(gm.rotate(v[:2], angle, center))
        new_angle = v[2] + angle
        vertices[k] = [new_pos[0], new_pos[1], new_angle, v[3]]


def scaffold_coords(rings):
    """ assign scaffold coordinate and angle
    {node: (coords, angle)}
    """
    rs = deque(sorted(rings, key=len, reverse=True))
    base = rs.popleft()
    first = base[0]
    coords = {first: [0, 0, 0, 1]}
    attach_spiro(coords, base, first)
    # display(coords)
    while rs:
        r = rs.popleft()
        isec = list(set(coords) & set(r))
        if not isec:
            rs.append(r)
        elif len(isec) == 1:  # spiro
            attach_spiro(coords, r, isec)
        elif len(isec) == 2:  # fuse
            attach_fused(coords, r, isec)
        else:  # bridge or append
            bridge(coords, r, isec)
        # display(coords)
    return coords


def bridge(coords, ring, isec):
    r = deque(ring)
    while r[0] not in isec or r[-1] in isec:
        r.rotate()
    v = r.popleft()
    u = None
    while True:
        p = r.popleft()
        if p in isec:
            u = p
        else:
            break
    r.appendleft(p)  # TODO: refactor
    ux, uy, ua, _ = coords[u]
    vx, vy, va, _ = coords[v]
    utov = gm.vector((ux, uy), (vx, vy))
    utovang = gm.rad(math.atan2(utov[1], utov[0]))
    dist = gm.distance((ux, uy), (vx, vy))
    step, angle = bridge_angle(dist, len(r))
    angle += utovang
    prev = u
    coords[prev][2] = gm.rad(utovang + math.pi)  # reverse
    while r:
        tail = r.popleft()
        x, y, a, _ = coords[prev]
        coords[tail] = [x + math.cos(angle), y + math.sin(angle),
                        gm.rad(a + step), 1]
        angle += step
        prev = tail


def bridge_angle(base, n):
    # TODO: There may be an elegant analytical solution
    length = 1
    res = 0
    dmin = 100
    for i in range(330):
        r = math.pi / 180 * i
        a = (math.pi - r / 2) / (n + 1)
        d = abs(length * math.sin(r / 2) - base * math.sin(a))
        if d < dmin:
            res = r
            dmin = d
        if d < 0.0001:
            break
    step = -(2 * math.pi - res) / (n + 1)
    angle = -step * n / 2
    return step, angle


def attach_spiro(coords, ring, isec):
    r = deque(ring)
    while r[0] != isec:
        r.rotate()
    prev = r.popleft()
    step = math.pi * 2 / len(ring)  # angle shift
    angle = gm.rad(coords[prev][2] + (math.pi - step) / 2)  # init angle
    coords[prev][2] += math.pi  # reverse
    while r:
        tail = r.popleft()
        x, y, a, _ = coords[prev]
        coords[tail] = [x + math.cos(angle), y + math.sin(angle),
                        gm.rad(a - step), 1]
        angle -= step
        prev = tail


def attach_fused(coords, ring, isec):
    r = deque(ring)
    while r[0] not in isec or r[-1] in isec:
        r.rotate()
    first = r.popleft()
    second = r.popleft()
    x, y = gm.vector(coords[first][:2], coords[second][:2])
    ftos = gm.rad(math.atan2(y, x))
    if ftos >= coords[second][2] or coords[second][2] > ftos + math.pi:
        # clockwise
        step = math.pi * -2 / len(ring)
        coords[second][2] = ftos + (math.pi + step) / 2
    else:  # counter-clockwise
        step = math.pi * 2 / len(ring)
        coords[second][2] = ftos - (math.pi - step) / 2
    angle = ftos + step
    prev = second
    while r:
        tail = r.popleft()
        x, y, a, _ = coords[prev]
        coords[tail] = [x + math.cos(angle), y + math.sin(angle),
                        gm.rad(a + step), 1]
        angle += step
        prev = tail
    coords[second][2] = ftos
    coords[first][2] = gm.rad(ftos + math.pi)


def display(coords):
    print("")
    for k, v in coords.items():
        print("{}: ({}, {}) {} {}".format(k,
                                          round(v[0], 2),
                                          round(v[1], 2),
                                          round(v[2], 2), v[3]))


def resolve_overlap(coords):
    while True:
        break
