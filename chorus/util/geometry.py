#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

from math import cos, acos, sin, hypot, pi
from itertools import cycle

from chorus.util import iterator


def vector(p1, p2):
    """ p1 to p2 vector
    Args:
      p1, p2: point (x, y)
    """
    return p2[0] - p1[0], p2[1] - p1[1]


def distance(p1, p2):
    """ distance between p1 and p2 """
    return hypot(*vector(p1, p2))


def translate(v, op):
    """ translate vector
    Args:
      v: vector (x, y)
      op: operator (x, y)
    """
    return v[0] + op[0], v[1] + op[1]


def scale(p, factor, o=(0, 0)):
    """ scale vector
    Args:
      p: point (x, y)
      factor: scaling factor
      o: origin (x, y)
    """
    v = vector(o, p)
    sv = v[0] * factor, v[1] * factor
    return translate(sv, o)


def unit(v, lg=1):
    """ unit vector
    Args:
        v: vector (x, y)
        lg: length
    Raises:
        ValueError: Null vector was given
    """
    try:
        res = scale(v, lg / distance((0, 0), v))
    except ZeroDivisionError:
        raise ValueError("Null vector was given")
    return res


def rotate(p, rad, o=(0, 0)):
    """ rotate vector
    Args:
      p: point (x, y)
      rad: angle(radian)
      o: origin (x, y)
    """
    v = vector(o, p)
    fx = lambda x, y, d: x * cos(d) - y * sin(d)
    fy = lambda x, y, d: x * sin(d) + y * cos(d)
    rv = fx(v[0], v[1], rad), fy(v[0], v[1], rad)
    return translate(rv, o)


def cross_product(p1, p2, o=(0, 0)):
    """ Returns cross product
    Args:
      p1, p2: point (x, y)
      o: origin
    """
    v1 = vector(o, p1)
    v2 = vector(o, p2)
    return v1[0] * v2[1] - v1[1] * v2[0]


def dot_product(p1, p2, o=(0, 0)):
    """ Returns dot product
    Args:
      p1, p2: point (x, y)
      o: origin
    """
    v1 = vector(o, p1)
    v2 = vector(o, p2)
    return v1[0] * v2[0] + v1[1] * v2[1]


def interior_angle(p1, p2, o=(0, 0)):
    """ Returns interior angle of two vector(0 <= θ <= pi)
    Args:
      p1, p2: point (x, y)
      o: origin
    Raises:
      ValueError: p1 or p2 is overlapped with origin
    """
    v1 = vector(o, p1)
    v2 = vector(o, p2)
    len1 = distance(o, p1)
    len2 = distance(o, p2)
    try:
        return acos(dot_product(v1, v2) / (len1 * len2))
    except ZeroDivisionError:
        raise ValueError("p1 or p2 is overlapped with origin")


def m_seg(p1, p2, rad, dist):
    """ move segment by distance
    Args:
      p1, p2: point(x, y)
      rad: relative direction angle(radian)
      dist: distance
    Return:
      translated segment(p1, p2)
    """
    v = vector(p1, p2)
    m = unit(rotate(v, rad), dist)
    return translate(p1, m), translate(p2, m)


def t_seg(p1, p2, t, align=0):
    """ trim segment
    Args:
      p1, p2: point(x, y)
      t: scaling factor (1 - trimed segment / original segment)
      align: 1: trim p2, 2: trim p1, 0: both side
    Return:
      trimmed segment(p1, p2)
    """
    v = vector(p1, p2)
    result = {
        1: lambda a, b: (a, translate(b, scale(v, -t))),
        2: lambda a, b: (translate(a, scale(v, t)), b),
        0: lambda a, b: (translate(a, scale(v, t / 2)),
                         translate(b, scale(v, -t / 2)))
    }
    return result[align](p1, p2)


def p_seg(p1, p2, cw, interval, trim=0, align=0):
    """ parallel segment
    Args:
      p1, p2: point(x, y)
      cw: m_seg rad True: -π/2, False: π/2
      interval: m_seg dist
      trim: t_seg trim
      align: t_seg align
    """
    case = {True: pi / -2, False: pi / 2}
    p1m, p2m = m_seg(p1, p2, case[cw], interval)
    return t_seg(p1m, p2m, trim, align)


def is_clockwise(vertices):
    """ Evaluate whether vertices are in clockwise order.
    Args:
      vertices: list of vertices (x, y) in polygon.
    Returns:
      True: clockwise, False: counter-clockwise
    Raises:
      ValueError: the polygon is complex or overlapped.
    """
    it = iterator.consecutive(cycle(vertices), 3)
    clockwise = 0
    counter = 0
    for _ in range(len(vertices)):
        p0, p1, p2 = next(it)
        cross = cross_product(p1, p2, p0)
        int_angle = interior_angle(p0, p2, p1)  # raises ValueError
        if cross < 0:
            clockwise += int_angle
            counter += 2 * pi - int_angle
        else:
            clockwise += 2 * pi - int_angle
            counter += int_angle
    if round(clockwise / pi) == len(vertices) - 2:
        return True
    elif round(counter / pi) == len(vertices) - 2:
        return False
    else:
        raise ValueError("the polygon is complex or overlapped")


def rad(radian):
    """ A radian in the range of 0 - 2pai """
    return (radian + 2 * pi) % (2 * pi)
