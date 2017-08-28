#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import chorus.util.geometry as gm

# Multiple bond interval / length
F_MULI = 0.15
# 1 - (inner / outer) double bond length
F_INNL = 0.2
# Factor to avoid overlap of bond segment with atom symbol
F_AOVL = 0.3


def draw(canvas, mol):
    """Draw molecule structure image.

    Args:
        canvas: draw.drawable.Drawable
        mol: model.graphmol.Compound
    """
    mol.require("ScaleAndCenter")
    mlb = mol.size2d[2]
    if not mol.atom_count():
        return
    bond_type_fn = {
        1: {
            0: single_bond,
            1: wedged_single,
            2: dashed_wedged_single,
            3: wave_single,
        }, 2: {
            0: cw_double,
            1: counter_cw_double,
            2: double_bond,
            3: cross_double
        }, 3: {
            0: triple_bond
        }
    }
    # Draw bonds
    for u, v, bond in mol.bonds_iter():
        if not bond.visible:
            continue
        if (u < v) == bond.is_lower_first:
            f, s = (u, v)
        else:
            s, f = (u, v)
        p1 = mol.atom(f).coords
        p2 = mol.atom(s).coords
        if p1 == p2:
            continue  # avoid zero division
        if mol.atom(f).visible:
            p1 = gm.t_seg(p1, p2, F_AOVL, 2)[0]
        if mol.atom(s).visible:
            p2 = gm.t_seg(p1, p2, F_AOVL, 1)[1]
        color1 = mol.atom(f).color
        color2 = mol.atom(s).color
        bond_type_fn[bond.order][bond.type](
            canvas, p1, p2, color1, color2, mlb)

    # Draw atoms
    for n, atom in mol.atoms_iter():
        if not atom.visible:
            continue
        p = atom.coords
        color = atom.color
        # Determine text direction
        if atom.H_count:
            cosnbrs = []
            hrzn = (p[0] + 1, p[1])
            for nbr in mol.graph.neighbors(n):
                pnbr = mol.atom(nbr).coords
                try:
                    cosnbrs.append(gm.dot_product(hrzn, pnbr, p) /
                                   gm.distance(p, pnbr))
                except ZeroDivisionError:
                    pass
            if not cosnbrs or min(cosnbrs) > 0:
                # [atom]< or isolated node(ex. H2O, HCl)
                text = atom.formula_html(True)
                canvas.draw_text(p, text, color, "right")
                continue
            elif max(cosnbrs) < 0:
                # >[atom]
                text = atom.formula_html()
                canvas.draw_text(p, text, color, "left")
                continue
        # -[atom]- or no hydrogens
        text = atom.formula_html()
        canvas.draw_text(p, text, color, "center")


def single_bond(canvas, p1, p2, color1, color2, mlb):
    canvas.draw_line(p1, p2, color1, color2)


def wedged_single(canvas, p1, p2, color1, color2, mlb):
    canvas.draw_wedge(p1, p2, color1)


def dashed_wedged_single(canvas, p1, p2, color1, color2, mlb):
    canvas.draw_dashed_wedge(p1, p2, color1)


def wave_single(canvas, p1, p2, color1, color2, mlb):
    canvas.draw_wave_line(p1, p2, color1)


def double_bond(canvas, p1, p2, color1, color2, mlb):
    ip1, ip2 = gm.p_seg(p1, p2, True, F_MULI / 2 * mlb)
    op1, op2 = gm.p_seg(p1, p2, False, F_MULI / 2 * mlb)
    canvas.draw_line(ip1, ip2, color1, color2)
    canvas.draw_line(op1, op2, color1, color2)


def cw_double(canvas, p1, p2, color1, color2, mlb):
    ip1, ip2 = gm.p_seg(p1, p2, True, F_MULI * mlb, F_INNL)
    canvas.draw_line(p1, p2, color1, color2)
    canvas.draw_line(ip1, ip2, color1, color2)


def counter_cw_double(canvas, p1, p2, color1, color2, mlb):
    ip1, ip2 = gm.p_seg(p1, p2, False, F_MULI * mlb, F_INNL)
    canvas.draw_line(p1, p2, color1, color2)
    canvas.draw_line(ip1, ip2, color1, color2)


def cross_double(canvas, p1, p2, color1, color2, mlb):
    ip1, ip2 = gm.p_seg(p1, p2, True, F_MULI / 2 * mlb)
    op1, op2 = gm.p_seg(p1, p2, False, F_MULI / 2 * mlb)
    canvas.draw_line(ip1, op2, color1, color2)
    canvas.draw_line(op1, ip2, color1, color2)


def triple_bond(canvas, p1, p2, color1, color2, mlb):
    ip1, ip2 = gm.p_seg(p1, p2, True, F_MULI * mlb)
    op1, op2 = gm.p_seg(p1, p2, False, F_MULI * mlb)
    canvas.draw_line(p1, p2, color1, color2)
    canvas.draw_line(ip1, ip2, color1, color2)
    canvas.draw_line(op1, op2, color1, color2)
