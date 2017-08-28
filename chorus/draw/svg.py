#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import base64
import functools
import math

from chorus.draw import drawable, drawer2d, helper
from chorus.util.geometry import scale, distance, m_seg, p_seg
from chorus.util.text import substitute

FONT_SIZE = 14   # Font size(pixel)


class SVG(drawable.Drawable):
    """Draw chemical structure in SVG format.

    Args:
        molecule: compound
        screen_size (width, height): pixel size to display

    Attributes:
        font_weight (str): font weight
        font_family (str): font family
        screen_size (width, height): pixel size to display
        logical_size (width, height): pixel size of actual viewport
    """
    def __init__(self, molecule, screen_size=(180, 180)):
        mol = helper.ready_to_draw(molecule)
        self.font_weight = "normal"
        self.font_family = "Helvetica"
        self.bgcolor = None  # background fill (if None, transparent)
        self.displaymlb = 30  # suitable for FONT_SIZE=14 and default linewidth
        self.margin = self.displaymlb
        w, h, mlb = mol.size2d
        self.scale_factor = self.displaymlb / mlb
        self.original_size = (w, h)  # for coords conversion
        self.screen_size = screen_size  # svg object display size
        self._header = ['<svg xmlns="http://www.w3.org/2000/svg"']
        self._header.append(' xmlns:xlink="http://www.w3.org/1999/xlink"')
        self._header.append(' version="1.2" baseProfile="tiny"')
        self._header.append(' text-rendering="geometricPrecision"')
        self._elems = []
        drawer2d.draw(self, mol)

    def contents(self):
        """Get svg string
        """
        c = self._header[:]
        c.append(' font-weight="{}"'.format(self.font_weight))
        c.append(' font-family="{}"'.format(self.font_family))
        c.append(' width="{}" height="{}"'.format(*self.screen_size))
        sclw = self.original_size[0] * self.scale_factor
        sclh = self.original_size[1] * self.scale_factor
        longside = max([sclw, sclh])
        width = round(longside + self.margin * 2, 2)
        height = round(longside + self.margin * 2, 2)
        xleft = round(-self.margin - (longside - sclw) / 2, 2)
        ytop = round(-self.margin - (longside - sclh) / 2, 2)
        c.append(' viewBox="{} {} {} {}">\n'.format(
            xleft, ytop, width, height))
        if self.bgcolor is not None:
            c.append('<rect x="{}", y="{}" width="{}" height="{}" fill="{}" \
                />\n'.format(xleft, ytop, width, height, self.bgcolor))
        c.extend(self._elems)
        c.append("</svg>")
        return "".join(c)

    def data_url_scheme(self):
        """Get svg in Data URL Scheme format.
        """
        # TODO: move to web.app or make it function
        # remove #svg from dataframe
        encoded = base64.b64encode(self.contents().encode())
        return "data:image/svg+xml;base64," + encoded.decode()

    def save(self, path):
        """Save svg as file(.svg)

        Args:
            path (str): destination to save file
        """
        with open(path, 'w') as f:
            f.write(self.contents())

    def _coords_conv(self, pos):
        """For Svg coordinate system, reflect over X axis and
        translate from center to top-left
        """
        px = (self.original_size[0] / 2 + pos[0]) * self.scale_factor
        py = (self.original_size[1] / 2 - pos[1]) * self.scale_factor
        return round(px, 2), round(py, 2)

    def line_conv(method):
        @functools.wraps(method)
        def wrapper(self, p1, p2, c1, c2=None):
            p1_ = self._coords_conv(p1)
            p2_ = self._coords_conv(p2)
            # rgb(255,255,255)
            c1_ = 'rgb({},{},{})'.format(*c1)
            if c2:
                c2_ = 'rgb({},{},{})'.format(*c2)
                method(self, p1_, p2_, c1_, c2_)
            else:
                method(self, p1_, p2_, c1_)
        return wrapper

    def text_conv(method):
        @functools.wraps(method)
        def wrapper(self, p, text, c, align):
            p_ = self._coords_conv(p)
            c_ = 'rgb({},{},{})'.format(*c)
            method(self, p_, text, c_, align)
        return wrapper

    def _draw_line_base(self, p1, p2, c1, c2=None, op=""):
        tmpl = '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}"{} />\n'
        if c2 and c2 != c1:
            pmx, pmy = (
                round((p1[0] + p2[0]) / 2, 2),
                round((p1[1] + p2[1]) / 2, 2)
            )
            self._elems.append(tmpl.format(p1[0], p1[1], pmx, pmy, c1, op))
            self._elems.append(tmpl.format(pmx, pmy, p2[0], p2[1], c2, op))
        else:
            self._elems.append(tmpl.format(p1[0], p1[1], p2[0], p2[1], c1, op))

    @line_conv
    def draw_line(self, p1, p2, c1, c2=None):
        self._draw_line_base(p1, p2, c1, c2)

    @line_conv
    def draw_dashed_line(self, p1, p2, c1, c2=None):
        op = ' stroke-dasharray="10,10"'
        self._draw_line_base(p1, p2, c1, c2, op)

    @line_conv
    def draw_wedge(self, head, tail, color):
        width = self.displaymlb * 0.2
        tmpl = '<polygon points="{}" fill="{}" />\n'
        t1 = m_seg(tail, head, math.pi / -2, width / 2)[0]
        t2 = m_seg(tail, head, math.pi / 2, width / 2)[0]
        ptext = ' '.join(['{},{}'.format(
            round(x, 2), round(y, 2)) for x, y in (head, t1, t2)])
        self._elems.append(tmpl.format(ptext, color))

    @line_conv
    def draw_dashed_wedge(self, head, tail, color):
        width = self.displaymlb * 0.2
        tmpl = '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" />\n'
        t1 = m_seg(tail, head, math.pi / -2, width / 2)[0]
        t2 = m_seg(tail, head, math.pi / 2, width / 2)[0]
        step_num = 7
        for i in range(step_num):
            step = distance(head, tail) / step_num * i
            trim = 1 / step_num * i
            tp1, tp2 = p_seg(t1, t2, True, step, trim)
            self._elems.append(
                tmpl.format(round(tp1[0], 2), round(tp1[1], 2),
                            round(tp2[0], 2), round(tp2[1], 2), color))

    @line_conv
    def draw_wave_line(self, p1, p2, color):
        width = self.displaymlb * 0.2
        tmpl = '<polyline points="{}" stroke="{}" fill="none" />\n'
        s1 = m_seg(p1, p2, math.pi / -2, width / 2)
        s2 = m_seg(p1, p2, math.pi / 2, width / 2)
        ps = [p1]
        for i, s in enumerate((s1, s2, s1, s2, s1, s2)):
            p = scale(s[1], (i + 1) / 7, s[0])
            ps.append(p)
        ps.append(p2)
        ptext = ' '.join(['{},{}'.format(
            round(x, 2), round(y, 2)) for x, y in ps])
        self._elems.append(tmpl.format(ptext, color))

    @text_conv
    def draw_text(self, pos, text, color, align="center"):
        tmpl = '<text x="{}" y="{}" font-size="{}" fill="{}"{}>{}</text>\n'
        stmpl = '<tspan baseline-shift="{}" font-size="{}">'
        rt = substitute({
            '<sub>': stmpl.format("-25%", FONT_SIZE * 0.7),
            '<sup>': stmpl.format("50%", FONT_SIZE * 0.7),
            '</sub>': '</tspan>',
            '</sup>': '</tspan>'
        }, text)
        op = {"center": ' text-anchor="middle"',
              "right": ' text-anchor="end"', "left": ''}
        xoffset = {"center": 0, "right": FONT_SIZE / 2,
                   "left": -FONT_SIZE / 2}
        px = pos[0] + xoffset[align]
        py = pos[1] + FONT_SIZE / 2
        self._elems.append(
            tmpl.format(round(px, 2), round(py, 2),
                        FONT_SIZE, color, op[align], rt))
