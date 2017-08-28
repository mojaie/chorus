#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import functools
import io
import math

import matplotlib as mpl
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.font_manager import FontProperties

from chorus.draw import drawable, drawer2d, helper
from chorus.util.geometry import scale, distance, m_seg, p_seg
from chorus.util.text import substitute

mpl.rcParams["mathtext.default"] = "regular"
mpl.rcParams["font.sans-serif"] = [
    "Arial", "Helvetica", "Liberation Sans", "Bitstream Vera Sans"]


def basefont(size):
    f = FontProperties()
    f.set_size(size)
    return f


class Matplotlib(drawable.Drawable):
    """Draw chemical structure with Matplotlib.

    Args:
        molecule: compound

    Attributes:
        font_weight (str): font weight
        font_family (str): font family
        screen_size (width, height): pixel size to display
        logical_size (width, height): pixel size of actual viewport
    """
    def __init__(self, molecule):
        mol = helper.ready_to_draw(molecule)
        self.dpi = 120
        self.fig = mpl.figure.Figure()
        FigureCanvas(self.fig)
        self.fig.subplots_adjust(left=0, right=1, top=1, bottom=0,
                                 wspace=0, hspace=0)
        self.ax = self.fig.add_subplot(111)
        w, h, mlb = mol.size2d
        self.displaymlb = 60
        self.margin = self.displaymlb
        self.scale_factor = self.displaymlb / mlb
        scaledw = w * self.scale_factor + self.margin * 2
        scaledh = h * self.scale_factor + self.margin * 2
        longside = max([scaledw, scaledh])
        self.font_size = 16
        self.fig.set_size_inches(longside / self.dpi, longside / self.dpi)
        self.ax.set_xlim(-scaledw / 2, scaledw / 2)
        self.ax.set_ylim(-scaledh / 2, scaledh / 2)
        self.ax.set_aspect("equal", "box")
        self.ax.axis('off')
        drawer2d.draw(self, mol)

    def to_bytesio(self):
        buf = io.BytesIO()
        self.fig.savefig(buf, dpi=self.dpi)  # resume scale
        return buf

    def save(self, path):
        self.fig.savefig(path, dpi=self.dpi)  # resume scale

    def get_size(self):
        return tuple(int(s * self.dpi) for s in self.fig.get_size_inches())

    def line_conv(method):
        @functools.wraps(method)
        def wrapper(self, p1, p2, c1, c2=None):
            p1_ = tuple(c * self.scale_factor for c in p1[:2])
            p2_ = tuple(c * self.scale_factor for c in p2[:2])
            # rgb in the range of 0 to 1
            c1_ = tuple(1 - ((255 - c) / 255) for c in c1)
            if c2:
                c2_ = tuple(1 - ((255 - c) / 255) for c in c2)
                method(self, p1_, p2_, c1_, c2_)
            else:
                method(self, p1_, p2_, c1_)
        return wrapper

    def text_conv(method):
        @functools.wraps(method)
        def wrapper(self, p, text, c, align):
            p_ = tuple(c * self.scale_factor for c in p[:2])
            c_ = tuple(1 - ((255 - c) / 255) for c in c)
            method(self, p_, text, c_, align)
        return wrapper

    def _draw_line_base(self, p1, p2, c1, c2=None, ls="solid"):
        if c2 and c2 != c1:
            pm = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
            px, py = [tuple(p) for p in zip(p1, pm)]
            self.ax.add_line(mpl.lines.Line2D(
                px, py, c=c1, lw=1.5, ls=ls, clip_on=False
            ))
            px, py = [tuple(p) for p in zip(pm, p2)]
            self.ax.add_line(mpl.lines.Line2D(
                px, py, c=c2, lw=1.5, ls=ls, clip_on=False
            ))
        else:
            px, py = [tuple(p) for p in zip(p1, p2)]
            self.ax.add_line(mpl.lines.Line2D(
                px, py, c=c1, lw=1.5, ls=ls, clip_on=False
            ))

    @line_conv
    def draw_line(self, p1, p2, c1, c2=None):
        self._draw_line_base(p1, p2, c1, c2)

    @line_conv
    def draw_dashed_line(self, p1, p2, c1, c2=None):
        self._draw_line_base(p1, p2, c1, c2, ls="dashed")

    @line_conv
    def draw_wedge(self, head, tail, color):
        width = self.displaymlb * 0.2
        self.ax.add_patch(mpl.patches.FancyArrowPatch(
            tail, head, ec=color, fc=color, clip_on=False,
            arrowstyle="wedge,tail_width={},shrink_factor=0.5".format(width)
        ))

    @line_conv
    def draw_dashed_wedge(self, head, tail, color):
        width = self.displaymlb * 0.2
        t1 = m_seg(tail, head, math.pi / -2, width / 2)[0]
        t2 = m_seg(tail, head, math.pi / 2, width / 2)[0]
        step_num = 7
        for i in range(step_num):
            step = distance(head, tail) / step_num * i
            trim = 1 / step_num * i
            tp1, tp2 = p_seg(t1, t2, True, step, trim)
            tpx, tpy = [tuple(p) for p in zip(tp1, tp2)]
            self.ax.add_line(mpl.lines.Line2D(
                tpx, tpy, c=color, lw=1.5, clip_on=False
            ))

    @line_conv
    def draw_wave_line(self, pos1, pos2, color):
        width = self.displaymlb * 0.2
        s1 = m_seg(pos1, pos2, math.pi / -2, width / 2)
        s2 = m_seg(pos1, pos2, math.pi / 2, width / 2)
        codes = [mpl.path.Path.MOVETO]
        vers = [pos1]
        for i, s in enumerate((s1, s2, s1, s2, s1, s2)):
            p = scale(s[1], (i + 1) / 7, s[0])
            codes.append(mpl.path.Path.LINETO)
            vers.append(p)
        codes.append(mpl.path.Path.LINETO)
        vers.append(pos2)
        path = mpl.path.Path(vers, codes)
        self.ax.add_patch(mpl.patches.PathPatch(
            path, fc="None", ec=color, lw=1.5, clip_on=False
        ))

    @text_conv
    def draw_text(self, pos, text, color, align="center"):
        replaced = substitute({
            '<sub>': '$_{',
            '<sup>': '$^{',
            '</sub>': '}$',
            '</sup>': '}$'
        }, text)
        xcorr = pos[0]
        if align == "left":
            xcorr -= self.font_size  # TODO: empirical
        elif align == "right":
            xcorr += self.font_size  # TODO: empirical
        self.ax.text(
            xcorr, pos[1], replaced, color=color, clip_on=False,
            fontproperties=basefont(self.font_size), ha=align, va="center"
        )
