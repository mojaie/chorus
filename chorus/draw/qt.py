#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

import io
import math

from PIL import Image
from PySide import QtCore, QtGui, QtSvg

from chorus.draw.drawermixin2d import DrawerMixin2D
from chorus.draw import helper
from chorus.util.geometry import distance, scale, m_seg, p_seg

MIN_LOGICAL_SIZE = 250
MARGIN = 40      # Output image margin
LINE_WEIGHT = 1  # Line weight(pixel)
RGB_TMPL = 'rgb({},{},{})'


class Qt(DrawerMixin2D):
    """Deprecated(no longer updated)"""
    def __init__(self):
        self.screen_size = (250, 250)

        self.contents = None
        self.logical_size = None
        self.painter = None
        # Get GUI instance
        app = QtGui.QApplication.instance()
        if app is None:
            QtGui.QApplication(())

    def draw(self, compound):
        helper.ready_to_draw(compound)
        size = max([compound.size2d[0] + MARGIN,
                    compound.size2d[1] + MARGIN,
                    MIN_LOGICAL_SIZE])
        self.logical_size = (size, size)
        self.contents = QtGui.QPixmap(*self.logical_size)
        self.painter = QtGui.QPainter(self.contents)
        self.painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
        self.painter.setRenderHint(QtGui.QPainter.SmoothPixmapTransform, True)
        self.painter.fillRect(0, 0, size, size, QtCore.Qt.white)
        super().draw(compound)
        self.painter.end()

    def load_svg(self, str_):
        xml = QtCore.QXmlStreamReader(str_)
        svg = QtSvg.QSvgRenderer(xml)
        w, h = svg.viewBoxF().size().toTuple()
        s = max((w, h,  MIN_LOGICAL_SIZE))
        self.contents = QtGui.QPixmap(s, s)
        painter = QtGui.QPainter(self.contents)
        painter.fillRect(0, 0, s, s, QtCore.Qt.white)
        svg.render(painter, QtCore.QRectF((s-w) / 2, (s-h) / 2, w, h))
        painter.end()
        # self.contents = self.contents.scaledToHeight(
        #    self.screen_size[1], QtCore.Qt.SmoothTransformation)

    def save_png_bytes_io(self):
        """ Save contents as PNG format file(BytesIO object)
        Note: DPI is platform dependent due to QPaintDevice DPI handling.
        Mac -> 72dpi, Ubuntu -> 96dpi
        """
        bytearr = QtCore.QByteArray()
        buf = QtCore.QBuffer(bytearr)
        buf.open(QtCore.QIODevice.WriteOnly)
        self.contents.save(buf, "PNG")
        bio = io.BytesIO(bytearr.data())
        # DPI correction
        img = Image.open(bio)
        img = img.resize((self.screen_size[0], self.screen_size[1]),
                         Image.ANTIALIAS)
        img.info["dpi"] = (96, 96)
        converted = io.BytesIO()
        img.save(converted, format='png')
        return converted

    def save_png(self, filename):
        self.contents.save(filename, "PNG")

    def widget(self):
        label = QtGui.QLabel()
        self.contents = self.contents.scaledToHeight(
            self.screen_size, QtCore.Qt.SmoothTransformation)
        label.setPixmap(self.contents)
        return label

    def _draw_line(self, p1, p2, c1=(0, 0, 0), c2=None,
                   op=QtCore.Qt.SolidLine):
        qp1 = QtCore.QPointF(*self._convert(p1))
        qp2 = QtCore.QPointF(*self._convert(p2))
        if c2 and c2 != c1:
            qpm = QtCore.QPointF(
                (qp1.x() + qp2.x()) / 2, (qp1.y() + qp2.y()) / 2)
            pen = QtGui.QPen(QtGui.QColor(*c1), LINE_WEIGHT, op)
            self.painter.setPen(pen)
            self.painter.drawLine(qp1, qpm)
            pen.setColor(QtGui.QColor(*c2))
            self.painter.setPen(pen)
            self.painter.drawLine(qpm, qp2)
        else:
            pen = QtGui.QPen(QtGui.QColor(*c1), LINE_WEIGHT, op)
            self.painter.setPen(pen)
            self.painter.drawLine(qp1, qp2)

    def _draw_dashed_line(self, p1, p2, c1=(0, 0, 0), c2=None):
        self._draw_line(p1, p2, c1, c2, QtCore.Qt.DashLine)

    def _draw_wave_line(self, pos1, pos2, width, color=(0, 0, 0)):
        cpos1 = self._convert(pos1)
        cpos2 = self._convert(pos2)
        b1 = m_seg(cpos1, cpos2, math.pi / -2, width / 2)
        b2 = m_seg(cpos1, cpos2, math.pi / 2, width / 2)
        path = QtGui.QPainterPath()
        path.moveTo(*cpos1)
        for i, b in enumerate([b1, b2, b1, b2, b1, b2]):
            p = scale(b[1], (i + 1) / 7, b[0])
            path.lineTo(*p)
        path.lineTo(*cpos2)
        pen = QtGui.QPen(QtGui.QColor(*color), 1, QtCore.Qt.SolidLine)
        self.painter.setPen(pen)
        self.painter.drawPath(path)

    def _draw_wedge(self, head, tail, width, color=(0, 0, 0)):
        chead = self._convert(head)
        ctail = self._convert(tail)
        b1 = m_seg(ctail, chead, math.pi / -2, width / 2)[0]
        b2 = m_seg(ctail, chead, math.pi / 2, width / 2)[0]
        ps = (chead, b1, b2)
        polygon = QtGui.QPolygonF()
        for ver in ps:
            polygon.append(QtCore.QPointF(*ver))
        pen = QtGui.QPen(QtGui.QColor(*color), 1, QtCore.Qt.SolidLine)
        self.painter.setPen(pen)
        self.painter.setBrush(QtGui.QColor(*color))
        self.painter.drawPolygon(polygon)

    def _draw_dashed_wedge(self, head, tail, width, color=(0, 0, 0)):
        chead = self._convert(head)
        ctail = self._convert(tail)
        b1 = m_seg(ctail, chead, math.pi / -2, width / 2)[0]
        b2 = m_seg(ctail, chead, math.pi / 2, width / 2)[0]
        pen = QtGui.QPen(QtGui.QColor(*color), 1, QtCore.Qt.SolidLine)
        self.painter.setPen(pen)
        step_num = 6
        # TODO: adjust step number
        for i in range(step_num):
            step = distance(chead, ctail) / step_num * i
            trim = 1 / step_num * i
            bp1, bp2 = p_seg(b1, b2, True, step, trim)
            qp1 = QtCore.QPointF(*bp1)
            qp2 = QtCore.QPointF(*bp2)
            self.painter.drawLine(qp1, qp2)

    def _draw_text(self, pos, text, size, color=(0, 0, 0), align="center"):
        qfont = QtGui.QFont("Helvetica", size)
        qtext = QtGui.QTextDocument()
        qtext.setDefaultFont(qfont)
        html_format = "<span style='color:rgb({},{},{})'>{}</span>"
        colored = list(color) + [text, ]
        qtext.setHtml(html_format.format(*colored))
        cpos = self._convert(pos)
        offset = size / 1.5
        if align == "center":
            qpos = QtCore.QPointF(
                cpos[0] - qtext.idealWidth() / 2, cpos[1] - offset)
        elif align == "right":
            qpos = QtCore.QPointF(
                cpos[0] - qtext.idealWidth() + offset, cpos[1] - offset)
        else:
            qpos = QtCore.QPointF(cpos[0] - offset, cpos[1] - offset)
        self.painter.save()
        self.painter.translate(qpos)
        qtext.drawContents(self.painter)
        self.painter.restore()

    def _convert(self, pos):
        """ For QPainter coordinate system, reflect over X axis and
        translate from center to top-left
        """
        px = pos[0] + self.logical_size.width() / 2
        py = self.logical_size.height() / 2 - pos[1]
        return px, py
