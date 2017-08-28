#
# (C) 2014-2017 Seiji Matsuoka
# Licensed under the MIT License (MIT)
# http://opensource.org/licenses/MIT
#

class Drawable(object):
    def draw_line(self, p1, p2, c1=(0, 0, 0), c2=None):
        """Draw line segment for bond.

        Args:
            p1 (x, y): the start point
            p2 (x, y): the end point
            c1 (R, G, B): RGB color(0-255) of p1 side half
            c2 (R, G, B): RGB color(0-255) of p2 side half.
                if c2 is undefined, c1 color will be used instead.
        """
        raise NotImplementedError()

    def draw_dashed_line(self, p1, p2, c1=(0, 0, 0), c2=None):
        """Draw dashed line segment for stereobond.

        Args:
            p1 (x, y): the start point
            p2 (x, y): the end point
            c1 (R, G, B): RGB color(0-255) of p1 side half
            c2 (R, G, B): RGB color(0-255) of p2 side half.
                if c2 is undefined, c1 color will be used instead.
        """
        raise NotImplementedError()

    def draw_wave_line(self, p1, p2, color=(0, 0, 0)):
        """Draw wave line segment for stereobond.

        Args:
            p1 (x, y): the start point
            p2 (x, y): the end point
            color (R, G, B): RGB color(0-255)
        """
        raise NotImplementedError()

    def draw_wedge(self, head, tail, color=(0, 0, 0)):
        """Draw wedge (isoscales triangle) for stereobond.

        Args:
            head (x, y): apex of the triangle
            tail (x, y): midpoint of the triangle base
            width (float): triangle base width
            color (R, G, B): RGB color(0-255)
        """
        raise NotImplementedError()

    def draw_dashed_wedge(self, head, tail, color=(0, 0, 0)):
        """Draw dashed wedge (isoscales triangle) for stereobond.

        Args:
            head (x, y): apex of the triangle
            tail (x, y): midpoint of the triangle base
            width (float): triangle base width
            color (R, G, B): RGB color(0-255)
        """
        raise NotImplementedError()

    def draw_text(self, pos, text, color=(0, 0, 0), align="center"):
        """Draw text for atom symbol.

        Args:
            pos (x, y): position of the text center
            text (str): contents
            color (R, G, B): RGB color(0-255)
            align ("right", "center" or "left"): text anchor position
        """
        raise NotImplementedError()
