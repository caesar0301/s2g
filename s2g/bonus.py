#!/usr/env/bin python
# encoding: utf-8

__all__ = ['perpend_to_line']


def perpend_to_line(p1, p2, p3):
    """Return the perpendicular line of a point to a line segment
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    px = x2 - x1
    py = y2 - y1
    sqr = float(px * px + py * py)
    if sqr == 0:
        return x1, y1
    u = ((x3 - x1) * px + (y3 - y1) * py) / sqr
    if u > 1:
        u = 1
    elif u < 0:
        u = 0
    x = x1 + u * px
    y = y1 + u * py
    return x, y
